// Copyright (C) 2015-2025 Severin Strobl <git@severin-strobl.de>
//
// SPDX-License-Identifier: MIT
//
// Exact calculation of the overlap volume of spheres and mesh elements.
// http://dx.doi.org/10.1016/j.jcp.2016.02.003

#ifndef OVERLAP_HPP
#define OVERLAP_HPP

// Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

// C++
#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>

#if !defined(overlap_assert)
#include <cassert>

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage,bugprone-macro-parentheses,readability-identifier-naming)
#define overlap_assert(expr, msg) assert((expr) && msg)
#endif

namespace overlap {

namespace detail {

// type aliases
using Scalar = double;

using Vector2 = Eigen::Matrix<Scalar, 2, 1, Eigen::DontAlign>;
using Vector3 = Eigen::Matrix<Scalar, 3, 1, Eigen::DontAlign>;
using Vector = Vector3;

// constants
const auto pi = Scalar{4} * std::atan(Scalar{1});

static constexpr auto tiny_epsilon =
    Scalar{4} * std::numeric_limits<Scalar>::epsilon();  // 4 ulp at 1.0
static constexpr auto medium_epsilon = Scalar{1e2} * tiny_epsilon;
static constexpr auto large_epsilon = Scalar{1e-10};

// Robust calculation of the normal vector of a polygon using Newell's method
// and a pre-calculated center.
// Ref: Christer Ericson - Real-Time Collision Detection (2005)
template<typename Iterator>
auto normal_newell(Iterator first, Iterator last, const Vector& center)
    -> Vector {
  const auto count = std::distance(first, last);
  auto normal = Vector::Zero().eval();
  for (auto i = decltype(count){0}; i < count; ++i) {
    normal +=
        (*(first + i) - center).cross(*(first + ((i + 1) % count)) - center);
  }

  const auto scale = normal.cwiseAbs().maxCoeff();
  if (const auto length = normal.stableNorm();
      length > scale * std::numeric_limits<Scalar>::epsilon()) {
    return normal / length;
  }

  return normal;
}

template<typename T>
class DoublePrecision {
  // This implementation of software double-precision is based on:
  // - T.J. Dekker, A floating-point technique for extending the available
  //   precision, https://doi.org/10.1007/BF01397083
  //
  // - J.-M. Muller, Elementary Functions - Algorithms and Implementation,
  //   https://doi.org/10.1007/978-1-4899-7983-4
  //
  // - X.S. Li et al., Design, implementation and testing of extended and mixed
  //   precision BLAS, https://doi.org/10.1145/567806.567808

  static_assert(std::is_floating_point_v<T>, "floating-point type required");

 public:
  constexpr DoublePrecision() = default;

  constexpr explicit DoublePrecision(const T value) {
    // std::tuple<Types...>::operator= not constexpr in C++17 -> no std::tie
    const auto [h, l] = split(value);
    high_ = h;
    low_ = l;
  }

  constexpr explicit DoublePrecision(const std::pair<T, T>& parts) :
      high_{parts.first}, low_{parts.second} {}

 private:
  constexpr DoublePrecision(const T h, const T l) : high_{h}, low_{l} {}

 public:
  // constant according to Veltkamp & Dekker: 2^(p - floor(p / 2)) + 1
  [[nodiscard]] static constexpr auto constant() -> T {
    constexpr auto digits =
        static_cast<std::uint32_t>(std::numeric_limits<T>::digits);

    return static_cast<T>((2u << (digits - digits / 2u - 1u)) + 1u);
  }

  [[nodiscard]] static constexpr auto split(const T value) -> std::pair<T, T> {
    const auto t = constant() * value;
    const auto h = t - (t - value);

    return std::make_pair(h, value - h);
  }

  [[nodiscard]] static constexpr auto fast_two_sum(const T x, const T y) {
    const auto s = x + y;
    const auto e = y - (s - x);

    return DoublePrecision{s, e};
  }

  [[nodiscard]] static constexpr auto two_sum(const T x, const T y) {
    const auto s = x + y;
    const auto v = s - x;
    const auto e = (x - (s - v)) + (y - v);

    return DoublePrecision{s, e};
  }

  [[nodiscard]] static constexpr auto two_product(const T x, const T y)
      -> DoublePrecision {
#if defined(__FMA__) || (defined(_MSC_VER) && defined(__AVX2__))
    const auto p = x * y;
    const auto e = x * y - p;
#else  // FMA not supported
    const auto p = x * y;
    const auto [x_h, x_l] = split(x);
    const auto [y_h, y_l] = split(y);
    const auto e = ((x_h * y_h - p) + x_h * y_l + x_l * y_h) + x_l * y_l;
#endif

    return DoublePrecision{p, e};
  }

  friend constexpr auto operator+(const DoublePrecision& lhs,
                                  const DoublePrecision& rhs)
      -> DoublePrecision {
    const auto s = two_sum(lhs.high(), rhs.high());
    const auto t = two_sum(lhs.low(), rhs.low());
    const auto v = fast_two_sum(s.high(), s.low() + t.high());

    return fast_two_sum(v.high(), v.low() + t.low());
  }

  friend constexpr auto operator-(const DoublePrecision& lhs,
                                  const DoublePrecision& rhs)
      -> DoublePrecision {
    return lhs + DoublePrecision{-rhs.high(), -rhs.low()};
  }

  friend constexpr auto operator*(const DoublePrecision& lhs,
                                  const DoublePrecision& rhs)
      -> DoublePrecision {
    const auto c = two_product(lhs.high(), rhs.high());
    const auto cc = (lhs.high() * rhs.low() + lhs.low() * rhs.high()) + c.low();

    return fast_two_sum(c.high(), cc);
  }

  [[nodiscard]] constexpr auto high() const -> T { return high_; }
  [[nodiscard]] constexpr auto low() const -> T { return low_; }
  [[nodiscard]] constexpr auto value() const -> T { return high_ + low_; }

  template<typename TOther>
  [[nodiscard]] constexpr auto as() const -> TOther {
    return TOther{high_} + TOther{low_};
  }

 private:
  T high_{0};  // Dekker: head/tail
  T low_{0};   // Dekker: tail
};

// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline auto orient2d(const Vector2& a, const Vector2& b, const Vector2& c)
    -> Scalar {
  using DoublePrecScalar = DoublePrecision<Vector2::Scalar>;

  const auto ax = DoublePrecScalar{a.x()};
  const auto ay = DoublePrecScalar{a.y()};
  const auto bx = DoublePrecScalar{b.x()};
  const auto by = DoublePrecScalar{b.y()};
  const auto cx = DoublePrecScalar{c.x()};
  const auto cy = DoublePrecScalar{c.y()};

  const auto result = (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);

  return result.template as<Vector2::Scalar>();
}

// Numerically robust calculation of the normal of the triangle defined by
// the points a, b, and c.
// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline auto triangle_normal(const Vector& a, const Vector& b, const Vector& c)
    -> Vector {
  const auto xy = orient2d({a.x(), a.y()}, {b.x(), b.y()}, {c.x(), c.y()});
  const auto yz = orient2d({a.y(), a.z()}, {b.y(), b.z()}, {c.y(), c.z()});
  const auto zx = orient2d({a.z(), a.x()}, {b.z(), b.x()}, {c.z(), c.x()});

  return Vector{yz, zx, xy}.normalized();
}

// Numerically robust routine to calculate the angle between normalized
// vectors.
// Ref: http://www.plunk.org/~hatch/rightway.html
inline auto angle(const Vector& u, const Vector& v) -> Scalar {
  if (u.dot(v) < Scalar{0}) {
    return pi - Scalar{2} * std::asin(Scalar{0.5} * (-v - u).stableNorm());
  }

  return Scalar{2} * std::asin(Scalar{0.5} * (v - u).stableNorm());
}

// Orthonormalize two unit vectors using the Gram–Schmidt process, returning
// two orthogonal unit vectors.
inline auto gram_schmidt(const Vector& v0, const Vector& v1)
    -> std::array<Vector, 2> {
  overlap_assert(std::abs(v0.norm() - Vector::Scalar{1}) < tiny_epsilon,
                 "vector v0 must be normalized");
  overlap_assert(std::abs(v1.norm() - Vector::Scalar{1}) < tiny_epsilon,
                 "vector v1 must be normalized");

  return {v0, (v1 - v1.dot(v0) * v0).normalized()};
}

inline auto clamp(Scalar value, Scalar min, Scalar max,
                  Scalar tolerance = Scalar{0}) -> Scalar {
  overlap_assert(min <= max && tolerance >= Scalar{0},
                 "invalid arguments for clamp()");

  value = (value < min && value > (min - tolerance)) ? min : value;
  value = (value > max && value < (max + tolerance)) ? max : value;

  return value;
}

struct Transformation {
  Transformation() = default;

  template<typename Derived>
  Transformation(const Eigen::MatrixBase<Derived>& t, Scalar s) :
      translation{t.eval()}, scaling{s} {}

  Vector translation = Vector::Zero();
  Scalar scaling = Scalar{1};
};

template<std::size_t VertexCount>
class Polygon {
  static_assert(VertexCount >= 3 && VertexCount <= 4,
                "only triangles and quadrilateral supported");

 public:
  static constexpr std::size_t vertex_count = VertexCount;

  Polygon() = default;

  template<typename... Types>
  explicit constexpr Polygon(const Vector& v0, Types... verts) :
      Polygon{std::array<Vector, VertexCount>{v0, verts...}} {}

  explicit Polygon(std::array<Vector, VertexCount> verts) :
      vertices(std::move(verts)) {
    center = (Scalar{1} / Scalar{VertexCount}) *
             std::accumulate(vertices.begin(), vertices.end(),
                             Vector::Zero().eval());

    // For a quadrilateral, Newell's method can be simplified significantly.
    // Ref: Christer Ericson - Real-Time Collision Detection (2005)
    if constexpr (VertexCount == 4) {
      normal = ((vertices[2] - vertices[0]).cross(vertices[3] - vertices[1]))
                   .normalized();

    } else {
      normal = detail::normal_newell(vertices.begin(), vertices.end(), center);
    }

    update_area();
  }

  void apply(const Transformation& t) {
    for (auto& v : vertices) {
      v = t.scaling * (v + t.translation);
    }

    center = t.scaling * (center + t.translation);

    update_area();
  }

  [[nodiscard]] auto is_planar(const Scalar tolerance = large_epsilon) const
      -> bool {
    if constexpr (VertexCount == 3U) {
      return true;
    }

    return std::all_of(std::begin(vertices), std::end(vertices),
                       [&](const Vector& v) {
                         return std::abs(normal.dot(v - center)) <= tolerance;
                       });
  }

 private:
  void update_area() {
    if constexpr (VertexCount == 4) {
      area = Scalar{0.5} *
             (((vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]))
                  .stableNorm() +
              ((vertices[2] - vertices[0]).cross(vertices[3] - vertices[0]))
                  .stableNorm());
    } else {
      area = Scalar{0.5} *
             ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]))
                 .stableNorm();
    }
  }

 public:
  std::array<Vector, VertexCount> vertices = {};
  Vector center = Vector::Zero();
  Vector normal = Vector::Identity();
  Scalar area = Scalar{0};
};

using Triangle = Polygon<3>;
using Quadrilateral = Polygon<4>;

// Forward declarations of the mesh elements.
class Tetrahedron;
class Wedge;
class Hexahedron;

// Some tricks are required to keep this code header-only.
template<typename T, typename Nil>
struct mappings;

template<typename Nil>
struct mappings<Tetrahedron, Nil> {
  // Map edges of a tetrahedron to vertices and faces.
  static const uint32_t edge_mapping[6][2][2];

  // Map vertices of a tetrahedron to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[4][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

template<typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::edge_mapping[6][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 0}, {0, 3}},
    {{0, 3}, {1, 3}}, {{1, 3}, {1, 2}}, {{2, 3}, {2, 3}}};

template<typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::vertex_mapping[4][3][3] = {
    {{0, 2, 3}, {0, 1, 0}, {0, 1, 3}},
    {{0, 1, 4}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 5}, {1, 0, 0}, {0, 2, 3}},
    {{3, 4, 5}, {1, 1, 1}, {1, 3, 2}}};

template<typename Nil>
const uint32_t mappings<Tetrahedron, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

using tet_mappings = mappings<Tetrahedron, void>;

class Tetrahedron : public tet_mappings {
 public:
  Tetrahedron() {
    std::fill(std::begin(vertices), std::end(vertices), Vector::Zero());
  }

  template<typename... Types>
  Tetrahedron(const Vector& v0, Types... verts) : vertices{{v0, verts...}} {
    // Make sure the ordering of the vertices is correct.
    overlap_assert((vertices[1] - vertices[0])
                           .cross(vertices[2] - vertices[0])
                           .dot(vertices[3] - vertices[0]) >= Scalar{0},
                   "invalid vertex order detected");

    init();
  }

  Tetrahedron(std::array<Vector, 4> verts) : vertices(std::move(verts)) {
    init();
  }

  void apply(const Transformation& t) {
    for (auto& v : vertices) {
      v = t.scaling * (v + t.translation);
    }

    for (auto& f : faces) {
      f.apply(t);
    }

    center = Scalar{0.25} * std::accumulate(vertices.begin(), vertices.end(),
                                            Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto surface_area() const -> Scalar {
    return std::accumulate(
        std::begin(faces), std::end(faces), Scalar{0},
        [](auto sum, const auto& face) { return sum + face.area; });
  }

 private:
  void init() {
    // 0: v2, v1, v0
    faces[0] = Triangle{vertices[2], vertices[1], vertices[0]};

    // 1: v0, v1, v3
    faces[1] = Triangle{vertices[0], vertices[1], vertices[3]};

    // 2: v1, v2, v3
    faces[2] = Triangle{vertices[1], vertices[2], vertices[3]};

    // 3: v2, v0, v3
    faces[3] = Triangle{vertices[2], vertices[0], vertices[3]};

    center = Scalar{0.25} * std::accumulate(vertices.begin(), vertices.end(),
                                            Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto calc_volume() const -> Scalar {
    return (Scalar{1} / Scalar{6}) *
           std::abs((vertices[0] - vertices[3])
                        .dot((vertices[1] - vertices[3])
                                 .cross(vertices[2] - vertices[3])));
  }

 public:
  std::array<Vector, 4> vertices;
  std::array<Triangle, 4> faces = {};
  Vector center = Vector::Zero();
  Scalar volume = Scalar{0};
};

template<typename Nil>
struct mappings<Wedge, Nil> {
  // Map edges of a wedge to vertices and faces.
  static const uint32_t edge_mapping[9][2][2];

  // Map vertices of a wedge to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[6][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

template<typename Nil>
const uint32_t mappings<Wedge, Nil>::edge_mapping[9][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}}, {{2, 0}, {0, 3}},
    {{0, 3}, {1, 3}}, {{1, 4}, {1, 2}}, {{2, 5}, {2, 3}},
    {{3, 4}, {1, 4}}, {{4, 5}, {2, 4}}, {{5, 3}, {3, 4}}};

// clang-format off
template<typename Nil>
const uint32_t mappings<Wedge, Nil>::vertex_mapping[6][3][3] = {
    {{0, 2, 3}, {0, 1, 0}, {0, 1, 3}},
    {{0, 1, 4}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 5}, {1, 0, 0}, {0, 2, 3}},

    {{3, 6, 8}, {1, 0, 1}, {1, 3, 4}},
    {{4, 6, 7}, {1, 1, 0}, {1, 2, 4}},
    {{5, 7, 8}, {1, 1, 0}, {2, 3, 4}}};
// clang-format on

template<typename Nil>
const uint32_t mappings<Wedge, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

using wedge_mappings = mappings<Wedge, void>;

class Wedge : public wedge_mappings {
 public:
  Wedge() {
    std::fill(std::begin(vertices), std::end(vertices), Vector::Zero());
  }

  template<typename... Types>
  Wedge(const Vector& v0, Types... verts) : vertices{{v0, verts...}} {
    init();
  }

  Wedge(std::array<Vector, 6> verts) : vertices{std::move(verts)} { init(); }

  void apply(const Transformation& t) {
    for (auto& v : vertices) {
      v = t.scaling * (v + t.translation);
    }

    for (auto& f : faces) {
      f.apply(t);
    }

    center = (Scalar{1} / Scalar{6}) * std::accumulate(vertices.begin(),
                                                       vertices.end(),
                                                       Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto surface_area() const -> Scalar {
    return std::accumulate(
        std::begin(faces), std::end(faces), Scalar{0},
        [](auto sum, const auto& face) { return sum + face.area; });
  }

 private:
  void init() {
    // All faces of the wedge are stored as quadrilaterals, so an
    // additional point is inserted between v0 and v1.
    // 0: v2, v1, v0, v02
    faces[0] = Quadrilateral(vertices[2], vertices[1], vertices[0],
                             Scalar{0.5} * (vertices[0] + vertices[2]));

    // 1: v0, v1, v4, v3
    faces[1] =
        Quadrilateral(vertices[0], vertices[1], vertices[4], vertices[3]);

    // 2: v1, v2, v5, v4
    faces[2] =
        Quadrilateral(vertices[1], vertices[2], vertices[5], vertices[4]);

    // 3: v2, v0, v3, v5
    faces[3] =
        Quadrilateral(vertices[2], vertices[0], vertices[3], vertices[5]);

    // All faces of the wedge are stored as quadrilaterals, so an
    // additional point is inserted between v3 and v5.
    // 4: v3, v4, v5, v53
    faces[4] = Quadrilateral(vertices[3], vertices[4], vertices[5],
                             Scalar{0.5} * (vertices[5] + vertices[3]));

    center = (Scalar{1} / Scalar{6}) * std::accumulate(vertices.begin(),
                                                       vertices.end(),
                                                       Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto calc_volume() const -> Scalar {
    // The wedge is treated as a degenerate hexahedron here by adding
    // two fake vertices v02 and v35.
    const Vector diagonal = vertices[5] - vertices[0];

    return (Scalar{1} / Scalar{6}) *
           (diagonal.dot(
               ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[4])) +
               ((vertices[3] - vertices[0])
                    .cross(vertices[4] -
                           Scalar{0.5} * (vertices[3] + vertices[5]))) +
               ((Scalar{0.5} * (vertices[0] + vertices[2]) - vertices[0])
                    .cross(Scalar{0.5} * (vertices[3] + vertices[5]) -
                           vertices[2]))));
  }

 public:
  std::array<Vector, 6> vertices;
  std::array<Quadrilateral, 5> faces = {};
  Vector center = Vector::Zero();
  Scalar volume = Scalar{0};
};

template<typename Nil>
struct mappings<Hexahedron, Nil> {
  // Map edges of a hexahedron to vertices and faces.
  static const uint32_t edge_mapping[12][2][2];

  // Map vertices of a hexahedron to edges and faces.
  // 0: local IDs of the edges intersecting at this vertex
  // 1: 0 if the edge is pointing away from the vertex, 1 otherwise
  // 2: faces joining at the vertex
  static const uint32_t vertex_mapping[8][3][3];

  // This mapping contains the three sets of the two edges for each of the
  // faces joining at a vertex. The indices are mapped to the local edge IDs
  // using the first value field of the 'vertex_mapping' table.
  static const uint32_t face_mapping[3][2];
};

// clang-format off
template<typename Nil>
const uint32_t mappings<Hexahedron, Nil>::edge_mapping[12][2][2] = {
    {{0, 1}, {0, 1}}, {{1, 2}, {0, 2}},
    {{2, 3}, {0, 3}}, {{3, 0}, {0, 4}},

    {{0, 4}, {1, 4}}, {{1, 5}, {1, 2}},
    {{2, 6}, {2, 3}}, {{3, 7}, {3, 4}},

    {{4, 5}, {1, 5}}, {{5, 6}, {2, 5}},
    {{6, 7}, {3, 5}}, {{7, 4}, {4, 5}}};

template<typename Nil>
const uint32_t mappings<Hexahedron, Nil>::vertex_mapping[8][3][3] = {
    {{0, 3, 4}, {0, 1, 0}, {0, 1, 4}},
    {{0, 1, 5}, {1, 0, 0}, {0, 1, 2}},
    {{1, 2, 6}, {1, 0, 0}, {0, 2, 3}},
    {{2, 3, 7}, {1, 0, 0}, {0, 3, 4}},

    {{4, 8, 11}, {1, 0, 1}, {1, 4, 5}},
    {{5, 8, 9}, {1, 1, 0}, {1, 2, 5}},
    {{6, 9, 10}, {1, 1, 0}, {2, 3, 5}},
    {{7, 10, 11}, {1, 1, 0}, {3, 4, 5}}};
// clang-format on

template<typename Nil>
const uint32_t mappings<Hexahedron, Nil>::face_mapping[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

using hex_mappings = mappings<Hexahedron, void>;

class Hexahedron : public hex_mappings {
 public:
  Hexahedron() {
    std::fill(std::begin(vertices), std::end(vertices), Vector::Zero());
  }

  template<typename... Types>
  Hexahedron(const Vector& v0, Types... verts) : vertices{{v0, verts...}} {
    init();
  }

  Hexahedron(std::array<Vector, 8> verts) : vertices{std::move(verts)} {
    init();
  }

  void apply(const Transformation& t) {
    for (auto& v : vertices) {
      v = t.scaling * (v + t.translation);
    }

    for (auto& f : faces) {
      f.apply(t);
    }

    center = (Scalar{1} / Scalar{8}) * std::accumulate(vertices.begin(),
                                                       vertices.end(),
                                                       Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto surface_area() const -> Scalar {
    return std::accumulate(
        std::begin(faces), std::end(faces), Scalar{0},
        [](auto sum, const auto& face) { return sum + face.area; });
  }

 private:
  void init() {
    // 0: v3, v2, v1, v0
    faces[0] =
        Quadrilateral(vertices[3], vertices[2], vertices[1], vertices[0]);

    // 1: v0, v1, v5, v4
    faces[1] =
        Quadrilateral(vertices[0], vertices[1], vertices[5], vertices[4]);

    // 2: v1, v2, v6, v5
    faces[2] =
        Quadrilateral(vertices[1], vertices[2], vertices[6], vertices[5]);

    // 3: v2, v3, v7, v6
    faces[3] =
        Quadrilateral(vertices[2], vertices[3], vertices[7], vertices[6]);

    // 4: v3, v0, v4, v7
    faces[4] =
        Quadrilateral(vertices[3], vertices[0], vertices[4], vertices[7]);

    // 5: v4, v5, v6, v7
    faces[5] =
        Quadrilateral(vertices[4], vertices[5], vertices[6], vertices[7]);

    center = (Scalar{1} / Scalar{8}) * std::accumulate(vertices.begin(),
                                                       vertices.end(),
                                                       Vector::Zero().eval());

    volume = calc_volume();
  }

  [[nodiscard]] auto calc_volume() const -> Scalar {
    const auto diagonal = vertices[6] - vertices[0];

    return (Scalar{1} / Scalar{6}) *
           diagonal.dot(
               ((vertices[1] - vertices[0]).cross(vertices[2] - vertices[5])) +
               ((vertices[4] - vertices[0]).cross(vertices[5] - vertices[7])) +
               ((vertices[3] - vertices[0]).cross(vertices[7] - vertices[2])));
  }

 public:
  std::array<Vector, 8> vertices = {};
  std::array<Quadrilateral, 6> faces = {};
  Vector center = Vector::Zero();
  Scalar volume = Scalar{0};
};

template<typename Element>
constexpr auto num_vertices() -> std::size_t {
  return std::tuple_size_v<decltype(Element::vertices)>;
}

template<typename Element>
constexpr auto num_faces() -> std::size_t {
  return std::tuple_size_v<decltype(Element::faces)>;
}

template<typename Element>
constexpr auto num_edges() -> std::size_t {
  if constexpr (std::is_same_v<Element, Hexahedron>) {
    return 12U;
  }

  if constexpr (std::is_same_v<Element, Wedge>) {
    return 9U;
  }

  if constexpr (std::is_same_v<Element, Tetrahedron>) {
    return 6U;
  }

  // older versions of GCC cannot handle exceptions in constexpr contexts
  return std::numeric_limits<std::size_t>::max();
}

template<typename T>
struct is_element
    : public std::integral_constant<bool, std::is_same_v<T, Tetrahedron> ||
                                              std::is_same_v<T, Wedge> ||
                                              std::is_same_v<T, Hexahedron>> {};

template<typename T>
inline constexpr bool is_element_v = is_element<T>::value;

class Sphere {
 public:
  Sphere() : Sphere{Vector::Zero(), Scalar{1}} {}

  Sphere(Vector c, Scalar r) :
      center{std::move(c)},
      radius{r},
      volume{((Scalar{4} / Scalar{3}) * pi) * r * r * r} {}

  [[nodiscard]] auto cap_volume(Scalar height) const -> Scalar {
    if (height <= Scalar{0}) {
      return Scalar{0};
    }

    if (height >= Scalar{2} * radius) {
      return volume;
    }

    return (pi / Scalar{3}) * height * height * (Scalar{3} * radius - height);
  }

  [[nodiscard]] auto cap_surface_area(Scalar height) const -> Scalar {
    if (height <= Scalar{0}) {
      return Scalar{0};
    }

    if (height >= Scalar{2} * radius) {
      return surface_area();
    }

    return (Scalar{2} * pi) * radius * height;
  }

  [[nodiscard]] auto disk_area(Scalar height) const -> Scalar {
    if (height <= Scalar{0} || height >= Scalar{2} * radius) {
      return Scalar{0};
    }

    return pi * height * (Scalar{2} * radius - height);
  }

  [[nodiscard]] auto surface_area() const -> Scalar {
    return (Scalar{4} * pi) * (radius * radius);
  }

  Vector center;
  Scalar radius;
  Scalar volume;
};

struct Plane {
  Vector center;
  Vector normal;
};

// Decomposition of a tetrahedron into 4 tetrahedra.
inline void decompose(const Tetrahedron& tet,
                      std::array<Tetrahedron, 4>& tets) {
  tets[0] = Tetrahedron(tet.vertices[0], tet.vertices[1], tet.vertices[2],
                        tet.center);

  tets[1] = Tetrahedron(tet.vertices[0], tet.vertices[1], tet.center,
                        tet.vertices[3]);

  tets[2] = Tetrahedron(tet.vertices[1], tet.vertices[2], tet.center,
                        tet.vertices[3]);

  tets[3] = Tetrahedron(tet.vertices[2], tet.vertices[0], tet.center,
                        tet.vertices[3]);
}

// Decomposition of a hexahedron into 2 wedges.
inline void decompose(const Hexahedron& hex, std::array<Wedge, 2>& wedges) {
  wedges[0] = Wedge(hex.vertices[0], hex.vertices[1], hex.vertices[2],
                    hex.vertices[4], hex.vertices[5], hex.vertices[6]);

  wedges[1] = Wedge(hex.vertices[0], hex.vertices[2], hex.vertices[3],
                    hex.vertices[4], hex.vertices[6], hex.vertices[7]);
}

// Decomposition of a hexahedron into 5 tetrahedra.
inline void decompose(const Hexahedron& hex, std::array<Tetrahedron, 5>& tets) {
  tets[0] = Tetrahedron(hex.vertices[0], hex.vertices[1], hex.vertices[2],
                        hex.vertices[5]);

  tets[1] = Tetrahedron(hex.vertices[0], hex.vertices[2], hex.vertices[7],
                        hex.vertices[5]);

  tets[2] = Tetrahedron(hex.vertices[0], hex.vertices[2], hex.vertices[3],
                        hex.vertices[7]);

  tets[3] = Tetrahedron(hex.vertices[0], hex.vertices[5], hex.vertices[7],
                        hex.vertices[4]);

  tets[4] = Tetrahedron(hex.vertices[2], hex.vertices[7], hex.vertices[5],
                        hex.vertices[6]);
}

// Decomposition of a hexahedron into 6 tetrahedra.
inline void decompose(const Hexahedron& hex, std::array<Tetrahedron, 6>& tets) {
  tets[0] = Tetrahedron(hex.vertices[0], hex.vertices[5], hex.vertices[7],
                        hex.vertices[4]);

  tets[1] = Tetrahedron(hex.vertices[0], hex.vertices[1], hex.vertices[7],
                        hex.vertices[5]);

  tets[2] = Tetrahedron(hex.vertices[1], hex.vertices[6], hex.vertices[7],
                        hex.vertices[5]);

  tets[3] = Tetrahedron(hex.vertices[0], hex.vertices[7], hex.vertices[2],
                        hex.vertices[3]);

  tets[4] = Tetrahedron(hex.vertices[0], hex.vertices[7], hex.vertices[1],
                        hex.vertices[2]);

  tets[5] = Tetrahedron(hex.vertices[1], hex.vertices[7], hex.vertices[6],
                        hex.vertices[2]);
}

inline auto contains(const Sphere& s, const Vector& p) -> bool {
  return (s.center - p).squaredNorm() <= s.radius * s.radius;
}

// The (convex!) polygon is assumed to be planar, making this a 2D problem.
// Check the projection of the point onto the plane of the polygon for
// containment within the polygon.
template<std::size_t VertexCount>
auto contains(const Polygon<VertexCount>& poly, const Vector& point) -> bool {
  const Vector proj =
      point - poly.normal.dot(point - poly.center) * poly.normal;

  for (std::size_t n = 0; n < poly.vertices.size(); ++n) {
    const auto& v0 = poly.vertices[n];
    const auto& v1 = poly.vertices[(n + 1) % poly.vertices.size()];

    // Note: Only the sign of the projection is of interest, so this vector
    // does not have to be normalized.
    const auto dir = (v1 - v0).cross(poly.normal);

    // Check whether the projection of the point lies inside of the
    // polygon.
    if (dir.dot(proj - Scalar{0.5} * (v0 + v1)) > Scalar{0}) {
      return false;
    }
  }

  return true;
}

template<typename Element, typename = std::enable_if_t<is_element_v<Element>>>
inline auto contains(const Element& element, const Vector& p) -> bool {
  return std::all_of(std::begin(element.faces), std::end(element.faces),
                     [&](const auto& face) -> bool {
                       return face.normal.dot(p - face.center) <= Scalar{0};
                     });
}

template<typename Element, typename = std::enable_if_t<is_element_v<Element>>>
inline auto contains(const Sphere& sphere, const Element& element) -> bool {
  return std::all_of(std::begin(element.vertices), std::end(element.vertices),
                     [&](const Vector& vertex) {
                       return (sphere.center - vertex).squaredNorm() <=
                              sphere.radius * sphere.radius;
                     });
}

inline auto intersects(const Sphere& s, const Plane& p) -> bool {
  const auto proj = p.normal.dot(s.center - p.center);

  return proj * proj - s.radius * s.radius < Scalar{0};
}

template<std::size_t VertexCount>
inline auto intersects(const Sphere& s, const Polygon<VertexCount>& poly)
    -> bool {
  return intersects(s, {poly.center, poly.normal}) && contains(poly, s.center);
}

template<typename Element, typename = std::enable_if_t<is_element_v<Element>>>
inline auto intersects_coarse(const Sphere& sphere, const Element& element)
    -> bool {
  using AABB = Eigen::AlignedBox<Scalar, 3>;

  const auto sphere_aabb =
      AABB{sphere.center - Vector::Constant(sphere.radius),
           sphere.center + Vector::Constant(sphere.radius)};

  auto element_aabb = AABB{};
  for (const auto& v : element.vertices) {
    element_aabb.extend(v);
  }

  return sphere_aabb.intersects(element_aabb);
}

inline auto line_sphere_intersection(const Vector& base,
                                     const Vector& direction,
                                     const Sphere& sphere)
    -> std::array<std::optional<Scalar>, 2> {
  const auto a = direction.squaredNorm();
  if (a == Scalar{0}) {
    return std::array<std::optional<Scalar>, 2>{};
  }

  const auto origin_relative = (base - sphere.center).eval();
  const auto b = Scalar{2} * direction.dot(origin_relative);
  const auto c = origin_relative.squaredNorm() - sphere.radius * sphere.radius;

  const auto discriminant = b * b - Scalar{4} * a * c;
  if (discriminant > Scalar{0}) {
    // two real roots
    const auto q =
        Scalar{-0.5} * (b + std::copysign(std::sqrt(discriminant), b));

    auto x1 = q / a;
    auto x2 = c / q;

    if (x1 > x2) {
      std::swap(x1, x2);
    }

    // if the midpoint of the two intersection points is not inside the sphere,
    // the intersection is spurious
    if (((base + (Scalar{0.5} * (x1 + x2)) * direction) - sphere.center)
            .squaredNorm() >= sphere.radius * sphere.radius) {
      return std::array<std::optional<Scalar>, 2>{};
    }

    return std::array{std::make_optional(x1), std::make_optional(x2)};
  }

  if (std::abs(discriminant) == Scalar{0}) {
    // double real root, line tangential to the sphere's surface
    return std::array{std::make_optional((Scalar{-0.5} * b) / a),
                      std::optional<Scalar>{}};
  }

  // no real roots
  return std::array<std::optional<Scalar>, 2>{};
}

// Calculate the volume of a regularized spherical wedge defined by the
// radius, the distance of the intersection point from the center of the
// sphere and the angle.
inline auto regularized_wedge(Scalar r, Scalar d, Scalar alpha) -> Scalar {
#ifndef NDEBUG
  // clamp slight deviations of the angle to valid range
  if (alpha < Scalar{0} && alpha > -detail::tiny_epsilon) {
    alpha = Scalar{0};
  }

  if (alpha > Scalar{0.5} * pi && alpha <= Scalar{0.5} * pi + tiny_epsilon) {
    alpha = Scalar{0.5} * pi;
  }
#endif

  overlap_assert(r > Scalar{0}, "invalid argument 'r' for regularized_wedge()");

  overlap_assert(d >= Scalar{0} && d <= r,
                 "invalid argument 'd' for regularized_wedge()");

  overlap_assert(alpha >= Scalar{0} && alpha <= Scalar{0.5} * pi,
                 "invalid argument 'alpha' for regularized_wedge()");

  const auto sin_alpha = std::sin(alpha);
  const auto cos_alpha = std::cos(alpha);

  const auto a = d * sin_alpha;
  const auto b = std::sqrt(std::abs(r * r - d * d));
  const auto c = d * cos_alpha;

  return (Scalar{1} / Scalar{3}) * a * b * c +
         a * ((Scalar{1} / Scalar{3}) * a * a - r * r) * std::atan2(b, c) +
         (Scalar{2} / Scalar{3}) * r * r * r *
             std::atan2(sin_alpha * b, cos_alpha * r);
}

// Wrapper around the above function handling correctly handling the case of
// alpha > pi/2 and negative z.
inline auto regularized_wedge(Scalar r, Scalar d, Scalar alpha, Scalar z)
    -> Scalar {
  if (z >= Scalar{0}) {
    if (alpha > Scalar{0.5} * pi) {
      const auto h = r - z;

      return (pi / Scalar{3}) * h * h * (Scalar{3} * r - h) -
             regularized_wedge(r, d, pi - alpha);
    }

    return regularized_wedge(r, d, alpha);
  }

  const auto hemisphere_volume = ((Scalar{2} / Scalar{3}) * pi) * r * r * r;
  if (alpha > Scalar{0.5} * pi) {
    return hemisphere_volume - regularized_wedge(r, d, pi - alpha);
  }

  const auto h = r + z;
  const auto cap_volume = (pi / Scalar{3}) * h * h * (Scalar{3} * r - h);

  return hemisphere_volume - (cap_volume - regularized_wedge(r, d, alpha));
}

// Calculate the surface area of a regularized spherical wedge defined by the
// radius, the distance of the intersection point from the center of the
// sphere and the angle.
// Ref: Gibson, K. D. & Scheraga, H. A.: Exact calculation of the volume and
//    surface area of fused hard-sphere molecules with unequal atomic radii,
//    Molecular Physics, 1987, 62, 1247-1265
inline auto regularized_wedge_area(Scalar r, Scalar z, Scalar alpha) -> Scalar {
#ifndef NDEBUG
  // clamp slight deviations of the angle to valid range
  if (alpha < Scalar{0} && alpha > -tiny_epsilon) {
    alpha = Scalar{0};
  }

  if (alpha > pi && alpha <= pi + tiny_epsilon) {
    alpha = pi;
  }
#endif

  overlap_assert(r > Scalar{0},
                 "invalid argument 'r' for regularized_wedge_area()");

  overlap_assert(z >= -r && z <= r,
                 "invalid argument 'z' for regularized_wedge_area()");

  overlap_assert(alpha >= Scalar{0} && alpha <= pi,
                 "invalid argument 'alpha' for regularized_wedge_area()");

  if (alpha < tiny_epsilon || std::abs(r * r - z * z) <= tiny_epsilon) {
    return Scalar{0};
  }

  const auto sin_alpha = std::sin(alpha);
  const auto cos_alpha = std::cos(alpha);
  const auto factor = Scalar{1} / std::sqrt(std::abs(r * r - z * z));

  // clamp slight deviations of the argument to acos() to valid range
  const auto arg0 =
      clamp(r * cos_alpha * factor, Scalar{-1}, Scalar{1}, tiny_epsilon);

  const auto arg1 = clamp((z * cos_alpha * factor) / sin_alpha, Scalar{-1},
                          Scalar{1}, tiny_epsilon);

  overlap_assert(Scalar{-1} <= arg0 && arg0 <= Scalar{1},
                 "invalid value for arg0 in regularized_wedge_area()");

  overlap_assert(Scalar{-1} <= arg1 && arg1 <= Scalar{1},
                 "invalid value for arg1 in regularized_wedge_area()");

  return Scalar{2} * r * (r * std::acos(arg0) - z * std::acos(arg1));
}

// calculate the volume of the spherical wedge or the area of the spherical
// lune, depending on the dimensionality
template<std::size_t Dim>
inline auto spherical_wedge(const Sphere& s, const Scalar angle) -> Scalar {
  static_assert(Dim == 2 || Dim == 3, "invalid dimensionality, must be 2 or 3");

  if constexpr (Dim == 2) {
    return Scalar{2} * s.radius * s.radius * angle;
  }

  return (Scalar{2} / Scalar{3}) * s.radius * s.radius * s.radius * angle;
}

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
template<std::size_t Dim>
inline auto general_wedge(const Sphere& s, const Plane& p0, const Plane& p1,
                          const Vector& d) -> Scalar {
  static_assert(Dim == 2 || Dim == 3, "invalid dimensionality, must be 2 or 3");

  const auto dist = d.stableNorm();
  if (dist < tiny_epsilon) {
    // the wedge (almost) touches the center, the volume/area depends only on
    // the angle
    return spherical_wedge<Dim>(s, pi - angle(p0.normal, p1.normal));
  }

  if (dist >= s.radius) {
    // intersection of the two planes (numerically) on the surface of the
    // sphere
    return Scalar{0};
  }

  const auto s0 = d.dot(p0.normal);
  const auto s1 = d.dot(p1.normal);

  // detect degenerated general spherical wedge that can be treated as
  // a regularized spherical wedge
  if (std::abs(s0) < tiny_epsilon || std::abs(s1) < tiny_epsilon) {
    const auto alpha = pi - angle(p0.normal, p1.normal);

    if constexpr (Dim == 2) {
      return regularized_wedge_area(
          s.radius, std::abs(s0) > std::abs(s1) ? s0 : s1, alpha);
    }

    return regularized_wedge(s.radius, dist, alpha,
                             std::abs(s0) > std::abs(s1) ? s0 : s1);
  }

  auto d_unit = Vector{d * (Scalar{1} / dist)};
  if (dist < large_epsilon) {
    d_unit =
        gram_schmidt(p0.normal.cross(p1.normal).stableNormalized(), d_unit)[1];
  }

  overlap_assert(p0.normal.dot(p1.center - p0.center) <= Scalar{0},
                 "invalid plane in general_wedge()");

  overlap_assert(p1.normal.dot(p0.center - p1.center) <= Scalar{0},
                 "invalid plane in general_wedge()");

  // calculate the angles between the vector from the sphere center
  // to the intersection line and the normal vectors of the two planes
  auto alpha0 = angle(p0.normal, d_unit);
  auto alpha1 = angle(p1.normal, d_unit);

  const auto pi_half = Scalar{0.5} * pi;
  const auto dir0 = d_unit.dot((s.center + d) - p0.center);
  const auto dir1 = d_unit.dot((s.center + d) - p1.center);

  if (s0 >= Scalar{0} && s1 >= Scalar{0}) {
    alpha0 = pi_half - std::copysign(alpha0, dir0);
    alpha1 = pi_half - std::copysign(alpha1, dir1);

    if constexpr (Dim == 2) {
      return regularized_wedge_area(s.radius, s0, alpha0) +
             regularized_wedge_area(s.radius, s1, alpha1);
    }

    return regularized_wedge(s.radius, dist, alpha0, s0) +
           regularized_wedge(s.radius, dist, alpha1, s1);
  }

  if (s0 < Scalar{0} && s1 < Scalar{0}) {
    alpha0 = pi_half + std::copysign(Scalar{1}, dir0) * (alpha0 - pi);
    alpha1 = pi_half + std::copysign(Scalar{1}, dir1) * (alpha1 - pi);

    if constexpr (Dim == 2) {
      return s.surface_area() - (regularized_wedge_area(s.radius, -s0, alpha0) +
                                 regularized_wedge_area(s.radius, -s1, alpha1));
    }

    return s.volume - (regularized_wedge(s.radius, dist, alpha0, -s0) +
                       regularized_wedge(s.radius, dist, alpha1, -s1));
  }

  alpha0 = pi_half - std::copysign(Scalar{1}, dir0 * s0) *
                         (alpha0 - (s0 < Scalar{0} ? pi : Scalar{0}));

  alpha1 = pi_half - std::copysign(Scalar{1}, dir1 * s1) *
                         (alpha1 - (s1 < Scalar{0} ? pi : Scalar{0}));

  if constexpr (Dim == 2) {
    const auto area0 = regularized_wedge_area(s.radius, std::abs(s0), alpha0);
    const auto area1 = regularized_wedge_area(s.radius, std::abs(s1), alpha1);

    return std::max(area0, area1) - std::min(area0, area1);
  }

  const auto volume0 = regularized_wedge(s.radius, dist, alpha0, std::abs(s0));
  const auto volume1 = regularized_wedge(s.radius, dist, alpha1, std::abs(s1));

  return std::max(volume0, volume1) - std::min(volume0, volume1);
}

template<typename Element>
using EdgeIntersections =
    std::array<std::optional<std::array<Vector, 2>>, num_edges<Element>()>;

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
template<std::size_t Dim, typename Element>
auto general_wedge(const Sphere& sphere, const Element& element,
                   std::size_t edge,
                   const EdgeIntersections<Element>& intersections) {
  static_assert(Dim == 2 || Dim == 3, "invalid dimensionality, must be 2 or 3");

  const auto& f0 = element.faces[Element::edge_mapping[edge][1][0]];
  const auto& f1 = element.faces[Element::edge_mapping[edge][1][1]];

  overlap_assert(intersections[edge].has_value(),
                 "inconsistent intersection detection for edge");

  const auto edge_midpoint = Vector{
      Scalar{0.5} * (((*intersections[edge])[0] +
                      element.vertices[Element::edge_mapping[edge][0][0]]) +
                     ((*intersections[edge])[1] +
                      element.vertices[Element::edge_mapping[edge][0][1]]))};

  const auto p0 = Plane{f0.center, f0.normal};
  const auto p1 = Plane{f1.center, f1.normal};

  return general_wedge<Dim>(sphere, p0, p1, edge_midpoint - sphere.center);
}

// if not all three edges intersecting at a vertex are marked, the
// sphere is only touching this vertex
template<typename Element, typename = std::enable_if_t<is_element_v<Element>>>
inline auto correct_marked_vertices(
    const std::bitset<num_vertices<Element>()>& marked_vertices,
    const std::bitset<num_edges<Element>()>& marked_edges)
    -> std::bitset<num_vertices<Element>()> {
  auto corrected_marked_vertices = marked_vertices;

  for (auto vertex_idx = 0u; vertex_idx < marked_vertices.size();
       ++vertex_idx) {
    if (!marked_vertices[vertex_idx]) {
      continue;
    }

    auto all_edges_marked = true;
    for (auto local_edge_idx = 0u; local_edge_idx < 3u; ++local_edge_idx) {
      const auto edge_idx =
          Element::vertex_mapping[vertex_idx][0][local_edge_idx];
      all_edges_marked &= static_cast<bool>(marked_edges[edge_idx]);
    }

    corrected_marked_vertices[vertex_idx] = all_edges_marked;
  }

  return corrected_marked_vertices;
}

template<typename Element>
struct EntityIntersections {
  std::bitset<num_vertices<Element>()> vertices;
  std::bitset<num_edges<Element>()> edges;
  std::bitset<num_faces<Element>()> faces;
};

template<typename Element>
auto unit_sphere_intersections(const Element& element)
    -> std::tuple<EntityIntersections<Element>, EdgeIntersections<Element>> {
  const auto unit_sphere = Sphere{};

  auto entity_intersections = EntityIntersections<Element>{};

  // the intersection points between the single edges and the sphere are cached
  auto edge_intersections = EdgeIntersections<Element>{};

  for (auto edge_idx = 0u; edge_idx < num_edges<Element>(); ++edge_idx) {
    const auto base = element.vertices[Element::edge_mapping[edge_idx][0][0]];
    const auto direction =
        (element.vertices[Element::edge_mapping[edge_idx][0][1]] - base).eval();

    const auto&& intersections =
        line_sphere_intersection(base, direction, unit_sphere);

    // no intersection between the edge and the sphere, where touching
    // (tangential) contacts are ignored
    if (!intersections[1].has_value() || *intersections[0] >= Scalar{1} ||
        *intersections[1] <= Scalar{0}) {
      continue;
    }

    entity_intersections.vertices[Element::edge_mapping[edge_idx][0][0]] =
        entity_intersections.vertices[Element::edge_mapping[edge_idx][0][0]] ||
        *intersections[0] < Scalar{0};

    entity_intersections.vertices[Element::edge_mapping[edge_idx][0][1]] =
        entity_intersections.vertices[Element::edge_mapping[edge_idx][0][1]] ||
        *intersections[1] > Scalar{1};

    // note: the intersection points are relative to the vertices
    edge_intersections[edge_idx] =
        std::array<Vector, 2>{*intersections[0] * direction,
                              (*intersections[1] - Scalar{1}) * direction};

    entity_intersections.edges[edge_idx] = true;

    // if the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well
    entity_intersections.faces[Element::edge_mapping[edge_idx][1][0]] = true;
    entity_intersections.faces[Element::edge_mapping[edge_idx][1][1]] = true;
  }

  // check whether the dependencies for a vertex intersection are fulfilled
  entity_intersections.vertices = correct_marked_vertices<Element>(
      entity_intersections.vertices, entity_intersections.edges);

  // check the interior of all faces for intersection with the unit sphere
  for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
    if (intersects(unit_sphere, element.faces[face_idx])) {
      entity_intersections.faces[face_idx] = true;
    }
  }

  return std::make_tuple(entity_intersections, edge_intersections);
}

template<std::size_t Dim, typename Element>
auto vertex_cone_correction(
    const Element& element,
    const EdgeIntersections<Element>& edge_intersections,
    const std::size_t vertex_idx) -> Scalar {
  static_assert(Dim == 2 || Dim == 3,
                "invalid dimension for computation of correction at vertex");

  // Collect the points where the three edges intersecting at this
  // vertex intersect the sphere.
  // Both the relative and the absolute positions are required.
  auto relative_intersection_points = std::array<Vector, 3>{};
  auto intersection_points = std::array<Vector, 3>{};
  for (auto local_edge_idx = 0u; local_edge_idx < 3; ++local_edge_idx) {
    const auto edge_idx =
        Element::vertex_mapping[vertex_idx][0][local_edge_idx];

    overlap_assert(edge_intersections[edge_idx].has_value(),
                   "inconsistent intersection detection for edge");

    relative_intersection_points[local_edge_idx] =
        (*edge_intersections[edge_idx])
            [Element::vertex_mapping[vertex_idx][1][local_edge_idx]];

    intersection_points[local_edge_idx] =
        relative_intersection_points[local_edge_idx] +
        element.vertices[vertex_idx];
  }

  // This triangle is constructed by hand to have more freedom of how
  // the normal vector is calculated.
  auto cone_triangle = Triangle{};
  cone_triangle.vertices = {
      {intersection_points[0], intersection_points[1], intersection_points[2]}};

  cone_triangle.center =
      (Scalar{1} / Scalar{3}) * std::accumulate(intersection_points.begin(),
                                                intersection_points.end(),
                                                Vector::Zero().eval());

  // Calculate the normal of the triangle defined by the intersection
  // points in relative coordinates to improve accuracy.
  // Also use double the normal precision to calculate this normal.
  cone_triangle.normal = triangle_normal(relative_intersection_points[0],
                                         relative_intersection_points[1],
                                         relative_intersection_points[2]);

  // The area of this triangle is never needed, so it is set to an
  // invalid value.
  cone_triangle.area = std::numeric_limits<Scalar>::infinity();

  auto distances = std::array<std::pair<std::size_t, Scalar>, 3>{};
  for (auto i = 0u; i < 3; ++i) {
    distances[i] =
        std::make_pair(i, relative_intersection_points[i].squaredNorm());
  }

  std::sort(distances.begin(), distances.end(),
            [](const std::pair<std::size_t, Scalar>& a,
               const std::pair<std::size_t, Scalar>& b) {
              return a.second < b.second;
            });

  const auto unit_sphere = Sphere{};

  if (distances[1].second < distances[2].second * detail::large_epsilon) {
    // Use the general spherical wedge defined by the edge with the
    // non-degenerated intersection point and the normals of the
    // two faces forming it.
    return general_wedge<Dim, Element>(
        unit_sphere, element,
        Element::vertex_mapping[vertex_idx][0][distances[2].first],
        edge_intersections);
  }

  // Make sure the normal points in the right direction i.e. away from
  // the center of the element.
  if (cone_triangle.normal.dot(element.center - cone_triangle.center) >
      Scalar{0}) {
    cone_triangle.normal = -cone_triangle.normal;
  }

  // Calculate the volume/surface area of the three spherical segments between
  // the faces joining at the vertex and the plane through the intersection
  // points.
  auto segment_correction = [&]() {
    const auto plane = Plane{cone_triangle.center, cone_triangle.normal};
    const auto local_face_indices = {0u, 1u, 2u};
    return std::accumulate(
        std::begin(local_face_indices), std::end(local_face_indices), Scalar{0},
        [&](const Scalar initial, const auto local_face_idx) {
          const auto& face =
              element.faces[Element::vertex_mapping[vertex_idx][2]
                                                   [local_face_idx]];

          const auto center =
              (Scalar{0.5} *
               (intersection_points[Element::face_mapping[local_face_idx][0]] +
                intersection_points[Element::face_mapping[local_face_idx][1]]))
                  .eval();

          return initial + general_wedge<Dim>(unit_sphere, plane,
                                              Plane{face.center, -face.normal},
                                              center);
        });
  };

  const auto dist = cone_triangle.normal.dot(-cone_triangle.center);

  if constexpr (Dim == 2) {
    const auto cap_surface =
        unit_sphere.cap_surface_area(unit_sphere.radius + dist);

    // If cap surface area is small, the corrections will be even smaller. There
    // is no way to actually calculate them with reasonable precision, so they
    // are just ignored.
    if (cap_surface < large_epsilon) {
      return Scalar{0};
    }

    // Calculate the surface area of the cone and clamp it to zero.
    return std::max(cap_surface - segment_correction(), Scalar{0});
  } else if constexpr (Dim == 3) {
    const auto tip_tet_volume =
        (Scalar{1} / Scalar{6}) *
        std::abs(-relative_intersection_points[2].dot(
            (relative_intersection_points[0] - relative_intersection_points[2])
                .cross(relative_intersection_points[1] -
                       relative_intersection_points[2])));

    const auto cap_volume = unit_sphere.cap_volume(unit_sphere.radius + dist);

    // The cap volume is tiny, so the corrections will be even smaller. There is
    // no way to actually calculate them with reasonable precision, so just the
    // volume of the tetrahedron at the tip is used.
    if (cap_volume < tiny_epsilon) {
      return tip_tet_volume;
    }

    // Calculate the volume of the cone and clamp it to zero.
    return std::max(tip_tet_volume + cap_volume - segment_correction(),
                    Scalar{0});
  }

  return Scalar{};
}

template<typename Element>
inline auto detect_non_planar_faces(const Element& element) -> void {
  for (const auto& face : element.faces) {
    if (!face.is_planar()) {
      throw std::invalid_argument{"non-planer face detected in element"};
    }
  }
}

// normalize the element w.r.t the unit sphere
template<typename Element>
inline auto normalize_element(const Sphere& sphere, const Element& element)
    -> Element {
  const auto transformation =
      Transformation{-sphere.center, Scalar{1} / sphere.radius};

  auto transformed_element = Element{element};
  transformed_element.apply(transformation);

  return transformed_element;
}

}  // namespace detail

// expose types required for public API
using Vector = detail::Vector;
using Scalar = detail::Scalar;

using Sphere = detail::Sphere;
using Tetrahedron = detail::Tetrahedron;
using Wedge = detail::Wedge;
using Hexahedron = detail::Hexahedron;

template<typename Element>
auto overlap_volume(const Sphere& sphere, const Element& element) -> Scalar {
  using namespace detail;

  static_assert(is_element_v<Element>, "invalid element type detected");

  if (!intersects_coarse(sphere, element)) {
    return Scalar{0};
  }

  // check for trivial case: element fully contained in sphere
  if (contains(sphere, element)) {
    return element.volume;
  }

  // sanity check: all faces of the mesh element have to be planar
  detect_non_planar_faces(element);

  // use unit sphere and transformed (scaled and shifted) version of the element
  const auto unit_sphere = Sphere{};
  auto transformed_element = normalize_element(sphere, element);

  const auto&& [entity_intersections, edge_intersections] =
      unit_sphere_intersections(transformed_element);

  // trivial case: the center of the sphere overlaps the element, but the sphere
  // does not intersect any of the faces of the element, meaning the sphere is
  // completely contained within the element
  if (!entity_intersections.faces.count() &&
      contains(transformed_element, unit_sphere.center)) {
    return sphere.volume;
  }

  // spurious intersection: The initial intersection test was positive, but the
  // detailed checks revealed no overlap
  if (!entity_intersections.vertices.count() &&
      !entity_intersections.edges.count() &&
      !entity_intersections.faces.count()) {
    return Scalar{0};
  }

  // initial value: volume of the full sphere
  auto result = unit_sphere.volume;

  // iterate over all the marked faces and subtract the volume of the cap cut
  // off by the plane
  for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
    if (!entity_intersections.faces[face_idx]) {
      continue;
    }

    const auto& face = transformed_element.faces[face_idx];
    const auto dist = face.normal.dot(-face.center);

    result -= unit_sphere.cap_volume(unit_sphere.radius + dist);
  }

  // handle the edges and add back the volume subtracted twice above in the
  // processing of the faces
  for (auto edge_idx = 0u; edge_idx < num_edges<Element>(); ++edge_idx) {
    if (!entity_intersections.edges[edge_idx]) {
      continue;
    }

    result += general_wedge<3, Element>(unit_sphere, transformed_element,
                                        edge_idx, edge_intersections);
  }

  // handle the vertices and subtract the volume added twice above in the
  // processing of the edges
  for (auto vertex_idx = 0u; vertex_idx < num_vertices<Element>();
       ++vertex_idx) {
    if (!entity_intersections.vertices[vertex_idx]) {
      continue;
    }

    result -= vertex_cone_correction<3>(transformed_element, edge_intersections,
                                        vertex_idx);

    // sanity check: detect negative intermediate result
    overlap_assert(result > -std::sqrt(detail::tiny_epsilon),
                   "negative intermediate result in overlap_volume()");
  }

  // in case of different sized objects, the error can become quite large, so a
  // relative limit is used
  const auto max_overlap_volume =
      std::min(unit_sphere.volume, transformed_element.volume);

  const auto limit =
      std::sqrt(std::numeric_limits<Scalar>::epsilon()) * max_overlap_volume;

  // clamp tiny negative volumes to zero
  if (result < Scalar{0} && result > -limit) {
    return Scalar{0};
  }

  // clamp results slightly too large
  if (result > max_overlap_volume && result - max_overlap_volume < limit) {
    return std::min(sphere.volume, element.volume);
  }

  // perform a final sanity check on the final result (debug version only)
  overlap_assert(result >= Scalar{0} && result <= max_overlap_volume,
                 "negative volume detected in overlap_volume()");

  // scale the overlap volume back for the original objects
  result = (result / unit_sphere.volume) * sphere.volume;

  return result;
}

template<typename Iterator>
auto overlap_volume(const Sphere& s, Iterator first, Iterator last) -> Scalar {
  static_assert(
      detail::is_element_v<typename std::iterator_traits<Iterator>::value_type>,
      "invalid element type detected");

  return std::accumulate(first, last, Scalar{0},
                         [&s](const Scalar partial, const auto& element) {
                           return partial + overlap_volume(s, element);
                         });
}

// Calculate the surface area of the sphere and the element that are contained
// within the common or intersecting part of the geometries, respectively.
// The returned array of size (N + 2), with N being the number of vertices,
// holds (in this order):
//   - surface area of the region of the sphere intersecting the element
//   - for each face of the element: area contained within the sphere
//   - total surface area of the element intersecting the sphere
template<typename Element>
auto overlap_area(const Sphere& sphere, const Element& element)
    -> std::array<Scalar, detail::num_faces<Element>() + 2> {
  using namespace detail;

  static_assert(is_element_v<Element>, "invalid element type detected");

  // initial value: zero overlap
  auto result = std::array<Scalar, num_faces<Element>() + 2>{};

  if (!intersects_coarse(sphere, element)) {
    return result;
  }

  // check for trivial case: element fully contained in sphere resulting in a
  // full coverage of all faces
  if (contains(sphere, element)) {
    for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
      result[face_idx + 1] = element.faces[face_idx].area;
      result.back() += element.faces[face_idx].area;
    }

    return result;
  }

  // sanity check: all faces of the mesh element have to be planar
  detect_non_planar_faces(element);

  // use unit sphere and transformed (scaled and shifted) version of the element
  const auto unit_sphere = Sphere{};
  auto transformed_element = normalize_element(sphere, element);

  const auto&& [entity_intersections, edge_intersections] =
      unit_sphere_intersections(transformed_element);

  // trivial case: the center of the sphere overlaps the element, but the sphere
  // does not intersect any of the faces of the element, meaning the sphere is
  // completely contained within the element
  if (!entity_intersections.faces.count() &&
      contains(transformed_element, unit_sphere.center)) {
    result[0] = sphere.surface_area();

    return result;
  }

  // spurious intersection: ihe initial intersection test was positive, but the
  // detailed checks revealed no overlap
  if (!entity_intersections.vertices.count() &&
      !entity_intersections.edges.count() &&
      !entity_intersections.faces.count()) {
    return result;
  }

  // initial value for the surface of the sphere: Surface area of the full
  // sphere
  result[0] = unit_sphere.surface_area();

  // iterate over all the marked faces and calculate the area of the disk
  // defined by the plane as well as the cap surfaces
  for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
    if (!entity_intersections.faces[face_idx]) {
      continue;
    }

    const auto& face = transformed_element.faces[face_idx];
    const auto dist = face.normal.dot(-face.center);

    result[0] -= unit_sphere.cap_surface_area(unit_sphere.radius + dist);
    result[face_idx + 1] = unit_sphere.disk_area(unit_sphere.radius + dist);
  }

  const auto circular_segment_area = [](const Scalar radius_sq,
                                        const Scalar chord_length) {
    const auto apothem = std::sqrt(std::max(
        Scalar{0}, radius_sq - Scalar{0.25} * chord_length * chord_length));

    const auto theta =
        Scalar{2} * std::atan2(chord_length, Scalar{2} * apothem);

    return Scalar{0.5} * radius_sq * (theta - std::sin(theta));
  };

  // cache the squared radius of the disk formed by the intersection between the
  // planes defined by each face and the sphere
  auto intersection_radius_sq = std::array<Scalar, num_faces<Element>()>{};

  // handle the edges and subtract the area of the respective disk cut off by
  // the edge and add back the surface area of the spherical wedge defined by
  // the edge
  for (auto edge_idx = 0u; edge_idx < num_edges<Element>(); ++edge_idx) {
    if (!entity_intersections.edges[edge_idx]) {
      continue;
    }

    // intersection area of the the sphere: add back the surface area of the
    // spherical wedge defined by the edge which was considered twice when
    // processing the two faces forming the edge
    result[0] += general_wedge<2, Element>(unit_sphere, transformed_element,
                                           edge_idx, edge_intersections);

    // intersection areas of the all faces: for each of the two faces forming
    // the edge, remove the part of the disk area beyond the edge
    //

    overlap_assert(edge_intersections[edge_idx].has_value(),
                   "inconsistent intersection detection for edge");

    // the chord length is given by the distance of the intersection points of
    // the sphere and the edge
    const auto chord =
        ((transformed_element.vertices[Element::edge_mapping[edge_idx][0][0]] +
          (*edge_intersections[edge_idx])[0]) -
         (transformed_element.vertices[Element::edge_mapping[edge_idx][0][1]] +
          (*edge_intersections[edge_idx])[1]))
            .eval();

    const auto chord_length = chord.stableNorm();

    // each edge belongs to two faces, indexed via
    // Element::edge_mapping[n][1][{0,1}]
    for (auto local_face_idx = 0u; local_face_idx < 2u; ++local_face_idx) {
      const auto face_idx = Element::edge_mapping[edge_idx][1][local_face_idx];
      const auto& face = transformed_element.faces[face_idx];

      // height of the spherical cap cut off by the plane containing the face
      const auto cap_height = unit_sphere.radius - face.normal.dot(face.center);
      const auto apothem = unit_sphere.radius - cap_height;
      intersection_radius_sq[face_idx] =
          cap_height * (unit_sphere.radius + apothem);

      // the part of the base of the spherical cap cut off by the edge
      auto segment_area =
          circular_segment_area(intersection_radius_sq[face_idx], chord_length);

      const auto chord_center =
          (Scalar{0.5} *
           ((transformed_element
                 .vertices[Element::edge_mapping[edge_idx][0][0]] +
             (*edge_intersections[edge_idx])[0]) +
            (transformed_element
                 .vertices[Element::edge_mapping[edge_idx][0][1]] +
             (*edge_intersections[edge_idx])[1])))
              .eval();

      // projection of the center of the sphere onto the face
      const auto proj =
          (unit_sphere.center -
           face.normal.dot(unit_sphere.center - face.center) * face.normal)
              .eval();

      // if the projected sphere center and the face center fall on opposite
      // sides of the edge, the area has to be inverted
      const auto invert_segment_area =
          chord.cross(proj - chord_center)
              .dot(chord.cross(face.center - chord_center)) < Scalar{0};

      if (invert_segment_area) {
        segment_area = intersection_radius_sq[face_idx] * pi - segment_area;
      }

      result[face_idx + 1] -= segment_area;
    }
  }

  for (auto vertex_idx = 0u; vertex_idx < num_vertices<Element>();
       ++vertex_idx) {
    if (!entity_intersections.vertices[vertex_idx]) {
      continue;
    }

    //
    // correct the intersection area of the the sphere
    //
    result[0] -= vertex_cone_correction<2>(transformed_element,
                                           edge_intersections, vertex_idx);

    // sanity checks: detect negative/excessively large intermediate result
    overlap_assert(result[0] > -std::sqrt(detail::tiny_epsilon),
                   "negative area as intermediate result in overlap_area()");

    overlap_assert(
        result[0] < unit_sphere.surface_area() + detail::tiny_epsilon,
        "invalid intermediate result in overlap_area()");

    //
    // correct the intersection areas of all facets
    //
    // iterate over all the faces joining at this vertex
    for (auto local_face_idx = 0u; local_face_idx < 3u; ++local_face_idx) {
      // determine the two edges of this face intersecting at the vertex
      const auto edge0 = Element::face_mapping[local_face_idx][0];
      const auto edge1 = Element::face_mapping[local_face_idx][1];
      const auto edge_indices =
          std::array{Element::vertex_mapping[vertex_idx][0][edge0],
                     Element::vertex_mapping[vertex_idx][0][edge1]};

      // extract the (relative) intersection points of these edges with the
      // sphere furthest from the vertex
      overlap_assert(edge_intersections[edge_indices[0]].has_value() &&
                         edge_intersections[edge_indices[1]].has_value(),
                     "inconsistent intersection detection for edge");

      const auto intersection_points =
          std::array{(*edge_intersections[edge_indices[0]])
                         [Element::vertex_mapping[vertex_idx][1][edge0]],

                     (*edge_intersections[edge_indices[1]])
                         [Element::vertex_mapping[vertex_idx][1][edge1]]};

      // together with the vertex, this determines the triangle representing one
      // part of the correction
      const auto triangle_area =
          Scalar{0.5} *
          (intersection_points[0].cross(intersection_points[1])).stableNorm();

      // the second component is the segment defined by the face and the
      // intersection points
      const auto chord_length =
          (intersection_points[0] - intersection_points[1]).stableNorm();

      const auto face_idx =
          Element::vertex_mapping[vertex_idx][2][local_face_idx];

      auto segment_area =
          circular_segment_area(intersection_radius_sq[face_idx], chord_length);

      // determine if the (projected) center of the sphere lies within the
      // triangle or not; if not, the segment area has to be corrected
      const auto chord_center =
          (Scalar{0.5} * (intersection_points[0] + intersection_points[1]))
              .eval();

      const auto& face = transformed_element.faces[face_idx];
      const auto proj = (-face.normal.dot(-face.center) * face.normal).eval();
      const auto invert_segment_area =
          chord_center.dot((proj - transformed_element.vertices[vertex_idx]) -
                           chord_center) > Scalar{0};

      if (invert_segment_area) {
        segment_area = intersection_radius_sq[face_idx] * pi - segment_area;
      }

      result[face_idx + 1] += triangle_area + segment_area;

      // sanity checks: detect excessively large intermediate result
      overlap_assert(
          result[face_idx + 1] < transformed_element.faces[face_idx].area +
                                     std::sqrt(detail::large_epsilon),
          "invalid intermediate result in overlap_area()");
    }
  }

  // scale the surface areas back for the original objects and clamp values
  // within reasonable limits
  const auto scaling = sphere.radius;
  const auto sphere_limit = std::sqrt(std::numeric_limits<Scalar>::epsilon()) *
                            unit_sphere.surface_area();

  // as the precision of the area calculation deteriorates quickly with a
  // increasing size ratio between the element and the sphere, the precision
  // limit applied to the sphere is used as the lower limit for the facets
  const auto face_limit =
      std::max(sphere_limit, std::sqrt(std::numeric_limits<Scalar>::epsilon()) *
                                 transformed_element.surface_area());

  // sanity checks: detect negative/excessively large results for the surface
  // area of the facets
  for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
    overlap_assert(result[face_idx + 1] > -face_limit,
                   "negative overlap area for face in overlap_area()");

    overlap_assert(result[face_idx + 1] <=
                       transformed_element.faces[face_idx].area + face_limit,
                   "invalid overlap area for face in overlap_area()");
  }

  // surface of the sphere within the element
  result[0] = (scaling * scaling) * detail::clamp(result[0], Scalar{0},
                                                  unit_sphere.surface_area(),
                                                  sphere_limit);

  // surfaces of the mesh element within the sphere
  for (auto face_idx = 0u; face_idx < num_faces<Element>(); ++face_idx) {
    result[face_idx + 1] =
        (scaling * scaling) *
        detail::clamp(result[face_idx + 1], Scalar{0},
                      transformed_element.faces[face_idx].area, face_limit);
  }

  result.back() =
      std::accumulate(result.begin() + 1, result.end() - 1, Scalar{0});

  // perform final sanity checks on the result (debug version only)
  overlap_assert(Scalar{0} <= result[0] && result[0] <= sphere.surface_area(),
                 "invalid overlap area for sphere surface in overlap_area()");

  overlap_assert(
      Scalar{0} <= result.back() && result.back() <= element.surface_area(),
      "invalid total overlap area for faces in overlap_area()");

  return result;
}

}  // namespace overlap

#endif  // OVERLAP_HPP

// Copyright (C) 2015-2024 Severin Strobl <git@severin-strobl.de>
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
#include <iterator>
#include <limits>
#include <numeric>
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
    Scalar{2} * std::numeric_limits<Scalar>::epsilon();
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

// This implementation of double_prec is based on:
// T.J. Dekker, A floating-point technique for extending the available
// precision, http://dx.doi.org/10.1007/BF01397083

template<typename T>
struct double_prec_constant;

template<>
struct double_prec_constant<float> {
  // Constant used to split double precision values:
  // 2^(24 - 24/2) + 1 = 2^12 + 1 = 4097
  static const uint32_t value = 4097;
};

template<>
struct double_prec_constant<double> {
  // Constant used to split double precision values:
  // 2^(53 - int(53/2)) + 1 = 2^27 + 1 = 134217729
  static const uint32_t value = 134217729;
};

// For GCC an attribute has to be used to control the FP precision...
#define ENFORCE_EXACT_FPMATH_ATTR

// ... whereas ICC requires a pragma.
#if defined(__ICC) || defined(__INTEL_COMPILER)
#define ENFORCE_EXACT_FPMATH_ATTR
#define USE_EXACT_FPMATH_PRAGMA 1
#endif

template<typename T>
class double_prec;

template<typename T>
inline double_prec<T> operator+(const double_prec<T>& lhs,
                                const double_prec<T>& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

template<typename T>
inline double_prec<T> operator-(const double_prec<T>& lhs,
                                const double_prec<T>& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

template<typename T>
inline double_prec<T> operator*(const double_prec<T>& lhs,
                                const double_prec<T>& rhs)
    ENFORCE_EXACT_FPMATH_ATTR;

template<typename T>
class double_prec {
 private:
  static const uint32_t c = detail::double_prec_constant<T>::value;

  template<typename TF>
  friend double_prec<TF> operator+(const double_prec<TF>&,
                                   const double_prec<TF>&);

  template<typename TF>
  friend double_prec<TF> operator-(const double_prec<TF>&,
                                   const double_prec<TF>&);

  template<typename TF>
  friend double_prec<TF> operator*(const double_prec<TF>&,
                                   const double_prec<TF>&);

 public:
  inline double_prec() : h_(0), l_(0) {}

  // This constructor requires floating point operations in accordance
  // with IEEE754 to perform the proper splitting. To allow full
  // optimization of all other parts of the code, precise floating point
  // ops are only requested here. Unfortunately the way to do this is
  // extremely compiler dependent.
  inline double_prec(const T& val) ENFORCE_EXACT_FPMATH_ATTR : h_(0),
                                                               l_(0) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

    T p = val * T(c);
    h_ = (val - p) + p;
    l_ = val - h_;
  }

 private:
  inline explicit double_prec(const T& h, const T& l) : h_(h), l_(l) {}

 public:
  inline const T& high() const { return h_; }

  inline const T& low() const { return l_; }

  inline T value() const { return h_ + l_; }

  template<typename TOther>
  inline TOther convert() const {
    return TOther(h_) + TOther(l_);
  }

 private:
  T h_;
  T l_;
};

template<typename T>
inline double_prec<T> operator+(const double_prec<T>& lhs,
                                const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  T h = lhs.h_ + rhs.h_;
  T l = std::abs(lhs.h_) >= std::abs(rhs.h_)
            ? ((((lhs.h_ - h) + rhs.h_) + lhs.l_) + rhs.l_)
            : ((((rhs.h_ - h) + lhs.h_) + rhs.l_) + lhs.l_);

  T c = h + l;

  return double_prec<T>(c, (h - c) + l);
}

template<typename T>
inline double_prec<T> operator-(const double_prec<T>& lhs,
                                const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  T h = lhs.h_ - rhs.h_;
  T l = std::abs(lhs.h_) >= std::abs(rhs.h_)
            ? ((((lhs.h_ - h) - rhs.h_) - rhs.l_) + lhs.l_)
            : ((((-rhs.h_ - h) + lhs.h_) + lhs.l_) - rhs.l_);

  T c = h + l;

  return double_prec<T>(c, (h - c) + l);
}

template<typename T>
inline double_prec<T> operator*(const double_prec<T>& lhs,
                                const double_prec<T>& rhs) {
#ifdef USE_EXACT_FPMATH_PRAGMA
#pragma float_control(precise, on)
#endif

  double_prec<T> l(lhs.h_);
  double_prec<T> r(rhs.h_);

  T p = l.h_ * r.h_;
  T q = l.h_ * r.l_ + l.l_ * r.h_;
  T v = p + q;

  double_prec<T> c(v, ((p - v) + q) + l.l_ * r.l_);
  c.l_ = ((lhs.h_ + lhs.l_) * rhs.l_ + lhs.l_ * rhs.h_) + c.l_;
  T z = c.value();

  return double_prec<T>(z, (c.h_ - z) + c.l_);
}

// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline auto orient2d(const Vector2& a, const Vector2& b, const Vector2& c)
    -> Scalar {
  using DoublePrecScalar = double_prec<Vector2::Scalar>;

  const auto a0 = DoublePrecScalar{a[0]};
  const auto a1 = DoublePrecScalar{a[1]};
  const auto b0 = DoublePrecScalar{b[0]};
  const auto b1 = DoublePrecScalar{b[1]};
  const auto c0 = DoublePrecScalar{c[0]};
  const auto c1 = DoublePrecScalar{c[1]};

  const auto result = (a0 - c0) * (b1 - c1) - (a1 - c1) * (b0 - c0);

  return result.template convert<Vector2::Scalar>();
}

// Numerically robust calculation of the normal of the triangle defined by
// the points a, b, and c.
// Ref: J.R. Shewchuk - Lecture Notes on Geometric Robustness
//      http://www.cs.berkeley.edu/~jrs/meshpapers/robnotes.pdf
inline auto triangle_normal(const Vector& a, const Vector& b, const Vector& c)
    -> Vector {
  const auto xy = orient2d({a[0], a[1]}, {b[0], b[1]}, {c[0], c[1]});
  const auto yz = orient2d({a[1], a[2]}, {b[1], b[2]}, {c[1], c[2]});
  const auto zx = orient2d({a[2], a[0]}, {b[2], b[0]}, {c[2], c[0]});

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

class Plane {
 public:
  Plane(Vector c, Vector n) : center{std::move(c)}, normal{std::move(n)} {}

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
    const auto base = (Scalar{1} / Scalar{2}) * (v0 + v1);
    if (dir.dot(proj - base) > Scalar{0}) {
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

inline auto line_sphere_intersection(const Vector& origin,
                                     const Vector& direction, const Sphere& s)
    -> std::pair<std::array<Scalar, 2>, std::size_t> {
  std::array<Scalar, 2U> solutions = {
      {std::numeric_limits<Scalar>::infinity(),
       std::numeric_limits<Scalar>::infinity()}};

  const auto a = direction.squaredNorm();
  if (a == Scalar{0}) {
    return std::make_pair(solutions, 0U);
  }

  const auto origin_relative = Vector{origin - s.center};
  const auto b = Scalar{2} * direction.dot(origin_relative);
  const auto c = origin_relative.squaredNorm() - s.radius * s.radius;

  const auto discriminant = b * b - Scalar{4} * a * c;
  if (discriminant > Scalar{0}) {
    // two real roots
    const auto q = (Scalar{-1} / Scalar{2}) *
                   (b + std::copysign(std::sqrt(discriminant), b));

    solutions[0] = q / a;
    solutions[1] = c / q;

    if (solutions[0] > solutions[1]) {
      std::swap(solutions[0], solutions[1]);
    }

    return std::make_pair(solutions, 2U);
  }

  if (std::abs(discriminant) == Scalar{0}) {
    // double real root
    solutions[0] = ((Scalar{-1} / Scalar{2}) * b) / a;
    solutions[1] = solutions[0];

    return std::make_pair(solutions, 1U);
  }

  // no real roots
  return std::make_pair(solutions, 0U);
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

  if (alpha > (Scalar{1} / Scalar{2}) * pi &&
      alpha <= (Scalar{1} / Scalar{2}) * pi + tiny_epsilon) {
    alpha = (Scalar{1} / Scalar{2}) * pi;
  }
#endif

  overlap_assert(r > Scalar{0}, "invalid argument 'r' for regularized_wedge()");

  overlap_assert(d >= Scalar{0} && d <= r,
                 "invalid argument 'd' for regularized_wedge()");

  overlap_assert(alpha >= Scalar{0} && alpha <= (Scalar{1} / Scalar{2}) * pi,
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
    if (alpha > (Scalar{1} / Scalar{2}) * pi) {
      const auto h = r - z;

      return (pi / Scalar{3}) * h * h * (Scalar{3} * r - h) -
             regularized_wedge(r, d, pi - alpha);
    }

    return regularized_wedge(r, d, alpha);
  }

  const auto hemisphere_volume = ((Scalar{2} / Scalar{3}) * pi) * r * r * r;
  if (alpha > (Scalar{1} / Scalar{2}) * pi) {
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
  if (alpha < Scalar{0} && alpha > -detail::tiny_epsilon) {
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
  const auto arg0 = clamp(r * cos_alpha * factor, Scalar{-1}, Scalar{1},
                          detail::tiny_epsilon);

  const auto arg1 = clamp((z * cos_alpha * factor) / sin_alpha, Scalar{-1},
                          Scalar{1}, detail::tiny_epsilon);

  overlap_assert(Scalar{-1} <= arg0 && arg0 <= Scalar(1),
                 "invalid value for arg0 in regularized_wedge_area()");

  overlap_assert(Scalar{-1} <= arg1 && arg1 <= Scalar(1),
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
  if (dist < detail::tiny_epsilon) {
    // the wedge (almost) touches the center, the volume/area depends only on
    // the angle
    return spherical_wedge<Dim>(s, pi - detail::angle(p0.normal, p1.normal));
  }

  const auto s0 = d.dot(p0.normal);
  const auto s1 = d.dot(p1.normal);

  // detect degenerated general spherical wedge that can be treated as
  // a regularized spherical wedge
  if (std::abs(s0) < detail::tiny_epsilon ||
      std::abs(s1) < detail::tiny_epsilon) {
    const auto angle = pi - detail::angle(p0.normal, p1.normal);

    if (Dim == 2) {
      return detail::regularized_wedge_area(
          s.radius, std::abs(s0) > std::abs(s1) ? s0 : s1, angle);
    }

    return detail::regularized_wedge(s.radius, dist, angle,
                                     std::abs(s0) > std::abs(s1) ? s0 : s1);
  }

  auto d_unit = Vector{d * (Scalar{1} / dist)};
  if (dist < detail::large_epsilon) {
    d_unit = detail::gram_schmidt(p0.normal.cross(p1.normal).stableNormalized(),
                                  d_unit)[1];
  }

  overlap_assert(p0.normal.dot(p1.center - p0.center) <= Scalar{0},
                 "invalid plane in general_wedge()");

  overlap_assert(p1.normal.dot(p0.center - p1.center) <= Scalar{0},
                 "invalid plane in general_wedge()");

  // calculate the angles between the vector from the sphere center
  // to the intersection line and the normal vectors of the two planes
  auto alpha0 = detail::angle(p0.normal, d_unit);
  auto alpha1 = detail::angle(p1.normal, d_unit);

  const auto pi_half = (Scalar{1} / Scalar{2}) * pi;
  const auto dir0 = d_unit.dot((s.center + d) - p0.center);
  const auto dir1 = d_unit.dot((s.center + d) - p1.center);

  if (s0 >= Scalar{0} && s1 >= Scalar{0}) {
    alpha0 = pi_half - std::copysign(alpha0, dir0);
    alpha1 = pi_half - std::copysign(alpha1, dir1);

    if (Dim == 2) {
      return detail::regularized_wedge_area(s.radius, s0, alpha0) +
             detail::regularized_wedge_area(s.radius, s1, alpha1);
    }

    return detail::regularized_wedge(s.radius, dist, alpha0, s0) +
           detail::regularized_wedge(s.radius, dist, alpha1, s1);
  }

  if (s0 < Scalar{0} && s1 < Scalar{0}) {
    alpha0 = pi_half + std::copysign(Scalar{1}, dir0) * (alpha0 - pi);
    alpha1 = pi_half + std::copysign(Scalar{1}, dir1) * (alpha1 - pi);

    if (Dim == 2) {
      return s.surface_area() -
             (detail::regularized_wedge_area(s.radius, -s0, alpha0) +
              detail::regularized_wedge_area(s.radius, -s1, alpha1));
    }

    return s.volume - (detail::regularized_wedge(s.radius, dist, alpha0, -s0) +
                       detail::regularized_wedge(s.radius, dist, alpha1, -s1));
  }

  alpha0 = pi_half - std::copysign(Scalar{1}, dir0 * s0) *
                         (alpha0 - (s0 < Scalar{0} ? pi : Scalar{0}));

  alpha1 = pi_half - std::copysign(Scalar{1}, dir1 * s1) *
                         (alpha1 - (s1 < Scalar{0} ? pi : Scalar{0}));

  if constexpr (Dim == 2) {
    const auto area0 =
        detail::regularized_wedge_area(s.radius, std::abs(s0), alpha0);

    const auto area1 =
        detail::regularized_wedge_area(s.radius, std::abs(s1), alpha1);

    return std::max(area0, area1) - std::min(area0, area1);
  }

  const auto volume0 =
      detail::regularized_wedge(s.radius, dist, alpha0, std::abs(s0));

  const auto volume1 =
      detail::regularized_wedge(s.radius, dist, alpha1, std::abs(s1));

  return std::max(volume0, volume1) - std::min(volume0, volume1);
}

// Depending on the dimensionality, either the volume or external surface area
// of the general wedge is computed.
template<std::size_t Dim, typename Element>
auto general_wedge(const Sphere& sphere, const Element& element,
                   std::size_t edge,
                   const std::array<std::array<Vector, 2>,
                                    num_edges<Element>()>& intersections) {
  static_assert(Dim == 2 || Dim == 3, "invalid dimensionality, must be 2 or 3");

  const auto& f0 = element.faces[Element::edge_mapping[edge][1][0]];
  const auto& f1 = element.faces[Element::edge_mapping[edge][1][1]];

  const auto edge_midpoint =
      Vector{(Scalar{1} / Scalar{2}) *
             ((intersections[edge][0] +
               element.vertices[Element::edge_mapping[edge][0][0]]) +
              (intersections[edge][1] +
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

  for (std::size_t vertex_idx = 0U; vertex_idx < marked_vertices.size();
       ++vertex_idx) {
    if (!marked_vertices[vertex_idx]) {
      continue;
    }

    auto all_edges_marked = true;
    for (std::size_t local_edge_idx = 0U; local_edge_idx < 3U;
         ++local_edge_idx) {
      auto edge_idx = Element::vertex_mapping[vertex_idx][0][local_edge_idx];
      all_edges_marked &= static_cast<bool>(marked_edges[edge_idx]);
    }

    corrected_marked_vertices[vertex_idx] = all_edges_marked;
  }

  return corrected_marked_vertices;
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
auto overlap_volume(const Sphere& sOrig, const Element& elementOrig) -> Scalar {
  using namespace detail;

  static_assert(is_element_v<Element>, "invalid element type detected");

  if (!intersects_coarse(sOrig, elementOrig)) {
    return Scalar{0};
  }

  // Check for trivial case: element fully contained in sphere.
  if (contains(sOrig, elementOrig)) {
    return elementOrig.volume;
  }

  // Sanity check: All faces of the mesh element have to be planar.
  for (const auto& face : elementOrig.faces) {
    if (!face.is_planar()) {
      throw std::invalid_argument{"non-planer face detected in element"};
    }
  }

  // Use scaled and shifted versions of the sphere and the element.
  Transformation transformation{-sOrig.center, Scalar{1} / sOrig.radius};

  const auto s = Sphere{};
  auto element = Element{elementOrig};
  element.apply(transformation);

  // Constants: Number of vertices and faces.
  constexpr auto nrVertices = num_vertices<Element>();
  constexpr auto nrFaces = num_faces<Element>();

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<num_edges<Element>()> eMarked;
  std::bitset<nrFaces> fMarked;

  // Initial value: Volume of the full sphere.
  auto result = s.volume;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  std::array<std::array<Vector, 2>, num_edges<Element>()> eIntersections;

  // Process all edges of the element.
  for (std::size_t n = 0; n < num_edges<Element>(); ++n) {
    Vector start(element.vertices[Element::edge_mapping[n][0][0]]);
    Vector direction(element.vertices[Element::edge_mapping[n][0][1]] - start);

    auto solutions = line_sphere_intersection(start, direction, s);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!solutions.second ||
        (solutions.first[0] >= Scalar{1} - detail::medium_epsilon) ||
        solutions.first[1] <= detail::medium_epsilon ||
        (solutions.first[0] > Scalar{0} && solutions.first[1] < Scalar{1} &&
         (solutions.first[1] - solutions.first[0] < detail::large_epsilon))) {
      continue;
    }

    vMarked[Element::edge_mapping[n][0][0]] = solutions.first[0] < Scalar{0};
    vMarked[Element::edge_mapping[n][0][1]] = solutions.first[1] > Scalar{1};

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    eIntersections[n][0] = solutions.first[0] * direction;
    eIntersections[n][1] = (solutions.first[1] - Scalar{1}) * direction;
    eMarked[n] = true;

    // If the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well.
    fMarked[Element::edge_mapping[n][1][0]] = true;
    fMarked[Element::edge_mapping[n][1][1]] = true;
  }

  // Check whether the dependencies for a vertex intersection are fulfilled.
  vMarked = correct_marked_vertices<Element>(vMarked, eMarked);

  // Process all faces of the element, ignoring the edges as those where
  // already checked above.
  for (std::size_t n = 0; n < nrFaces; ++n) {
    if (intersects(s, element.faces[n])) {
      fMarked[n] = true;
    }
  }

  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(element, s.center)) {
    return sOrig.volume;
  }

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count()) {
    return Scalar{0};
  }

  // Iterate over all the marked faces and subtract the volume of the cap cut
  // off by the plane.
  for (std::size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) {
      continue;
    }

    const auto& f = element.faces[n];
    const auto dist = f.normal.dot(s.center - f.center);

    result -= s.cap_volume(s.radius + dist);
  }

  // Handle the edges and add back the volume subtracted twice above in the
  // processing of the faces.
  for (std::size_t n = 0; n < num_edges<Element>(); ++n) {
    if (!eMarked[n]) {
      continue;
    }

    result += general_wedge<3, Element>(s, element, n, eIntersections);
  }

  // Handle the vertices and subtract the volume added twice above in the
  // processing of the edges.
  for (std::size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) {
      continue;
    }

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    std::array<Vector, 3> intersectionPointsRelative;
    std::array<Vector, 3> intersectionPoints;
    for (std::size_t e = 0; e < 3; ++e) {
      const auto edgeIdx = Element::vertex_mapping[n][0][e];
      intersectionPointsRelative[e] =
          eIntersections[edgeIdx][Element::vertex_mapping[n][1][e]];

      intersectionPoints[e] =
          intersectionPointsRelative[e] + element.vertices[n];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.
    Triangle coneTria;
    coneTria.vertices = {
        {intersectionPoints[0], intersectionPoints[1], intersectionPoints[2]}};

    coneTria.center =
        (Scalar{1} / Scalar{3}) * std::accumulate(intersectionPoints.begin(),
                                                  intersectionPoints.end(),
                                                  Vector::Zero().eval());

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use double the normal precision to calculate this normal.
    coneTria.normal = detail::triangle_normal(intersectionPointsRelative[0],
                                              intersectionPointsRelative[1],
                                              intersectionPointsRelative[2]);

    // The area of this triangle is never needed, so it is set to an
    // invalid value.
    coneTria.area = std::numeric_limits<Scalar>::infinity();

    std::array<std::pair<std::size_t, Scalar>, 3> distances;
    for (std::size_t i = 0; i < 3; ++i) {
      distances[i] =
          std::make_pair(i, intersectionPointsRelative[i].squaredNorm());
    }

    std::sort(distances.begin(), distances.end(),
              [](const std::pair<std::size_t, Scalar>& a,
                 const std::pair<std::size_t, Scalar>& b) -> bool {
                return a.second < b.second;
              });

    if (distances[1].second < distances[2].second * detail::large_epsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      const auto correction = general_wedge<3, Element>(
          s, element, Element::vertex_mapping[n][0][distances[2].first],
          eIntersections);

      result -= correction;

      continue;
    }

    const auto tipTetVolume =
        (Scalar{1} / Scalar{6}) *
        std::abs(-intersectionPointsRelative[2].dot(
            (intersectionPointsRelative[0] - intersectionPointsRelative[2])
                .cross(intersectionPointsRelative[1] -
                       intersectionPointsRelative[2])));

    // Make sure the normal points in the right direction i.e. away from
    // the center of the element.
    if (coneTria.normal.dot(element.center - coneTria.center) > Scalar{0}) {
      coneTria.normal = -coneTria.normal;
    }

    const auto dist = coneTria.normal.dot(s.center - coneTria.center);
    const auto cap_volume = s.cap_volume(s.radius + dist);

    // The cap volume is tiny, so the corrections will be even smaller.
    // There is no way to actually calculate them with reasonable
    // precision, so just the volume of the tetrahedron at the tip is
    // used.
    if (cap_volume < detail::tiny_epsilon) {
      result -= tipTetVolume;
      continue;
    }

    // Calculate the volume of the three spherical segments between
    // the faces joining at the vertex and the plane through the
    // intersection points.
    const auto plane = Plane{coneTria.center, coneTria.normal};
    auto segmentVolume = Scalar{0};

    for (std::size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[Element::vertex_mapping[n][2][e]];
      const auto e0 = Element::face_mapping[e][0];
      const auto e1 = Element::face_mapping[e][1];
      const Vector center = (Scalar{1} / Scalar{2}) *
                            (intersectionPoints[e0] + intersectionPoints[e1]);

      const auto wedgeVolume = general_wedge<3>(
          s, plane, Plane(f.center, -f.normal), center - s.center);

      segmentVolume += wedgeVolume;
    }

    // Calculate the volume of the cone and clamp it to zero.
    const auto coneVolume =
        std::max(tipTetVolume + cap_volume - segmentVolume, Scalar{0});

    // Sanity check: detect negative cone volume.
    overlap_assert(coneVolume > -std::sqrt(detail::tiny_epsilon),
                   "negative cone volume in overlap_volume()");

    result -= coneVolume;

    // Sanity check: detect negative intermediate result.
    overlap_assert(result > -std::sqrt(detail::tiny_epsilon),
                   "negative intermediate result in overlap_volume()");
  }

  // In case of different sized objects the error can become quite large,
  // so a relative limit is used.
  const auto maxOverlap = std::min(s.volume, element.volume);
  const auto limit =
      std::sqrt(std::numeric_limits<Scalar>::epsilon()) * maxOverlap;

  // Clamp tiny negative volumes to zero.
  if (result < Scalar{0} && result > -limit) {
    return Scalar{0};
  }

  // Clamp results slightly too large.
  if (result > maxOverlap && result - maxOverlap < limit) {
    return std::min(sOrig.volume, elementOrig.volume);
  }

  // Perform a sanity check on the final result (debug version only).
  overlap_assert(result >= Scalar{0} && result <= maxOverlap,
                 "negative volume detected in overlap_volume()");

  // Scale the overlap volume back for the original objects.
  result = (result / s.volume) * sOrig.volume;

  return result;
}

template<typename Iterator>
auto overlap_volume(const Sphere& s, Iterator eBegin, Iterator eEnd) -> Scalar {
  auto sum = Scalar{0};
  for (auto it = eBegin; it != eEnd; ++it) {
    sum += overlap_volume(s, *it);
  }

  return sum;
}

// Calculate the surface area of the sphere and the element that are contained
// within the common or intersecting part of the geometries, respectively.
// The returned array of size (N + 2), with N being the number of vertices,
// holds (in this order):
//   - surface area of the region of the sphere intersecting the element
//   - for each face of the element: area contained within the sphere
//   - total surface area of the element intersecting the sphere
template<typename Element,
         std::size_t NrFaces = detail::num_faces<Element>() + 2>
auto overlap_area(const Sphere& sOrig, const Element& elementOrig)
    -> std::array<Scalar, NrFaces> {
  using namespace detail;

  static_assert(is_element_v<Element>, "invalid element type detected");
  static_assert(NrFaces == num_faces<Element>() + 2,
                "invalid number of faces for the element provided");

  constexpr auto nrVertices = num_vertices<Element>();
  constexpr auto nrFaces = num_faces<Element>();

  // Initial value: Zero overlap.
  std::array<Scalar, nrFaces + 2> result{};

  if (!intersects_coarse(sOrig, elementOrig)) {
    return result;
  }

  // Check for trivial case: element fully contained in sphere resulting in
  // a full coverage of all faces.
  if (contains(sOrig, elementOrig)) {
    for (std::size_t n = 0; n < nrFaces; ++n) {
      result[n + 1] = elementOrig.faces[n].area;
      result[nrFaces + 1] += elementOrig.faces[n].area;
    }

    return result;
  }

  // Sanity check: All faces of the mesh element have to be planar.
  for (const auto& face : elementOrig.faces) {
    if (!face.is_planar()) {
      throw std::runtime_error{"non-planer face detected in element"};
    }
  }

  // Use scaled and shifted versions of the sphere and the element.
  Transformation transformation(-sOrig.center, Scalar(1) / sOrig.radius);

  Sphere s(Vector::Zero(), Scalar(1));

  Element element(elementOrig);
  element.apply(transformation);

  // Sets of overlapping primitives.
  std::bitset<nrVertices> vMarked;
  std::bitset<num_edges<Element>()> eMarked;
  std::bitset<nrFaces> fMarked;

  // The intersection points between the single edges and the sphere, this
  // is needed later on.
  std::array<std::array<Vector, 2>, num_edges<Element>()> eIntersections;

  // Cache the squared radius of the disk formed by the intersection between
  // the planes defined by each face and the sphere.
  std::array<Scalar, nrFaces> intersectionRadiusSq;

  // Process all edges of the element.
  for (std::size_t n = 0; n < num_edges<Element>(); ++n) {
    Vector start(element.vertices[Element::edge_mapping[n][0][0]]);
    Vector direction(element.vertices[Element::edge_mapping[n][0][1]] - start);

    auto solutions = line_sphere_intersection(start, direction, s);

    // No intersection between the edge and the sphere, where intersection
    // points close to the surface of the sphere are ignored.
    // Or:
    // The sphere cuts the edge twice, no vertex is inside of the
    // sphere, but the case of the edge only touching the sphere has to
    // be avoided.
    if (!solutions.second ||
        solutions.first[0] >= Scalar(1) - detail::medium_epsilon ||
        solutions.first[1] <= detail::medium_epsilon ||
        (solutions.first[0] > Scalar(0) && solutions.first[1] < Scalar(1) &&
         solutions.first[1] - solutions.first[0] < detail::large_epsilon)) {
      continue;
    }

    vMarked[Element::edge_mapping[n][0][0]] = solutions.first[0] < Scalar(0);

    vMarked[Element::edge_mapping[n][0][1]] = solutions.first[1] > Scalar(1);

    // Store the two intersection points of the edge with the sphere for
    // later usage.
    eIntersections[n][0] = solutions.first[0] * direction;
    eIntersections[n][1] = (solutions.first[1] - Scalar(1)) * direction;

    eMarked[n] = true;

    // If the edge is marked as having an overlap, the two faces forming it
    // have to be marked as well.
    fMarked[Element::edge_mapping[n][1][0]] = true;
    fMarked[Element::edge_mapping[n][1][1]] = true;
  }

  // Check whether the dependencies for a vertex intersection are fulfilled.
  for (std::size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) {
      continue;
    }

    bool edgesValid = true;
    for (std::size_t eN = 0; eN < 3; ++eN) {
      const auto edgeId = Element::vertex_mapping[n][0][eN];
      edgesValid &= eMarked[edgeId];
    }

    // If not all three edges intersecting at this vertex where marked, the
    // sphere is only touching.
    if (!edgesValid) {
      vMarked[n] = false;
    }
  }

  // Process all faces of the element, ignoring the edges as those where
  // already checked above.
  for (std::size_t n = 0; n < nrFaces; ++n) {
    if (intersects(s, element.faces[n])) {
      fMarked[n] = true;
    }
  }
  // Trivial case: The center of the sphere overlaps the element, but the
  // sphere does not intersect any of the faces of the element, meaning the
  // sphere is completely contained within the element.
  if (!fMarked.count() && contains(element, s.center)) {
    result[0] = sOrig.surface_area();

    return result;
  }

  // Spurious intersection: The initial intersection test was positive, but
  // the detailed checks revealed no overlap.
  if (!vMarked.count() && !eMarked.count() && !fMarked.count()) return result;

  // Initial value for the surface of the sphere: Surface area of the full
  // sphere.
  result[0] = s.surface_area();

  // Iterate over all the marked faces and calculate the area of the disk
  // defined by the plane as well as the cap surfaces.
  for (std::size_t n = 0; n < nrFaces; ++n) {
    if (!fMarked[n]) {
      continue;
    }
    const auto& f = element.faces[n];
    Scalar dist = f.normal.dot(s.center - f.center);
    result[0] -= s.cap_surface_area(s.radius + dist);
    result[n + 1] = s.disk_area(s.radius + dist);
  }

  // Handle the edges and subtract the area of the respective disk cut off by
  // the edge and add back the surface area of the spherical wedge defined
  // by the edge.
  for (std::size_t n = 0; n < num_edges<Element>(); ++n) {
    if (!eMarked[n]) {
      continue;
    }

    result[0] += general_wedge<2, Element>(s, element, n, eIntersections);

    // The intersection points are relative to the vertices forming the
    // edge.
    const Vector chord = ((element.vertices[Element::edge_mapping[n][0][0]] +
                           eIntersections[n][0]) -
                          (element.vertices[Element::edge_mapping[n][0][1]] +
                           eIntersections[n][1]));

    const Scalar chordLength = chord.stableNorm();

    // Each edge belongs to two faces, indexed via
    // Element::edge_mapping[n][1][{0,1}].
    for (std::size_t e = 0; e < 2; ++e) {
      const auto faceIdx = Element::edge_mapping[n][1][e];
      const auto& f = element.faces[faceIdx];

      // Height of the spherical cap cut off by the plane containing the
      // face.
      const Scalar dist = f.normal.dot(s.center - f.center) + s.radius;
      intersectionRadiusSq[faceIdx] = dist * (Scalar(2) * s.radius - dist);

      // Calculate the height of the triangular segment in the plane of
      // the base.
      const Scalar factor = std::sqrt(
          std::max(Scalar(0), intersectionRadiusSq[faceIdx] -
                                  Scalar(0.25) * chordLength * chordLength));

      const Scalar theta =
          Scalar(2) * std::atan2(chordLength, Scalar(2) * factor);

      Scalar area = Scalar(0.5) * intersectionRadiusSq[faceIdx] *
                    (theta - std::sin(theta));

      // FIXME: Might not be necessary to use the center of the chord.
      const Vector chordCenter =
          Scalar(0.5) * ((element.vertices[Element::edge_mapping[n][0][0]] +
                          eIntersections[n][0]) +
                         (element.vertices[Element::edge_mapping[n][0][1]] +
                          eIntersections[n][1]));

      const Vector proj(s.center -
                        f.normal.dot(s.center - f.center) * f.normal);

      // If the projected sphere center and the face center fall on
      // opposite sides of the edge, the area has to be inverted.
      if (chord.cross(proj - chordCenter)
              .dot(chord.cross(f.center - chordCenter)) < Scalar(0)) {
        area = intersectionRadiusSq[faceIdx] * pi - area;
      }

      result[faceIdx + 1] -= area;
    }
  }

  // Handle the vertices and add the area subtracted twice above in the
  // processing of the edges.

  // First, handle the spherical surface area of the intersection.
  // This is to a large part code duplicated from the volume calculation.
  // TODO: Unify the area and volume calculation to remove duplicate code.
  for (std::size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) {
      continue;
    }

    // Collect the points where the three edges intersecting at this
    // vertex intersect the sphere.
    // Both the relative and the absolute positions are required.
    std::array<Vector, 3> intersectionPointsRelative;
    std::array<Vector, 3> intersectionPoints;
    for (std::size_t e = 0; e < 3; ++e) {
      auto edgeIdx = Element::vertex_mapping[n][0][e];
      intersectionPointsRelative[e] =
          eIntersections[edgeIdx][Element::vertex_mapping[n][1][e]];

      intersectionPoints[e] =
          intersectionPointsRelative[e] + element.vertices[n];
    }

    // This triangle is constructed by hand to have more freedom of how
    // the normal vector is calculated.
    Triangle coneTria;
    coneTria.vertices = {
        {intersectionPoints[0], intersectionPoints[1], intersectionPoints[2]}};

    coneTria.center =
        (Scalar{1} / Scalar{3}) * std::accumulate(intersectionPoints.begin(),
                                                  intersectionPoints.end(),
                                                  Vector::Zero().eval());

    // Calculate the normal of the triangle defined by the intersection
    // points in relative coordinates to improve accuracy.
    // Also use double the normal precision to calculate this normal.
    coneTria.normal = detail::triangle_normal(intersectionPointsRelative[0],
                                              intersectionPointsRelative[1],
                                              intersectionPointsRelative[2]);

    // The area of this triangle is never needed, so it is set to an
    // invalid value.
    coneTria.area = std::numeric_limits<Scalar>::infinity();

    std::array<std::pair<std::size_t, Scalar>, 3> distances;
    for (std::size_t i = 0; i < 3; ++i) {
      distances[i] =
          std::make_pair(i, intersectionPointsRelative[i].squaredNorm());
    }

    std::sort(distances.begin(), distances.end(),
              [](const std::pair<std::size_t, Scalar>& a,
                 const std::pair<std::size_t, Scalar>& b) -> bool {
                return a.second < b.second;
              });

    if (distances[1].second < distances[2].second * detail::large_epsilon) {
      // Use the general spherical wedge defined by the edge with the
      // non-degenerated intersection point and the normals of the
      // two faces forming it.
      Scalar correction = general_wedge<2, Element>(
          s, element, Element::vertex_mapping[n][0][distances[2].first],
          eIntersections);

      result[0] -= correction;

      continue;
    }

    // Make sure the normal points in the right direction, i.e., away from
    // the center of the element.
    if (coneTria.normal.dot(element.center - coneTria.center) > Scalar{0}) {
      coneTria.normal = -coneTria.normal;
    }

    Plane plane(coneTria.center, coneTria.normal);

    Scalar dist = coneTria.normal.dot(s.center - coneTria.center);
    Scalar capSurface = s.cap_surface_area(s.radius + dist);

    // If cap surface area is small, the corrections will be even smaller.
    // There is no way to actually calculate them with reasonable
    // precision, so they are just ignored.
    if (capSurface < detail::large_epsilon) {
      continue;
    }

    // Calculate the surface area of the three spherical segments between
    // the faces joining at the vertex and the plane through the
    // intersection points.
    Scalar segmentSurface = 0;
    for (std::size_t e = 0; e < 3; ++e) {
      const auto& f = element.faces[Element::vertex_mapping[n][2][e]];
      uint32_t e0 = Element::face_mapping[e][0];
      uint32_t e1 = Element::face_mapping[e][1];

      Vector center(Scalar{0.5} *
                    (intersectionPoints[e0] + intersectionPoints[e1]));

      segmentSurface += general_wedge<2>(s, plane, Plane(f.center, -f.normal),
                                         center - s.center);
    }

    // Calculate the surface area of the cone and clamp it to zero.
    Scalar coneSurface = std::max(capSurface - segmentSurface, Scalar{0});

    result[0] -= coneSurface;

    // Sanity checks: detect negative/excessively large intermediate
    // result.
    overlap_assert(result[0] > -std::sqrt(detail::tiny_epsilon),
                   "negative area as intermediate result in overlap_area()");

    overlap_assert(result[0] < s.surface_area() + detail::tiny_epsilon,
                   "invalid intermediate result in overlap_area()");
  }

  // Second, correct the intersection area of the facets.
  for (std::size_t n = 0; n < nrVertices; ++n) {
    if (!vMarked[n]) {
      continue;
    }

    // Iterate over all the faces joining at this vertex.
    for (std::size_t f = 0; f < 3; ++f) {
      // Determine the two edges of this face intersecting at the
      // vertex.
      uint32_t e0 = Element::face_mapping[f][0];
      uint32_t e1 = Element::face_mapping[f][1];
      std::array<uint32_t, 2> edgeIndices = {
          {Element::vertex_mapping[n][0][e0],
           Element::vertex_mapping[n][0][e1]}};

      // Extract the (relative) intersection points of these edges with
      // the sphere furthest from the vertex.
      std::array<Vector, 2> intersectionPoints = {
          {eIntersections[edgeIndices[0]][Element::vertex_mapping[n][1][e0]],

           eIntersections[edgeIndices[1]][Element::vertex_mapping[n][1][e1]]}};

      // Together with the vertex, this determines the triangle
      // representing one part of the correction.
      const Scalar triaArea =
          Scalar{0.5} *
          (intersectionPoints[0].cross(intersectionPoints[1])).stableNorm();

      // The second component is the segment defined by the face and the
      // intersection points.
      const Scalar chordLength =
          (intersectionPoints[0] - intersectionPoints[1]).stableNorm();

      const auto faceIdx = Element::vertex_mapping[n][2][f];

      // TODO: Cache theta for each edge.
      const Scalar theta =
          Scalar{2} *
          std::atan2(chordLength,
                     Scalar{2} * std::sqrt(std::max(
                                     Scalar{0}, intersectionRadiusSq[faceIdx] -
                                                    Scalar{0.25} * chordLength *
                                                        chordLength)));

      Scalar segmentArea = Scalar{0.5} * intersectionRadiusSq[faceIdx] *
                           (theta - std::sin(theta));

      // Determine if the (projected) center of the sphere lies within
      // the triangle or not. If not, the segment area has to be
      // corrected.
      const Vector d(Scalar{0.5} *
                     (intersectionPoints[0] + intersectionPoints[1]));

      const auto& face = element.faces[faceIdx];
      const Vector proj(s.center -
                        face.normal.dot(s.center - face.center) * face.normal);

      if (d.dot((proj - element.vertices[n]) - d) > Scalar{0}) {
        segmentArea = intersectionRadiusSq[faceIdx] * pi - segmentArea;
      }

      result[faceIdx + 1] += triaArea + segmentArea;

      // Sanity checks: detect excessively large intermediate result.
      overlap_assert(result[faceIdx + 1] < element.faces[faceIdx].area +
                                               std::sqrt(detail::large_epsilon),
                     "invalid intermediate result in overlap_area()");
    }
  }

  // Scale the surface areas back for the original objects and clamp
  // values within reasonable limits.
  const Scalar scaling = sOrig.radius / s.radius;
  const Scalar sLimit(std::sqrt(std::numeric_limits<Scalar>::epsilon()) *
                      s.surface_area());

  // As the precision of the area calculation deteriorates quickly with a
  // increasing size ratio between the element and the sphere, the precision
  // limit applied to the sphere is used as the lower limit for the facets.
  const Scalar fLimit(
      std::max(sLimit, std::sqrt(std::numeric_limits<Scalar>::epsilon()) *
                           element.surface_area()));

  // Sanity checks: detect negative/excessively large results for the
  // surface area of the facets.
  for (std::size_t n = 0; n < nrFaces; ++n) {
    overlap_assert(result[n + 1] > -fLimit,
                   "negative overlap area for face in overlap_area()");

    overlap_assert(result[n + 1] <= element.faces[n].area + fLimit,
                   "invalid overlap area for face in overlap_area()");
  }

  // Surface of the sphere.
  result[0] = detail::clamp(result[0], Scalar{0}, s.surface_area(), sLimit);
  result[0] *= (scaling * scaling);

  // Surface of the mesh element.
  for (std::size_t f = 0; f < nrFaces; ++f) {
    auto& value = result[f + 1];
    value = detail::clamp(value, Scalar{0}, element.faces[f].area, fLimit);

    value = value * (scaling * scaling);
  }

  result.back() =
      std::accumulate(result.begin() + 1, result.end() - 1, Scalar{0});

  // Perform some more sanity checks on the final result (debug version
  // only).
  overlap_assert(Scalar(0) <= result[0] && result[0] <= sOrig.surface_area(),
                 "invalid overlap area for sphere surface in overlap_area()");

  overlap_assert(
      Scalar(0) <= result.back() && result.back() <= elementOrig.surface_area(),
      "invalid total overlap area for faces in overlap_area()");

  return result;
}

}  // namespace overlap

#endif  // OVERLAP_HPP

#pragma once

namespace steps::util {

namespace detail {
/// type traits to figure out if the container given to \a strong_random_access
/// should be copied (for gls::span) or treated as reference
/// (for std::vector, Omega_h::array, ...)
template <typename Container>
struct strong_ra_traits {
    using member_type = Container;
    static constexpr bool container_owned = true;
};

template <typename Container>
struct strong_ra_traits<const Container&> {
    static constexpr bool specialization = 0;
    using member_type = const Container&;
    static constexpr bool container_owned = false;
};

template <typename Container>
struct strong_ra_traits<Container&> {
    using member_type = Container&;
    static constexpr bool container_owned = false;
};

template <typename Container>
struct strong_ra_traits<Container&&> {
    using member_type = Container;
    static constexpr bool container_owned = true;
};

template <typename T>
struct strong_ra_traits<gsl::span<T>&> {
    using member_type = gsl::span<T>;
    static constexpr bool container_owned = true;
};

#ifdef STEPS_USE_DIST_MESH

template <typename T>
struct strong_ra_traits<Omega_h::Read<T>&> {
    using member_type = const Omega_h::Read<T>;
    static constexpr bool container_owned = true;
};

#endif  // STEPS_USE_DIST_MESH

}  // namespace detail

/**
 * Wrapper class for random-access containers like std::vector, Omega_h arrays, or gsl::span
 * where indices to access data are strong identifier, not integral values.
 */
template <typename StrongId, typename Container>
class strong_random_access {
  private:
    using decay_container_type = std::decay_t<Container>;

  public:
    using container_type = typename detail::strong_ra_traits<Container>::member_type;
    using value_type = typename decay_container_type::value_type;
    using reference = value_type&;
    using const_reference = const value_type&;

    static constexpr bool container_owned = detail::strong_ra_traits<Container>::container_owned;

  private:
    template <typename T>
    using enable_if_non_const =
        typename std::enable_if<!std::is_const_v<std::remove_reference_t<T>>, bool>::type;

  public:
    strong_random_access(strong_random_access<StrongId, Container>& other)
        : strong_random_access(
              const_cast<const strong_random_access<StrongId, Container>&>(other)) {}

    strong_random_access(const strong_random_access<StrongId, Container>& other)
        : container_(other.container_) {}

    strong_random_access(strong_random_access<StrongId, Container>&& other)
        : container_(std::move(other.container_)) {}

    strong_random_access<StrongId, Container>& operator=(
        const strong_random_access<StrongId, Container>& other) {
        container_ = other.container_;
        return *this;
    }

    template <typename... Args>
    constexpr strong_random_access(Args&&... args)
        : container_(std::forward<Args>(args)...) {}

    auto range() const noexcept {
        return StrongId::range(size());
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    reference at(StrongId pos) {
        return container_.at(pos.get());
    }

    const_reference at(StrongId pos) const {
        return container_.at(pos.get());
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    reference operator[](StrongId pos) {
        return container_[pos.get()];
    }

    const_reference operator[](StrongId pos) const {
        return container_[pos.get()];
    }

    bool operator==(const strong_random_access<StrongId, Container>& other) const noexcept {
        return container() == other.container();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    reference front() {
        return container_.front();
    }

    const_reference front() const {
        return container_.front();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    reference back() {
        return container_.back();
    }

    const_reference back() const {
        return container_.back();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    auto begin() noexcept {
        return container_.begin();
    }

    auto begin() const noexcept {
        return container_.begin();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    auto end() noexcept {
        return container_.end();
    }

    auto end() const noexcept {
        return container_.end();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    value_type* data() noexcept {
        return container_.data();
    }

    const value_type* data() const noexcept {
        return container_.data();
    }

    typename StrongId::value_type size() const noexcept {
        return container_.size();
    }

    template <typename Member = container_type, enable_if_non_const<Member> = true>
    container_type& container() noexcept {
        return container_;
    }

    const container_type& container() const noexcept {
        return container_;
    }

  private:
    container_type container_;
};

/// utility functions to ease instantiation of a \a strong_random_access,
///
///   std::vector<int> ra(3);
///   auto sra = make_strong_random_accessor<solver::local_specie_id>(ra);
///
template <typename StrongId, typename Container>
auto make_strong_random_accessor(Container& container) {
    return strong_random_access<StrongId, Container&>(container);
}

/// utility function to ease instantiation of a \a strong_random_access,
/// where only the StrongId type has to be specified i.e. This overloading
/// transfer the ownership of the container to the \a strong_random_access instance
///
///   std::vector<int> ra(3);
///   auto sra = make_strong_random_accessor<solver::local_specie_id>(std::move(ra));
///
template <typename StrongId, typename Container>
auto make_strong_random_accessor(Container&& container) {
    return strong_random_access<StrongId, Container&&>(std::forward<Container>(container));
}

/// utility function to ease instantiation of a \a strong_random_access,
/// With this overloading, the \a strong_random_access instance has ownership of the
/// random-access container. The arguments specified in this function are forwarded
/// to the container constructor.
///
///   auto sra = make_strong_random_accessor<solver::local_specie_id, std::vector<int>>(3, 5);
///
template <typename StrongId, typename Container, typename... Args>
auto make_strong_random_accessor(Args&&... args) {
    return strong_random_access<StrongId, Container&&>(std::forward<Args>(args)...);
}

template <typename StrongId, typename... Args>
using strongid_vector = strong_random_access<StrongId, std::vector<Args...>>;

}  // namespace steps::util

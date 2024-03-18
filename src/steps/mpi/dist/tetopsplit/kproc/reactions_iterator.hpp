/**
 * \file reactions_iterator.hpp
 * Provide iterator constructs for \a ReactionsType classes
 */

#pragma once

#include <iterator>

namespace steps::dist::kproc {

template <typename ReactionClass>
struct reactions_iterator;

template <typename ReactionClass>
struct reaction_partial {
    inline reaction_partial(ReactionClass& reactions, size_t pos) noexcept
        : reactions_(reactions)
        , pos_(pos) {}

    mesh::tetrahedron_id_t getOwnerPoint() const noexcept {
        return reactions_.getOwnerPoint(pos_);
    }

    osh::Real computeRate(const MolState& mol_state) const noexcept {
        return reactions_.computeRate(mol_state, pos_);
    }

    size_t getIndex() const noexcept {
        return pos_;
    }

    const std::vector<MolStateElementID>& getMolStateElementsUpdates() const noexcept {
        return reactions_.getMolStateElementsUpdates(pos_);
    }

    const std::vector<MolStateElementID>& getPropensityDependency() const noexcept {
        return reactions_.getPropensityDependency(pos_);
    }

  private:
    ReactionClass& reactions_;
    size_t pos_;
    friend struct reactions_iterator<ReactionClass>;
};

/**
 * Iterator class for \a ReactionsType classes
 */
template <typename ReactionClass>
struct reactions_iterator {
    inline reactions_iterator(ReactionClass& reactions, size_t pos = 0) noexcept;

    inline reactions_iterator& operator++() noexcept;
    inline const reactions_iterator operator++(int) noexcept;

    inline const reaction_partial<ReactionClass>& operator*() const noexcept;

    inline bool operator==(const reactions_iterator& other) const noexcept;
    inline bool operator!=(const reactions_iterator& other) const noexcept;

  private:
    reaction_partial<ReactionClass> partial_;
};

// class reactions_iterator
//------------------------------------------------------------------

template <class ReactionClass>
inline reactions_iterator<ReactionClass>::reactions_iterator(ReactionClass& reactions,
                                                             size_t pos) noexcept
    : partial_(reactions, pos) {}

//-------------------------------------------------------------

template <class ReactionClass>
inline reactions_iterator<ReactionClass>& reactions_iterator<ReactionClass>::operator++() noexcept {
    ++partial_.pos_;
    return *this;
}

//------------------------------------------------------------------

template <class ReactionClass>
inline const reactions_iterator<ReactionClass> reactions_iterator<ReactionClass>::operator++(
    int) noexcept {
    reactions_iterator copy(*this);
    ++partial_.pos_;
    return copy;
}

//------------------------------------------------------------------

template <class ReactionClass>
inline const reaction_partial<ReactionClass>& reactions_iterator<ReactionClass>::operator*()
    const noexcept {
    return partial_;
}

//------------------------------------------------------------------

template <class ReactionClass>
inline bool reactions_iterator<ReactionClass>::operator==(
    const reactions_iterator& other) const noexcept {
    return &partial_.reactions_ == &other.partial_.reactions_ &&
           partial_.pos_ == other.partial_.pos_;
}

//------------------------------------------------------------------

template <class ReactionClass>
inline bool reactions_iterator<ReactionClass>::operator!=(
    const reactions_iterator& other) const noexcept {
    return partial_.pos_ != other.partial_.pos_ ||
           &partial_.reactions_ != &other.partial_.reactions_;
}

template <typename ReactionClass>
inline size_t operator-(const reactions_iterator<ReactionClass>& lhs,
                        const reactions_iterator<ReactionClass>& rhs) {
    return (*lhs).getIndex() - (*rhs).getIndex();
}

}  // namespace steps::dist::kproc
namespace std {
template <typename ReactionClass>
struct iterator_traits<steps::dist::kproc::reactions_iterator<ReactionClass>> {
    using difference_type = std::ptrdiff_t;
    using value_type = typename steps::dist::kproc::reaction_partial<ReactionClass>;
    using pointer = typename steps::dist::kproc::reaction_partial<ReactionClass>*;
    using reference = typename steps::dist::kproc::reaction_partial<ReactionClass>&;
    using iterator_category = random_access_iterator_tag;
};
}  // namespace std

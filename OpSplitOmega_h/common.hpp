#pragma once

#include <initializer_list>
#include <memory>

#include <Omega_h_defines.hpp>
#include <petscsys.h>

#include "opsplit/fwd.hpp"
#include "opsplit/vocabulary.hpp"

#if PETSC_VERSION_LT(3, 10, 0)
#define PETSC_DM_SET_SECTION DMSetDefaultSection
#else
#define PETSC_DM_SET_SECTION DMSetSection
#endif

namespace zee {

static const PetscScalar AVOGADRO = 6.02214076e23;
static const PetscInt UNKNOWN_IDX = -1;

// Forward declarations
enum class SSAMethod { SSA, RSSA };

template <typename RNG>
class SimulationInput;

template <osh::Int Dim, SSAMethod SSA, typename RNG>
class SimulationData;

template <osh::Int Dim, SSAMethod SSA, typename RNG>
using SimulationDataPtr = std::unique_ptr<SimulationData<Dim, SSA, RNG>>;

template <osh::Int Dim>
class OmegaHMesh;

template <osh::Int Dim, typename RNG>
class DiffusionOperator;

template <typename RNG>
class SSAOperator;

template <typename RNG>
class RSSAOperator;

template <osh::Int Dim>
class PoolsIncrements;

template <osh::Int Dim>
class DiffusionDiscretizedRates;

class MolState;
template <typename RNG>
class LeavingMolecules;
class DiffusionVariables;

template <osh::LO Dim, typename RNG>
class Diffusions;

namespace details {

// sink to consume expanded arguments
struct sink {
    template <typename T>
    sink(T const& /*v*/) {}
};

template <int... Values>
struct integer_sequence {};

template <std::size_t Size, int... Accu>
struct ones_traits {
    using type = typename ones_traits<Size - 1, 1, Accu...>::type;
};

template <int... Accu>
struct ones_traits<0, Accu...> {
    using type = integer_sequence<Accu...>;
};

template <typename T, int... Ones>
std::initializer_list<T> initializer_list_impl(T value, const integer_sequence<Ones...>&) {
    static auto array = {(sink{Ones}, value)...};
    return array;
}

}  // namespace details

template <std::size_t Size, typename T>
std::initializer_list<T> initializer_list(T value) {
    return details::initializer_list_impl<T>(value, typename details::ones_traits<Size>::type());
}

}  // namespace zee

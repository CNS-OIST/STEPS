#pragma once

#include "kproc/kproc_id.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

template <typename RNG, typename NumMolecules> class SimulationInput;

enum class NextEventSearchMethod { Direct, GibsonBruck };

enum class SSAMethod { SSA, RSSA };

template <SSAMethod SSA, typename RNG, typename NumMolecules,
          NextEventSearchMethod SearchMethod>
class SimulationData;

template <SSAMethod SSA, typename RNG, typename NumMolecules,
          NextEventSearchMethod SearchMethod>
using SimulationDataPtr =
    std::unique_ptr<SimulationData<SSA, RNG, NumMolecules, SearchMethod>>;

template <typename NumMolecules> class MolState;

template <typename NumMolecules> struct propensity_function_traits {
  using value =
      std::function<osh::Real(kproc::KProcID, const MolState<NumMolecules> &)>;
};
} // namespace dist
} // namespace steps

namespace std {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
static inline std::string to_string(steps::dist::SSAMethod method) {
  switch (method) {
  case steps::dist::SSAMethod::SSA:
    return "DIRECT";
  case steps::dist::SSAMethod::RSSA:
    return "RSSA";
  }
}
#pragma GCC diagnostic pop
} // namespace std

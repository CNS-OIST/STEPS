#pragma once

#include <functional>

#include "kproc/kproc_id.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

class SimulationInput;

enum class NextEventSearchMethod { Direct = 0, GibsonBruck, RLeaping };

enum class SSAMethod { SSA = 0, RSSA, RLeaping };

enum class DiffusionMethod { ConstantDiffDt = 0, TauLeapingDiffDt };

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
class SimulationData;

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
using SimulationDataPtr = std::unique_ptr<SimulationData<SSA, SearchMethod, DiffMethod>>;

class MolState;

struct propensity_function_traits {
    using value = std::function<osh::Real(kproc::KProcID, const MolState&)>;
};

}  // namespace steps::dist

namespace std {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
static inline std::string to_string(steps::dist::SSAMethod method) {
    switch (method) {
    case steps::dist::SSAMethod::SSA:
        return "DIRECT";
    case steps::dist::SSAMethod::RLeaping:
        return "RLEAPING";
    case steps::dist::SSAMethod::RSSA:
        return "RSSA";
    }
}
#pragma GCC diagnostic pop
}  // namespace std

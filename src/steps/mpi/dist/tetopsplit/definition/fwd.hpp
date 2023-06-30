#pragma once

#include <functional>
#include <optional>
#include <string>

#include <Omega_h_defines.hpp>

#include "util/vocabulary.hpp"

namespace steps::dist {

namespace osh = Omega_h;

// Forward declaration
class Compdef;
class Diffdef;
template <typename P>
class SReacdefBase;
struct SReacInfo {
    osh::Real kCst;
};
using SReacdef = SReacdefBase<SReacInfo>;
using vdep_propensity_fun_t = std::function<osh::Real(osh::Real)>;
struct VDepInfo {
    vdep_propensity_fun_t kCstFun;
};

/// Parameters of the particular GHK reaction. Note that, conversely from
/// parallel steps, there is no concept of "global temperature" in distributed
/// steps.
using VDepSReacdef = SReacdefBase<VDepInfo>;
struct GHKInfo {
    /// GHK current identifier
    model::ghk_current_id curr_id;
    /// if true, the surface reaction involves an ion transfer from the inner to outer compartment,
    /// and conversely otherwise.
    bool in2out;
    /// permeability per ion channel
    osh::Real permeability;
    /// ion valence
    osh::I64 valence;
    /// optional locked-in inner compartment concentration of the ion
    std::optional<osh::Real> inner_conc;
    /// optional locked-in outer compartment concentration of the ion
    std::optional<osh::Real> outer_conc;
};
using GHKSReacdef = SReacdefBase<GHKInfo>;
class Reacdef;
class Patchdef;

class Simulation;
class Statedef;

}  // namespace steps::dist

#pragma once

#include <functional>
#include <optional>
#include <string>

#include <Omega_h_defines.hpp>
#include <type_traits>

#include "model/complexsreac.hpp"
#include "model/fwd.hpp"
#include "model/ghkcurr.hpp"
#include "model/sreac.hpp"
#include "model/vdepsreac.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

namespace osh = Omega_h;

// Forward declaration
class Compdef;
class Diffdef;
class Complexdef;
template <typename P>
class ReacdefBase;
using Reacdef = ReacdefBase<steps::model::Reac>;
class ComplexReacdef;

template <typename Entity>
struct EntityMolecules;
struct ComplexFilterDescr;
template <typename Entity>
class ComplexFilter;
struct FilterHash;
template <typename Entity>
class ComplexState;

template <typename P>
class SReacdefBase;
template <typename T>
class ModelSReacdef;
template <typename T>
class ModelComplexSReacdef;
struct SReacInfo {
    osh::Real kCst;
    osh::I64 charge;
};
using SReacdef = ModelSReacdef<steps::model::SReac>;
using ComplexSReacdef = ModelComplexSReacdef<steps::model::ComplexSReac>;
using VDepComplexSReacdef = ModelComplexSReacdef<steps::model::VDepComplexSReac>;

using vdep_propensity_fun_t = std::function<osh::Real(osh::Real)>;
struct VDepInfo {
    vdep_propensity_fun_t kCstFun;
    osh::I64 charge;
};
using VDepSReacdef = ModelSReacdef<steps::model::VDepSReac>;

/// Parameters of the particular GHK reaction.
struct GHKInfo {
    /// GHK current identifier
    model::ghk_current_id curr_id;
    /// if true, the surface reaction involves an ion transfer from the inner to outer compartment,
    /// and conversely otherwise.
    bool in2out;
    /// permeability per ion channel
    osh::Real permeability;
    /// ion valence
    osh::I64 charge;
    /// optional locked-in inner compartment concentration of the ion
    std::optional<osh::Real> inner_conc;
    /// optional locked-in outer compartment concentration of the ion
    std::optional<osh::Real> outer_conc;
};
class GHKSReacdef;
class GHKCurrdef;
class ComplexGHKCurrdef;
class OhmicCurrdef;
class ComplexOhmicCurrdef;
class ComplexGHKSReacdef;
class Patchdef;

class Simulation;
class Statedef;

template <typename SReacT>
struct Model2Def;

template <>
struct Model2Def<steps::model::Reac> {
    using def_type = Reacdef;
    using model_id_type = model::reaction_id;
};

template <>
struct Model2Def<steps::model::ComplexReac> {
    using def_type = ComplexReacdef;
    using model_id_type = model::complex_reaction_id;
};

template <>
struct Model2Def<steps::model::SReac> {
    using def_type = SReacdef;
    using model_id_type = model::surface_reaction_id;
};

template <>
struct Model2Def<steps::model::ComplexSReac> {
    using def_type = ComplexSReacdef;
    using model_id_type = model::complex_surface_reaction_id;
};

template <>
struct Model2Def<steps::model::VDepComplexSReac> {
    using def_type = VDepComplexSReacdef;
    using model_id_type = model::vdep_complex_surface_reaction_id;
};

template <>
struct Model2Def<steps::model::VDepSReac> {
    using def_type = VDepSReacdef;
    using model_id_type = model::vdep_surface_reaction_id;
};

template <>
struct Model2Def<steps::model::GHKcurr> {
    using def_type = GHKSReacdef;
    using curr_def_type = GHKCurrdef;
    using model_id_type = model::ghk_current_id;
};

template <>
struct Model2Def<steps::model::ComplexGHKcurr> {
    using def_type = ComplexGHKSReacdef;
    using curr_def_type = ComplexGHKCurrdef;
    using model_id_type = model::complex_ghk_current_id;
};

template <>
struct Model2Def<steps::model::OhmicCurr> {
    using curr_def_type = OhmicCurrdef;
    using model_id_type = model::ohmic_current_id;
};

template <>
struct Model2Def<steps::model::ComplexOhmicCurr> {
    using curr_def_type = ComplexOhmicCurrdef;
    using model_id_type = model::complex_ohmic_current_id;
};

}  // namespace steps::dist

#include "simulation.hpp"
#include "model/ohmiccurr.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "util/common.hpp"
#include "util/vocabulary.hpp"

#include <Omega_h_defines.hpp>
#include <limits>
#include <memory>
#include <numeric>

#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>
#include <stdexcept>
#include <string>

#if USE_PETSC
#include "mpi/dist/tetopsplit/operator/efield_operator.hpp"
#endif  // USE_PETSC

#include "geom/dist/distcomp.hpp"
#include "geom/dist/distmemb.hpp"
#include "geom/dist/distpatch.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "rng/rng.hpp"
#include "util/error.hpp"
#include "util/mesh.hpp"
#include "util/mpitools.hpp"
#include "util/profile/profiler_interface.hpp"
#include "util/tracker/time_tracker.hpp"

#undef MPI_Scatter

namespace steps::dist {

//-----------------------------------------------

Simulation::Simulation(DistMesh& t_mesh, rng::RNG& t_rng)
    : comm_rank(util::mpi_comm_rank(t_mesh.comm_impl()))
    , comm_size(util::mpi_comm_size(t_mesh.comm_impl()))
    , mesh(t_mesh)
    , rng(t_rng) {}

//-----------------------------------------------

Simulation::~Simulation() noexcept = default;

//-----------------------------------------------

///////////////////////////////////////////////////////////////////////////////////////////////
// Convenience methods
///////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////
// Location: Tetrahedron //
///////////////////////////

//-----------------------------------------------

double Simulation::getTetSpecCount(osh::GO tet, const model::species_name& s, bool local) const {
    double res;
    getBatchTetSpecCountsNP(&tet, 1, s, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setTetSpecCount(osh::GO tet,
                                 const model::species_name& s,
                                 double count,
                                 bool local) {
    setBatchTetSpecCountsNP(&tet, 1, s, &count, 1, local);
}

//-----------------------------------------------

double Simulation::getTetSpecConc(osh::GO tet, const model::species_name& s, bool local) const {
    double res;
    getBatchTetSpecConcsNP(&tet, 1, s, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setTetSpecConc(osh::GO tet,
                                const model::species_name& s,
                                double conc,
                                bool local) {
    setBatchTetSpecConcsNP(&tet, 1, s, &conc, 1, local);
}

//-----------------------------------------------

#if USE_PETSC
// E-field value

//-----------------------------------------------

double Simulation::getTetV(osh::GO tet, bool local) const {
    double res;
    getBatchTetVsNP(&tet, 1, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setTetV(osh::GO tet, double v, bool local) {
    setBatchTetVsNP(&tet, 1, &v, 1, local);
}

//-----------------------------------------------

#endif  // USE_PETSC

//-----------------------------------------------

////////////////////////
// Location: Triangle //
////////////////////////

//-----------------------------------------------

double Simulation::getTriSpecCount(osh::GO tri, const model::species_name& s, bool local) const {
    double res;
    getBatchTriSpecCountsNP(&tri, 1, s, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setTriSpecCount(osh::GO tri,
                                 const model::species_name& s,
                                 double count,
                                 bool local) {
    setBatchTriSpecCountsNP(&tri, 1, s, &count, 1, local);
}

//-----------------------------------------------

#if USE_PETSC
// E-field value

//-----------------------------------------------

double Simulation::getTriV(osh::GO tri, bool local) const {
    double res;
    getBatchTriVsNP(&tri, 1, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setTriV(osh::GO tri, double v, bool local) {
    setBatchTriVsNP(&tri, 1, &v, 1, local);
}

//-----------------------------------------------

double Simulation::getTriOhmicErev(osh::GO tri,
                                   const model::ohmic_current_id& curr,
                                   bool local) const {
    double res;
    getBatchTriOhmicErevsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriComplexOhmicErev(osh::GO tri,
                                          const model::complex_ohmic_current_id& curr,
                                          bool local) const {
    double res;
    getBatchTriComplexOhmicErevsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriSReacI(osh::GO tri,
                                const model::surface_reaction_id& reac,
                                bool local) const {
    double res;
    getBatchTriSReacIsNP(&tri, 1, reac, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriComplexSReacI(osh::GO tri,
                                       const model::complex_surface_reaction_id& reac,
                                       bool local) const {
    double res;
    getBatchTriComplexSReacIsNP(&tri, 1, reac, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriVDepSReacI(osh::GO tri,
                                    const model::vdep_surface_reaction_id& reac,
                                    bool local) const {
    double res;
    getBatchTriVDepSReacIsNP(&tri, 1, reac, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriVDepComplexSReacI(osh::GO tri,
                                           const model::vdep_complex_surface_reaction_id& reac,
                                           bool local) const {
    double res;
    getBatchTriVDepComplexSReacIsNP(&tri, 1, reac, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriOhmicI(osh::GO tri,
                                const model::ohmic_current_id& curr,
                                bool local) const {
    double res;
    getBatchTriOhmicIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriComplexOhmicI(osh::GO tri,
                                       const model::complex_ohmic_current_id& curr,
                                       bool local) const {
    double res;
    getBatchTriComplexOhmicIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriGHKI(osh::GO tri, const model::ghk_current_id& curr, bool local) const {
    double res;
    getBatchTriGHKIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriComplexGHKI(osh::GO tri,
                                     const model::complex_ghk_current_id& curr,
                                     bool local) const {
    double res;
    getBatchTriComplexGHKIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

//-----------------------------------------------

double Simulation::getTriI(osh::GO tri, bool local) const {
    double res;
    getBatchTriIsNP(&tri, 1, &res, 1, local);
    return res;
}

//-----------------------------------------------

#endif  // USE_PETSC

//////////////////////
// Location: Vertex //
//////////////////////

#if USE_PETSC
// E-field value

//-----------------------------------------------

double Simulation::getVertV(osh::GO vert, bool local) const {
    double res;
    getBatchVertVsNP(&vert, 1, &res, 1, local);
    return res;
}

//-----------------------------------------------

void Simulation::setVertV(osh::GO vert, double v, bool local) {
    setBatchVertVsNP(&vert, 1, &v, 1, local);
}

//-----------------------------------------------

#endif  // USE_PETSC


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::OmegaHSimulation(steps::model::Model& model,
                                                                  DistMesh& t_mesh,
                                                                  const rng::RNGptr& r,
                                                                  bool t_indepKProcs,
                                                                  bool isEfield)
    : super_type(t_mesh, *r)
    , mesh(t_mesh)
    , indepKProcs(t_indepKProcs) {
    auto stateDefPtr = std::make_unique<Statedef>(model, mesh);
    if (!isEfield) {
        stateDefPtr->disableEField();
    }

    init(std::move(stateDefPtr));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::~OmegaHSimulation() = default;

//-----------------------------------------------

///////////////////////////////////////////////////////////////////////////////////////////////
// Required methods
///////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////
// General methods //
/////////////////////

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::reset() {
    this->num_iterations = 0;
    this->state_time = 0;
    statedef->reset();
    mesh.resetDiffBoundaries();
    input->reset(*statedef, mesh);
    data->reset(this->state_time);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::run(osh::Real end_time) {
    Instrumentor::phase p("OmegaHSimulation::run()");

    if (allReduce(outdated_diffusions, MPI_LOR)) {
        data->initialize_diffusions();
        outdated_diffusions = false;
    }

    assert(end_time >= 0.0);

#if USE_PETSC
    const osh::Real ef_dt_std = data->efield ? data->efield->getDt()
                                             : std::numeric_limits<double>::infinity();
#else
    const osh::Real ef_dt_std = std::numeric_limits<double>::infinity();
#endif


    // number of standard steps -1. It can be negative
    const int n_steps_std = std::floor((end_time - state_time) / ef_dt_std) - 1;
    // std steps loop -1. n_steps_std can be negative so int is the correct type
    for (int i_ef = 0; i_ef < n_steps_std; ++i_ef) {
        evolve(ef_dt_std);
    }

    // sync time steps
    // we compare always with end_time because comparing dts can fail due to numerical error
    const auto next_std_state_time = state_time + ef_dt_std;
    if (!steps::util::almost_equal(next_std_state_time, end_time) &&
        next_std_state_time < end_time) {
        evolve(ef_dt_std);
        assert(end_time > state_time);
        assert(!steps::util::almost_equal(end_time, state_time));
        evolve(end_time - state_time);
    } else if (!steps::util::almost_equal(end_time, state_time)) {
        evolve(end_time - state_time);
    }

    assert(steps::util::almost_equal(end_time, state_time));
}

//-----------------------------------------------

// Data getting / setting

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionTolerance() const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getTolerance();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getDiffusionTolerance.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionTolerance(osh::Real tolerance) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setTolerance(tolerance);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setDiffusionTolerance.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionNormalApproximationThreshold() const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getNormalApproximationThreshold();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling "
            "getDiffusionNormalApproximationThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionNormalApproximationThreshold(
    osh::Real threshold) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setNormalApproximationThreshold(threshold);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling "
            "setDiffusionNormalApproximationThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionCrankNicolsonThreshold()
    const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getCrankNicolsonThreshold();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling "
            "getDiffusionCrankNicolsonThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionCrankNicolsonThreshold(
    osh::Real threshold) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setCrankNicolsonThreshold(threshold);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling "
            "setDiffusionCrankNicolsonThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
uint OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionLeapThreshold() const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getLeapThreshold();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getDiffusionLeapThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionLeapThreshold(uint leap_thresh) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setLeapThreshold(leap_thresh);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setDiffusionLeapThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
uint OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionMaxDtSkips() const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getMaxDtSkips();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getDiffusionMaxDtSkips.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionMaxDtSkips(uint max_skips) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setMaxDtSkips(max_skips);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setDiffusionMaxDtSkips.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionMinDtFactor() const {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        return data->diffOp.getMinDtFactor();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getDiffusionMinDtFactor.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffusionMinDtFactor(osh::Real factor) {
    if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        data->diffOp.setMinDtFactor(factor);
        outdated_diffusions = true;
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setDiffusionMinDtFactor.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
uint OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionSSAThreshold() const {
    if constexpr (SSA == SSAMethod::RLeaping) {
        return data->ssaOp.getSSAThreshold();
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling getReactionSSAThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setReactionSSAThreshold(uint thresh) {
    if constexpr (SSA == SSAMethod::RLeaping) {
        data->ssaOp.setSSAThreshold(thresh);
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling setReactionSSAThreshold.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
uint OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionSSASteps() const {
    if constexpr (SSA == SSAMethod::RLeaping) {
        return data->ssaOp.getSSASteps();
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling getReactionSSASteps.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setReactionSSASteps(uint steps) {
    if constexpr (SSA == SSAMethod::RLeaping) {
        data->ssaOp.setSSASteps(steps);
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling setReactionSSASteps.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
uint OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionLComputePeriod() const {
    if constexpr (SSA == SSAMethod::RLeaping) {
        return data->ssaOp.getLComputePeriod();
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling getReactionLComputePeriod.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setReactionLComputePeriod(uint period) {
    if constexpr (SSA == SSAMethod::RLeaping) {
        data->ssaOp.setLComputePeriod(period);
    } else {
        throw std::logic_error(
            "The chosen reaction operator does not allow calling setReactionLComputePeriod.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionTolerance() const {
    if constexpr (SSA == SSAMethod::RLeaping) {
        return data->ssaOp.getTolerance();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getReactionTolerance.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setReactionTolerance(osh::Real tolerance) {
    if constexpr (SSA == SSAMethod::RLeaping) {
        data->ssaOp.setTolerance(tolerance);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setReactionTolerance.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionTheta() const {
    if constexpr (SSA == SSAMethod::RLeaping) {
        return data->ssaOp.getTheta();
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling getReactionTheta.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setReactionTheta(osh::Real theta) {
    if constexpr (SSA == SSAMethod::RLeaping) {
        data->ssaOp.setTheta(theta);
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setReactionTheta.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
std::string OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getSolverName() const {
    return "disttetopsplit";
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffApplyThreshold(osh::Real threshold) {
    if constexpr (DiffMethod == DiffusionMethod::ConstantDiffDt) {
        data->diffOp.setBinomialThreshold(static_cast<osh::GO>(threshold));
    } else {
        throw std::logic_error(
            "The chosen diffusion method does not allow calling setDiffApplyThreshold.");
    }
}

//-----------------------------------------------

#if USE_PETSC
// E-field specific

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getEfieldDT() const {
    if (data->efield) {
        return data->efield->getDt();
    }
    return 0;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setEfieldDT(const osh::Real dt) const {
    if (data->efield) {
        data->efield->setDt(dt);
    } else {
        throw std::logic_error("E-Field is not in use.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setPetscOptions(const std::string& s) {
    auto err = PetscOptionsInsertString(nullptr, s.c_str());
    CHKERRABORT(mesh.comm_impl(), err);

    // the check to see if we have the efield active is done in simulation
    if (data->efield) {
        data->efield->setPetscOptions();
    } else {
        throw std::logic_error("E-Field is not in use.");
    }
}

//-----------------------------------------------

#endif  // USE_PETSC

// Debugging / monitoring

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::dumpDepGraphToFile(
    const std::string& path) const {
    std::ofstream ostr(path);
    data->kproc_state.write_dependency_graph(ostr);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
std::string OmegaHSimulation<SSA, SearchMethod, DiffMethod>::createStateReport() const {
    // TODO(TCL) FIXME
    return "";
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffExtent(bool local) const {
    osh::I64 extent = data->diffOp.getExtent();
    if (local) {
        return extent;
    }
    return allReduce(extent);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReacExtent(bool local) const {
    const osh::I64 extent = data->ssaOp.getExtent();
    if (local) {
        return extent;
    }
    return allReduce(extent);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
std::map<std::string, double> OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getReactionDebugInfo(
    bool local) const {
    const auto info = data->ssaOp.getDebugInfo();
    if (local) {
        return info;
    }
    std::map<std::string, double> total_info;
    for (auto& [name, val]: info) {
        total_info[name] = allReduce(val);
    }
    return total_info;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
std::map<std::string, double>
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffusionDebugInfo(bool local) const {
    const auto info = data->diffOp.getDebugInfo();
    if (local) {
        return info;
    }
    std::map<std::string, double> total_info;
    for (auto& [name, val]: info) {
        total_info[name] = allReduce(val);
    }
    return total_info;
}

//-----------------------------------------------

///////////////////////////
// Location: Compartment //
///////////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    return allReduce(getOwnedCompSpecCount(compartment, species));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& spec,
    osh::Real n,
    const math::DistributionMethod distribution) {
    const auto& [elems, volumes, rank_volume] = mesh.measure(compartment);
    std::vector<osh::Real> rank_volumes;
    if (this->comm_rank == 0) {
        rank_volumes.resize(static_cast<size_t>(this->comm_size));
    }
    int err = MPI_Gather(
        &rank_volume, 1, MPI_DOUBLE, rank_volumes.data(), 1, MPI_DOUBLE, 0, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    if (n >= static_cast<double>(INT64_MAX)) {
        std::ostringstream oss;
        oss << "Unsupported number of molecules: " << std::setprecision(20) << n
            << " but maximum value is " << std::numeric_limits<osh::GO>::max()
            << " (max 64 bits integral value)";
        ArgErrLog(oss.str());
    }
    auto dist = math::make_dist(static_cast<osh::GO>(n), rank_volumes);
    // only rank 0 generates non-zero values
    const std::vector<osh::GO>& num_molecules_on_ranks = dist.distribute(this->rng, distribution);


    // send `num_molecules_on_ranks[i]` to rank `i` and store it in
    // `num_molecules_on_rank`.
    osh::GO num_molecules_on_rank;
    err = MPI_Scatter(num_molecules_on_ranks.data(),
                      1,
                      MPI_INT64_T,
                      &num_molecules_on_rank,
                      1,
                      MPI_INT64_T,
                      0,
                      this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    setOwnedCompSpecCount(compartment, spec, num_molecules_on_rank, distribution);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    const auto spec_count = getCompSpecCount(compartment, species);
    return spec_count / (1.0e3 * mesh.total_measure(compartment) * math::AVOGADRO);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& spec,
    osh::Real conc,
    const math::DistributionMethod distribution) {
    const auto factor = mesh.total_measure(compartment) * 1.0e3 * math::AVOGADRO;
    setCompSpecCount(compartment, spec, conc * factor, distribution);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompSpecClamped(
    const model::compartment_id& compartment,
    const model::species_name& spec) const {
    auto& compdef = statedef->getCompdef(compartment);
    const auto spec_id = statedef->getCompSpecContainerIdx(compartment, spec);
    return compdef.getSpecClamped(spec_id);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompSpecClamped(
    const model::compartment_id& compartment,
    const model::species_name& spec,
    bool clamped) {
    auto& compdef = statedef->getCompdef(compartment);
    const auto spec_id = statedef->getCompSpecContainerIdx(compartment, spec);
    compdef.setSpecClamped(spec_id, clamped);
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompReacK(
    const model::compartment_id& compartment,
    const model::reaction_id& reac) const {
    auto& compdef = statedef->getCompdef(compartment);
    auto& reacdef = compdef.template getReacdef<Reacdef>(reac);
    return reacdef.getKcst();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompReacK(
    const model::compartment_id& compartment,
    const model::reaction_id& reac,
    osh::Real kcst) {
    auto& compdef = statedef->getCompdef(compartment);
    auto& reacdef = compdef.template getReacdef<Reacdef>(reac);
    data->kproc_state.reactions().clearKcst(compdef.getIdx());
    reacdef.setKcst(kcst);
    data->kproc_state.reactions().updateKcst();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompReacExtent(
    const model::compartment_id& compartment,
    const model::reaction_id& reac) const {
    auto& compdef = statedef->getCompdef(compartment);
    auto& reacdef = compdef.template getReacdef<Reacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompComplexReacExtent(
    const model::compartment_id& compartment,
    const model::complex_reaction_id& reac) const {
    auto& compdef = statedef->getCompdef(compartment);
    auto& reacdef = compdef.template getReacdef<ComplexReacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompDiffD(
    const model::compartment_id& compartment,
    const model::diffusion_id& diff) const {
    auto& diffusion = statedef->getCompdef(compartment).getDiffdef(diff);
    return diffusion.getDcst();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompDiffD(
    const model::compartment_id& compartment,
    const model::diffusion_id& diff,
    osh::Real dcst) {
    auto& diffusion = statedef->getCompdef(compartment).getDiffdef(diff);
    diffusion.setDcst(dcst);
    // Clear individually set values (with setTetDiffD)
    auto diffId = diffusion.getDiffContainerIdx();
    for (auto tet: mesh.getOwnedEntities(compartment)) {
        data->diffusions.clear_tet_dcst(tet, diffId, -1);
        for (auto& d: mesh.tet_neighbors_int_data()[tet.get()]) {
            data->diffusions.clear_tet_dcst(tet, diffId, d[1]);
        }
    }
    data->initialize_diffusions();
}

//-----------------------------------------------

// Complexes
template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompComplexCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const {
    return allReduce(getOwnedCompComplexCount(compartment, complex, _convertComplexFilters(f)));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setCompComplexCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
    osh::Real num_molecules,
    math::DistributionMethod distribution) {
    const auto& [elems, volumes, rank_volume] = mesh.measure(compartment);
    std::vector<osh::Real> rank_volumes;
    if (this->comm_rank == 0) {
        rank_volumes.resize(static_cast<size_t>(this->comm_size));
    }
    int err = MPI_Gather(
        &rank_volume, 1, MPI_DOUBLE, rank_volumes.data(), 1, MPI_DOUBLE, 0, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    auto dist = math::make_dist(static_cast<osh::GO>(num_molecules), rank_volumes);
    // only rank 0 generates non-zero values
    const std::vector<osh::GO>& num_molecules_on_ranks = dist.distribute(this->rng, distribution);

    // send `num_molecules_on_ranks[i]` to rank `i` and store it in
    // `num_molecules_on_rank`.
    osh::GO num_molecules_on_rank;
    err = MPI_Scatter(num_molecules_on_ranks.data(),
                      1,
                      MPI_INT64_T,
                      &num_molecules_on_rank,
                      1,
                      MPI_INT64_T,
                      0,
                      this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    setOwnedCompComplexCount(
        compartment, complex, _convertComplexState(i), num_molecules_on_rank, distribution);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getCompComplexSUSCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
    model::complex_substate_id m) const {
    return allReduce(
        getOwnedCompComplexSUSCount(compartment, complex, _convertComplexFilters(f), m));
}

//-----------------------------------------------

/////////////////////
// Location: Patch //
/////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchSpecCount(
    const model::patch_id& patch,
    const model::species_name& species) const {
    const auto& boundaries = mesh.getOwnedEntities(patch);
    osh::Write<osh::LO> mols_counts(boundaries.size());
    const container::species_id spec_id{statedef->getPatchdef(patch).getSpecContainerIdx(species)};
    const auto& molecules = data->pools.moleculesOnPatchBoundaries();

    std::transform(boundaries.begin(), boundaries.end(), mols_counts.begin(), [&](auto bound) {
        return molecules(bound, spec_id);
    });

    return mesh.get_MPI_sum(osh::LOs(mols_counts));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setPatchSpecCount(
    const model::patch_id& patch,
    const model::species_name& species,
    osh::Real num_molecules,
    const math::DistributionMethod distribution) {
    const auto& [elems, areas, rank_area] = mesh.measure(patch);
    osh::GO num_molecules_on_rank;

    {
        std::vector<osh::Real> rank_areas;
        if (this->comm_rank == 0) {
            rank_areas.resize(static_cast<size_t>(this->comm_size));
        }
        int err = MPI_Gather(
            &rank_area, 1, MPI_DOUBLE, rank_areas.data(), 1, MPI_DOUBLE, 0, this->comm());

        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }

        std::vector<osh::GO> num_molecules_on_ranks;
        auto dist = math::make_dist(static_cast<osh::GO>(num_molecules), rank_areas);
        num_molecules_on_ranks = dist.distribute(this->rng, distribution);

        // send `num_molecules_on_ranks[i]` to rank `i` and store it in
        // `num_molecules_on_rank`.
        err = MPI_Scatter(num_molecules_on_ranks.data(),
                          1,
                          MPI_UNSIGNED_LONG,
                          &num_molecules_on_rank,
                          1,
                          MPI_UNSIGNED_LONG,
                          0,
                          this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }

    osh::Write<molecules_t> mols_on_elements;
    auto dist = math::make_dist(num_molecules_on_rank, areas);
    mols_on_elements = dist.distribute(this->rng, distribution);

    container::species_id cont_spec_id = statedef->getPatchdef(patch).getSpecContainerIdx(species);
    for (auto k = 0; k < mols_on_elements.size(); ++k) {
        const mesh::triangle_id_t boundary(elems[k]);
        const auto molecules = mols_on_elements[k];
        this->data->pools.assign(boundary, cont_spec_id, molecules);
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchSpecClamped(
    const model::patch_id& patch,
    const model::species_name& spec) const {
    auto& patchdef = statedef->getPatchdef(patch);
    const auto spec_id = patchdef.getSpecContainerIdx(spec);
    return patchdef.getSpecClamped(spec_id);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setPatchSpecClamped(
    const model::patch_id& patch,
    const model::species_name& spec,
    bool clamped) {
    auto& patchdef = statedef->getPatchdef(patch);
    const auto spec_id = patchdef.getSpecContainerIdx(spec);
    patchdef.setSpecClamped(spec_id, clamped);
}

//-----------------------------------------------

// Complexes

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchComplexCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) const {
    return allReduce(getOwnedPatchComplexCount(patch, complex, _convertComplexFilters(f)));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setPatchComplexCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& i,
    osh::Real num_molecules,
    math::DistributionMethod distribution) {
    const auto& [elems, areas, rank_area] = mesh.measure(patch);
    std::vector<osh::Real> rank_areas;
    if (this->comm_rank == 0) {
        rank_areas.resize(static_cast<size_t>(this->comm_size));
    }
    int err =
        MPI_Gather(&rank_area, 1, MPI_DOUBLE, rank_areas.data(), 1, MPI_DOUBLE, 0, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    auto dist = math::make_dist(static_cast<osh::GO>(num_molecules), rank_areas);
    // only rank 0 generates non-zero values
    const std::vector<osh::GO>& num_molecules_on_ranks = dist.distribute(this->rng, distribution);

    // send `num_molecules_on_ranks[i]` to rank `i` and store it in
    // `num_molecules_on_rank`.
    osh::GO num_molecules_on_rank;
    err = MPI_Scatter(num_molecules_on_ranks.data(),
                      1,
                      MPI_INT64_T,
                      &num_molecules_on_rank,
                      1,
                      MPI_INT64_T,
                      0,
                      this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    setOwnedPatchComplexCount(
        patch, complex, _convertComplexState(i), num_molecules_on_rank, distribution);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchComplexSUSCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f,
    model::complex_substate_id m) const {
    return allReduce(getOwnedPatchComplexSUSCount(patch, complex, _convertComplexFilters(f), m));
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchSReacK(
    const model::patch_id& patchId,
    const model::surface_reaction_id& reactionId) const {
    auto& patchdef = statedef->getPatchdef(patchId);
    auto reacId = patchdef.getReacIdx(reactionId);
    auto& reacdef = *patchdef.template reacdefs<SReacdef>().at(reacId.get());
    return reacdef.getInfo().kCst;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setPatchSReacK(
    const model::patch_id& patchId,
    const model::surface_reaction_id& reactionId,
    osh::Real kCst) {
    auto& patchdef = statedef->getPatchdef(patchId);
    auto reacId = patchdef.getReacIdx(reactionId);
    auto& reacdef = *patchdef.template reacdefs<SReacdef>().at(reacId.get());
    // Clear values set per triangle (with setTriSReacK)
    data->kproc_state.surfaceReactions().clear_Kcst(patchdef.getIdx());

    reacdef.getInfo().kCst = kCst;
    data->kproc_state.surfaceReactions().updateCcst();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchSReacExtent(
    const model::patch_id& patch,
    const model::surface_reaction_id& reac) const {
    auto& patchdef = statedef->getPatchdef(patch);
    auto& reacdef = patchdef.template getReacdef<SReacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchComplexSReacExtent(
    const model::patch_id& patch,
    const model::complex_surface_reaction_id& reac) const {
    auto& patchdef = statedef->getPatchdef(patch);
    auto& reacdef = patchdef.template getReacdef<ComplexSReacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchVDepSReacExtent(
    const model::patch_id& patch,
    const model::vdep_surface_reaction_id& reac) const {
    auto& patchdef = statedef->getPatchdef(patch);
    auto& reacdef = patchdef.template getReacdef<VDepSReacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getPatchVDepComplexSReacExtent(
    const model::patch_id& patch,
    const model::vdep_complex_surface_reaction_id& reac) const {
    auto& patchdef = statedef->getPatchdef(patch);
    auto& reacdef = patchdef.template getReacdef<VDepComplexSReacdef>(reac);
    return allReduce(reacdef.getExtent());
}

//-----------------------------------------------

////////////////////////
// Location: Membrane //
////////////////////////

#if USE_PETSC
// E-field values

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
MembraneResistivity OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getMembRes(
    const model::membrane_id& membrane) const {
    return {statedef->getResistivity(membrane), statedef->getReversalPotential(membrane)};
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setMembRes(const model::membrane_id& membrane,
                                                                 osh::Real resistivity,
                                                                 osh::Real reversal_potential) {
    statedef->setResistivity(membrane, resistivity);
    statedef->setReversalPotential(membrane, reversal_potential);
    for (const auto& patchId: statedef->getMembrane(membrane).getPatchesIds()) {
        for (const auto tri: mesh.getOwnedEntities(patchId)) {
            input->conductivity_on_triangles_w[tri.get()] = 1.0 / resistivity;
            input->reversal_potential_on_triangles_w[tri.get()] = reversal_potential;
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setMembVolRes(const model::membrane_id& memb,
                                                                    osh::Real ro) {
    auto& membrane = statedef->getMembrane(memb);
    membrane.setConductivity(1.0 / ro);
    for (const auto& patch: membrane.getPatchesIds()) {
        auto compid = statedef->getPatchdef(patch).getInnerCompId();
        auto& comp = statedef->getCompdef(compid);
        comp.setConductivity(1.0 / ro);
    }
    if (data->efield) {
        data->efield->resetStiffnessMatrix();
    } else {
        throw std::logic_error("E-Field is not in use, cannot set VolRes.");
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setMembCapac(const model::membrane_id& memb,
                                                                   osh::Real capacitance) {
    auto& membrane = statedef->getMembrane(memb);
    membrane.setCapacitance(capacitance);
    for (const auto& patch: membrane.getPatchesIds()) {
        for (const auto tri: mesh.getOwnedEntities(patch)) {
            input->capacitance_on_triangles_w[tri.get()] = capacitance;
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setMembPotential(
    const model::membrane_id& memb,
    osh::Real value) {
    auto membit = mesh.membranes().find(memb);
    if (membit == mesh.membranes().end()) {
        throw std::invalid_argument("Invalid membrane " + memb);
    }
    const auto& allPatches = mesh.getAllPatches();
    for (const auto& patchid: membit->second->patches()) {
        auto pmeshid = mesh.getPatchID(patchid);
        const auto* patch = allPatches[pmeshid];
        const auto* icomp = dynamic_cast<const DistComp*>(&patch->getIComp());
        if (icomp == nullptr) {
            continue;
        }
        for (const auto tet: icomp->getLocalTetIndices(false)) {
            const auto verts = osh::gather_verts<4>(mesh.ask_elem_verts(), tet.get());
            for (const auto& vert: verts) {
                input->potential_on_vertices_w[vert] = value;
            }
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setMembIClamp(
    const model::membrane_id& membrane,
    osh::Real current) {
    statedef->setStimulus(membrane, current);
}

//-----------------------------------------------

#endif  // USE_PETSC

//////////////////////////////////
// Location: Diffusion boundary //
//////////////////////////////////

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getDiffBoundarySpecDiffusionActive(
    const mesh::diffusion_boundary_name& diffusion_boundary_name,
    const model::species_name& spec_id) const {
    model::species_id mdl_spec_id = statedef->getSpecModelIdx(spec_id);
    DistMesh::DiffusionBoundary& db = mesh.getDiffusionBoundary(diffusion_boundary_name);
    Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
    Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
    container::species_id sp1 = comp1.getSpecContainerIdx(mdl_spec_id);
    container::species_id sp2 = comp2.getSpecContainerIdx(mdl_spec_id);
    return db.comp1_spec2dcst[sp1.get()] != 0.0 && db.comp2_spec2dcst[sp2.get()] != 0.0;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffBoundarySpecDiffusionActive(
    const mesh::diffusion_boundary_name& diffusion_boundary_name,
    const model::species_name& spec_id,
    bool set_active) {
    model::species_id mdl_spec_id = statedef->getSpecModelIdx(spec_id);
    DistMesh::DiffusionBoundary& db = mesh.getDiffusionBoundary(diffusion_boundary_name);
    Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
    Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
    container::species_id sp1 = comp1.getSpecContainerIdx(mdl_spec_id);
    container::species_id sp2 = comp2.getSpecContainerIdx(mdl_spec_id);
    db.comp1_spec2dcst[sp1.get()] = set_active ? -1.0 : 0.0;
    db.comp2_spec2dcst[sp2.get()] = set_active ? -1.0 : 0.0;
    data->initialize_diffusions();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setDiffBoundarySpecDcst(
    const mesh::diffusion_boundary_name& diffb,
    const model::species_name& spec,
    osh::Real dcst) {
    model::species_id mdl_spec_id = statedef->getSpecModelIdx(spec);
    DistMesh::DiffusionBoundary& db = mesh.getDiffusionBoundary(diffb);
    Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
    Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
    container::species_id sp1 = comp1.getSpecContainerIdx(mdl_spec_id);
    container::species_id sp2 = comp2.getSpecContainerIdx(mdl_spec_id);
    db.comp1_spec2dcst[sp1.get()] = dcst;
    db.comp2_spec2dcst[sp2.get()] = dcst;
    data->initialize_diffusions();
}

//-----------------------------------------------

///////////////////////////
// Location: Tetrahedron //
///////////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTetSpecCountsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* counts,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchElemValsNP(indices, input_size, s, counts, false, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchTetSpecCountsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* counts,
    size_t output_size,
    bool local) {
    assert(input_size == output_size);
    (void) output_size;
    setBatchElemValsNP(indices, input_size, s, counts, false, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTetSpecConcsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* counts,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchElemValsNP(indices, input_size, s, counts, true, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchTetSpecConcsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* concs,
    size_t output_size,
    bool local) {
    assert(input_size == output_size);
    (void) output_size;
    setBatchElemValsNP(indices, input_size, s, concs, true, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetSpecClamped(
    osh::GO tet,
    const model::species_name& s,
    bool local) const {
    bool clamped = false;
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    if (localInd.valid()) {
        const auto compartment_id = mesh.getCompartment(localInd);
        const auto spec_id = statedef->getCompSpecContainerIdx(compartment_id, s);
        clamped = data->pools.get_clamped(localInd, spec_id);
    }
    return allReduce(clamped, MPI_LOR);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetSpecClamped(
    osh::GO tet,
    const model::species_name& s,
    bool clamped,
    bool local) {
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    if (localInd.valid()) {
        const auto compartment_id = mesh.getCompartment(localInd);
        const auto spec_id = statedef->getCompSpecContainerIdx(compartment_id, s);
        data->pools.set_clamped(localInd, spec_id, clamped);
    }
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetReacK(
    osh::GO tet,
    const model::reaction_id reac,
    bool local) const {
    return getTetReacK(data->kproc_state.reactions(), tet, reac, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetReacK(osh::GO tet,
                                                                  const model::reaction_id reac,
                                                                  osh::Real kcst,
                                                                  bool local) {
    setTetReacK(data->kproc_state.reactions(), tet, reac, kcst, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetComplexReacK(
    osh::GO tet,
    const model::complex_reaction_id reac,
    bool local) const {
    return getTetReacK(data->kproc_state.complexReactions(), tet, reac, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetComplexReacK(
    osh::GO tet,
    const model::complex_reaction_id reac,
    osh::Real kcst,
    bool local) {
    setTetReacK(data->kproc_state.complexReactions(), tet, reac, kcst, local);
}

//-----------------------------------------------

// Diffusions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetDiffD(
    osh::GO tet,
    const model::diffusion_id diff,
    osh::GO direc_tet,
    bool local) const {
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    osh::Real dcst = 0.0;
    if (localInd.valid()) {
        const auto compartment_id = mesh.getCompartment(localInd);
        const auto& compartment = statedef->getCompdef(compartment_id);
        const auto& diffusion = compartment.getDiffdef(diff);
        auto direc_localInd = getLocalInd<mesh::tetrahedron_local_id_t>(direc_tet, local, true);
        if (direc_localInd.valid()) {
            int face_idx = -1;
            for (auto& d: mesh.tet_neighbors_int_data()[localInd.get()]) {
                if (d[0] == direc_localInd.get()) {
                    face_idx = d[1];
                    break;
                }
            }
            if (face_idx == -1) {
                throw std::invalid_argument("Direction tetrahedron " + std::to_string(direc_tet) +
                                            " is not a neighbor of " + std::to_string(tet));
            }
            dcst =
                data->diffusions.get_tet_dcst(localInd, diffusion.getDiffContainerIdx(), face_idx);
        } else {
            dcst = data->diffusions.get_tet_dcst(localInd, diffusion.getDiffContainerIdx(), -1);
        }
    }
    return allReduce(dcst);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetDiffD(osh::GO tet,
                                                                  const model::diffusion_id diff,
                                                                  double dcst,
                                                                  osh::GO direc_tet,
                                                                  bool local) {
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    if (localInd.valid()) {
        const auto compartment_id = mesh.getCompartment(localInd);
        const auto& compartment = statedef->getCompdef(compartment_id);
        const auto& diffusion = compartment.getDiffdef(diff);
        auto direc_localInd = getLocalInd<mesh::tetrahedron_local_id_t>(direc_tet, local, true);
        if (direc_localInd.valid()) {
            bool found = false;
            for (auto& d: mesh.tet_neighbors_int_data()[localInd.get()]) {
                if (d[0] == direc_localInd.get()) {
                    data->diffusions.set_tet_dcst(localInd,
                                                  diffusion.getDiffContainerIdx(),
                                                  d[1],
                                                  dcst);
                    found = true;
                    break;
                }
            }
            if (not found) {
                throw std::invalid_argument("Direction tetrahedron " + std::to_string(direc_tet) +
                                            " is not a neighbor of " + std::to_string(tet));
            }
        } else {
            data->diffusions.set_tet_dcst(localInd, diffusion.getDiffContainerIdx(), -1, dcst);
        }
    }
    outdated_diffusions = true;
}

//-----------------------------------------------

#if USE_PETSC
// E-field value

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTetVsNP(const osh::GO* indices,
                                                                      size_t input_size,
                                                                      osh::Real* voltages,
                                                                      size_t output_size,
                                                                      bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    std::fill(voltages, voltages + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto tet2verts = osh::gather_verts<4>(mesh.ask_elem_verts(), localInd.get());
            for (auto vert: tet2verts) {
                voltages[i] += input->potential_on_vertices_w[vert] / 4.0;
            }
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, voltages, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchTetVsNP(const osh::GO* indices,
                                                                      size_t input_size,
                                                                      osh::Real* voltages,
                                                                      size_t output_size,
                                                                      bool local) {
    assert(input_size == output_size);
    (void) output_size;
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto tet2verts = osh::gather_verts<4>(mesh.ask_elem_verts(), localInd.get());
            for (auto vert: tet2verts) {
                input->potential_on_vertices_w[vert] = voltages[i];
            }
        }
    }
    const auto& syncedv = mesh.sync_array(osh::VERT, osh::Reals(input->potential_on_vertices_w), 1);
    std::copy(syncedv.begin(), syncedv.end(), input->potential_on_vertices_w.begin());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetVClamped(osh::GO vertex,
                                                                     bool local) const {
    auto tet = getLocalInd<mesh::tetrahedron_local_id_t>(vertex, local, true);
    bool clamped = false;
    if (tet.valid()) {
        if (data->efield) {
            const auto& tets2verts = mesh.ask_elem_verts();
            const auto verts = osh::gather_verts<4>(tets2verts, tet.get());
            clamped = true;
            for (auto v: verts) {
                clamped &= data->efield->getVertVClamped(mesh::vertex_local_id_t(v));
            }
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
    return allReduce(clamped, MPI_LOR);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetVClamped(osh::GO vertex,
                                                                     bool clamped,
                                                                     bool local) {
    auto tet = getLocalInd<mesh::tetrahedron_local_id_t>(vertex, local, true);
    if (tet.valid()) {
        if (data->efield) {
            const auto& tets2verts = mesh.ask_elem_verts();
            const auto verts = osh::gather_verts<4>(tets2verts, tet.get());
            for (auto v: verts) {
                data->efield->setVertVClamped(mesh::vertex_local_id_t(v), clamped);
            }
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
}

//-----------------------------------------------

#endif  // USE_PETSC

////////////////////////
// Location: Triangle //
////////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriSpecCountsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* counts,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchBoundSpecCountNP(indices, input_size, s, counts, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchTriSpecCountsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& s,
    double* counts,
    size_t output_size,
    bool local) {
    assert(input_size == output_size);
    (void) output_size;
    setBatchBoundSpecCountNP(indices, input_size, s, counts, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriSpecClamped(
    osh::GO tri,
    const model::species_name& s,
    bool local) const {
    bool clamped = false;
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, true);
    if (localInd.valid()) {
        const auto patch_id = model::patch_id(mesh.getTriPatch(localInd)->getID());
        auto spec_id = statedef->getPatchdef(patch_id).getSpecContainerIdx(s);
        clamped = data->pools.get_clamped(localInd, spec_id);
    }
    return allReduce(clamped, MPI_LOR);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriSpecClamped(
    osh::GO tri,
    const model::species_name& s,
    bool clamped,
    bool local) {
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, true);
    if (localInd.valid()) {
        const auto patch_id = model::patch_id(mesh.getTriPatch(localInd)->getID());
        auto spec_id = statedef->getPatchdef(patch_id).getSpecContainerIdx(s);
        data->pools.set_clamped(localInd, spec_id, clamped);
    }
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriSReacK(
    osh::GO triangle,
    const model::surface_reaction_id& reactionId,
    bool local) const {
    return getTriSReacK(data->kproc_state.surfaceReactions(), triangle, reactionId, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriSReacK(
    osh::GO triangle,
    const model::surface_reaction_id& reactionId,
    osh::Real kCst,
    bool local) {
    setTriSReacK(data->kproc_state.surfaceReactions(), triangle, reactionId, kCst, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriComplexSReacK(
    osh::GO triangle,
    const model::complex_surface_reaction_id& reactionId,
    bool local) const {
    return getTriSReacK(data->kproc_state.complexSurfaceReactions(), triangle, reactionId, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriComplexSReacK(
    osh::GO triangle,
    const model::complex_surface_reaction_id& reactionId,
    osh::Real kCst,
    bool local) {
    setTriSReacK(data->kproc_state.complexSurfaceReactions(), triangle, reactionId, kCst, local);
}

//-----------------------------------------------

#if USE_PETSC
// E-field value

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriVsNP(const osh::GO* indices,
                                                                      size_t input_size,
                                                                      osh::Real* voltages,
                                                                      size_t output_size,
                                                                      bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    std::fill(voltages, voltages + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto tri2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                        localInd.get());
            for (auto vert: tri2verts) {
                voltages[i] += input->potential_on_vertices_w[vert] / 3.0;
            }
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, voltages, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchTriVsNP(const osh::GO* indices,
                                                                      size_t input_size,
                                                                      osh::Real* voltages,
                                                                      size_t output_size,
                                                                      bool local) {
    assert(input_size == output_size);
    (void) output_size;
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto tri2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                        localInd.get());
            for (auto vert: tri2verts) {
                input->potential_on_vertices_w[vert] = voltages[i];
            }
        }
    }
    const auto& syncedv = mesh.sync_array(osh::VERT, osh::Reals(input->potential_on_vertices_w), 1);
    std::copy(syncedv.begin(), syncedv.end(), input->potential_on_vertices_w.begin());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriSReacIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::surface_reaction_id reac,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    getBatchTriSReacIsNP(
        data->kproc_state.surfaceReactions(), indices, input_size, reac, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriComplexSReacIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::complex_surface_reaction_id reac,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    getBatchTriSReacIsNP(
        data->kproc_state.complexSurfaceReactions(), indices, input_size, reac, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriVDepSReacIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::vdep_surface_reaction_id reac,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    getBatchTriSReacIsNP(
        data->kproc_state.vDepSurfaceReactions(), indices, input_size, reac, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriVDepComplexSReacIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::vdep_complex_surface_reaction_id reac,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    getBatchTriSReacIsNP(data->kproc_state.vDepComplexSurfaceReactions(),
                         indices,
                         input_size,
                         reac,
                         currents,
                         local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriOhmicIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::ohmic_current_id curr,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    (void) output_size;
    assert(input_size == output_size);
    getBatchTriOhmicIsNP<OhmicCurrdef>(indices, input_size, curr, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriComplexOhmicIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::complex_ohmic_current_id curr,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchTriOhmicIsNP<ComplexOhmicCurrdef>(indices, input_size, curr, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriGHKIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::ghk_current_id curr,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchTriGHKIsNP(
        data->kproc_state.ghkSurfaceReactions(), indices, input_size, curr, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriComplexGHKIsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::complex_ghk_current_id curr,
    osh::Real* currents,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchTriGHKIsNP(
        data->kproc_state.complexGhkSurfaceReactions(), indices, input_size, curr, currents, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriIsNP(const osh::GO* indices,
                                                                      size_t input_size,
                                                                      osh::Real* currents,
                                                                      size_t /*output_size*/,
                                                                      bool local) const {
    std::fill(currents, currents + input_size, 0);
    auto addSReacCurrents = [&](const auto& surfReacs) {
        for (size_t i = 0; i < input_size; ++i) {
            auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
            if (localInd.valid()) {
                auto& info = mesh.getTri(localInd);
                auto patch = info.patchPtr;
                if (patch != nullptr and info.cont_id.valid()) {
                    auto& patchdef = statedef->getPatchdef(model::patch_id(patch->getID()));
                    currents[i] += surfReacs.getCurrent(patchdef.getIdx(), info.cont_id);
                }
            }
        }
    };
    auto addOhmicCurrents = [&](const auto& getCurrs) {
        for (size_t i = 0; i < input_size; ++i) {
            auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
            if (localInd.valid()) {
                auto* patch = mesh.getTriPatch(localInd);
                if (patch != nullptr) {
                    container::patch_id patch_id(patch->getMeshID().get());
                    const auto& patchdef = statedef->getPatchdef(patch_id);
                    for (const auto& currPtr: getCurrs(patchdef)) {
                        const auto& face_bf2verts =
                            osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE), localInd.get());
                        for (const auto& vert_id: face_bf2verts) {
                            currents[i] += currPtr->getTriCurrentOnVertex(
                                input->potential_on_vertices_w[vert_id],
                                localInd,
                                input->pools,
                                mesh,
                                state_time);
                        }
                    }
                }
            }
        }
    };
    addSReacCurrents(data->kproc_state.surfaceReactions());
    addSReacCurrents(data->kproc_state.vDepSurfaceReactions());
    addSReacCurrents(data->kproc_state.complexSurfaceReactions());
    addSReacCurrents(data->kproc_state.vDepComplexSurfaceReactions());
    addSReacCurrents(data->kproc_state.ghkSurfaceReactions());
    addSReacCurrents(data->kproc_state.complexGhkSurfaceReactions());
    addOhmicCurrents(
        [](const Patchdef& pd) -> const auto& { return pd.template currents<OhmicCurrdef>(); });
    addOhmicCurrents([](const Patchdef& pd) -> const auto& {
        return pd.template currents<ComplexOhmicCurrdef>();
    });

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, currents, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriOhmicErevsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::ohmic_current_id& ohmic_current,
    double* rv,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchTriOhmicErevsNP<OhmicCurrdef>({indices, input_size},
                                          ohmic_current,
                                          {rv, output_size},
                                          local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriComplexOhmicErevsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::complex_ohmic_current_id& ohmic_current,
    double* rv,
    size_t output_size,
    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    getBatchTriOhmicErevsNP<ComplexOhmicCurrdef>({indices, input_size},
                                                 ohmic_current,
                                                 {rv, output_size},
                                                 local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriOhmicErev(
    osh::GO triangle,
    const model::ohmic_current_id& ohmic_current,
    double reversal_potential,
    bool local) {
    setTriOhmicErev<OhmicCurrdef>(triangle, ohmic_current, reversal_potential, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriComplexOhmicErev(
    osh::GO triangle,
    const model::complex_ohmic_current_id& ohmic_current,
    double reversal_potential,
    bool local) {
    setTriOhmicErev<ComplexOhmicCurrdef>(triangle, ohmic_current, reversal_potential, local);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriVClamped(osh::GO vertex,
                                                                     bool local) const {
    auto tri = getLocalInd<mesh::triangle_local_id_t>(vertex, local, true);
    bool clamped = false;
    if (tri.valid()) {
        if (data->efield) {
            clamped = true;
            const auto& tris2verts = mesh.ask_verts_of(Omega_h::FACE);
            const auto verts = osh::gather_verts<3>(tris2verts, tri.get());
            for (auto v: verts) {
                clamped &= data->efield->getVertVClamped(mesh::vertex_local_id_t(v));
            }
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
    return allReduce(clamped, MPI_LOR);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriVClamped(osh::GO vertex,
                                                                     bool clamped,
                                                                     bool local) {
    auto tri = getLocalInd<mesh::triangle_local_id_t>(vertex, local, true);
    if (tri.valid()) {
        if (data->efield) {
            const auto& tris2verts = mesh.ask_verts_of(Omega_h::FACE);
            const auto verts = osh::gather_verts<3>(tris2verts, tri.get());
            for (auto v: verts) {
                data->efield->setVertVClamped(mesh::vertex_local_id_t(v), clamped);
            }
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriIClamp(osh::GO tri,
                                                                        bool local) const {
    osh::Real local_val(0.0);
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, true);
    if (localInd.valid()) {
        local_val = input->current_on_triangles_w[localInd.get()];
    }
    if (local) {
        return local_val;
    } else {
        return allReduce(local_val);
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriIClamp(osh::GO tri,
                                                                   osh::Real current,
                                                                   bool local) {
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, false);
    if (localInd.valid()) {
        auto& info = mesh.getTri(localInd);
        auto patch = info.patchPtr;
        if (patch == nullptr) {
            throw std::invalid_argument("Triangle " + std::to_string(localInd) +
                                        " is not part of a patch, cannot set IClamp on it.");
        }
        input->current_on_triangles_w[localInd.get()] = current;
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
MembraneResistivity OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriRes(osh::GO tri,
                                                                               bool local) const {
    osh::Real local_res(0.0);
    osh::Real local_erev(0.0);

    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, true);
    if (localInd.valid()) {
        local_res = 1.0 / input->conductivity_on_triangles_w[localInd.get()];
        local_erev = input->reversal_potential_on_triangles_w[localInd.get()];
    }
    if (local) {
        return {local_res, local_erev};
    } else {
        osh::Real res(0.0);
        osh::Real erev(0.0);
        auto err = MPI_Allreduce(&local_res, &res, 1, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        err = MPI_Allreduce(&local_erev, &erev, 1, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        return {res, erev};
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriRes(const osh::GO tri,
                                                                osh::Real res,
                                                                osh::Real erev,
                                                                bool local) {
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, false);
    if (localInd.valid()) {
        input->conductivity_on_triangles_w[localInd.get()] = 1.0 / res;
        input->reversal_potential_on_triangles_w[localInd.get()] = erev;
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriCapac(osh::GO tri,
                                                                       bool local) const {
    osh::Real local_val(0.0);

    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, true);
    if (localInd.valid()) {
        local_val = input->capacitance_on_triangles_w[localInd.get()];
    }
    if (local) {
        return local_val;
    } else {
        osh::Real res(0.0);
        auto err = MPI_Allreduce(&local_val, &res, 1, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        return res;
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriCapac(const osh::GO tri,
                                                                  osh::Real c,
                                                                  bool local) {
    auto localInd = getLocalInd<mesh::triangle_local_id_t>(tri, local, false);
    if (localInd.valid()) {
        input->capacitance_on_triangles_w[localInd.get()] = c;
    }
}

//-----------------------------------------------

#endif  // USE_PETSC


//////////////////////
// Location: Vertex //
//////////////////////

#if USE_PETSC
// E-field value

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchVertVsNP(const osh::GO* indices,
                                                                       size_t input_size,
                                                                       osh::Real* voltages,
                                                                       size_t output_size,
                                                                       bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    if (not local) {
        std::fill(voltages, voltages + input_size, 0);
    }
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::vertex_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            voltages[i] = input->potential_on_vertices_w[localInd.get()];
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, voltages, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchVertVsNP(const osh::GO* indices,
                                                                       size_t input_size,
                                                                       osh::Real* voltages,
                                                                       size_t output_size,
                                                                       bool local) {
    assert(input_size == output_size);
    (void) output_size;
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::vertex_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            input->potential_on_vertices_w[localInd.get()] = voltages[i];
        }
    }
    const auto& syncedv = mesh.sync_array(osh::VERT, osh::Reals(input->potential_on_vertices_w), 1);
    std::copy(syncedv.begin(), syncedv.end(), input->potential_on_vertices_w.begin());
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
bool OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getVertVClamped(osh::GO vertex,
                                                                      bool local) const {
    auto vert = getLocalInd<mesh::vertex_local_id_t>(vertex, local);
    bool clamped = false;
    if (vert.valid()) {
        if (data->efield) {
            clamped = data->efield->getVertVClamped(vert);
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
    return allReduce(clamped, MPI_LOR);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setVertVClamped(osh::GO vertex,
                                                                      bool clamped,
                                                                      bool local) {
    auto vert = getLocalInd<mesh::vertex_local_id_t>(vertex, local);
    if (vert.valid()) {
        if (data->efield) {
            data->efield->setVertVClamped(vert, clamped);
        } else {
            throw std::logic_error("Efield is not enabled, cannot clamp vertex potential.");
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getVertIClamp(const osh::GO vertex,
                                                                         bool local) const {
    osh::Real local_val(0.0);

    auto localInd = getLocalInd<mesh::vertex_local_id_t>(vertex, local, true);
    if (localInd.valid()) {
        local_val = input->current_on_vertices_w[localInd.get()];
    }
    if (local) {
        return local_val;
    } else {
        osh::Real res(0.0);
        auto err = MPI_Allreduce(&local_val, &res, 1, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        return res;
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setVertIClamp(const osh::GO vertex,
                                                                    const osh::Real current,
                                                                    bool local) {
    auto localInd = getLocalInd<mesh::vertex_local_id_t>(vertex, local, false);
    if (localInd.valid()) {
        input->current_on_vertices_w[localInd.get()] = current;
    }
}

//-----------------------------------------------

#endif  // USE_PETSC

///////////////////////////////////////////////////////////////////////////////////////////////
// Internal methods
///////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////
// General methods //
/////////////////////

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::compute_num_species_per_elements(
    DistMesh& t_mesh,
    const Statedef& statedef,
    osh::LOs& num_species_per_owned_elems,
    osh::LOs& num_species_per_elems,
    std::optional<osh::LOs>& num_species_per_bounds) {
    const auto& owned_elems_mask = t_mesh.owned_elems_mask();

    {
        osh::Write<osh::LO> num_species_per_owned_elems_w(owned_elems_mask.size(), 0);
        osh::Write<osh::LO> num_species_per_elems_w(owned_elems_mask.size(), 0);
        for (const auto& compartment: statedef.compdefs()) {
            const auto num_species = compartment->getNSpecs();
            for (auto elem: t_mesh.getEntities(compartment->getID())) {
                if (owned_elems_mask[elem.get()] != 0) {
                    num_species_per_owned_elems_w[elem.get()] = num_species;
                }
                num_species_per_elems_w[elem.get()] = num_species;
            }
        }
        num_species_per_owned_elems = num_species_per_owned_elems_w;
        num_species_per_elems = num_species_per_elems_w;
    }

    if (!statedef.patchdefs().container().empty()) {
        // initialize a vector to record the number of species owned by a patch
        // element and owned by the process
        osh::Write<osh::LO> num_species_per_bounds_w(t_mesh.owned_bounds_mask().size(), 0);
        for (const auto& patch: statedef.patchdefs()) {
            for (const auto boundary: t_mesh.getOwnedEntities(patch->getID())) {
                num_species_per_bounds_w[boundary.get()] = patch->getNSpecs();
            }
        }
        num_species_per_bounds = num_species_per_bounds_w;
    } else {
        num_species_per_bounds = std::nullopt;
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::init(std::unique_ptr<Statedef>&& t_statedef) {
    this->statedef.swap(t_statedef);
    assert(statedef != nullptr);
    this->mesh.init();

    // Initialize diffusion boundaries
    for (auto& db: mesh.diffusionBoundaries()) {
        Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
        Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
        db.comp1_spec2dcst.resize(comp1.getNSpecs(), 0.0);
        db.comp2_spec2dcst.resize(comp2.getNSpecs(), 0.0);
        db.conv_12.resize(comp1.getNSpecs());
        for (auto sp1: container::species_id::range(comp1.getNSpecs())) {
            db.conv_12[sp1.get()] = comp2.getSpecContainerIdx(comp1.getSpecModelIdx(sp1));
        }
        db.conv_21.resize(comp2.getNSpecs());
        for (auto sp2: container::species_id::range(comp2.getNSpecs())) {
            db.conv_21[sp2.get()] = comp1.getSpecContainerIdx(comp2.getSpecModelIdx(sp2));
        }
    }

    osh::LOs num_species_per_owned_elems;
    osh::LOs num_species_per_elems;
    std::optional<osh::LOs> num_species_per_bounds;
    compute_num_species_per_elements(mesh,
                                     *statedef,
                                     num_species_per_owned_elems,
                                     num_species_per_elems,
                                     num_species_per_bounds);
    this->input = std::make_unique<SimulationInput>(num_species_per_owned_elems,
                                                    num_species_per_bounds,
                                                    num_species_per_elems,
                                                    statedef->substates_per_complexes(),
                                                    0 /*Unused in context*/,
                                                    this->rng,
                                                    mesh.owned_verts_mask().size(),
                                                    mesh,
                                                    *statedef);

    // Initialize triangle capacitance on membranes
    for (auto& memb: statedef->membranes()) {
        auto capac = memb->capacitance();
        for (const auto& patch: memb->getPatchesIds()) {
            for (const auto tri: mesh.getOwnedEntities(patch)) {
                input->capacitance_on_triangles_w[tri.get()] = capac;
            }
        }
    }

    data = std::make_unique<SimulationData<SSA, SearchMethod, DiffMethod>>(
        mesh, *this->statedef, *input, this->rng, indepKProcs);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::evolve_rd(const osh::Real rd_dt) {
    util::TimeTracker t;
    t.start();

    data->pools.reset_occupancy_rd(state_time);
    data->ssaOp.run(rd_dt, state_time);
    t.stop();
    this->reactions_timer += t.diff();

    if (data->diffOp.has_active_diffusions()) {
        t.start();
        data->diffOp(rd_dt, state_time);
        t.stop();
        this->diffusions_timer += t.diff();
    }
    state_time += rd_dt;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::run_rd(const osh::Real end_time) {
    data->ssaOp.resetAndUpdateAll(state_time, end_time);

    osh::Real rd_dt;
    if constexpr (DiffMethod == DiffusionMethod::ConstantDiffDt) {
        rd_dt = data->diffOp.getDt();
        // number of standard steps (this can be negative)
        const int n_steps_std = std::floor((end_time - state_time) / rd_dt) - 1;
        // std steps loop -1. n_steps_std can be negative so int is the correct type
        for (int i_rd = 0; i_rd < n_steps_std; ++i_rd) {
            evolve_rd(std::min(rd_dt, end_time - state_time));
        }

    } else if constexpr (DiffMethod == DiffusionMethod::TauLeapingDiffDt) {
        rd_dt = allReduce(data->diffOp.getDt(), MPI_MIN);
        bool recompute = rd_dt > data->diffOp.getDefaultDt();
        unsigned int skip = 0;
        while (state_time + rd_dt < end_time) {
            evolve_rd(rd_dt);
            if (recompute) {
                rd_dt = allReduce(data->diffOp.getDt(), MPI_MIN);
                recompute = rd_dt > data->diffOp.getDefaultDt();
            } else {
                // If the rd_dt is not bigger than the default value, we run several steps of normal
                // diffusion
                if (++skip > data->diffOp.getMaxDtSkips()) {
                    skip = 0;
                    recompute = true;
                }
            }
        }

    } else {
        static_assert(steps::util::always_false_v<decltype(DiffMethod)>,
                      "Unknown diffusion method");
    }

    // sync time steps
    // we compare always with end_time because comparing dts can fail due to numerical error
    const auto next_std_state_time = state_time + rd_dt;
    if (!steps::util::almost_equal(next_std_state_time, end_time) &&
        next_std_state_time < end_time) {
        evolve_rd(rd_dt);
        assert(end_time > state_time);
        assert(!steps::util::almost_equal(end_time, state_time));
        evolve_rd(end_time - state_time);
    } else if (!steps::util::almost_equal(end_time, state_time)) {
        evolve_rd(end_time - state_time);
    }

    assert(steps::util::almost_equal(end_time, state_time));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::evolve(const osh::Real ef_dt) {
    data->pools.reset_occupancy_ef(state_time);

    ++this->num_iterations;
    data->kproc_state.resetCurrents();
    osh::Reals potential_on_vertices(input->potential_on_vertices_w);
    data->kproc_state.updateVDepSReacs(potential_on_vertices);

    run_rd(state_time + ef_dt);

    // We divide the charge_flows in currents_ by ef_dt so we really get the currents
    data->kproc_state.finalizeCurrents(ef_dt);
#if USE_PETSC
    if (data->efield) {
        util::TimeTracker t;
        t.start();
        data->efield->evolve(input->potential_on_vertices_w,
                             input->current_on_vertices_w,
                             input->current_on_triangles_w,
                             input->capacitance_on_triangles_w,
                             input->conductivity_on_triangles_w,
                             input->reversal_potential_on_triangles_w,
                             input->pools,
                             data->kproc_state,
                             state_time,
                             ef_dt);
        t.stop();
        this->efield_timer += t.diff();
    }
#endif  // USE_PETSC
}

//-----------------------------------------------

///////////////////////////
// Location: Compartment //
///////////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& spec_id) const {
    const auto species = statedef->getCompSpecContainerIdx(compartment, spec_id);
    const auto lambda = [=](osh::GO accu, mesh::tetrahedron_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, species));
    };
    const auto& elements = mesh.getOwnedEntities(compartment);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setOwnedCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& species,
    osh::Real num_molecules,
    const math::DistributionMethod distribution) {
    const auto& [elems, volumes, rank_volume] = mesh.measure(compartment);
    const container::species_id species_id = statedef->getCompSpecContainerIdx(compartment,
                                                                               species);

    osh::Write<osh::GO> mols_on_elements;

    auto dist = math::make_dist(static_cast<osh::I64>(num_molecules), volumes);
    mols_on_elements = dist.distribute(this->rng, distribution);

    for (auto k = 0; k < mols_on_elements.size(); ++k) {
        if (mols_on_elements[k] >= static_cast<osh::GO>(INT32_MAX)) {
            std::ostringstream oss;
            oss << "Unsupported number of molecules per tetrahedron: " << std::setprecision(20)
                << mols_on_elements[k] << " but maximum value is "
                << std::numeric_limits<osh::LO>::max() << " (max 32 bits integral value)";
            ArgErrLog(oss.str());
        }
        data->pools.assign(mesh::tetrahedron_id_t(elems[k]),
                           species_id,
                           static_cast<osh::LO>(mols_on_elements[k]));
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    const auto spec_count = getOwnedCompSpecCount(compartment, species);
    return spec_count / (1.0e3 * mesh.getMeasure().rank_measure() * math::AVOGADRO);
}

//-----------------------------------------------

// Complexes

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedCompComplexCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const std::vector<util::strongid_vector<model::complex_substate_id,
                                            steps::model::SubunitStateFilter>>& f) const {
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);
    const auto& filt = data->pools.moleculesOnElements().updatedComplexFilter(cplxIdx, f);

    const auto lambda = [=](osh::GO accu, mesh::tetrahedron_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, cplxIdx, filt));
    };
    const auto& elements = mesh.getOwnedEntities(compartment);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setOwnedCompComplexCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const util::strongid_vector<model::complex_substate_id, uint>& i,
    osh::Real num_molecules,
    const math::DistributionMethod distribution) {
    const auto& [elems, volumes, rank_volume] = mesh.measure(compartment);
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);

    osh::Write<osh::GO> mols_on_elements;

    auto dist = math::make_dist(static_cast<osh::I64>(num_molecules), volumes);
    mols_on_elements = dist.distribute(this->rng, distribution);

    for (auto k = 0; k < mols_on_elements.size(); ++k) {
        if (mols_on_elements[k] >= static_cast<osh::GO>(INT32_MAX)) {
            std::ostringstream oss;
            oss << "Unsupported number of molecules per tetrahedron: " << std::setprecision(20)
                << mols_on_elements[k] << " but maximum value is "
                << std::numeric_limits<osh::LO>::max() << " (max 32 bits integral value)";
            ArgErrLog(oss.str());
        }
        data->pools.assign(mesh::tetrahedron_id_t(elems[k]),
                           cplxIdx,
                           i,
                           static_cast<osh::LO>(mols_on_elements[k]));
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedCompComplexSUSCount(
    const model::compartment_id& compartment,
    const model::complex_name& complex,
    const std::vector<
        util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>& f,
    model::complex_substate_id m) const {
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);
    const auto& filt = data->pools.moleculesOnElements().updatedComplexFilter(cplxIdx, f);

    const auto lambda = [=](osh::GO accu, mesh::tetrahedron_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, cplxIdx, filt, m));
    };
    const auto& elements = mesh.getOwnedEntities(compartment);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

//-----------------------------------------------

/////////////////////
// Location: Patch //
/////////////////////

// Species

//-----------------------------------------------

//-----------------------------------------------

// Complexes

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedPatchComplexCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const std::vector<util::strongid_vector<model::complex_substate_id,
                                            steps::model::SubunitStateFilter>>& f) const {
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);
    const auto& filt = data->pools.moleculesOnPatchBoundaries().updatedComplexFilter(cplxIdx, f);

    const auto lambda = [=](osh::GO accu, mesh::triangle_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, cplxIdx, filt));
    };
    const auto& elements = mesh.getOwnedEntities(patch);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setOwnedPatchComplexCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const util::strongid_vector<model::complex_substate_id, uint>& i,
    osh::Real num_molecules,
    const math::DistributionMethod distribution) {
    const auto& [elems, areas, rank_area] = mesh.measure(patch);
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);

    osh::Write<osh::GO> mols_on_elements;

    auto dist = math::make_dist(static_cast<osh::I64>(num_molecules), areas);
    mols_on_elements = dist.distribute(this->rng, distribution);

    for (auto k = 0; k < mols_on_elements.size(); ++k) {
        if (mols_on_elements[k] >= static_cast<osh::GO>(INT32_MAX)) {
            std::ostringstream oss;
            oss << "Unsupported number of molecules per triangle: " << std::setprecision(20)
                << mols_on_elements[k] << " but maximum value is "
                << std::numeric_limits<osh::LO>::max() << " (max 32 bits integral value)";
            ArgErrLog(oss.str());
        }
        data->pools.assign(mesh::triangle_id_t(elems[k]),
                           cplxIdx,
                           i,
                           static_cast<osh::LO>(mols_on_elements[k]));
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getOwnedPatchComplexSUSCount(
    const model::patch_id& patch,
    const model::complex_name& complex,
    const std::vector<
        util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>& f,
    model::complex_substate_id m) const {
    const model::complex_id cplxIdx = statedef->getComplexModelIdx(complex);
    const auto& filt = data->pools.moleculesOnPatchBoundaries().updatedComplexFilter(cplxIdx, f);

    const auto lambda = [=](osh::GO accu, mesh::triangle_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, cplxIdx, filt, m));
    };
    const auto& elements = mesh.getOwnedEntities(patch);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

//-----------------------------------------------

///////////////////////////
// Location: Tetrahedron //
///////////////////////////

// Species
template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setOwnedElementSpecCount(
    const model::compartment_id& compartment,
    mesh::tetrahedron_id_t element,
    const model::species_name& species,
    osh::Real num_molecules) {
    const auto species_id = statedef->getCompSpecContainerIdx(compartment, species);
    if (num_molecules >= static_cast<osh::Real>(INT32_MAX)) {
        std::ostringstream oss;
        oss << "Unsupported number of molecules per tetrahedron: " << std::setprecision(20)
            << num_molecules << " but maximum value is " << std::numeric_limits<osh::LO>::max()
            << " (max 32 bits integral value)";
        ArgErrLog(oss.str());
    }

    data->pools.assign(element, species_id, static_cast<osh::LO>(num_molecules));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchElemValsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& species,
    osh::Real* vals,
    bool useConc,
    bool local) const {
    if (not local) {
        std::fill(vals, vals + input_size, 0);
    }
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto compartment_id = mesh.getCompartment(localInd);
            const auto spec_id = statedef->getCompSpecContainerIdx(compartment_id, species);
            vals[i] = data->pools(localInd, spec_id);
            if (useConc) {
                vals[i] /= mesh.getTet(localInd).vol * 1.0e3 * math::AVOGADRO;
            }
        }
    }
    if (not local) {
        auto err = MPI_Allreduce(MPI_IN_PLACE, vals, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchElemValsNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& species,
    osh::Real* vals,
    bool useConc,
    bool local) const {
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto compartment_id = mesh.getCompartment(localInd);
            const auto spec_id = statedef->getCompSpecContainerIdx(compartment_id, species);
            osh::LO nb;
            if (useConc) {
                auto v = vals[i] * mesh.getTet(localInd).vol * 1e3 * math::AVOGADRO;
                auto v_trunc = static_cast<osh::LO>(v);
                if (v_trunc < v) {
                    std::uniform_real_distribution<double> uniform(0.0, 1.0);
                    if (uniform(this->rng) < v - v_trunc) {
                        ++v_trunc;
                    }
                }
                nb = v_trunc;
            } else {
                nb = static_cast<osh::LO>(vals[i]);
            }

            data->pools.assign(mesh::tetrahedron_id_t(localInd), spec_id, nb);
        }
    }
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename Reacs, typename MReacID>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTetReacK(const Reacs& reacs,
                                                                       osh::GO tet,
                                                                       const MReacID reac,
                                                                       bool local) const {
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    osh::Real kcst = 0.0;
    if (localInd.valid()) {
        auto& info = mesh.getTet(localInd);
        auto comp = info.compPtr;
        if (comp != nullptr) {
            auto& compdef = statedef->getCompdef(model::compartment_id(comp->getID()));
            kcst = reacs.getKcst(compdef.getIdx(), info.cont_id, compdef.getReacIdx(reac));
        } else {
            throw std::invalid_argument("Tetrahedron " + std::to_string(localInd.get()) +
                                        " is not associated to a compartment.");
        }
    }
    return allReduce(kcst);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename Reacs, typename MReacID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTetReacK(Reacs& reacs,
                                                                  osh::GO tet,
                                                                  const MReacID reac,
                                                                  osh::Real kcst,
                                                                  bool local) {
    auto localInd = getLocalInd<mesh::tetrahedron_local_id_t>(tet, local, true);
    if (localInd.valid()) {
        auto& info = mesh.getTet(localInd);
        auto comp = info.compPtr;
        if (comp != nullptr) {
            auto& compdef = statedef->getCompdef(model::compartment_id(comp->getID()));
            reacs.setKcst(compdef.getIdx(), info.cont_id, compdef.getReacIdx(reac), kcst);
        } else {
            throw std::invalid_argument("Tetrahedron " + std::to_string(localInd.get()) +
                                        " is not associated to a compartment.");
        }
    }
}

//-----------------------------------------------

////////////////////////
// Location: Triangle //
////////////////////////

// Species

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchBoundSpecCountNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& species,
    osh::Real* counts,
    bool local) const {
    const auto spec_model_idx = statedef->getSpecModelIdx(species);
    const auto& molecules = data->pools.moleculesOnPatchBoundaries();

    if (not local) {
        std::fill(counts, counts + input_size, 0);
    }
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto patch_id = model::patch_id(
                mesh.getTriPatch(mesh::triangle_id_t(localInd))->getID());
            auto spec_id = statedef->getPatchdef(patch_id).getSpecContainerIdx(spec_model_idx);
            counts[i] = molecules(localInd, spec_id);
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, counts, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setBatchBoundSpecCountNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& species,
    osh::Real* counts,
    bool local) const {
    const auto spec_model_idx = statedef->getSpecModelIdx(species);

    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            const auto patch_id = model::patch_id(
                mesh.getTriPatch(mesh::triangle_id_t(localInd))->getID());
            auto spec_id = statedef->getPatchdef(patch_id).getSpecContainerIdx(spec_model_idx);
            data->pools.assign(localInd, spec_id, counts[i]);
        }
    }
}

//-----------------------------------------------

// Reactions

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename Reacs, typename MReacID>
osh::Real OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getTriSReacK(const Reacs& reacs,
                                                                        osh::GO triangle,
                                                                        const MReacID& reactionId,
                                                                        bool local) const {
    osh::Real kcst = 0.0;
    auto tri = getLocalInd<mesh::triangle_local_id_t>(triangle, local, true);
    if (tri.valid()) {
        auto& info = mesh.getTri(tri);
        auto patch = info.patchPtr;
        if (patch != nullptr) {
            auto& patchdef = statedef->getPatchdef(model::patch_id(patch->getID()));
            kcst = reacs.get_Kcst(patchdef.getIdx(), info.cont_id, patchdef.getReacIdx(reactionId));
        } else {
            throw std::invalid_argument("Triangle " + std::to_string(triangle) +
                                        " is not associated to a patch.");
        }
    }
    return allReduce(kcst);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename Reacs, typename MReacID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriSReacK(Reacs& reacs,
                                                                   osh::GO triangle,
                                                                   const MReacID& reactionId,
                                                                   osh::Real kCst,
                                                                   bool local) {
    auto tri = getLocalInd<mesh::triangle_local_id_t>(triangle, local, true);
    if (tri.valid()) {
        auto& info = mesh.getTri(tri);
        auto patch = info.patchPtr;
        if (patch != nullptr) {
            auto& patchdef = statedef->getPatchdef(model::patch_id(patch->getID()));
            reacs.set_Kcst(patchdef.getIdx(), info.cont_id, patchdef.getReacIdx(reactionId), kCst);
        } else {
            throw std::invalid_argument("Triangle " + std::to_string(triangle) +
                                        " is not associated to a patch.");
        }
    }
}

//-----------------------------------------------

#if USE_PETSC
// E-field value

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename SReacT, typename MReacID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriSReacIsNP(const SReacT& reacs,
                                                                           const osh::GO* indices,
                                                                           size_t input_size,
                                                                           const MReacID reac,
                                                                           osh::Real* currents,
                                                                           bool local) const {
    std::fill(currents, currents + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            auto& info = mesh.getTri(localInd);
            auto patch = info.patchPtr;
            if (patch != nullptr and info.cont_id.valid()) {
                auto& patchdef = statedef->getPatchdef(model::patch_id(patch->getID()));
                currents[i] +=
                    reacs.getCurrent(patchdef.getIdx(), info.cont_id, patchdef.getReacIdx(reac));
            }
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, currents, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename CurrdefT, typename MCurrID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriOhmicIsNP(const osh::GO* indices,
                                                                           size_t input_size,
                                                                           const MCurrID curr,
                                                                           osh::Real* currents,
                                                                           bool local) const {
    std::fill(currents, currents + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            auto* patch = mesh.getTriPatch(localInd);
            if (patch != nullptr) {
                container::patch_id patch_id(patch->getMeshID().get());
                const auto& patchdef = statedef->getPatchdef(patch_id);
                const auto& currId = patchdef.getCurrIdx(curr);
                const auto& currPtr = patchdef.template currents<CurrdefT>().at(currId.get());
                const auto& face_bf2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                                 localInd.get());
                for (const auto& vert_id: face_bf2verts) {
                    currents[i] +=
                        currPtr->getTriCurrentOnVertex(input->potential_on_vertices_w[vert_id],
                                                       localInd,
                                                       input->pools,
                                                       mesh,
                                                       state_time);
                }
            }
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, currents, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename SReacT, typename MReacID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriGHKIsNP(const SReacT& surfReacs,
                                                                         const osh::GO* indices,
                                                                         size_t input_size,
                                                                         const MReacID curr,
                                                                         osh::Real* currents,
                                                                         bool local) const {
    std::fill(currents, currents + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        auto localInd = getLocalInd<mesh::triangle_local_id_t>(indices[i], local, true);
        if (localInd.valid()) {
            auto& info = mesh.getTri(localInd);
            auto patch = info.patchPtr;
            if (patch != nullptr and info.cont_id.valid()) {
                auto& patchdef = statedef->getPatchdef(model::patch_id(patch->getID()));
                currents[i] += surfReacs.getCurrent(patchdef.getIdx(),
                                                    info.cont_id,
                                                    patchdef.getReacIdx(curr));
            }
        }
    }

    if (not local) {
        auto err =
            MPI_Allreduce(MPI_IN_PLACE, currents, input_size, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename CurrdefT, typename MCurrID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::getBatchTriOhmicErevsNP(
    const gsl::span<const osh::GO>& triangles,
    const MCurrID& ohmic_current,
    const gsl::span<double>& erev,
    bool local) const {
    for (size_t i = 0; i < triangles.size(); ++i) {
        auto triangle = getLocalInd<mesh::triangle_local_id_t>(triangles[i], local, true);
        if (triangle.valid()) {
            auto* patch = mesh.getTriPatch(triangle);
            if (patch != nullptr) {
                container::patch_id patch_id(patch->getMeshID().get());
                const auto& patchdef = statedef->getPatchdef(patch_id);
                const auto& currId = patchdef.getCurrIdx(ohmic_current);
                const auto& currPtr = patchdef.template currents<CurrdefT>().at(currId.get());
                erev[i] = currPtr->getReversalPotential(triangle);
            }
        }
    }

    if (local) {
        return;
    } else {
        auto err = MPI_Allreduce(
            MPI_IN_PLACE, erev.data(), erev.size(), MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
template <typename CurrdefT, typename MCurrID>
void OmegaHSimulation<SSA, SearchMethod, DiffMethod>::setTriOhmicErev(osh::GO triangle,
                                                                      const MCurrID& ohmic_current,
                                                                      double reversal_potential,
                                                                      bool local) {
    auto local_index = getLocalInd<mesh::triangle_local_id_t>(triangle, local, true);
    if (local_index.valid()) {
        auto* patch = mesh.getTriPatch(local_index);
        if (patch != nullptr) {
            container::patch_id patch_id(patch->getMeshID().get());
            const auto& patchdef = statedef->getPatchdef(patch_id);
            const auto& currId = patchdef.getCurrIdx(ohmic_current);
            const auto& currPtr = patchdef.template currents<CurrdefT>().at(currId.get());
            currPtr->setReversalPotential(local_index, reversal_potential);
        } else {
            throw std::invalid_argument("Triangle " + std::to_string(local_index.get()) +
                                        " not in a patch.");
        }
    }
}

//-----------------------------------------------

#endif  // USE_PETSC

/////////////////////
// Utility methods //
/////////////////////

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
util::strongid_vector<model::complex_substate_id, uint>
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::_convertComplexState(
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) {
    util::strongid_vector<model::complex_substate_id, uint> state;
    AssertLog(f.size() == 1);
    state.container().reserve(f[0].size());
    for (const auto& susfilt: f[0]) {
        AssertLog(susfilt.min == susfilt.max);
        state.container().push_back(susfilt.min);
    }
    return state;
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod, DiffusionMethod DiffMethod>
std::vector<util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
OmegaHSimulation<SSA, SearchMethod, DiffMethod>::_convertComplexFilters(
    const std::vector<std::vector<steps::model::SubunitStateFilter>>& f) {
    std::vector<util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
        filts;
    filts.reserve(f.size());
    for (auto& filt: f) {
        filts.emplace_back(filt);
    }
    return filts;
}

//-----------------------------------------------

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<Simulation> GetSimulation(steps::model::Model& model,
                                          DistMesh& mesh,
                                          const rng::RNGptr& r,
                                          int _ssaMethod,
                                          int _searchMethod,
                                          int _diffMethod,
                                          bool indepKProcs,
                                          bool isEfield) {
    auto ssaMethod = static_cast<SSAMethod>(_ssaMethod);
    auto searchMethod = static_cast<NextEventSearchMethod>(_searchMethod);
    auto diffMethod = static_cast<DiffusionMethod>(_diffMethod);

    switch (ssaMethod) {
    case SSAMethod::SSA:
        return _getSimulation<SSAMethod::SSA>(
            model, mesh, r, searchMethod, diffMethod, indepKProcs, isEfield);
    case SSAMethod::RSSA:
        return _getSimulation<SSAMethod::RSSA>(
            model, mesh, r, searchMethod, diffMethod, indepKProcs, isEfield);
    case SSAMethod::RLeaping:
        if (searchMethod == NextEventSearchMethod::RLeaping) {
            return _getSimulation<SSAMethod::RLeaping>(
                model, mesh, r, searchMethod, diffMethod, indepKProcs, isEfield);
        } else {
            throw std::invalid_argument(
                "The RLeaping SSA method can only be used with the RLeaping next event "
                "search method");
        }
    default:
        throw std::invalid_argument("Invalid SSA method");
    }
}

template <SSAMethod SSA>
std::unique_ptr<Simulation> _getSimulation(steps::model::Model& model,
                                           DistMesh& mesh,
                                           const rng::RNGptr& r,
                                           NextEventSearchMethod searchMethod,
                                           DiffusionMethod diffMethod,
                                           bool indepKProcs,
                                           bool isEfield) {
    switch (searchMethod) {
    case NextEventSearchMethod::Direct:
        return _getSimulation<SSA, NextEventSearchMethod::Direct>(
            model, mesh, r, diffMethod, indepKProcs, isEfield);
    case NextEventSearchMethod::GibsonBruck:
        if constexpr (SSA == SSAMethod::RSSA) {
            throw std::invalid_argument(
                "Cannot use GibsonBruck next event search method with RSSA method");
        } else {
            return _getSimulation<SSA, NextEventSearchMethod::GibsonBruck>(
                model, mesh, r, diffMethod, indepKProcs, isEfield);
        }
    case NextEventSearchMethod::RLeaping:
        if constexpr (SSA == SSAMethod::RLeaping) {
            return _getSimulation<SSA, NextEventSearchMethod::RLeaping>(
                model, mesh, r, diffMethod, indepKProcs, isEfield);
        } else {
            throw std::invalid_argument(
                "The RLeaping SSA method can only be used in conjunction with the RLeaping "
                "next event search method.");
        }
    default:
        throw std::invalid_argument("Invalid next event search method");
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::unique_ptr<Simulation> _getSimulation(steps::model::Model& model,
                                           DistMesh& mesh,
                                           const rng::RNGptr& r,
                                           DiffusionMethod diffMethod,
                                           bool indepKProcs,
                                           bool isEfield) {
    switch (diffMethod) {
    case DiffusionMethod::ConstantDiffDt:
        return std::make_unique<
            OmegaHSimulation<SSA, SearchMethod, DiffusionMethod::ConstantDiffDt>>(
            model, mesh, r, indepKProcs, isEfield);
    case DiffusionMethod::TauLeapingDiffDt:
        return std::make_unique<
            OmegaHSimulation<SSA, SearchMethod, DiffusionMethod::TauLeapingDiffDt>>(
            model, mesh, r, indepKProcs, isEfield);
    default:
        throw std::invalid_argument("Invalid diffusion method");
    }
}

// explicit template instantiation definitions

template class OmegaHSimulation<SSAMethod::SSA,
                                NextEventSearchMethod::GibsonBruck,
                                DiffusionMethod::ConstantDiffDt>;
template class OmegaHSimulation<SSAMethod::SSA,
                                NextEventSearchMethod::Direct,
                                DiffusionMethod::ConstantDiffDt>;
template class OmegaHSimulation<SSAMethod::RSSA,
                                NextEventSearchMethod::Direct,
                                DiffusionMethod::ConstantDiffDt>;
template class OmegaHSimulation<SSAMethod::RLeaping,
                                NextEventSearchMethod::RLeaping,
                                DiffusionMethod::ConstantDiffDt>;
template class OmegaHSimulation<SSAMethod::SSA,
                                NextEventSearchMethod::GibsonBruck,
                                DiffusionMethod::TauLeapingDiffDt>;
template class OmegaHSimulation<SSAMethod::SSA,
                                NextEventSearchMethod::Direct,
                                DiffusionMethod::TauLeapingDiffDt>;
template class OmegaHSimulation<SSAMethod::RSSA,
                                NextEventSearchMethod::Direct,
                                DiffusionMethod::TauLeapingDiffDt>;
template class OmegaHSimulation<SSAMethod::RLeaping,
                                NextEventSearchMethod::RLeaping,
                                DiffusionMethod::TauLeapingDiffDt>;

}  // namespace steps::dist

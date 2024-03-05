#include "tetopsplit.hpp"

#include <iostream>
#ifdef USE_PETSC
#include <petscksp.h>
#endif  // USE_PETSC

#include "definition/compdef.hpp"
#include "definition/diffdef.hpp"
#include "definition/patchdef.hpp"
#include "definition/reacdef.hpp"
#include "definition/sreacdef.hpp"

#include "util/mpitools.hpp"

namespace steps::dist {

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
TetOpSplit<SSA, SearchMethod>::TetOpSplit(steps::model::Model& model,
                                          steps::dist::DistMesh& mesh,
                                          const rng::RNGptr& r,
                                          bool indepKProcs,
                                          bool isEfield)
    : meshref(mesh) {
#if USE_PETSC
    PetscErrorCode ierr = PetscInitialize(nullptr, nullptr, nullptr, "OpSplit solver");
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif  // USE_PETSC

    sim =
        std::make_unique<OmegaHSimulation<SSA, SearchMethod>>(meshref, *r, std::clog, indepKProcs);

    auto stateDefPtr = std::make_unique<Statedef>(model, mesh);
    if (!isEfield) {
        stateDefPtr->disableEField();
    }

    sim->init(std::move(stateDefPtr));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
TetOpSplit<SSA, SearchMethod>::~TetOpSplit() {
#if USE_PETSC
    PetscErrorCode ierr = PetscFinalize();
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif  // USE_PETSC
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverName() const {
    return "disttetopsplit";
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverDesc() const {
    return "";
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverAuthors() const {
    return "Stefan Wils, Iain Hepburn, Weiliang Chen";
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::reset() {
    sim->reset();
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::run(double seconds) {
    sim->run(seconds);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTime() const {
    return sim->getTime();
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getCompSpecCount(std::string const& c,
                                                       std::string const& s) const {
    return sim->getCompSpecCount(steps::dist::model::compartment_id(c),
                                 steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setCompSpecCount(std::string const& c,
                                                     std::string const& s,
                                                     double n,
                                                     const math::DistributionMethod distribution) {
    sim->setCompSpecCount(steps::dist::model::compartment_id(c),
                          steps::dist::model::species_name(s),
                          n,
                          distribution);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getCompSpecConc(std::string const& c,
                                                      std::string const& s) const {
    return sim->getCompSpecConc(steps::dist::model::compartment_id(c),
                                steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setCompSpecConc(std::string const& c,
                                                    std::string const& s,
                                                    double conc,
                                                    const math::DistributionMethod distribution) {
    sim->setCompSpecConc(steps::dist::model::compartment_id(c),
                         steps::dist::model::species_name(s),
                         conc,
                         distribution);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getPatchSpecCount(std::string const& p,
                                                        std::string const& s) const {
    return sim->getPatchSpecCount(steps::dist::model::patch_id(p),
                                  steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setPatchSpecCount(std::string const& p,
                                                      std::string const& s,
                                                      double n,
                                                      const math::DistributionMethod distribution) {
    return sim->setPatchSpecCount(steps::dist::model::patch_id(p),
                                  steps::dist::model::species_name(s),
                                  n,
                                  distribution);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setPatchSReacK(std::string const& p,
                                                   std::string const& r,
                                                   double kf) {
    sim->setPatchSReacK(steps::dist::model::patch_id(p),
                        steps::dist::model::surface_reaction_id(r),
                        kf);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getPatchMaxV(std::string const& p) const {
    return sim->getMaxPotentialOnTriangles(steps::dist::model::patch_id(p));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setMembPotential(const std::string& memb, double pot) {
    sim->setMembPotential(steps::dist::model::membrane_id(memb), pot);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setMembRes(const std::string& membrane,
                                               osh::Real resistivity,
                                               osh::Real reversal_potential) {
    sim->setMembRes(steps::dist::model::membrane_id(membrane), resistivity, reversal_potential);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
MembraneResistivity TetOpSplit<SSA, SearchMethod>::getMembRes(const std::string& membrane) const {
    return sim->getMembRes(steps::dist::model::membrane_id(membrane));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
MembraneResistivity TetOpSplit<SSA, SearchMethod>::getTriRes(osh::GO tri, bool local) const {
    return sim->getTriRes(tri, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTriRes(const osh::GO tri,
                                              osh::Real res,
                                              osh::Real erev,
                                              bool local) {
    sim->setTriRes(tri, res, erev, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriCapac(osh::GO tri, bool local) const {
    return sim->getTriCapac(tri, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTriCapac(const osh::GO tri, osh::Real c, bool local) {
    sim->setTriCapac(tri, c, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getVertIClamp(const osh::GO vert, bool local) const {
    return sim->getVertIClamp(vert, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setVertIClamp(const osh::GO vert,
                                                  const double current,
                                                  bool local) {
    sim->setVertIClamp(vert, current, local);
}

#if USE_PETSC

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriOhmicErev(osh::GO triangle,
                                                      const std::string& ohmic_current,
                                                      bool local) const {
    return sim->getTriOhmicErev(triangle, model::ohmic_current_id(ohmic_current), local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriOhmicErevsNP(const osh::GO* indices,
                                                            size_t input_size,
                                                            const std::string& ohmic_current,
                                                            double* rv,
                                                            size_t output_size,
                                                            bool local) const {
    assert(input_size == output_size);
    sim->getBatchTriOhmicErevsNP(gsl::span(indices, input_size),
                                 model::ohmic_current_id(ohmic_current),
                                 gsl::span(rv, output_size),
                                 local);
}


template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTriOhmicErev(osh::GO triangle,
                                                    const std::string& ohmic_current,
                                                    double reversal_potential,
                                                    bool local) {
    sim->setTriOhmicErev(triangle,
                         steps::dist::model::ohmic_current_id(ohmic_current),
                         reversal_potential,
                         local);
}

#endif  // USE_PETSC

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetSpecCount(osh::GO tet,
                                                      const std::string& s,
                                                      bool local) const {
    double res;
    getBatchTetSpecCountsNP(&tet, 1, s, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetSpecConc(osh::GO tet,
                                                     const std::string& s,
                                                     bool local) const {
    double res;
    getBatchTetSpecConcsNP(&tet, 1, s, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTetSpecCount(osh::GO tet,
                                                    const std::string& s,
                                                    double count,
                                                    bool local) {
    setBatchTetSpecCountsNP(&tet, 1, s, &count, 1, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTetSpecConc(osh::GO tet,
                                                   const std::string& s,
                                                   double conc,
                                                   bool local) {
    setBatchTetSpecConcsNP(&tet, 1, s, &conc, 1, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriSpecCount(osh::GO tri,
                                                      const std::string& s,
                                                      bool local) const {
    double res;
    getBatchTriSpecCountsNP(&tri, 1, s, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTriSpecCount(osh::GO tri,
                                                    const std::string& s,
                                                    double count,
                                                    bool local) {
    setBatchTriSpecCountsNP(&tri, 1, s, &count, 1, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getVertV(osh::GO vert, bool local) const {
    double res;
    getBatchVertVsNP(&vert, 1, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriV(osh::GO tri, bool local) const {
    double res;
    getBatchTriVsNP(&tri, 1, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetV(osh::GO tet, bool local) const {
    double res;
    getBatchTetVsNP(&tet, 1, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriOhmicI(osh::GO tri,
                                                   const std::string& curr,
                                                   bool local) const {
    double res;
    getBatchTriOhmicIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriGHKI(osh::GO tri,
                                                 const std::string& curr,
                                                 bool local) const {
    double res;
    getBatchTriGHKIsNP(&tri, 1, curr, &res, 1, local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTetSpecCounts(
    const std::vector<osh::GO>& tets,
    std::string const& s,
    bool local) const {
    std::vector<double> res(tets.size());
    getBatchTetSpecCountsNP(tets.data(), tets.size(), s, res.data(), res.size(), local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTetSpecConcs(
    const std::vector<osh::GO>& tets,
    std::string const& s,
    bool local) const {
    std::vector<double> res(tets.size());
    getBatchTetSpecConcsNP(tets.data(), tets.size(), s, res.data(), res.size(), local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetSpecCounts(const std::vector<osh::GO>& tets,
                                                          std::string const& s,
                                                          std::vector<double>& counts,
                                                          bool local) {
    setBatchTetSpecCountsNP(tets.data(), tets.size(), s, counts.data(), counts.size(), local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetSpecConcs(const std::vector<osh::GO>& tets,
                                                         std::string const& s,
                                                         std::vector<double>& concs,
                                                         bool local) {
    setBatchTetSpecConcsNP(tets.data(), tets.size(), s, concs.data(), concs.size(), local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTriSpecCounts(
    const std::vector<osh::GO>& tris,
    std::string const& s,
    bool local) const {
    std::vector<double> res(tris.size());
    getBatchTriSpecCountsNP(tris.data(), tris.size(), s, res.data(), res.size(), local);
    return res;
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTriSpecCounts(const std::vector<osh::GO>& tris,
                                                          std::string const& s,
                                                          std::vector<double>& counts,
                                                          bool local) {
    setBatchTriSpecCountsNP(tris.data(), tris.size(), s, counts.data(), counts.size(), local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetSpecCountsNP(const osh::GO* indices,
                                                            size_t input_size,
                                                            std::string const& s,
                                                            double* counts,
                                                            size_t output_size,
                                                            bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchElemValsNP(
        indices, input_size, steps::dist::model::species_name(s), counts, false, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetSpecConcsNP(const osh::GO* indices,
                                                           size_t input_size,
                                                           std::string const& s,
                                                           double* concs,
                                                           size_t output_size,
                                                           bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchElemValsNP(
        indices, input_size, steps::dist::model::species_name(s), concs, true, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetSpecCountsNP(const osh::GO* indices,
                                                            size_t input_size,
                                                            std::string const& s,
                                                            double* counts,
                                                            size_t output_size,
                                                            bool local) {
    assert(input_size == output_size);
    (void) output_size;
    sim->setBatchElemValsNP(
        indices, input_size, steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetSpecConcsNP(const osh::GO* indices,
                                                           size_t input_size,
                                                           std::string const& s,
                                                           double* concs,
                                                           size_t output_size,
                                                           bool local) {
    assert(input_size == output_size);
    (void) output_size;
    sim->setBatchElemValsNP(
        indices, input_size, steps::dist::model::species_name(s), concs, true, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriSpecCountsNP(const osh::GO* indices,
                                                            size_t input_size,
                                                            std::string const& s,
                                                            double* counts,
                                                            size_t output_size,
                                                            bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchBoundSpecCountNP(
        indices, input_size, steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTriSpecCountsNP(const osh::GO* indices,
                                                            size_t input_size,
                                                            std::string const& s,
                                                            double* counts,
                                                            size_t output_size,
                                                            bool local) {
    assert(input_size == output_size);
    (void) output_size;
    sim->setBatchBoundSpecCountNP(
        indices, input_size, steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchVertVsNP(const osh::GO* indices,
                                                     size_t input_size,
                                                     double* voltages,
                                                     size_t output_size,
                                                     bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchVertVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriVsNP(const osh::GO* indices,
                                                    size_t input_size,
                                                    double* voltages,
                                                    size_t output_size,
                                                    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchTriVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetVsNP(const osh::GO* indices,
                                                    size_t input_size,
                                                    double* voltages,
                                                    size_t output_size,
                                                    bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchTetVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriOhmicIsNP(const osh::GO* indices,
                                                         size_t input_size,
                                                         const std::string& curr,
                                                         double* currents,
                                                         size_t output_size,
                                                         bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchTriOhmicIsNP(
        indices, input_size, steps::dist::model::ohmic_current_id(curr), currents, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriGHKIsNP(const osh::GO* indices,
                                                       size_t input_size,
                                                       const std::string& curr,
                                                       double* currents,
                                                       size_t output_size,
                                                       bool local) const {
    assert(input_size == output_size);
    (void) output_size;
    sim->getBatchTriGHKIsNP(
        indices, input_size, steps::dist::model::ghk_current_id(curr), currents, local);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setDiffBoundarySpecDiffusionActive(const std::string& name,
                                                                       const std::string& spec,
                                                                       bool active) {
    sim->setDiffusionBoundaryActive(steps::dist::mesh::diffusion_boundary_name(name),
                                    steps::dist::model::species_name(spec),
                                    active);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
bool TetOpSplit<SSA, SearchMethod>::getDiffBoundarySpecDiffusionActive(const std::string& name,
                                                                       const std::string& spec) {
    return sim->getDiffusionBoundaryActive(steps::dist::mesh::diffusion_boundary_name(name),
                                           steps::dist::model::species_name(spec));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setDiffApplyThreshold(int threshold) {
    sim->setDiffOpBinomialThreshold(threshold);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setMembIClamp(const std::string& memb, double stim) {
    sim->setMembIClamp(steps::dist::model::membrane_id(memb), stim);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::dumpDepGraphToFile(const std::string& path) const {
    std::ofstream ostr(path);
    sim->dumpDepGraphToFile(path);
}

#ifdef USE_PETSC

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getEfieldDT() const {
    return *sim->getEfieldDt();
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setEfieldDT(double dt) {
    sim->setEfieldDt(dt);
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setPetscOptions(const std::string& s) {
    auto err = PetscOptionsInsertString(nullptr, s.c_str());
    CHKERRABORT(meshref.comm_impl(), err);

    // the check to see if we have the efield active is done in simulation
    sim->setPetscOptions();
}

#endif  // USE_PETSC

// explicit template instantiation definitions

template class TetOpSplit<steps::dist::SSAMethod::SSA, steps::dist::NextEventSearchMethod::Direct>;
template class TetOpSplit<steps::dist::SSAMethod::SSA,
                          steps::dist::NextEventSearchMethod::GibsonBruck>;
template class TetOpSplit<steps::dist::SSAMethod::RSSA, steps::dist::NextEventSearchMethod::Direct>;

}  // namespace steps::dist

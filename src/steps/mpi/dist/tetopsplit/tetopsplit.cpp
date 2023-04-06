#include "tetopsplit.hpp"

#include <iostream>

#include "definition/compdef.hpp"
#include "definition/diffdef.hpp"
#include "definition/patchdef.hpp"
#include "definition/reacdef.hpp"
#include "definition/sreacdef.hpp"

#include "util/mpitools.hpp"

namespace steps {
namespace dist {

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
TetOpSplit<SSA, SearchMethod>::TetOpSplit(steps::model::Model &model,
                                          steps::dist::DistMesh &mesh,
                                          const rng::RNGptr &r,
                                          bool indepKProcs)
    : meshref(mesh) {
#if USE_PETSC
    PetscErrorCode ierr =
        PetscInitialize(nullptr, nullptr, nullptr, "OpSplit solver");
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif // USE_PETSC

    steps::dist::ScenarioInput input{""};
    input.sim_indep_k_procs = indepKProcs;

    sim =
        std::make_unique<steps::dist::OmegaHSimulation<SSA, steps::rng::RNG,
                                                       osh::I64, SearchMethod>>(
            input, meshref, *r, std::clog);

    sim->init(std::make_unique<steps::dist::Statedef>(model, mesh));
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
TetOpSplit<SSA, SearchMethod>::~TetOpSplit() {
#if USE_PETSC
    PetscErrorCode ierr = PetscFinalize();
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif // USE_PETSC
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverName() const {
    return "disttetopsplit";
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverDesc() const {
    return "";
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverAuthors() const {
    return "Stefan Wils, Iain Hepburn, Weiliang Chen";
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::string TetOpSplit<SSA, SearchMethod>::getSolverEmail() const {
    return "steps.dev@gmail.com";
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::reset() { sim->reset(); }

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::run(double seconds) {
    sim->run(seconds);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTime() const {
    return sim->getTime();
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getCompCount(std::string const &c,
                                                   std::string const &s) const {
    return sim->getCompCount(steps::dist::model::compartment_id(c),
                             steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setCompCount(std::string const& c,
                                                 std::string const& s,
                                                 double n,
                                                 const math::DistributionMethod distribution) {
    sim->setCompCount(steps::dist::model::compartment_id(c),
                      steps::dist::model::species_name(s),
                      n,
                      distribution);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getCompConc(std::string const &c,
                                                  std::string const &s) const {
    return sim->getCompConc(steps::dist::model::compartment_id(c),
                            steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA, steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setCompConc(std::string const& c,
                                                std::string const& s,
                                                double conc,
                                                const math::DistributionMethod distribution) {
    sim->setCompConc(steps::dist::model::compartment_id(c),
                     steps::dist::model::species_name(s),
                     conc,
                     distribution);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double
TetOpSplit<SSA, SearchMethod>::getPatchCount(std::string const &p,
                                             std::string const &s) const {
  return sim->getPatchCount(steps::dist::model::patch_id(p),
                            steps::dist::model::species_name(s));
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setPatchCount(
    std::string const &p, std::string const &s, double n,
    const math::DistributionMethod distribution) {
  return sim->setPatchCount(steps::dist::model::patch_id(p),
                            steps::dist::model::species_name(s), n,
                            distribution);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setPatchSReacK(std::string const &p,
                                                   std::string const &r,
                                                   double kf) {
  sim->setPatchSReacK(steps::dist::model::patch_id(p),
                      steps::dist::model::surface_reaction_id(r), kf);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getPatchMaxV(std::string const &p) const {
  return sim->getMaxPotentialOnTriangles(steps::dist::model::patch_id(p));
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setMembPotential(const std::string &memb,
                                                     double pot) {
  sim->setMembPotential(steps::dist::model::membrane_id(memb), pot);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getVertIClamp(const osh::GO vert,
                                                    bool local) const {
  return sim->getVertIClamp(vert, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setVertIClamp(const osh::GO vert,
                                                  const double current,
                                                  bool local) {
  sim->setVertIClamp(vert, current, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetCount(osh::GO tet,
                                                  const std::string &s,
                                                  bool local) const {
  double res;
  getBatchTetCountsNP(&tet, 1, s, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetConc(osh::GO tet,
                                                 const std::string &s,
                                                 bool local) const {
  double res;
  getBatchTetConcsNP(&tet, 1, s, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTetCount(osh::GO tet,
                                                const std::string &s,
                                                double count, bool local) {
  setBatchTetCountsNP(&tet, 1, s, &count, 1, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTetConc(osh::GO tet,
                                               const std::string &s,
                                               double conc, bool local) {
  setBatchTetConcsNP(&tet, 1, s, &conc, 1, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriCount(osh::GO tri,
                                                  const std::string &s,
                                                  bool local) const {
  double res;
  getBatchTriCountsNP(&tri, 1, s, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setTriCount(osh::GO tri,
                                                const std::string &s,
                                                double count, bool local) {
  setBatchTriCountsNP(&tri, 1, s, &count, 1, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getVertV(osh::GO vert, bool local) const {
  double res;
  getBatchVertVsNP(&vert, 1, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriV(osh::GO tri, bool local) const {
  double res;
  getBatchTriVsNP(&tri, 1, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTetV(osh::GO tet, bool local) const {
  double res;
  getBatchTetVsNP(&tet, 1, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriOhmicI(osh::GO tri,
                                                   const std::string &curr,
                                                   bool local) const {
  double res;
  getBatchTriOhmicIsNP(&tri, 1, curr, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getTriGHKI(osh::GO tri,
                                                 const std::string &curr,
                                                 bool local) const {
  double res;
  getBatchTriGHKIsNP(&tri, 1, curr, &res, 1, local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTetCounts(
    const std::vector<osh::GO> &tets, std::string const &s, bool local) const {
  std::vector<double> res(tets.size());
  getBatchTetCountsNP(tets.data(), tets.size(), s, res.data(), res.size(),
                      local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTetConcs(
    const std::vector<osh::GO> &tets, std::string const &s, bool local) const {
  std::vector<double> res(tets.size());
  getBatchTetConcsNP(tets.data(), tets.size(), s, res.data(), res.size(),
                     local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetCounts(
    const std::vector<osh::GO> &tets, std::string const &s,
    std::vector<double> &counts, bool local) {
  setBatchTetCountsNP(tets.data(), tets.size(), s, counts.data(), counts.size(),
                      local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetConcs(
    const std::vector<osh::GO> &tets, std::string const &s,
    std::vector<double> &concs, bool local) {
  setBatchTetConcsNP(tets.data(), tets.size(), s, concs.data(), concs.size(),
                     local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
std::vector<double> TetOpSplit<SSA, SearchMethod>::getBatchTriCounts(
    const std::vector<osh::GO> &tris, std::string const &s, bool local) const {
  std::vector<double> res(tris.size());
  getBatchTriCountsNP(tris.data(), tris.size(), s, res.data(), res.size(),
                      local);
  return res;
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTriCounts(
    const std::vector<osh::GO> &tris, std::string const &s,
    std::vector<double> &counts, bool local) {
  setBatchTriCountsNP(tris.data(), tris.size(), s, counts.data(), counts.size(),
                      local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetCountsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *counts, size_t output_size, bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchElemValsNP(indices, input_size,
                          steps::dist::model::species_name(s), counts, false,
                          local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetConcsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *concs, size_t output_size, bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchElemValsNP(indices, input_size,
                          steps::dist::model::species_name(s), concs, true,
                          local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetCountsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *counts, size_t output_size, bool local) {
  assert(input_size == output_size);
  (void)output_size;
  sim->setBatchElemValsNP(indices, input_size,
                          steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTetConcsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *concs, size_t output_size, bool local) {
  assert(input_size == output_size);
  (void)output_size;
  sim->setBatchElemValsNP(indices, input_size,
                          steps::dist::model::species_name(s), concs, true,
                          local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriCountsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *counts, size_t output_size, bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchBoundCountNP(indices, input_size,
                            steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setBatchTriCountsNP(
    const osh::GO *indices, size_t input_size, std::string const &s,
    double *counts, size_t output_size, bool local) {
  assert(input_size == output_size);
  (void)output_size;
  sim->setBatchBoundCountNP(indices, input_size,
                            steps::dist::model::species_name(s), counts, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchVertVsNP(const osh::GO *indices,
                                                     size_t input_size,
                                                     double *voltages,
                                                     size_t output_size,
                                                     bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchVertVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriVsNP(const osh::GO *indices,
                                                    size_t input_size,
                                                    double *voltages,
                                                    size_t output_size,
                                                    bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchTriVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTetVsNP(const osh::GO *indices,
                                                    size_t input_size,
                                                    double *voltages,
                                                    size_t output_size,
                                                    bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchTetVsNP(indices, input_size, voltages, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriOhmicIsNP(
    const osh::GO *indices, size_t input_size, const std::string &curr,
    double *currents, size_t output_size, bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchTriOhmicIsNP(indices, input_size,
                            steps::dist::model::ohmic_current_id(curr),
                            currents, local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::getBatchTriGHKIsNP(
    const osh::GO *indices, size_t input_size, const std::string &curr,
    double *currents, size_t output_size, bool local) const {
  assert(input_size == output_size);
  (void)output_size;
  sim->getBatchTriGHKIsNP(indices, input_size,
                          steps::dist::model::ghk_current_id(curr), currents,
                          local);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setDiffBoundaryDiffusionActive(
    const std::string &name, const std::string &spec, bool active) {
  sim->setDiffusionBoundaryActive(
      steps::dist::mesh::diffusion_boundary_name(name),
      steps::dist::model::species_name(spec), active);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
bool TetOpSplit<SSA, SearchMethod>::getDiffBoundaryDiffusionActive(
    const std::string &name, const std::string &spec) {
  return sim->getDiffusionBoundaryActive(
      steps::dist::mesh::diffusion_boundary_name(name),
      steps::dist::model::species_name(spec));
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setDiffApplyThreshold(int threshold) {
  sim->setDiffOpBinomialThreshold(threshold);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setMembIClamp(const std::string &memb,
                                                  double stim) {
  sim->setMembIClamp(steps::dist::model::membrane_id(memb), stim);
}

#ifdef USE_PETSC

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
double TetOpSplit<SSA, SearchMethod>::getEfieldDT() const {
    return sim->getEfieldDt().get();
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setEfieldDT(double dt) {
    sim->setEfieldDt(dt);
}

template <steps::dist::SSAMethod SSA,
          steps::dist::NextEventSearchMethod SearchMethod>
void TetOpSplit<SSA, SearchMethod>::setEfieldTolerances(double atol,
                                                        double rtol,
                                                        KSPNormType norm_type) {
  sim->setEfieldTolerances(atol, rtol, norm_type);
}

#endif // USE_PETSC

// explicit template instantiation definitions

template class TetOpSplit<steps::dist::SSAMethod::SSA,
                          steps::dist::NextEventSearchMethod::Direct>;
template class TetOpSplit<steps::dist::SSAMethod::SSA,
                          steps::dist::NextEventSearchMethod::GibsonBruck>;
template class TetOpSplit<steps::dist::SSAMethod::RSSA,
                          steps::dist::NextEventSearchMethod::Direct>;

}  // namespace dist
}  // namespace steps

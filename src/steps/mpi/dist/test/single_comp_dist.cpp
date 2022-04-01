#include "single_comp_dist.hpp"

#include <complex>
#include <map>

#include "geom/dist/distmesh.hpp"
#include "model/model.hpp"
#include "mpi/dist/test/simulation.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/diffdef.hpp"
#include "mpi/dist/tetopsplit/definition/patchdef.hpp"
#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/collections.hpp"

#undef MPI_Scatter

namespace steps {
namespace dist {

SingleCompDist::SingleCompDist(const ScenarioInput &t_input)
    : Scenario("SingleCompDist", "Single compartment distribution", t_input) {}

std::unique_ptr<Statedef>
SingleCompDist::createStatedef(const simulation_t &simulation) const {
  steps::model::Model model;
  SingleCompDistSimdef simdef(model, simulation.getMesh());
  return std::move(simdef.getStatedef());
}

void SingleCompDist::register_compartments(DistMesh &mesh) const {
  mesh.addComp("comp1", model::compartment_label(1));
}

void SingleCompDist::fill_compartments(simulation_t &simulation) const {
  steps::model::Model model;
  simulation.setCompCount(SingleCompDistSimdef(model, simulation.getMesh()).getCompartementCounts(),
                          math::DistributionMethod::DIST_MULTINOMIAL);
}


int SingleCompDist::check_and_log_results_impl(simulation_t &simulation) const {
  bool status = EXIT_SUCCESS;

  auto tot_elem_counts = simulation.getElemCount("C");

  for (std::size_t i = 0; i < tot_elem_counts.first.size(); ++i) {
    simulation.log_once(
        "Elem global id: " + std::to_string(tot_elem_counts.first[i]) +
        " mol count: " + std::to_string(tot_elem_counts.second[i]));
  }

  const auto comp1 = model::compartment_id("comp1");

  osh::LO n_mols = std::accumulate(tot_elem_counts.second.begin(),
                                   tot_elem_counts.second.end(), 0);

  simulation.log_once("total number of mols: " + std::to_string(n_mols));

  // I did not find any other way to get simdef in a member function that is not
  // const. This should be equal to what is written in simdef.cpp
  osh::LO expected_n_mols = 100000;
  if (!tot_elem_counts.first.empty() &&
      n_mols != expected_n_mols) { // compCount is non-zero only for rank 0
    status = EXIT_FAILURE;
  }

  // check that the distribution is multinomial with the chisquared test
  const osh::Real tot_vol = simulation.getMesh().total_measure(comp1);

  osh::Reals volumes;
  osh::LOs elem_local_IDs;
  osh::Real local_volume;
  std::tie(elem_local_IDs, volumes, local_volume) =
      simulation.getMesh().measure(comp1);
  auto local_elem_counts = simulation.getOwnedElemCount("C");

  // this is \sum _ID (Observed - Expected)^2 / Expected
  osh::Real local_chi_frac(0.0);
  for (osh::LO i = 0; i < elem_local_IDs.size(); ++i) {
    const osh::Real expected_val = expected_n_mols * volumes[i] / tot_vol;
    const osh::Real diff = local_elem_counts.second[i] - expected_val;
    local_chi_frac += diff * diff / expected_val;
  }

  osh::Real tot_chi_frac(0.0);
  int err = MPI_Allreduce(&local_chi_frac, &tot_chi_frac, 1, MPI_DOUBLE,
                          MPI_SUM, simulation.comm());

  if (err != MPI_SUCCESS) {
    MPI_Abort(simulation.comm(), err);
  }

  // chi square distribution. Confidence level 0.01
  const osh::Real threshold = 18.48;

  simulation.log_once("total chi ratio: " + std::to_string(tot_chi_frac) +
                      " (< " + std::to_string(threshold) + ")");

  if (tot_chi_frac > threshold) {
    status = EXIT_FAILURE;
  }

  //////////////// Uniform distribution
  simulation.log_once("Uniform distribution");
  steps::model::Model model;
  simulation.setCompCount(SingleCompDistSimdef(model, simulation.getMesh()).getCompartementCounts(),
                          math::DistributionMethod::DIST_UNIFORM);

  const osh::Real count_in_rank = simulation.getOwnedCompCount("comp1", "C");

  simulation.log_all("Count fraction: " + std::to_string(count_in_rank / expected_n_mols) +
                     " volume fraction: " + std::to_string(local_volume / tot_vol) + '\n');

  if (static_cast<int>(1000 * count_in_rank / expected_n_mols) !=
      static_cast<int>(1000 * local_volume / tot_vol)) {
      status = EXIT_FAILURE;
  }

  local_elem_counts = simulation.getOwnedElemCount("C");
  simulation.log_once("Elem fractions approximated at the 3rd digit");
  for (osh::LO i = 0; i < elem_local_IDs.size(); ++i) {
      const int expected_val = std::round(1000 * volumes[i] / tot_vol);

      const int val = std::round(1000 * static_cast<double>(local_elem_counts.second[i]) /
                                 static_cast<double>(expected_n_mols));

      if (expected_val != val) {
          status = EXIT_FAILURE;
      }

      simulation.log_all("Elem local id: " + std::to_string(i) +
                         " mol fraction*1000: " + std::to_string(val) +
                         " "
                         "expected*1000: " +
                         std::to_string(expected_val));
  }

  return status;
}

} // namespace dist
} // namespace steps

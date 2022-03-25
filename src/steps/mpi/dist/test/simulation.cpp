#include "simulation.hpp"

#include <iostream>

#include "geom/dist/distmesh.hpp"
#include "util/mpitools.hpp"

namespace steps {
namespace dist {

template <typename RNG>
Simulation<RNG>::Simulation(const ScenarioInput &t_scenario, DistMesh &t_mesh,
                            RNG &t_rng, std::ostream &t_outstream)
    : comm_rank(util::mpi_comm_rank(t_mesh.comm_impl())),
      comm_size(util::mpi_comm_size(t_mesh.comm_impl())), scenario(t_scenario),
      mesh(t_mesh), rng(t_rng), outstream(t_outstream) {}

template <typename RNG> Simulation<RNG>::~Simulation() noexcept = default;

template <typename RNG>
void Simulation<RNG>::setCompCount(const Simdef::compartment_counts_t& counts,
                                   const math::DistributionMethod distribution) {
    for (const auto& comp_counts: counts) {
        setCompCount(comp_counts.first, comp_counts.second, distribution);
    }
}

template <typename RNG>
void Simulation<RNG>::log_all(const std::string &message) const {
  this->outstream << '[' << this->comm_rank << "] " << message << '\n';
}

template <typename RNG>
void Simulation<RNG>::log_once(const std::string &message,
                               bool force_stdout) const {
  if (this->comm_rank == 0) {
    if (force_stdout) {
      std::cout << message << '\n';
    } else {
      this->outstream << message << '\n';
    }
  }
}

template <typename RNG>
void Simulation<RNG>::setCompConc(const Simdef::compartment_concs_t& concentrations,
                                  const math::DistributionMethod distribution) {
    for (const auto& comp_concs: concentrations) {
        setCompConc(comp_concs.first, comp_concs.second, distribution);
    }
}

template <typename RNG> void Simulation<RNG>::log_diffusion_exchanges() const {
  const auto &local_exchanges = get_diffusion_rank_exchanges();
  std::vector<unsigned int> global_exchanges(local_exchanges.size());
  MPI_Gather(local_exchanges.data() + comm_rank * comm_size,  // NOLINT
             comm_size,
             MPI_UNSIGNED,
             global_exchanges.data(),
             comm_size,
             MPI_UNSIGNED,
             0,
             mesh.comm_impl());
  std::ostringstream oss;
  oss << "Diffusion exchanges rates:\n";
  for (int i = 0; i < comm_size; ++i) {
    for (int j = 0; j < comm_size; ++j) {
      oss << global_exchanges[static_cast<size_t>(i * comm_size + j)] << ';';
    }
    oss << '\n';
  }
  log_once(oss.str());
}

template <typename RNG>
void Simulation<RNG>::log_progress(const double i, const double tot,
                                   const std::string &name) const {
  std::stringstream s;
  s << name << " progress: " << std::round(1000 * i / tot) / 10 << "%";
  log_once(s.str(), true);
}

// explicit template instantiation definition
template class Simulation<std::mt19937>;

template class Simulation<steps::rng::RNG>;

} // namespace dist
} // namespace steps

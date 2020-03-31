#include <mpitools.hpp>

#include "opsplit/test/simulation.hpp"

namespace zee {

template <typename RNG>
Simulation<RNG>::Simulation(const ScenarioInput& t_scenario, DistMesh& t_mesh, RNG& t_rng)
    : comm_rank(mpi_comm_rank())
    , comm_size(mpi_comm_size())
    , scenario(t_scenario)
    , mesh(t_mesh)
    , rng(t_rng) {}

template <typename RNG>
Simulation<RNG>::~Simulation() noexcept = default;

template <typename RNG>
void Simulation<RNG>::log_diffusion_exchanges() const {
    MPI_Barrier(MPI_COMM_WORLD);
    const auto& local_exchanges = get_diffusion_rank_exchanges();
    std::vector<unsigned int> global_exchanges(local_exchanges.size());
    MPI_Gather(local_exchanges.data() + comm_rank * comm_size,  // NOLINT
               comm_size,
               MPI_UNSIGNED,
               global_exchanges.data(),
               comm_size,
               MPI_UNSIGNED,
               0,
               MPI_COMM_WORLD);
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

// explicit template instantiation definition
template class Simulation<std::mt19937>;

}  // namespace zee

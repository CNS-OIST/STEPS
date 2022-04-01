#include "measure.hpp"

#include "mpi/mpi_init.hpp"

namespace steps {
namespace dist {

Measure::Measure(MPI_Comm comm,
                 osh::Int num_compartments,
                 const element_measure_func& t_element_measure_func)
    : comm_(comm)
    , rank_(steps::mpi::getRank(comm_))
    , num_measures_per_rank_(1 + num_compartments)
    , element_measure_func_(t_element_measure_func) {}

void Measure::init(const mesh::tetrahedron_ids &t_owned_elements,
                   const osh::LOs &t_elem2compid) {
  const auto comm_size = steps::mpi::getNHosts(comm_);
  osh::Write<osh::Real> rank2measures(num_measures_per_rank_ * comm_size, 0);
  auto rankmeasures = rank2measures.data() + rank_ * num_measures_per_rank_;
  for (auto element : t_owned_elements) {
    const auto measure = element_measure_func_(element);
    const auto compid = t_elem2compid[element.get()];
    rankmeasures[0] += measure;
    if (compid != INITIAL_COMPARTMENT_ID) {
      rankmeasures[1 + compid] += measure;
    }
  }

  {
    int err = MPI_Allgather(MPI_IN_PLACE, num_measures_per_rank_, MPI_DOUBLE,
                            rank2measures.data(), num_measures_per_rank_,
                            MPI_DOUBLE, comm_);
    if (err != 0) {
      MPI_Abort(comm_, err);
    }
    rank2measures_ = rank2measures;
  }

  {
    osh::Write<osh::Real> mesh_measures(num_measures_per_rank_, 0);
    for (auto rank = 0; rank < comm_size; ++rank) {
      for (auto i = 0; i < num_measures_per_rank_; ++i) {
        mesh_measures[i] += rank2measures[rank * num_measures_per_rank_ + i];
      }
    }
    mesh_measures_ = mesh_measures;
  }
}

}  // namespace dist
}  // namespace steps

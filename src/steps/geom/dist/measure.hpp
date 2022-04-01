#pragma once

#include "util/vocabulary.hpp"

#include <Omega_h_array.hpp>
#include <mpi.h>

namespace steps {
namespace dist {

/**
 * Provide data and operations dealing with measure of mesh elements i.e
 * surface in 2D and volume in 3D.
 */
class Measure {
  public:
    /// type of a function taking an element identifier in parameter
    /// and return its measure.
    using element_measure_func =
        std::function<osh::Real(mesh::tetrahedron_local_id_t)>;

    /**
     *
     * \param t_num_elements Number of elements in the current rank
     * \param t_elem2compid get compartment identifier of an element
     * \param t_element_measure_func Functor to get measure of an element
     */
    Measure(MPI_Comm comm,
            osh::Int num_compartments,
            const element_measure_func& t_element_measure_func);

    void init(const mesh::tetrahedron_ids &t_owned_elements,
              const osh::LOs &t_elem2compid);

    /**
     * \return the sum of all element measures assigned to this MPI rank
     */
    inline osh::Real rank_measure() const noexcept {
        return rank_measure(rank_);
    }

    /**
     * \return the sum of all element measures in the given MPI rank
     */
    inline osh::Real rank_measure(osh::LO rank) const noexcept {
        return rank2measures_[rank * num_measures_per_rank_];
    }

    /**
     * \return the sum of all element measures that belong to the given compartment
     * in the current rank
     */
    inline osh::Real rank_measure(mesh::compartment_id compartment) const noexcept {
        return rank_measure(rank_, compartment);
    }

    /**
     * \return the sum of all element measures part of the given rank and compartment
     */
    inline osh::Real rank_measure(osh::Int rank, mesh::compartment_id compartment) const noexcept {
      return rank2measures_[rank * num_measures_per_rank_ + compartment.get() +
                            1];
    }

    /**
     * \return measure of the entire mesh
     */
    inline osh::Real mesh_measure() const noexcept {
        return mesh_measures_[0];
    }

    /**
     * \return measure of the given compartment in the entire mesh
     */
    inline osh::Real mesh_measure(mesh::compartment_id compartment) const noexcept {
        return mesh_measures_[1 + compartment.get()];
    }

    /**
     * \return measure of the given element
     */
    inline osh::Real element_measure(mesh::tetrahedron_id_t element) const
        noexcept {
      return element_measure_func_(element);
    }

    /**
     * Given a total number of molecules to spread in the current rank, this function
     * returns, for a given element and compartment, the number of molecules in regards of its
     * measure compared to the other elements. \param compartment target compartment identifier
     * \param element a mesh element
     * \param num_molecules number of molecules to spread in the current rank
     */
    inline osh::Real molecules_in_element(mesh::compartment_id compartment,
                                          mesh::tetrahedron_id_t element,
                                          osh::Real num_molecules) const {
      const auto element_measure = element_measure_func_(element);
      const auto rank_measure_p = rank_measure(rank_, compartment);
      return num_molecules * element_measure / rank_measure_p;
    }

    /**
     * \return function to get measure of a given element
     */
    inline const element_measure_func& measure_func() const noexcept {
        return element_measure_func_;
    }

  private:
    /// MPI Communicator of the mesh
    const MPI_Comm comm_;
    /// current rank
    const osh::Int rank_;
    /// 1 + number of compartments
    const osh::Int num_measures_per_rank_;
    /// functor to get measure of an element
    const element_measure_func& element_measure_func_;
    /// sum of element measure in the entire mesh
    osh::Reals mesh_measures_;
    /// measures of every compartment for every rank
    osh::Reals rank2measures_;
};

}  // namespace dist
}  // namespace steps

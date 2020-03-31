#include "simulation.hpp"

#include <iostream>
#include <memory>

#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>
#include <petsctime.h>

#include "mesh_utils.hpp"
#include "opsplit/compdef.hpp"
#include "opsplit/diffdef.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/reacdef.hpp"
#include "opsplit/sreacdef.hpp"

#include "math.hpp"

namespace zee {

template <osh::Int Dim, SSAMethod SSA, typename RNG>
OmegaHSimulation<Dim, SSA, RNG>::OmegaHSimulation(const ScenarioInput& t_scenario,
                                                  OmegaHMesh<Dim>& t_mesh,
                                                  RNG& t_rng)
    : super_type(t_scenario, t_mesh, t_rng)
    , mesh(t_mesh)
    , elems2verts(mesh.getMesh().ask_elem_verts())
    , coords(mesh.getMesh().coords()) {}


template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::distribute_num_molecules_on_ranks(
    const model::compartment_id& compartmentId,
    PetscScalar num_molecules,
    std::vector<PetscScalar>& result) {
    const auto compartment_id = mesh.getCompID(compartmentId);
    const auto& measureInfo = mesh.getMeasureInfo();
    auto remain_count = stochastic_round<osh::GO>(num_molecules, this->rng);
    result.resize(static_cast<size_t>(this->comm_size));
    for (auto rank = 0; rank < static_cast<int>(result.size()); ++rank) {
        const auto approximate_num_molecules =
            measureInfo.molecules_in_rank(rank, compartment_id, num_molecules);
        const auto part_num_molecules =
            stochastic_round<osh::GO>(approximate_num_molecules, this->rng, remain_count);
        result[static_cast<size_t>(rank)] = static_cast<PetscScalar>(part_num_molecules);
        remain_count -= part_num_molecules;
    }
    const auto mesh_measure = measureInfo.mesh_measure(compartment_id);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    while (remain_count > 0) {
        PetscScalar accum{};
        PetscScalar selector = dist(this->rng) * mesh_measure;
        for (auto rank = 0; rank < static_cast<osh::Int>(this->comm_size); ++rank) {
            accum += measureInfo.rank_measure(rank, compartment_id);
            if (selector < accum) {
                result[static_cast<size_t>(rank)] += 1;
                remain_count -= 1;
                break;
            }
        }
    }
}


//-----------------------------------------------

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getPatchCount(const model::patch_id& patch_id,
                                                           const model::specie_name& specie) const {
    const auto& elems = mesh.getOwnedEntities(patch_id);
    osh::Write<osh::LO> mols_counts(elems.size());
    container::specie_id spec_id{
        statedef->getPatchdef(patch_id).getSpecPatchIdx(statedef->getSpecModelIdx(specie))};
    for (auto i = 0; i < elems.size(); ++i) {
        mols_counts[i] = data->pools.moleculesCountOnPatchBoundaries()(elems[i].get(), spec_id);
    }
    return osh::get_sum(mesh.getMesh().comm(), osh::LOs(mols_counts));
}

//-----------------------------------------------

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setPatchCount(const model::patch_id& patch_id,
                                                    const model::specie_name& specie,
                                                    PetscScalar num_molecules) {
    osh::Read<osh::Real> areas;
    osh::LOs elems;
    osh::Real rank_area;
    std::tie(elems, areas, rank_area) = mesh.measure(patch_id);
    size_t num_molecules_on_rank;

    {
        std::vector<osh::Real> rank_areas;
        if (this->comm_rank == 0) {
            rank_areas.resize(static_cast<size_t>(this->comm_size));
        }
        int err = MPI_Gather(
            &rank_area, 1, MPIU_REAL, rank_areas.data(), 1, MPIU_REAL, 0, MPI_COMM_WORLD);

        if (err != 0) {
            MPI_Abort(MPI_COMM_WORLD, err);
        }
        MultinomialDistribution<size_t, PetscScalar> dist(static_cast<size_t>(num_molecules),
                                                          rank_areas);
        const auto& num_molecules_on_ranks = dist(this->rng);
        // send `num_molecules_on_ranks[i]` to rank `i` and store it in `num_molecules_on_rank`.
        err = MPI_Scatter(num_molecules_on_ranks.data(),
                          1,
                          MPI_UNSIGNED_LONG,
                          &num_molecules_on_rank,
                          1,
                          MPI_UNSIGNED_LONG,
                          0,
                          MPI_COMM_WORLD);
        if (err != 0) {
            MPI_Abort(MPI_COMM_WORLD, err);
        }
    }

    std::vector<osh::Real> areas_v(static_cast<size_t>(areas.size()));
    std::copy(areas.begin(), areas.end(), areas_v.begin());
    MultinomialDistribution<size_t, PetscScalar> dist(num_molecules_on_rank, areas_v);
    const auto& mols_on_elements = dist(this->rng);
    container::specie_id cont_spec_id = statedef->getPatchdef(patch_id).getSpecPatchIdx(
        statedef->getSpecModelIdx(specie));
    for (size_t k = 0; k < mols_on_elements.size(); ++k) {
        this->data->pools.moleculesCountOnPatchBoundaries()(elems[static_cast<osh::LO>(k)],
                                                            cont_spec_id) =
            static_cast<osh::LO>(mols_on_elements[k]);
    }
    osh::Write<osh::LO> vec(static_cast<osh::LO>(mols_on_elements.size()));
    for (osh::LO k = 0; k < vec.size(); ++k) {
        vec[k] = static_cast<osh::LO>(mols_on_elements[static_cast<size_t>(k)]);
    }
    assert(get_sum(mesh.getMesh().comm(), osh::LOs(vec)) == static_cast<osh::GO>(num_molecules));
}

//-----------------------------------------------

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setPatchCount(const Simdef::patch_counts_t& counts) {
    for (const auto& count: counts) {
        this->setPatchCount(count.patch, count.specie, count.num_mols);
    }
}

//-----------------------------------------------

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setCompCount(const model::compartment_id& comp_id,
                                                   const model::specie_name& specie,
                                                   PetscScalar num_molecules) {
    std::vector<PetscScalar> num_molecules_on_ranks;

    if (this->comm_rank == 0) {
        distribute_num_molecules_on_ranks(comp_id, num_molecules, num_molecules_on_ranks);
    }
    // send `num_molecules_on_ranks[i]` to rank `i` and store it in `spec_count`.
    int err = MPI_Scatter(num_molecules_on_ranks.data(),
                          1,
                          MPIU_REAL,
                          &num_molecules,
                          1,
                          MPIU_REAL,
                          0,
                          MPI_COMM_WORLD);
    if (err != 0) {
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    this->setOwnedCompCount(comp_id, specie, num_molecules);
    data->kproc_state.updateAllPropensities(data->pools);
}


template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setCompCount(const Simdef::compartment_counts_t& counts) {
    for (const auto& count: counts) {
        this->setCompCount(count.compartment, count.specie, count.num_mols);
    }
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setCompConc(const model::compartment_id& compartment,
                                                  const model::specie_name& specie,
                                                  PetscScalar concentration) {
    const auto spec_count = concentration * mesh.getTotalCompVol(mesh.getCompLabel(compartment)) *
                            1.0e3 * AVOGADRO;
    setCompCount(compartment, specie, spec_count);
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setCompConc(const Simdef::compartment_concs_t& concs) {
    for (const auto& conc: concs) {
        this->setCompConc(conc.compartment, conc.specie, conc.concentration);
    }
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setOwnedElementCount(const model::compartment_id& compartment,
                                                           PetscInt element,
                                                           const model::specie_name& specie,
                                                           PetscScalar num_molecules) {
    const auto specie_id = statedef->getCompSpecContainerIdx(compartment, specie);
    data->pools(mesh::element_id(static_cast<osh::LO>(element)),
                specie_id) = static_cast<osh::LO>(num_molecules);
    data->kproc_state.updateAllPropensities(data->pools);
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setOwnedCompCount(const model::compartment_id& compartment,
                                                        const model::specie_name& specie,
                                                        PetscScalar num_molecules) {
    auto& molecules = data->pools;
    const auto comp_id = mesh.getCompID(compartment);
    const container::specie_id specie_id = statedef->getCompSpecContainerIdx(compartment, specie);
    if (static_cast<PetscInt>(num_molecules) == 0) {
        for (auto elem: mesh.getOwnedElems(comp_id)) {
            molecules(elem, specie_id) = 0;
        }
        return;
    }

    std::uniform_real_distribution<double> dist;
    const MeasureInfo& measureInfo = mesh.getMeasureInfo();
    auto remaining_spec_count = num_molecules;
    for (const auto elem: mesh.getOwnedElems(comp_id)) {
        auto num_mols_in_elem = measureInfo.molecules_in_element(comp_id, elem, num_molecules);
        {
            auto num_mols_in_elem_actual = std::floor(num_mols_in_elem);
            auto num_mols_in_elem_fraction = num_mols_in_elem_actual - num_mols_in_elem;
            if (dist(this->rng) < num_mols_in_elem_fraction && num_mols_in_elem_fraction > 0) {
                num_mols_in_elem_actual += 1;
            }
            num_mols_in_elem = num_mols_in_elem_actual;
        }
        molecules(elem, specie_id) = static_cast<osh::LO>(num_mols_in_elem);
        remaining_spec_count -= static_cast<osh::LO>(num_mols_in_elem);
    }
    while (remaining_spec_count > 0) {
        PetscScalar accum = 0;
        PetscScalar selector = dist(this->rng) * measureInfo.rank_measure(comp_id);
        for (const auto elem: mesh.getOwnedElems(comp_id)) {
            accum += measureInfo.measure_func()(elem);
            if (selector < accum) {
                molecules(elem, specie_id) += 1;
                remaining_spec_count -= 1;
                break;
            }
        }
    }
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getOwnedCompConc(
    const model::compartment_id& compartment,
    const model::specie_name& specie) const {
    const auto spec_count = getOwnedCompCount(compartment, specie);
    return spec_count / (1.0e3 * mesh.getMeasureInfo().rank_measure() * AVOGADRO);
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getCompConc(const model::compartment_id& compartment,
                                                         const model::specie_name& specie) const {
    const auto spec_count = getCompCount(compartment, specie);
    return spec_count / (1.0e3 * mesh.getTotalCompVol(mesh.getCompLabel(compartment)) * AVOGADRO);
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getOwnedCompCount(
    const model::compartment_id& compartment,
    const model::specie_name& spec_id) const {
    const auto specie = statedef->getCompSpecContainerIdx(compartment, spec_id);
    const auto lambda = [=](osh::LO accu, mesh::element_id elem) {
        return accu + data->pools(elem, specie);
    };
    return static_cast<PetscScalar>(
        std::accumulate(mesh.getOwnedElems(mesh.getCompID(compartment)).begin(),
                        mesh.getOwnedElems(mesh.getCompID(compartment)).end(),
                        0,
                        lambda));
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getCompCount(const model::compartment_id& compartment,
                                                          const model::specie_name& specie) const {
    auto local_num_molecules = getOwnedCompCount(compartment, specie);
    PetscScalar global_num_molecules{};
    auto err = MPI_Allreduce(
        &local_num_molecules, &global_num_molecules, 1, MPIU_REAL, MPI_SUM, MPI_COMM_WORLD);
    if (err != 0) {
        MPI_Abort(MPI_COMM_WORLD, err);
    }
    return global_num_molecules;
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::setDiffOpBinomialThreshold(PetscScalar threshold) {
    data->diffOp.setBinomialThreshold(static_cast<osh::GO>(threshold));
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscScalar OmegaHSimulation<Dim, SSA, RNG>::getIterationTimeStep() const noexcept {
    return data->time_delta;
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::exportMolStateToVTK(const std::string& /* filename */) {
    // TODO(TCL) FIXME
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscInt64 OmegaHSimulation<Dim, SSA, RNG>::getDiffOpExtent() const {
    PetscInt64 global_ext, extent = data->diffOp.getExtent();
    MPI_Reduce(&extent, &global_ext, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_ext;
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscInt64 OmegaHSimulation<Dim, SSA, RNG>::getSSAOpExtent() const {
    return data->ssaOp.getExtent();
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
PetscInt64 OmegaHSimulation<Dim, SSA, RNG>::getNIterations() const noexcept {
    return num_iterations;
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
std::string OmegaHSimulation<Dim, SSA, RNG>::createStateReport() {
    // TODO(TCL) FIXME
    return "";
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::init(std::unique_ptr<Statedef>&& t_statedef) {
    this->statedef.swap(t_statedef);
    assert(statedef != nullptr);
    this->mesh.init();

    // fill a vector providing the number of species per element

    KProcState kproc_state_dry_run(*statedef, mesh, true);
    osh::LOs num_species_per_elements;
    boost::optional<osh::LOs> num_species_per_boundary_elements;
    std::tie(num_species_per_elements, num_species_per_boundary_elements) =
        kproc_state_dry_run.getNumberOfSpeciesPerOwnedElement();

    this->input =
        std::make_unique<SimulationInput<RNG>>(num_species_per_elements,
                                               num_species_per_boundary_elements,
                                               kproc_state_dry_run.getNumberOfSpeciesPerElement(),
                                               0 /*Unused in context*/,
                                               this->rng);

    data = std::make_unique<SimulationData<Dim, SSA, RNG>>(
        mesh,
        *this->statedef,
        *input,
        this->scenario.molecules_pools_force_dist_for_variable_sized,
        this->rng);
    initialize_discretized_rates();
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::reset() {
    this->num_iterations = 0;
    this->state_time = 0;
}


template <osh::LO Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::initialize_discretized_rates() {
    const auto& measure_info = mesh.getMeasureInfo();

    osh::parallel_for(
        this->mesh.getMesh().nelems(),
        [&measure_info, this](osh::LO e) __attribute__((always_inline, flatten)) {
            const mesh::element_id elem(e);
            const auto elem_measure = measure_info.element_measure(elem);
            const auto compartment_id = this->mesh.getCompartment(elem);
            const auto& compartment = this->statedef->getCompdef(compartment_id);
            for (const auto& diffusion: compartment.diffdefs()) {
                const auto specie = diffusion->getSpecContainerIdx();
                const auto dcst = diffusion->getDcst();
                data->diffusions.rates_sum(elem, specie) = 0;
                const auto num_neighbors = mesh.neighbors_per_element()[elem.get()];
                for (auto neighbor = 0; neighbor < num_neighbors; ++neighbor) {
                    const auto neighbor_boundary_distance =
                        mesh.tet_neighbors_real_data()(elem.get(), neighbor)[0];
                    const auto neighbor_boundary_measure =
                        mesh.tet_neighbors_real_data()(elem.get(), neighbor)[1];
                    const auto propensity = dcst * neighbor_boundary_measure / elem_measure /
                                            neighbor_boundary_distance;
                    data->diffusions.ith_rate(elem, specie, neighbor) = propensity;
                    data->diffusions.rates_sum(elem, specie) += propensity;
                }
            }
        });
    data->updateIterationTimeStep();
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::run(PetscScalar end_time) {
    while (state_time < end_time) {
        ++this->num_iterations;

        PetscScalar opsplit_period = std::min(data->time_delta, end_time - state_time);
        this->reactions_timer.resume();
        data->ssaOp.run(opsplit_period);
        this->reactions_timer.stop();

        this->diffusion_timer.resume();
        // TODO: Diff op should take as argument opsplit_period
        data->diffOp();
        this->diffusion_timer.stop();

        data->reset();
        state_time += opsplit_period;
    }
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::log_all(const std::string& message) const {
    std::clog << '[' << this->comm_rank << "] " << message << '\n';
}

template <osh::Int Dim, SSAMethod SSA, typename RNG>
void OmegaHSimulation<Dim, SSA, RNG>::log_once(const std::string& message) const {
    if (this->comm_rank == 0) {
        std::clog << message << '\n';
    }
}

// explicit template instantiation definitions
template class OmegaHSimulation<2, SSAMethod::SSA, std::mt19937>;
template class OmegaHSimulation<2, SSAMethod::RSSA, std::mt19937>;
template class OmegaHSimulation<3, SSAMethod::SSA, std::mt19937>;
template class OmegaHSimulation<3, SSAMethod::RSSA, std::mt19937>;

}  // namespace zee

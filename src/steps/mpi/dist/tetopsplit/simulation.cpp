#include "simulation.hpp"

#include <limits>
#include <memory>
#include <numeric>

#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

#if USE_PETSC
#include "mpi/dist/tetopsplit/operator/efield_operator.hpp"
#endif  // USE_PETSC

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


Simulation::Simulation(DistMesh& t_mesh, rng::RNG& t_rng, std::ostream& t_outstream)
    : comm_rank(util::mpi_comm_rank(t_mesh.comm_impl()))
    , comm_size(util::mpi_comm_size(t_mesh.comm_impl()))
    , mesh(t_mesh)
    , rng(t_rng)
    , outstream(t_outstream) {}

Simulation::~Simulation() noexcept = default;

void Simulation::setCompSpecCount(const compartment_counts_t& counts,
                                  const math::DistributionMethod distribution) {
    for (const auto& comp_counts: counts) {
        setCompSpecCount(comp_counts.first, comp_counts.second, distribution);
    }
}

void Simulation::log_all(const std::string& message) const {
    this->outstream << '[' << this->comm_rank << "] " << message << '\n';
}

void Simulation::log_once(const std::string& message, bool force_stdout) const {
    if (this->comm_rank == 0) {
        if (force_stdout) {
            std::cout << message << '\n';
        } else {
            this->outstream << message << '\n';
        }
    }
}

void Simulation::setCompSpecConc(const compartment_concs_t& concentrations,
                                 const math::DistributionMethod distribution) {
    for (const auto& comp_concs: concentrations) {
        setCompSpecConc(comp_concs.first, comp_concs.second, distribution);
    }
}

void Simulation::log_diffusion_exchanges() const {
    const auto& local_exchanges = get_diffusion_rank_exchanges();
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

void Simulation::log_progress(const double i, const double tot, const std::string& name) const {
    std::stringstream s;
    s << name << " progress: " << std::round(1000 * i / tot) / 10 << "%";
    log_once(s.str(), true);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
OmegaHSimulation<SSA, SearchMethod>::OmegaHSimulation(DistMesh& t_mesh,
                                                      rng::RNG& t_rng,
                                                      std::ostream& t_outstream,
                                                      bool t_indepKProcs)
    : super_type(t_mesh, t_rng, t_outstream)
    , mesh(t_mesh)
    , elems2verts(mesh.ask_elem_verts())
    , coords(mesh.coords())
    , indepKProcs(t_indepKProcs) {}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getPatchSpecCount(
    const model::patch_id& patch,
    const model::species_name& species) const {
    const auto& boundaries = mesh.getOwnedEntities(patch);
    osh::Write<osh::LO> mols_counts(boundaries.size());
    const container::species_id spec_id{
        statedef->getPatchdef(patch).getSpecPatchIdx(statedef->getSpecModelIdx(species))};
    const auto& molecules = data->pools.moleculesOnPatchBoundaries();

    std::transform(boundaries.begin(), boundaries.end(), mols_counts.begin(), [&](auto bound) {
        return molecules(bound, spec_id);
    });

    return mesh.get_MPI_sum(osh::LOs(mols_counts));
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setPatchSpecCount(
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

    container::species_id cont_spec_id = statedef->getPatchdef(patch).getSpecPatchIdx(
        statedef->getSpecModelIdx(species));
    for (auto k = 0; k < mols_on_elements.size(); ++k) {
        const mesh::triangle_id_t boundary(elems[k]);
        const auto molecules = mols_on_elements[k];
        this->data->pools.assign(boundary, cont_spec_id, molecules);
    }
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setPatchSpecCount(
    const patch_counts_t& counts,
    const math::DistributionMethod distribution) {
    for (const auto& count: counts) {
        this->setPatchSpecCount(count.patch, count.species, count.num_mols, distribution);
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& species,
    osh::Real num_molecules,
    const math::DistributionMethod distribution) {
    setCompSpecCount(compartment, {{species, num_molecules}}, distribution);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setCompSpecCount(
    const model::compartment_id& compartment,
    const std::vector<CompartmentCount>& counts,
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

    for (auto mol_idx = 0u; mol_idx < counts.size(); ++mol_idx) {
        if (counts[mol_idx].num_mols >= static_cast<double>(INT64_MAX)) {
            std::ostringstream oss;
            oss << "Unsupported number of molecules: " << std::setprecision(20)
                << counts[mol_idx].num_mols << " but maximum value is "
                << std::numeric_limits<osh::GO>::max() << " (max 64 bits integral value)";
            ArgErrLog(oss.str());
        }
        auto dist = math::make_dist(static_cast<osh::GO>(counts[mol_idx].num_mols), rank_volumes);
        // only rank 0 generates non-zero values
        const std::vector<osh::GO>& num_molecules_on_ranks = dist.distribute(this->rng,
                                                                             distribution);


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

        setOwnedCompSpecCount(compartment,
                              counts[mol_idx].species,
                              num_molecules_on_rank,
                              distribution);
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setOwnedCompSpecCount(
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& species,
    osh::Real concentration,
    const math::DistributionMethod distribution) {
    setCompSpecConc(compartment, {{species, concentration}}, distribution);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setCompSpecConc(
    const model::compartment_id& compartment,
    const std::vector<CompartmentConc>& concs,
    const math::DistributionMethod distribution) {
    const auto factor = mesh.total_measure(compartment) * 1.0e3 * math::AVOGADRO;
    std::vector<CompartmentCount> counts;
    counts.reserve(concs.size());
    for (const auto& conc: concs) {
        counts.emplace_back(conc.species, conc.concentration * factor);
    }
    setCompSpecCount(compartment, counts, distribution);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setOwnedElementSpecCount(
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


template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getOwnedCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    const auto spec_count = getOwnedCompSpecCount(compartment, species);
    return spec_count / (1.0e3 * mesh.getMeasure().rank_measure() * math::AVOGADRO);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
const std::vector<mesh::triangle_id_t>& OmegaHSimulation<SSA, SearchMethod>::getGHKBoundaries()
    const {
    return data->kproc_state.ghkCurrentsBoundaries();
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Reals OmegaHSimulation<SSA, SearchMethod>::getGHKCurrents() const {
    return data->kproc_state.ghkSurfaceReactions().currents();
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getTotalGHKCurrent() const {
    return mesh.get_MPI_sum(getGHKCurrents());
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setPatchSReacK(
    const model::patch_id& patchId,
    const model::surface_reaction_id& reactionId,
    osh::Real kCst) {
    auto reacId = statedef->getSReacIdx(reactionId);
    data->kproc_state.surfaceReactions().updateKCst(patchId, reacId, kCst, mesh);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Reals OmegaHSimulation<SSA, SearchMethod>::getPotentialOnVertices(
    const model::patch_id& patch) {
    const auto& ents = mesh.getOwnedEntities(patch);
    osh::Write<osh::Real> vals(3 * ents.size());
    const auto& all_verts = mesh.ask_verts_of(osh::FACE);
    const auto fill_vals = OMEGA_H_LAMBDA(osh::LO entity_idx) {
        const auto& verts = osh::gather_verts<DistMesh::dim()>(all_verts, ents[entity_idx].get());
        for (osh::LO l = 0; l < DistMesh::dim(); ++l) {
            vals[DistMesh::dim() * entity_idx + l] = input->potential_on_vertices_w[verts[l]];
        }
    };
    osh::parallel_for(ents.size(), fill_vals, "OmegaHSimulation::getPotentialOnVertices");

    return osh::Reals(vals);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Reals OmegaHSimulation<SSA, SearchMethod>::getPotentialOnTriangles(
    const model::patch_id& patch) {
    const auto& ents = mesh.getOwnedEntities(patch);
    osh::Write<osh::Real> vals = osh::Write<osh::Real>(ents.size(), 0.0);

    const auto fill_vals = OMEGA_H_LAMBDA(osh::LO entity_idx) {
        const mesh::triangle_id_t triangle_id{ents[entity_idx]};
        const auto& face_bf2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                         triangle_id.get());
        for (const auto& vert_id: face_bf2verts) {
            vals[entity_idx] += input->potential_on_vertices_w[vert_id] / 3.0;
        }
    };
    osh::parallel_for(ents.size(), fill_vals, "OmegaHSimulation::getPotentialOnTriangles");
    return osh::Reals(vals);
}

//-----------------------------------------------

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getCompSpecConc(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    const auto spec_count = getCompSpecCount(compartment, species);
    return spec_count / (1.0e3 * mesh.total_measure(compartment) * math::AVOGADRO);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getOwnedCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& spec_id) const {
    const auto species = statedef->getCompSpecContainerIdx(compartment, spec_id);
    const auto lambda = [=](osh::GO accu, mesh::tetrahedron_id_t elem) -> osh::GO {
        return accu + static_cast<osh::GO>(data->pools(elem, species));
    };
    const auto& elements = mesh.getOwnedEntities(compartment);
    return static_cast<osh::Real>(std::accumulate(elements.begin(), elements.end(), 0, lambda));
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::pair<std::reference_wrapper<const mesh::tetrahedron_ids>, std::vector<osh::LO>>
OmegaHSimulation<SSA, SearchMethod>::getOwnedElemSpecCount(
    const model::species_name& species) const {
    std::vector<osh::LO> counts;
    counts.reserve(mesh.owned_elems().size());
    for (const auto elem: mesh.owned_elems()) {
        const auto spec_model_idx = statedef->getSpecModelIdx(species);
        if (spec_model_idx.unknown()) {
            counts.push_back(0);
            continue;
        }
        const auto compartment_id = this->mesh.getCompartment(elem);
        const auto comp_model_idx = statedef->getCompModelIdx(compartment_id);
        const auto spec_id =
            statedef->compdefs()[static_cast<size_t>(comp_model_idx.get())]->getSpecContainerIdx(
                spec_model_idx);
        counts.push_back(data->pools(elem, spec_id));
    }

    return {mesh.owned_elems(), counts};
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::pair<std::vector<mesh::tetrahedron_global_id_t>, std::vector<osh::LO>>
OmegaHSimulation<SSA, SearchMethod>::getElemSpecCount(const model::species_name& species) const {
    const auto local_ID_and_counts = getOwnedElemSpecCount(species);
    const mesh::tetrahedron_ids& local_IDs = local_ID_and_counts.first;
    const std::vector<osh::LO>& local_counts = local_ID_and_counts.second;

    const int local_counts_size = local_ID_and_counts.second.size();

    std::vector<mesh::tetrahedron_global_id_t> owned_global_ids;
    owned_global_ids.reserve(local_counts_size);

    std::transform(local_IDs.begin(),
                   local_IDs.end(),
                   std::back_inserter(owned_global_ids),
                   [this](const mesh::tetrahedron_id_t& local_id) {
                       return mesh.getGlobalIndex(local_id);
                   });

    std::vector<int> count_sizes;
    if (this->comm_rank == 0) {
        count_sizes.resize(static_cast<size_t>(this->comm_size));
    }

    int err =
        MPI_Gather(&local_counts_size, 1, MPI_INT, count_sizes.data(), 1, MPI_INT, 0, this->comm());

    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    std::vector<int> offsets(count_sizes.size() + 1);
    std::partial_sum(count_sizes.begin(), count_sizes.end(), offsets.begin() + 1);

    std::vector<osh::LO> counts;
    if (this->comm_rank == 0) {
        counts.resize(offsets.back());
    }

    err = MPI_Gatherv(local_counts.data(),
                      local_counts.size(),
                      MPI_INT32_T,
                      counts.data(),
                      count_sizes.data(),
                      offsets.data(),
                      MPI_INT32_T,
                      0,
                      this->comm());

    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    std::vector<mesh::tetrahedron_global_id_t> global_ids;
    if (this->comm_rank == 0) {
        global_ids.resize(offsets.back());
    }

    err = MPI_Gatherv(owned_global_ids.data(),
                      owned_global_ids.size(),
                      MPI_INT64_T,
                      global_ids.data(),
                      count_sizes.data(),
                      offsets.data(),
                      MPI_INT64_T,
                      0,
                      this->comm());

    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }

    return {global_ids, counts};
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchElemValsNP(const osh::GO* indices,
                                                             size_t input_size,
                                                             const model::species_name& species,
                                                             osh::Real* vals,
                                                             bool useConc,
                                                             bool local) const {
    const auto spec_model_idx = statedef->getSpecModelIdx(species);
    if (not local) {
        std::fill(vals, vals + input_size, 0);
    }
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::tetrahedron_global_id_t(indices[i]), true).get();
        }
        mesh::tetrahedron_id_t localInd(ind);
        if (localInd.valid()) {
            // TODO Maybe getting the spec_id could be faster
            const auto compartment_id = mesh.getCompartment(localInd);
            const auto comp_model_idx = statedef->getCompModelIdx(compartment_id);
            const auto spec_id = statedef->compdefs()[static_cast<size_t>(comp_model_idx.get())]
                                     ->getSpecContainerIdx(spec_model_idx);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setBatchElemValsNP(const osh::GO* indices,
                                                             size_t input_size,
                                                             const model::species_name& species,
                                                             osh::Real* vals,
                                                             bool useConc,
                                                             bool local) const {
    const auto spec_model_idx = statedef->getSpecModelIdx(species);

    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::tetrahedron_global_id_t(indices[i]), true).get();
        }
        mesh::tetrahedron_id_t localInd(ind);
        if (localInd.valid()) {
            // TODO Maybe getting the spec_id could be faster
            const auto compartment_id = mesh.getCompartment(localInd);
            const auto comp_model_idx = statedef->getCompModelIdx(compartment_id);
            const auto spec_id = statedef->compdefs()[static_cast<size_t>(comp_model_idx.get())]
                                     ->getSpecContainerIdx(spec_model_idx);

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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchBoundSpecCountNP(
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
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::triangle_global_id_t(indices[i]), true).get();
        }
        mesh::triangle_id_t localInd(ind);
        if (localInd.valid()) {
            const auto patch_id = model::patch_id(
                mesh.getTriPatch(mesh::triangle_id_t(localInd))->getID());
            auto spec_id = statedef->getPatchdef(patch_id).getSpecPatchIdx(spec_model_idx);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setBatchBoundSpecCountNP(
    const osh::GO* indices,
    size_t input_size,
    const model::species_name& species,
    osh::Real* counts,
    bool local) const {
    const auto spec_model_idx = statedef->getSpecModelIdx(species);

    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::triangle_global_id_t(indices[i]), true).get();
        }
        mesh::triangle_id_t localInd(ind);
        if (localInd.valid()) {
            const auto patch_id = model::patch_id(
                mesh.getTriPatch(mesh::triangle_id_t(localInd))->getID());
            auto spec_id = statedef->getPatchdef(patch_id).getSpecPatchIdx(spec_model_idx);
            data->pools.assign(localInd, spec_id, counts[i]);
        }
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchVertVsNP(const osh::GO* indices,
                                                           size_t input_size,
                                                           osh::Real* voltages,
                                                           bool local) const {
    if (not local) {
        std::fill(voltages, voltages + input_size, 0);
    }
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::vertex_global_id_t(indices[i]), true).get();
        }
        mesh::vertex_id_t localInd(ind);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTriVsNP(const osh::GO* indices,
                                                          size_t input_size,
                                                          osh::Real* voltages,
                                                          bool local) const {
    std::fill(voltages, voltages + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::triangle_global_id_t(indices[i]), true).get();
        }
        mesh::triangle_id_t localInd(ind);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTetVsNP(const osh::GO* indices,
                                                          size_t input_size,
                                                          osh::Real* voltages,
                                                          bool local) const {
    std::fill(voltages, voltages + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::tetrahedron_global_id_t(indices[i]), true).get();
        }
        mesh::tetrahedron_id_t localInd(ind);
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

#ifdef USE_PETSC

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTriOhmicIsNP(const osh::GO* indices,
                                                               size_t input_size,
                                                               const model::ohmic_current_id curr,
                                                               osh::Real* currents,
                                                               bool local) const {
    const auto& currs = statedef->ohmicCurrents();
    const auto curr_it = currs.find(curr);
    if (curr_it == currs.end()) {
        throw std::logic_error("No ohmic current : " + curr);
    }
    const auto& h = *curr_it->second;

    std::fill(currents, currents + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::triangle_global_id_t(indices[i]), true).get();
        }
        mesh::triangle_id_t localInd(ind);
        if (localInd.valid()) {
            const auto& face_bf2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                             localInd.get());
            for (const auto& vert_id: face_bf2verts) {
                currents[i] += h.getTriCurrentOnVertex(input->potential_on_vertices_w[vert_id],
                                                       localInd,
                                                       input->pools,
                                                       mesh,
                                                       state_time);
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

#else

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTriOhmicIsNP(
    const osh::GO* /*indices*/,
    size_t /*input_size*/,
    const model::ohmic_current_id /*curr*/,
    osh::Real* /*currents*/,
    bool /*local*/) const {
    throw std::logic_error("PETSc is required to compute ohmic current");
}

#endif  // !USE_PETSC

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTriGHKIsNP(const osh::GO* indices,
                                                             size_t input_size,
                                                             const model::ghk_current_id curr,
                                                             osh::Real* currents,
                                                             bool local) const {
    const auto& surfReacs = data->kproc_state.ghkSurfaceReactions();

    std::fill(currents, currents + input_size, 0);
    for (size_t i = 0; i < input_size; ++i) {
        osh::LO ind{static_cast<osh::LO>(indices[i])};
        if (not local) {
            ind = mesh.getLocalIndex(mesh::triangle_global_id_t(indices[i]), true).get();
        }
        mesh::triangle_id_t localInd(ind);
        if (localInd.valid()) {
            const osh::Write<osh::GO>& tri2Curr = surfReacs.getTri2Curr(curr);
            for (uint k = 0; k < surfReacs.rpt(); ++k) {
                const auto& ridx = tri2Curr[localInd.get() * surfReacs.rpt() + k];
                if (ridx != -1) {
                    currents[i] += surfReacs.currents()[ridx];
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getCompSpecCount(
    const model::compartment_id& compartment,
    const model::species_name& species) const {
    auto local_num_molecules = getOwnedCompSpecCount(compartment, species);
    osh::Real global_num_molecules{};
    auto err = MPI_Allreduce(
        &local_num_molecules, &global_num_molecules, 1, MPI_DOUBLE, MPI_SUM, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }
    return global_num_molecules;
}

#if USE_PETSC
template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::pair<mesh::triangle_ids, osh::Reals> OmegaHSimulation<SSA, SearchMethod>::getOhmicCurrents(
    const model::membrane_id& mem_id,
    const model::channel_id& chan_id) const {
    // mesh::triangle_ids tri_ids does not create a fully functional object
    mesh::triangle_ids tri_ids(0);
    osh::Write<osh::Real> ohm_currs(0);

    auto membrane_ptr = statedef->membranes().find(mem_id);
    if (membrane_ptr != statedef->membranes().end()) {
        const auto* membrane = membrane_ptr->second.get();

        auto channel_ptr = membrane->channels().find(chan_id);
        if (channel_ptr != membrane->channels().end()) {
            const auto& chan = channel_ptr->second;
            const auto& patch_id = membrane->getPatch();
            const auto& patch_tris = mesh.getOwnedEntities(patch_id);

            tri_ids = patch_tris;
            ohm_currs = osh::Write<osh::Real>(patch_tris.size(), 0.0);

            const auto collect_currents = OMEGA_H_LAMBDA(osh::LO patch_tris_idx) {
                const mesh::triangle_id_t triangle_id{patch_tris[patch_tris_idx]};
                const auto& face_bf2verts = osh::gather_verts<3>(mesh.ask_verts_of(osh::FACE),
                                                                 triangle_id.get());

                for (const auto& h: chan.ohmic_currents) {
                    for (const auto& vert_id: face_bf2verts) {
                        ohm_currs[patch_tris_idx] +=
                            h.get().getTriCurrentOnVertex(input->potential_on_vertices_w[vert_id],
                                                          triangle_id,
                                                          input->pools,
                                                          mesh,
                                                          state_time);
                    }
                }
            };
            osh::parallel_for(patch_tris.size(), collect_currents);
        } else {
            throw std::logic_error("No channel: " + chan_id);
        }
    } else {
        throw std::logic_error("No membrane: " + mem_id);
    }

    return {tri_ids, ohm_currs};
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getTotalOhmicCurrent(
    const model::membrane_id& mem_id,
    const model::channel_id& chan_id) const {
    return mesh.get_MPI_sum(getOhmicCurrents(mem_id, chan_id).second);
}
#endif  // USE_PETSC

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setMembPotential(const model::membrane_id& memb,
                                                           osh::Real value) {
    auto membit = mesh.membranes().find(memb);
    if (membit == mesh.membranes().end()) {
        throw std::invalid_argument("Invalid membrane " + memb);
    }
    const auto& allPatches = mesh.getAllPatches();
    for (const auto& patchid: membit->second->patches()) {
        auto pmeshid = mesh.getPatchID(patchid);
        const auto* patch = allPatches[pmeshid.get()];
        const auto* icomp = dynamic_cast<const DistComp*>(&patch->getIComp());
        if (icomp == nullptr) {
            continue;
        }
        for (const auto& tet: icomp->getLocalTetIndices(false)) {
            const auto verts = osh::gather_verts<4>(mesh.ask_elem_verts(), tet.get());
            for (const auto& vert: verts) {
                input->potential_on_vertices_w[vert] = value;
            }
        }
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setMembRes(const model::membrane_id& membrane,
                                                     osh::Real resistivity,
                                                     osh::Real reversal_potential) {
    statedef->setResistivity(membrane, resistivity);
    statedef->setReversalPotential(membrane, reversal_potential);
    auto it = statedef->membranes().find(membrane);
    // Initialize triangle resistivity on membranes
    if (it != statedef->membranes().end()) {
        for (const auto tri: mesh.getOwnedEntities(it->second->getPatch())) {
            input->conductivity_on_triangles_w[tri.get()] = 1.0 / resistivity;
            input->reversal_potential_on_triangles_w[tri.get()] = reversal_potential;
        }
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
MembraneResistivity OmegaHSimulation<SSA, SearchMethod>::getMembRes(
    const model::membrane_id& membrane) {
    return {statedef->getResistivity(membrane), statedef->getReversalPotential(membrane)};
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
MembraneResistivity OmegaHSimulation<SSA, SearchMethod>::getTriRes(osh::GO tri, bool local) const {
    osh::Real local_res(0.0);
    osh::Real local_erev(0.0);

    auto ind{static_cast<osh::LO>(tri)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::triangle_global_id_t(tri), true).get();
    }
    mesh::triangle_local_id_t localInd(ind);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setTriRes(const osh::GO tri,
                                                    osh::Real res,
                                                    osh::Real erev,
                                                    bool local) {
    osh::LO ind{static_cast<osh::LO>(tri)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::triangle_global_id_t(tri), false).get();
    }
    mesh::vertex_id_t localInd(ind);
    if (localInd.valid()) {
        input->conductivity_on_triangles_w[localInd.get()] = 1.0 / res;
        input->reversal_potential_on_triangles_w[localInd.get()] = erev;
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getTriCapac(osh::GO tri, bool local) const {
    osh::Real local_val(0.0);

    osh::LO ind{static_cast<osh::LO>(tri)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::triangle_global_id_t(tri), true).get();
    }
    mesh::triangle_local_id_t localInd(ind);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setTriCapac(const osh::GO tri, osh::Real c, bool local) {
    osh::LO ind{static_cast<osh::LO>(tri)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::triangle_global_id_t(tri), false).get();
    }
    mesh::vertex_id_t localInd(ind);
    if (localInd.valid()) {
        input->capacitance_on_triangles_w[localInd.get()] = c;
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getVertIClamp(const osh::GO vertex,
                                                             bool local) const {
    osh::Real local_val(0.0);

    osh::LO ind{static_cast<osh::LO>(vertex)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::vertex_global_id_t(vertex), true).get();
    }
    mesh::vertex_id_t localInd(ind);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setVertIClamp(const osh::GO vertex,
                                                        const osh::Real current,
                                                        bool local) {
    osh::LO ind{static_cast<osh::LO>(vertex)};
    if (not local) {
        ind = mesh.getLocalIndex(mesh::vertex_global_id_t(vertex), false).get();
    }
    mesh::vertex_id_t localInd(ind);
    if (localInd.valid()) {
        input->current_on_vertices_w[localInd.get()] = current;
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::getBatchTriOhmicErevsNP(
    const gsl::span<const osh::GO>& triangles,
    const model::ohmic_current_id& ohmic_current,
    const gsl::span<double>& erev,
    bool local) const {
    const auto& oc = statedef->ohmicCurrents().at(ohmic_current);
    for (size_t i = 0; i < triangles.size(); ++i) {
        mesh::triangle_local_id_t triangle(static_cast<osh::LO>(triangles[i]));
        if (!local) {
            triangle = mesh.getLocalIndex(mesh::triangle_global_id_t(triangles[i]), true);
        }
        if (triangle.valid()) {
            erev[i] = oc->getReversalPotential(triangle);
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
double OmegaHSimulation<SSA, SearchMethod>::getTriOhmicErev(
    osh::GO triangle,
    const model::ohmic_current_id& ohmic_current,
    bool local) const {
    mesh::triangle_id_t local_index(static_cast<osh::LO>(triangle));
    if (not local) {
        local_index = mesh.getLocalIndex(mesh::triangle_global_id_t(triangle), true);
    }
    double rp{};
    if (local_index.valid()) {
        const auto& oc = statedef->ohmicCurrents().at(ohmic_current);
        rp = oc->getReversalPotential(local_index);
    }

    if (local) {
        return rp;
    } else {
        double result{};
        auto err = MPI_Allreduce(&rp, &result, 1, MPI_DOUBLE, MPI_SUM, this->comm());
        if (err != MPI_SUCCESS) {
            MPI_Abort(this->comm(), err);
        }
        return result;
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setTriOhmicErev(
    osh::GO triangle,
    const model::ohmic_current_id& ohmic_current,
    double reversal_potential,
    bool local) {
    mesh::triangle_id_t local_index(static_cast<osh::LO>(triangle));
    if (not local) {
        local_index = mesh.getLocalIndex(mesh::triangle_global_id_t(triangle), true);
    }
    if (local_index.valid()) {
        auto& oc = statedef->ohmicCurrents().at(ohmic_current);
        oc->setReversalPotential(local_index, reversal_potential);
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setDiffOpBinomialThreshold(osh::Real threshold) {
    data->diffOp.setBinomialThreshold(static_cast<osh::GO>(threshold));
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::Real OmegaHSimulation<SSA, SearchMethod>::getIterationTimeStep() const noexcept {
    return data->time_delta;
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::exportMolStateToVTK(const std::string& /* filename */) {
    // TODO(TCL) FIXME
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod>::getDiffOpExtent(bool local) const {
    osh::I64 extent = data->diffOp.getExtent();
    if (local) {
        return extent;
    }

    osh::I64 global_ext{};
    auto err = MPI_Allreduce(&extent, &global_ext, 1, MPI_INT64_T, MPI_SUM, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }
    return global_ext;
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod>::getSSAOpExtent(bool local) const {
    const osh::I64 extent = data->ssaOp.getExtent();
    if (local) {
        return extent;
    }

    osh::I64 global_extent{};
    auto err = MPI_Allreduce(&extent, &global_extent, 1, MPI_INT64_T, MPI_SUM, this->comm());
    if (err != MPI_SUCCESS) {
        MPI_Abort(this->comm(), err);
    }
    return global_extent;
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
osh::I64 OmegaHSimulation<SSA, SearchMethod>::getNIterations() const noexcept {
    return num_iterations;
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
std::string OmegaHSimulation<SSA, SearchMethod>::createStateReport() {
    // TODO(TCL) FIXME
    return "";
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::compute_num_species_per_elements(
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

    if (!statedef.patchdefs().empty()) {
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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::init(std::unique_ptr<Statedef>&& t_statedef) {
    this->statedef.swap(t_statedef);
    assert(statedef != nullptr);
    this->mesh.init();

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
                                                    0 /*Unused in context*/,
                                                    this->rng,
                                                    mesh.owned_verts_mask().size());

    // Initialize triangle capacitance on membranes
    for (auto& memb: statedef->membranes()) {
        auto capac = memb.second->capacitance();
        for (const auto tri: mesh.getOwnedEntities(memb.second->getPatch())) {
            input->capacitance_on_triangles_w[tri.get()] = capac;
        }
    }

    data = std::make_unique<SimulationData<SSA, SearchMethod>>(
        mesh, *this->statedef, *input, this->rng, indepKProcs);
    initialize_discretized_rates();
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::reset() {
    this->num_iterations = 0;
    this->state_time = 0;
    setPotential(DEFAULT_MEMB_POT);
    data->reset(this->state_time);
    {
        // reset ohmic currents manually set with method \a setTriOhmicErev
        for (auto& [_, current]: statedef->ohmicCurrents()) {
            current->reset();
        }
    }
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::initialize_discretized_rates() {
    const auto& measure_info = mesh.getMeasure();
    data->diffusions.reset();
    osh::parallel_for(
        this->mesh.owned_elems().size(),
        [&measure_info, this](osh::LO e) __attribute__((always_inline, flatten)) {
            const mesh::tetrahedron_id_t elem(mesh.owned_elems()[e]);
            const auto elem_measure = measure_info.element_measure(elem);
            const auto compartment_id = this->mesh.getCompartment(elem);
            const auto compartment_mid = this->mesh.getCompartmentMeshID(elem);
            const auto& compartment = this->statedef->getCompdef(compartment_id);
            for (const auto& diffusion: compartment.diffdefs()) {
                const auto species = diffusion->getSpecContainerIdx();
                const auto dcst = diffusion->getDcst();
                data->diffusions.rates_sum(elem, species) = 0;
                const auto num_neighbors = mesh.tet_neighbors_int_data().size(elem.get());
                for (auto face_idx = 0; face_idx < num_neighbors; ++face_idx) {
                    auto d = mesh.tet_neighbors_int_data()(elem.get(), face_idx);

                    for (auto& db: mesh.diffusionBoundaries()) {
                        if (static_cast<size_t>(species.get()) >=
                            db.comp2_diffusing_species.size()) {
                            db.comp2_diffusing_species.resize(species.get() + 1, false);
                        }
                        if (static_cast<size_t>(species.get()) >=
                            db.comp1_diffusing_species.size()) {
                            db.comp1_diffusing_species.resize(species.get() + 1, false);
                        }
                    }

                    if (mesh.getCompartmentMeshID(mesh::tetrahedron_id_t(d[0])) ==
                            compartment_mid ||
                        mesh.isActiveDiffusionBoundary(mesh::triangle_id_t(d[2]),
                                                       compartment_mid,
                                                       species)) {
                        const auto neighbor_boundary_distance =
                            mesh.tet_neighbors_real_data()(elem.get(), face_idx)[0];
                        const auto neighbor_boundary_measure =
                            mesh.tet_neighbors_real_data()(elem.get(), face_idx)[1];
                        const auto propensity = dcst * neighbor_boundary_measure / elem_measure /
                                                neighbor_boundary_distance;
                        data->diffusions.ith_rate(elem, species, face_idx) = propensity;
                        data->diffusions.rates_sum(elem, species) += propensity;
                    }
                }
            }
        });
    data->updateIterationTimeStep();
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::evolve_rd(const osh::Real rd_dt) {
    util::TimeTracker t;
    t.start();

    data->pools.reset_occupancy_rd(state_time);
    data->ssaOp.run(rd_dt, state_time);
    t.stop();
    this->reactions_timer += t.diff();

    if (data->active_diffusions) {
        t.start();
        data->diffOp(rd_dt, state_time);
        t.stop();
        this->diffusions_timer += t.diff();
    }
    state_time += rd_dt;
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::run_rd(const osh::Real end_time) {
    data->ssaOp.resetAndUpdateAll(state_time, end_time);
    // number of standard steps
    const auto rd_dt_std = data->time_delta;
    // this can be negative
    const int n_steps_std = std::floor((end_time - state_time) / rd_dt_std) - 1;
    // std steps loop -1. n_steps_std can be negative so int is the correct type
    for (int i_rd = 0; i_rd < n_steps_std; ++i_rd) {
        evolve_rd(std::min(rd_dt_std, end_time - state_time));
    }

    // sync time steps
    // we compare always with end_time because comparing dts can fail due to numerical error
    const auto next_std_state_time = state_time + rd_dt_std;
    if (!steps::util::almost_equal(next_std_state_time, end_time) &&
        next_std_state_time < end_time) {
        evolve_rd(rd_dt_std);
        assert(end_time > state_time);
        assert(!steps::util::almost_equal(end_time, state_time));
        evolve_rd(end_time - state_time);
    } else if (!steps::util::almost_equal(end_time, state_time)) {
        evolve_rd(end_time - state_time);
    }

    assert(steps::util::almost_equal(end_time, state_time));
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::evolve(const osh::Real ef_dt) {
    data->pools.reset_occupancy_ef(state_time);

    ++this->num_iterations;
    data->kproc_state.resetCurrents();
    osh::Reals potential_on_vertices(input->potential_on_vertices_w);
    data->kproc_state.updateVDepSReacs(potential_on_vertices);

    run_rd(state_time + ef_dt);

    // We divide the charge_flows in currents_ by ef_dt so we really get the currents
    data->kproc_state.ghkSurfaceReactions().finalizeCurrents(ef_dt);
#if USE_PETSC
    if (data->efield) {
        util::TimeTracker t;
        t.start();
        data->efield->evolve(input->potential_on_vertices_w,
                             input->current_on_vertices_w,
                             input->capacitance_on_triangles_w,
                             input->conductivity_on_triangles_w,
                             input->reversal_potential_on_triangles_w,
                             input->pools,
                             data->kproc_state.ghkSurfaceReactions().currents(),
                             state_time,
                             ef_dt);
        t.stop();
        this->efield_timer += t.diff();
    }
#endif  // USE_PETSC
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::run(osh::Real end_time) {
    Instrumentor::phase p("OmegaHSimulation::run()");

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

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setDiffusionBoundaryActive(
    const mesh::diffusion_boundary_name& diffusion_boundary_name,
    const model::species_name& spec_id,
    bool set_active) {
    auto diffusion_boundary_id = mesh.getDiffusionBoundaryIndex(diffusion_boundary_name);
    model::species_id mdl_spec_id = statedef->getSpecModelIdx(spec_id);
    if (diffusion_boundary_id >= mesh.diffusionBoundaries().size()) {
        throw std::invalid_argument("Invalid diffusion boundary " +
                                    std::to_string(diffusion_boundary_id));
    }
    DistMesh::DiffusionBoundary& db = mesh.diffusionBoundaries()[diffusion_boundary_id];
    Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
    Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
    container::species_id sp1 = comp1.getSpecContainerIdx(mdl_spec_id);
    container::species_id sp2 = comp2.getSpecContainerIdx(mdl_spec_id);
    auto& comp1_specs = db.comp1_diffusing_species;
    auto& comp2_specs = db.comp2_diffusing_species;
    if (comp1_specs.size() <= static_cast<size_t>(sp1.get())) {
        comp1_specs.resize(sp1.get() + 1, false);
    }
    if (comp2_specs.size() <= static_cast<size_t>(sp2.get())) {
        comp2_specs.resize(sp2.get() + 1, false);
    }
    if (db.conv_12.size() <= static_cast<size_t>(sp1.get())) {
        db.conv_12.resize(sp1.get() + 1);
    }
    if (db.conv_21.size() <= static_cast<size_t>(sp2.get())) {
        db.conv_21.resize(sp2.get() + 1);
    }
    db.conv_12[sp1.get()] = sp2;
    db.conv_21[sp2.get()] = sp1;
    comp1_specs[sp1.get()] = set_active;
    comp2_specs[sp2.get()] = set_active;
    initialize_discretized_rates();
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
bool OmegaHSimulation<SSA, SearchMethod>::getDiffusionBoundaryActive(
    const mesh::diffusion_boundary_name& diffusion_boundary_name,
    const model::species_name& spec_id) {
    auto diffusion_boundary_id = mesh.getDiffusionBoundaryIndex(diffusion_boundary_name);
    model::species_id mdl_spec_id = statedef->getSpecModelIdx(spec_id);
    if (diffusion_boundary_id >= mesh.diffusionBoundaries().size()) {
        throw std::invalid_argument("Invalid diffusion boundary " +
                                    std::to_string(diffusion_boundary_id));
    }
    DistMesh::DiffusionBoundary& db = mesh.diffusionBoundaries()[diffusion_boundary_id];
    Compdef& comp1 = statedef->getCompdef(db.mdl_comp1);
    Compdef& comp2 = statedef->getCompdef(db.mdl_comp2);
    container::species_id sp1 = comp1.getSpecContainerIdx(mdl_spec_id);
    container::species_id sp2 = comp2.getSpecContainerIdx(mdl_spec_id);
    auto& comp1_specs = db.comp1_diffusing_species;
    auto& comp2_specs = db.comp2_diffusing_species;
    if (comp1_specs.size() <= static_cast<size_t>(sp1.get())) {
        comp1_specs.resize(sp1.get() + 1, false);
    }
    if (comp2_specs.size() <= static_cast<size_t>(sp2.get())) {
        comp2_specs.resize(sp2.get() + 1, false);
    }
    if (db.conv_12.size() <= static_cast<size_t>(sp1.get())) {
        db.conv_12.resize(sp1.get() + 1);
    }
    if (db.conv_21.size() <= static_cast<size_t>(sp2.get())) {
        db.conv_21.resize(sp2.get() + 1);
    }
    db.conv_12[sp1.get()] = sp2;
    db.conv_21[sp2.get()] = sp1;
    return comp1_specs[sp1.get()] && comp2_specs[sp2.get()];
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::setMembIClamp(const model::membrane_id& membrane,
                                                        osh::Real current) {
    statedef->setStimulus(membrane, current);
}

template <SSAMethod SSA, NextEventSearchMethod SearchMethod>
void OmegaHSimulation<SSA, SearchMethod>::dumpDepGraphToFile(const std::string& path) const {
    std::ofstream ostr(path);
    data->kproc_state.write_dependency_graph(ostr);
}

// explicit template instantiation definitions

template class OmegaHSimulation<SSAMethod::SSA, NextEventSearchMethod::GibsonBruck>;
template class OmegaHSimulation<SSAMethod::SSA, NextEventSearchMethod::Direct>;
template class OmegaHSimulation<SSAMethod::RSSA, NextEventSearchMethod::Direct>;

}  // namespace steps::dist

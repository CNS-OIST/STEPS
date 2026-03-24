#include "diffusion_operator.hpp"

#include <Omega_h_defines.hpp>
#include <algorithm>
#include <limits>
#include <numeric>
#include <random>

#include <Omega_h_for.hpp>
#include <stdexcept>

#include "geom/dist/distmesh.hpp"
#include "rng/rng.hpp"
#include "util/flat_multimap.hpp"
#include "util/profile/profiler_interface.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

DiffusionOperator::DiffusionOperator(DistMesh& t_mesh,
                                     rng::RNG& t_rng,
                                     MolState& t_pools,
                                     kproc::Diffusions& t_diffusions,
                                     kproc::KProcState& t_kproc_state)
    : mesh(t_mesh)
    , rng(t_rng)
    , pools(t_pools)
    , diffusions_(t_diffusions)
    , kproc_state_(t_kproc_state)
    , ur_distribution(0, 1) {}

void DiffusionOperator::initialize() {
    auto global_max_sums = diffusions_.global_rates_max_sum();
    active_diffusions = (global_max_sums > std::numeric_limits<osh::Real>::epsilon());

    // if no diffusions are present. time delta can be set to infinity.
    time_delta = active_diffusions ? 1.0 / global_max_sums
                                   : std::numeric_limits<osh::Real>::infinity();
}

void DiffusionOperator::operator()(const osh::Real opsplit_period, const osh::Real state_time) {
    Instrumentor::phase p("DiffusionOperator::operator()");

    diffusions_.leaving_molecules_reset();
    species_leaving_elements(opsplit_period, state_time);

    num_diffusions_ += diffusions_.leaving_molecules().sum();

    Instrumentor::phase_begin("sync_delta_pools()");
    diffusions_.leaving_molecules().sync_values();
    Instrumentor::phase_end("sync_delta_pools()");

    species_entering_elements();
}

std::map<std::string, double> DiffusionOperator::getDebugInfo() const {
    std::map<std::string, double> info;
    info["total_diff_steps"] = num_diffusions_;
    return info;
}

void DiffusionOperator::reset() {
    num_diffusions_ = 0;
}

void DiffusionOperator::species_leaving_elements(const osh::Real opsplit_period,
                                                 const osh::Real state_time) {
    Instrumentor::phase p("DiffusionOperator::species_leaving_elements()");

    const auto lambda = [this, opsplit_period, state_time](osh::LO elemIdx)
        __attribute__((always_inline, flatten)) {
        const auto element = mesh.owned_elems()[elemIdx];
        for (auto species: pools.species(element)) {
            const auto num_molecules = pools(element, species);
            if (num_molecules > 0) {
                species_leaving_element(
                    element, species, num_molecules, opsplit_period, state_time);
            }
        }
    };
    osh::parallel_for(mesh.owned_elems().size(), lambda, "species_leaving_elements");
}

void DiffusionOperator::species_leaving_element(mesh::tetrahedron_id_t element,
                                                container::species_id species,
                                                molecules_t num_molecules,
                                                const osh::Real opsplit_period,
                                                const osh::Real state_time) {
    const osh::Real rates_sum = diffusions_.rates_sum(element, species);
    const auto delta_pool_total = this->get_leaving_molecules(
        element, species, num_molecules, rates_sum, opsplit_period, state_time);
    if (delta_pool_total > 0) {
        if (delta_pool_total <= diffusion_threshold_) {
            species_leaving_element_standard(element, species, delta_pool_total, rates_sum);
        } else {
            species_leaving_element_binomial(element, species, delta_pool_total, rates_sum);
        }
    }
}


// Possible improvement:
// the STL function "std::discrete_distribution" generates a categorical random variables.
// This function is more efficient when the number of categories (i.e. the size of rates vector) is
// big. For a small vector the current implementation is slightly faster.

void DiffusionOperator::species_leaving_element_standard(mesh::tetrahedron_id_t element,
                                                         container::species_id species,
                                                         molecules_t delta_pool_total,
                                                         osh::Real scaled_dcst) {
    while (delta_pool_total > 0) {
        const auto selector = ur_distribution(rng) * scaled_dcst;
        osh::Real partial_sum_scaled_dcst{0};
        const auto& rates = diffusions_.rates().rates(element, species);
        const auto num_rates = static_cast<osh::LO>(rates.size());
        for (auto direction = 0; direction < num_rates; ++direction) {
            auto current_scaled_dcst = rates[static_cast<size_t>(direction)];
            partial_sum_scaled_dcst += current_scaled_dcst;
            if (selector < partial_sum_scaled_dcst) {
                diffusions_.leaving_molecules().value(element, species, direction) += 1;
                delta_pool_total -= 1;
                break;
            }
        }
    }
}

void DiffusionOperator::species_leaving_element_binomial(mesh::tetrahedron_id_t element,
                                                         container::species_id species,
                                                         molecules_t delta_pool_total,
                                                         osh::Real scaled_dcst) {
    const auto& rates = diffusions_.rates().rates(element, species);
    for (auto e = 0; e < static_cast<osh::LO>(rates.size()); ++e) {  // loop over boundary/faces
        const auto probability_e = rates[static_cast<size_t>(e)] / scaled_dcst;
        auto delta_pool_e = delta_pool_total;
        if (probability_e < 1) {
            delta_pool_e = static_cast<molecules_t>(
                this->rng.getBinom(static_cast<uint>(delta_pool_total),
                                   static_cast<double>(probability_e)));
            delta_pool_total -= delta_pool_e;
            scaled_dcst -= rates[static_cast<size_t>(e)];
        }
        diffusions_.leaving_molecules().value(element, species, e) += delta_pool_e;
    }
}

/**
 * Functor to update molecules pools according to diffusion increment vector
 */

struct SpeciesEnteringElement {
    SpeciesEnteringElement(const DistMesh& t_mesh,
                           MolState& t_pools,
                           kproc::PoolsIncrements& t_leaving_molecules,
                           kproc::KProcState& t_kproc_state)
        : mesh(t_mesh)
        , pools(t_pools)
        , leaving_molecules(t_leaving_molecules)
        , kproc_state_(t_kproc_state) {}

    inline void operator()(mesh::tetrahedron_id_t elem) const {
        const auto& dep_map = kproc_state_.get_dependency_map_elems();
        auto num_elements = mesh.tet_neighbors_int_data().size(elem.get());
        auto elem_comp_id = mesh.getRegionMeshID(elem);
        for (auto face_id = 0; face_id < num_elements; face_id++) {
            const auto& neighbor = mesh.tet_neighbors_int_data()(elem.get(), face_id);
            const auto neighbor_id = mesh::tetrahedron_id_t(neighbor[0]);
            auto comp_id = mesh.getRegionMeshID(neighbor_id);
            const auto neighbour_face_id = neighbor[1];
            const auto triangle_id = mesh::triangle_id_t(neighbor[2]);
            if (comp_id == elem_comp_id) {
                // Careful, pools.species is only defined for an owned element.
                for (osh::LO sp = 0; sp < leaving_molecules.size(neighbor_id.get()); ++sp) {
                    container::species_id species(sp);
                    auto num_entering = leaving_molecules.synced_value(
                        mesh::tetrahedron_id_t(neighbor_id), species, neighbour_face_id);
                    // Remove the leaving molecules from the current tetrahedron
                    auto num_leaving = leaving_molecules.synced_value(elem, species, face_id);
                    auto change = num_entering - num_leaving;
                    // no need to update occupancy here because channels cannot diffuse and rd
                    // occupancy should not track diffusion changes
                    if (change != 0) {
                        pools.add<false>(elem, species, change);
                        auto ab = pools.moleculesOnElements().ab(elem, species);
                        kproc_state_.add_outdated_kprocs(dep_map[ab]);
                    }
                }
            } else {
                for (osh::LO sp = 0; sp < leaving_molecules.size(neighbor_id.get()); ++sp) {
                    container::species_id species(sp);
                    auto num_entering = leaving_molecules.synced_value(
                        mesh::tetrahedron_id_t(neighbor_id), species, neighbour_face_id);
                    if (mesh.getDiffusionBoundaryDcst(triangle_id, comp_id, species) != 0.0) {
                        auto sID = mesh.convertSpeciesID(triangle_id, comp_id, species);
                        // Remove the leaving molecules from the current tetrahedron
                        auto num_leaving = leaving_molecules.synced_value(elem, sID, face_id);
                        auto change = num_entering - num_leaving;
                        // no need to update occupancy here because channels cannot diffuse and
                        // rd occupancy should not track diffusion changes
                        if (change != 0) {
                            pools.add<true>(elem, sID, change);
                            auto ab = pools.moleculesOnElements().ab(elem, sID);
                            kproc_state_.add_outdated_kprocs(dep_map[ab]);
                        }
                    }
                }
            }
        }
    }

  private:
    const DistMesh& mesh;
    MolState& pools;
    const kproc::PoolsIncrements& leaving_molecules;
    kproc::KProcState& kproc_state_;
};

void DiffusionOperator::species_entering_elements() {
    Instrumentor::phase p("DiffusionOperator::species_entering_elements()");

    const SpeciesEnteringElement func_per_element(mesh,
                                                  pools,
                                                  diffusions_.leaving_molecules(),
                                                  kproc_state_);

    const auto lambda = [this, &func_per_element](osh::LO ownedElemIdx)
        __attribute__((always_inline, flatten)) {
        const auto elem = mesh.owned_elems()[ownedElemIdx];
        func_per_element(mesh::tetrahedron_id_t(elem));
    };

    osh::parallel_for(mesh.owned_elems().size(), lambda, "species_entering_elements");
}

molecules_t DiffusionOperator::get_leaving_molecules(mesh::tetrahedron_id_t elem,
                                                     container::species_id species,
                                                     molecules_t num_molecules,
                                                     osh::Real sum_rates,
                                                     const osh::Real opsplit_period,
                                                     const osh::Real state_time) {
    // if probability is 0 we do not diffuse
    if (sum_rates == 0) {
        return 0;
    }

    // diffusions happen at the end of the rd_dt time step. Event_time = state_time + opsplit_period
    const auto occupancy = pools.get_occupancy_rd(elem, species, state_time + opsplit_period);

    const auto mean_population = rng.stochastic_round<molecules_t>(occupancy, num_molecules);
    return diffusions_.total_leaving()(mean_population, sum_rates, opsplit_period);
}

//-----------------------------------------------

TauLeapingDiffusionOperator::TauLeapingDiffusionOperator(DistMesh& t_mesh,
                                                         rng::RNG& t_rng,
                                                         MolState& t_pools,
                                                         kproc::Diffusions& t_diffusions,
                                                         kproc::KProcState& t_kproc_state)
    : mesh(t_mesh)
    , rng(t_rng)
    , pools(t_pools)
    , diffusions_(t_diffusions)
    , kproc_state_(t_kproc_state)
    , out_fluxes(t_mesh, diffusions_.leaving_molecules().species_per_elements())
    , tet_outbound(out_fluxes.species_per_elements(), 0)
    , tri2face(neighbs_per_triangle(t_mesh))
    , border_tris_info(border_triangles_info(t_mesh))
    , border_tris(border_triangles_ids(border_tris_info))
    , border_tets(border_tetrahedrons(border_tris_info))
    , border_tet_pop(t_mesh,
                     diff_species_per_border_tetrahedron(border_tets,
                                                         out_fluxes.species_per_elements()),
                     border_tets)
    , border_tri_diff(t_mesh,
                      diff_species_per_border_triangle(border_tris_info,
                                                       out_fluxes.species_per_elements()),
                      border_tris)
    , internal_tris(internal_triangles(t_mesh)) {
    // Fill tri2face
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());
    for (auto tri: mesh.getTriInfo().range()) {
        for (auto idx = 0; idx < tri2face.size(tri.get()); ++idx) {
            mesh::tetrahedron_local_id_t tet(tri2tets_ab2b[tri2tets_a2ab[tri.get()] + idx]);

            auto num_elements = mesh.tet_neighbors_int_data().size(tet.get());
            for (auto face_id = 0; face_id < num_elements; ++face_id) {
                const auto& neighbor = mesh.tet_neighbors_int_data()(tet.get(), face_id);
                if (neighbor[2] == tri) {
                    tri2face(tri.get(), idx)[0] = face_id;
                    tri2face(tri.get(), idx)[1] = mesh.getRegionMeshID(tet).get();
                    break;
                }
            }
        }
    }

    // Precompute the border faces (owned or not) for all border tetrahedrons
    // The value is kept as a binary mask composed of 4 bits. The LSB corresponds to the face with
    // index 0 of the tetrahedron. If the bit is 1, the corresponding face is a border face.
    osh::Write<osh::LO> owned_border_faces(border_tets.size(), 0);
    for (const auto& info: border_tris_info) {
        if (info.owned) {
            owned_border_faces[info.itet0.get()] |= 1 << info.face0;
            owned_border_faces[info.itet1.get()] |= 1 << info.face1;
        }
    }

    osh::Write<osh::LO> border_tets_w(border_tets.size());
    osh::parallel_for(border_tets.size(), [&](osh::LO i) {
        border_tets_w[i] = border_tets[mesh::tetrahedron_internal_id_t(i)].get();
    });
    auto dist = mesh.create_subset_dist(mesh.dim(), border_tets_w);
    // In that case, since the masks in each rank are non-overlaping, reducing with SUM is
    // equivalent to reducing with OR, but the latter is not available in Omega_h
    auto sum_border_faces =
        dist.exch_reduce(osh::Read<osh::LO>(owned_border_faces), 1, OMEGA_H_SUM);
    auto total_border_faces =
        mesh.sync_subset_array(mesh.dim(), sum_border_faces, border_tets_w, 0, 1);

    auto bitcount = [](unsigned short v) {
        int nb = 0;
        while (v > 0) {
            nb += v & 1;
            v >>= 1;
        }
        return nb;
    };
    border_tet_info.container().resize(border_tets.size());
    osh::parallel_for(border_tets.size(), [&](osh::LO i) {
        mesh::tetrahedron_internal_id_t itet(i);

        border_tet_info[itet].owned_border_faces = bitcount(owned_border_faces[itet.get()]);
        border_tet_info[itet].total_border_faces = bitcount(total_border_faces[itet.get()]);
        border_tet_info[itet].border_faces = total_border_faces[itet.get()];
    });
}

void TauLeapingDiffusionOperator::initialize() {
    auto max_sums = diffusions_.global_rates_max_sum();
    auto mean_sums = diffusions_.global_rates_mean_sum();
    auto min_sums = diffusions_.global_rates_min_sum();
    active_diffusions = (max_sums > std::numeric_limits<osh::Real>::epsilon());
    osh::Real diff_propensity = 0.0;
    if (min_dt_factor >= 0 and min_dt_factor <= 0.5) {
        diff_propensity = min_sums + (mean_sums - min_sums) * 2.0 * min_dt_factor;
    } else if (min_dt_factor <= 1.0) {
        diff_propensity = mean_sums + (max_sums - mean_sums) * 2.0 * (min_dt_factor - 0.5);
    } else {
        throw std::invalid_argument("MinDtFactor should be >= 0 and <= 1.");
    }
    const_dt = active_diffusions ? 1.0 / diff_propensity
                                 : std::numeric_limits<osh::Real>::infinity();
}

void TauLeapingDiffusionOperator::reset() {
    num_diffusions_ = 0;
    num_steps = 0;
    relative_step_sizes.clear();
}

void TauLeapingDiffusionOperator::operator()(const osh::Real opsplit_period,
                                             const osh::Real state_time) {
    Instrumentor::phase p("DiffusionOperator::operator()");

    ++num_steps;
    relative_step_sizes.push_back(opsplit_period / const_dt);

    sync_out_fluxes(state_time + opsplit_period);
    compute_available_border_pop();

    tet_outbound.assign(0.0);

    // Border diffusions
    border_tri_diff.reset();
    sample_triangle_fluxes<true>(opsplit_period);
    border_tri_diff.sync_values();
    apply_border_diffusion();

    // Internal diffusions
    sample_triangle_fluxes<false>(opsplit_period);
}

void TauLeapingDiffusionOperator::sync_out_fluxes(osh::Real end_time) {
    Instrumentor::phase p("TauLeapingDiffusionOperator::sync_out_fluxes()");

    out_fluxes.reset();

    for (auto elem: mesh.owned_elems()) {
        for (auto spec: pools.species(elem)) {
            auto num_elements = mesh.tet_neighbors_int_data().size(elem.get());
            for (auto face_id = 0; face_id < num_elements; face_id++) {
                osh::Real out = diffusions_.ith_rate(elem, spec, face_id) *
                                pools.get_occupancy_rd(elem, spec, end_time);
                out_fluxes.value(elem, spec, face_id) += out;
            }
        }
    }

    out_fluxes.sync_values();
}

void TauLeapingDiffusionOperator::compute_available_border_pop() {
    static int internal_counter = 0;
    ++internal_counter;

    const auto tet2tris = mesh.ask_down(mesh.dim(), mesh.dim() - 1);

    // First compute available populations on border tetrahedrons
    border_tet_pop.assign(0);
    for (const mesh::tetrahedron_internal_id_t itet: border_tets.range()) {
        mesh::tetrahedron_local_id_t tet(border_tets[itet]);
        if (mesh.isOwned(tet)) {
            for (auto spec: pools.species(tet)) {
                // Cannot use occupancy here because we care about actually available pop
                border_tet_pop.value(itet, spec) = pools(tet, spec);
            }
        }
    }
    border_tet_pop.sync_values();

    // Allocate available population depending on how many owned border faces each border
    // tetrahedron has.
    for (const mesh::tetrahedron_internal_id_t itet: border_tets.range()) {
        mesh::tetrahedron_local_id_t tet(border_tets[itet]);
        const auto tris = osh::gather_down<DistMesh::dim() + 1>(tet2tris.ab2b, tet.get());
        for (auto spec: border_tet_pop.species(itet)) {
            const auto& info = border_tet_info[itet];
            // Available population for each border face of the tetrahedron
            auto avg = border_tet_pop.synced_value(itet, spec) / info.total_border_faces;
            // Remaning population that is not assigned to any triangle yet
            auto rem = border_tet_pop.synced_value(itet, spec) % info.total_border_faces;

            border_tet_pop.value(itet, spec) = avg * info.owned_border_faces;

            // Assign the remaining population to triangles owned in that rank
            auto face_id = internal_counter % (mesh.dim() + 1);
            // Since internal_counter is the same for all ranks, the faces are processed in the same
            // order for all ranks
            while (rem > 0) {
                // If the face is a border face
                if ((info.border_faces & (1 << face_id)) != 0) {
                    --rem;
                    // TODO: Improvement, maybe precompute which faces are owned by the current
                    // rank?
                    mesh::triangle_local_id_t tri(tris[face_id]);
                    if (mesh.isOwned(tri)) {
                        border_tet_pop.value(itet, spec)++;
                    }
                }
                face_id = (face_id + 1) % (mesh.dim() + 1);
            }
        }
    }
}

template <bool Border>
void TauLeapingDiffusionOperator::sample_triangle_fluxes(osh::Real dt) {
    const auto& dep_map = kproc_state_.get_dependency_map_elems();
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());

    const util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t>&
        triangles = [this]() constexpr -> const auto& {
        if constexpr (Border) {
            return border_tris;
        } else {
            return internal_tris;
        }
    }();

    for (const auto itri: triangles.range()) {
        const auto& tri = triangles[itri];

        if constexpr (Border) {
            if (not mesh.isOwned(tri)) {
                continue;
            }
        }

        mesh::tetrahedron_local_id_t tet0{tri2tets_ab2b[tri2tets_a2ab[tri.get()]]};
        mesh::tetrahedron_local_id_t tet1{tri2tets_ab2b[tri2tets_a2ab[tri.get()] + 1]};
        auto face_id0 = tri2face(tri.get(), 0)[0];
        auto face_id1 = tri2face(tri.get(), 1)[0];
        auto comp0 = mesh::compartment_id(tri2face(tri.get(), 0)[1]);
        auto comp1 = mesh::compartment_id(tri2face(tri.get(), 1)[1]);

        auto species = [this, tet0, itri]() constexpr -> auto {
            if constexpr (Border) {
                return border_tet_pop.species(border_tris_info[itri].itet0);
            } else {
                return pools.species(tet0);
            }
        }();
        for (auto spec: species) {
            auto spec_ind = container::species_id::unknown_value();

            // Get in and out fluxes for that triangle
            osh::Real Dout = diffusions_.ith_rate(tet0, spec, face_id0);
            osh::Real out = out_fluxes.synced_value(tet0, spec, face_id0);
            osh::Real in = 0.0;
            osh::Real Din = 0.0;
            if (comp0 == comp1) {
                Din = diffusions_.ith_rate(tet1, spec, face_id1);
                spec_ind = spec.get();
                in = out_fluxes.synced_value(tet1, spec, face_id1);

            } else if (mesh.getDiffusionBoundaryDcst(tri, comp1, spec) != 0.0) {
                auto ls = mesh.convertSpeciesID(tri, comp1, spec);
                Din = diffusions_.ith_rate(tet1, ls, face_id1);
                in = out_fluxes.synced_value(tet1, ls, face_id1);
                spec_ind = ls.get();
            }
            container::species_id spec1(spec_ind);

            if (in == 0.0 and out == 0.0) {
                continue;
            }

            // Sample net flux
            molecules_t diff = 0;
            if (std::max(in, out) * dt > normal_approx_thresh) {
                osh::Real net_flux = (out - in) * dt;
                osh::Real var = (out + in) * dt;
                if (dt > crank_nicolson_thresh * const_dt) {
                    // Crank-Nicolson scheme
                    // This prevents oscillations for larger periods, and, in contrast with the
                    // backward euler method, it leads to less underestimation of diffusion.
                    diff = std::round(
                        2.0 *
                        (net_flux + std::sqrt(var) * static_cast<osh::Real>(rng.getStdNrm())) /
                        (2.0 + dt * (Dout + Din)));
                } else {
                    diff = std::round(net_flux +
                                      std::sqrt(var) * static_cast<osh::Real>(rng.getStdNrm()));
                }
            } else {
                // If the normal approximation does not hold, sample using correlated binomial
                // distributions. For this, we need to keep track of the number of molecules
                // that are still available for diffusion in the source tetrahedron (nbLeft),
                // and of the cumulated diffusion rates that were already sampled from (Dcum).
                if (out != 0.0) {
                    osh::Real Dtot0 = diffusions_.rates_sum(tet0, spec);
                    molecules_t nbLeft = std::max(0.0,
                                                  std::round(out / Dout) -
                                                      tet_outbound(tet0.get(), spec.get())[0]);
                    osh::Real Dcum = tet_outbound(tet0.get(), spec.get())[1];
                    if (Dtot0 * dt >= 1) {
                        diff += molecules_t(rng.getBinom(nbLeft, Dout / (Dtot0 - Dcum)));
                    } else {
                        diff += molecules_t(rng.getBinom(nbLeft, Dout * dt / (1.0 - Dcum * dt)));
                    }
                }
                if (in != 0.0) {
                    osh::Real Dtot1 = diffusions_.rates_sum(tet1, spec1);
                    molecules_t nbLeft = std::max(0.0,
                                                  std::round(in / Din) -
                                                      tet_outbound(tet1.get(), spec1.get())[0]);
                    osh::Real Dcum = tet_outbound(tet1.get(), spec1.get())[1];
                    if (Dtot1 * dt >= 1) {
                        diff -= molecules_t(rng.getBinom(nbLeft, Din / (Dtot1 - Dcum)));
                    } else {
                        diff -= molecules_t(rng.getBinom(nbLeft, Din * dt / (1.0 - Dcum * dt)));
                    }
                }
            }

            molecules_t pop0, pop1;
            mesh::tetrahedron_internal_id_t itet0, itet1;
            if constexpr (Border) {
                const auto& t = border_tris_info[itri];
                itet0 = t.itet0;
                itet1 = t.itet1;
                pop0 = border_tet_pop.value(itet0, spec);
                pop1 = border_tet_pop.value(itet1, spec1);
            } else {
                pop0 = pools(tet0, spec);
                pop1 = pools(tet1, spec1);
            }

            // Restrict the flux to what is available in the source tet
            if (diff > 0 and diff > pop0) {
                diff = pop0;
            } else if (diff < 0 and -diff > pop1) {
                diff = -pop1;
            }

            if (diff != 0) {
                if constexpr (Border) {
                    // Set the sampled diffusion
                    border_tri_diff.value(itri, spec) = diff;
                    // Update available border tetrahedron populations after the sampled diffusion
                    border_tet_pop.value(itet0, spec) -= diff;
                    border_tet_pop.value(itet1, spec1) += diff;
                } else {
                    // Apply the sampled diffusion
                    if (comp0 == comp1) {
                        pools.add<false>(tet0, spec, -diff);
                        pools.add<false>(tet1, spec1, diff);
                    } else {
                        pools.add<true>(tet0, spec, -diff);
                        pools.add<true>(tet1, spec1, diff);
                    }
                    auto ab = pools.moleculesOnElements().ab(tet0, spec);
                    auto ab1 = pools.moleculesOnElements().ab(tet1, spec1);
                    kproc_state_.add_outdated_kprocs(dep_map[ab]);
                    kproc_state_.add_outdated_kprocs(dep_map[ab1]);
                }

                // Update nbLeft in tet_outbound
                if (diff > 0) {
                    tet_outbound(tet0.get(), spec.get())[0] += diff;
                } else {
                    tet_outbound(tet1.get(), spec1.get())[0] += -diff;
                }
            }

            num_diffusions_ += std::abs(diff);
            // Update Dcum in tet_outbound
            tet_outbound(tet0.get(), spec.get())[1] += Dout;
            tet_outbound(tet1.get(), spec1.get())[1] += Din;
        }
    }
}

void TauLeapingDiffusionOperator::apply_border_diffusion() {
    // Apply the sampled and synced border diffusions
    const auto& dep_map = kproc_state_.get_dependency_map_elems();

    for (const auto itri: border_tris_info.range()) {
        const auto& tri_info = border_tris_info[itri];
        auto tri = tri_info.tri;
        auto tet0 = tri_info.tet0;
        auto tet1 = tri_info.tet1;
        auto comp0 = mesh.getRegionMeshID(tet0);
        auto comp1 = mesh.getRegionMeshID(tet1);
        if (mesh.isOwned(tet0)) {
            for (auto spec0: pools.species(tet0)) {
                auto num_mols = -border_tri_diff.synced_value(itri, spec0);
                if (num_mols != 0) {
                    if (comp0 == comp1) {
                        pools.add<false>(tet0, spec0, num_mols);
                    } else {
                        pools.add<true>(tet0, spec0, num_mols);
                    }
                    auto ab = pools.moleculesOnElements().ab(tet0, spec0);
                    kproc_state_.add_outdated_kprocs(dep_map[ab]);
                }
            }
        } else {
            for (auto spec1: pools.species(tet1)) {
                if (comp0 == comp1) {
                    auto num_mols = border_tri_diff.synced_value(itri, spec1);
                    if (num_mols != 0) {
                        pools.add<false>(tet1, spec1, num_mols);
                        auto ab = pools.moleculesOnElements().ab(tet1, spec1);
                        kproc_state_.add_outdated_kprocs(dep_map[ab]);
                    }
                } else if (mesh.getDiffusionBoundaryDcst(tri, comp1, spec1) != 0.0) {
                    auto spec0 = mesh.convertSpeciesID(tri, comp0, spec1);
                    auto num_mols = border_tri_diff.synced_value(itri, spec0);
                    if (num_mols != 0) {
                        pools.add<true>(tet1, spec1, num_mols);
                        auto ab = pools.moleculesOnElements().ab(tet1, spec1);
                        kproc_state_.add_outdated_kprocs(dep_map[ab]);
                    }
                }
            }
        }
    }
}

osh::Real TauLeapingDiffusionOperator::getDt() {
    // TODO Add code for surface diffusion once this is implemented
    osh::Real rd_dt = std::numeric_limits<osh::Real>::infinity();
    for (auto tet0: mesh.owned_elems()) {
        for (auto spec0: pools.species(tet0)) {
            osh::Real cnt = pools(tet0, spec0);

            if (cnt > 0 and cnt < leap_threshold) {
                return const_dt * (1 - std::numeric_limits<osh::Real>::epsilon());
            }

            osh::Real out = diffusions_.rates_sum(tet0, spec0) * cnt;
            osh::Real mu = -out;
            osh::Real sigma = out;

            auto num_elements = mesh.tet_neighbors_int_data().size(tet0.get());
            auto comp0 = mesh.getRegionMeshID(tet0);
            for (auto face_id0 = 0; face_id0 < num_elements; face_id0++) {
                const auto& neighbor = mesh.tet_neighbors_int_data()(tet0.get(), face_id0);
                const mesh::tetrahedron_local_id_t tet1(neighbor[0]);
                const auto face_id1 = neighbor[1];
                const mesh::triangle_local_id_t tri(neighbor[2]);

                auto comp1 = mesh.getRegionMeshID(tet1);

                osh::Real in = 0.0;
                if (comp0 == comp1) {
                    if (mesh.isOwned(tet1)) {
                        in = diffusions_.ith_rate(tet1, spec0, face_id1) *
                             static_cast<osh::Real>(pools(tet1, spec0));
                    } else {
                        // If the tetrahedron is not owned, use the previously synchronized flux
                        // as an approximation for the actual flux, which would be too costly to
                        // synchronize
                        in = out_fluxes.synced_value(tet1, spec0, face_id1);
                    }
                } else if (mesh.getDiffusionBoundaryDcst(mesh::triangle_local_id_t(tri),
                                                         comp1,
                                                         spec0) != 0.0) {
                    auto ls = mesh.convertSpeciesID(tri, comp1, spec0);
                    if (mesh.isOwned(tet1)) {
                        in = diffusions_.ith_rate(tet1, ls, face_id1) *
                             static_cast<osh::Real>(pools(tet1, ls));
                    } else {
                        in = out_fluxes.synced_value(tet1, ls, face_id1);
                    }
                }

                mu += in;
                sigma += in;
            }

            osh::Real thresh = std::max(tolerance_ * cnt, 1.0);
            if (mu != 0) {
                rd_dt = std::min(rd_dt, thresh / std::abs(mu));
            }
            if (sigma > 0) {
                rd_dt = std::min(rd_dt, thresh * thresh / sigma);
            }
        }
        if (rd_dt <= const_dt) {
            rd_dt = const_dt * (1 - std::numeric_limits<osh::Real>::epsilon());
            break;
        }
    }
    return rd_dt;
}

std::map<std::string, double> TauLeapingDiffusionOperator::getDebugInfo() const {
    std::map<std::string, double> info;
    info["total_diff_steps"] = num_diffusions_;
    info["total_diffleaping_steps"] = num_steps;

    auto copy = relative_step_sizes;
    std::sort(copy.begin(), copy.end());
    info["diffleaping_cum_rel_step_sizes"] = std::accumulate(copy.begin(), copy.end(), 0.0);

    info["diffleaping_cum_squared_rel_step_sizes"] =
        std::transform_reduce(copy.begin(), copy.end(), 0.0, std::plus<>(), [](osh::Real& v) {
            return v * v;
        });

    return info;
}

osh::LOs TauLeapingDiffusionOperator::diff_species_per_triangle(DistMesh& mesh,
                                                                const osh::LOs spec_per_elems) {
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());
    auto ntris = mesh.getTriInfo().size();
    osh::Write<osh::LO> spec_per_tri(ntris);
    for (auto tri: mesh.getTriInfo().range()) {
        mesh::tetrahedron_local_id_t tet0{tri2tets_ab2b[tri2tets_a2ab[tri.get()]]};
        spec_per_tri[tri.get()] = spec_per_elems[tet0.get()];
    }
    return osh::LOs(spec_per_tri);
}

osh::LOs TauLeapingDiffusionOperator::neighbs_per_triangle(DistMesh& mesh) {
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    auto ntris = mesh.getTriInfo().size();
    osh::Write<osh::LO> neighbs_per_tri(ntris);
    std::adjacent_difference(tri2tets_a2ab.begin() + 1,
                             tri2tets_a2ab.end(),
                             neighbs_per_tri.begin());
    return osh::LOs(neighbs_per_tri);
}

util::strongid_vector<mesh::triangle_internal_id_t, TauLeapingDiffusionOperator::BorderTriangle>
TauLeapingDiffusionOperator::border_triangles_info(DistMesh& mesh) {
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());
    const auto tet2tris = mesh.ask_down(mesh.dim(), mesh.dim() - 1);

    util::strongid_vector<mesh::triangle_internal_id_t, TauLeapingDiffusionOperator::BorderTriangle>
        triangles;
    std::map<mesh::tetrahedron_local_id_t, mesh::tetrahedron_internal_id_t> tets;
    for (auto tri: mesh::triangle_local_id_t::range(mesh.getTriInfo().size())) {
        osh::LO noffsets = tri2tets_a2ab[tri.get() + 1] - tri2tets_a2ab[tri.get()];
        assert(noffsets <= 2);
        if (noffsets == 2) {
            mesh::tetrahedron_local_id_t tet0{tri2tets_ab2b[tri2tets_a2ab[tri.get()]]};
            mesh::tetrahedron_local_id_t tet1{tri2tets_ab2b[tri2tets_a2ab[tri.get()] + 1]};
            if (mesh.isOwned(tet0) != mesh.isOwned(tet1)) {
                auto [it0, _] = tets.emplace(tet0, tets.size());
                auto [it1, __] = tets.emplace(tet1, tets.size());
                short face0 = 0;
                short face1 = 0;
                const auto tris0 = osh::gather_down<DistMesh::dim() + 1>(tet2tris.ab2b, tet0.get());
                for (auto tri0: tris0) {
                    if (tri0 == tri.get()) {
                        break;
                    }
                    ++face0;
                }
                assert(face0 < mesh.dim() + 1);
                const auto tris1 = osh::gather_down<DistMesh::dim() + 1>(tet2tris.ab2b, tet1.get());
                for (auto tri1: tris1) {
                    if (tri1 == tri.get()) {
                        break;
                    }
                    ++face1;
                }
                assert(face1 < mesh.dim() + 1);

                triangles.container().emplace_back(
                    tri, mesh.isOwned(tri), tet0, tet1, it0->second, it1->second, face0, face1);
            }
        }
    }
    return triangles;
}

util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t>
TauLeapingDiffusionOperator::border_triangles_ids(
    const util::strongid_vector<mesh::triangle_internal_id_t, BorderTriangle>& border_tris_info) {
    util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t> ids(
        border_tris_info.size());
    for (const auto itri: border_tris_info.range()) {
        const auto& info = border_tris_info[itri];
        ids[itri] = info.tri;
    }
    return ids;
}

util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t>
TauLeapingDiffusionOperator::internal_triangles(DistMesh& mesh) {
    const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
    const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());

    util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t> tris;
    tris.container().reserve(mesh.owned_bounds().size());

    for (auto tri: mesh.owned_bounds()) {
        osh::LO noffsets = tri2tets_a2ab[tri + 1] - tri2tets_a2ab[tri];
        assert(noffsets <= 2);
        if (noffsets == 2) {
            mesh::tetrahedron_local_id_t tet0{tri2tets_ab2b[tri2tets_a2ab[tri]]};
            mesh::tetrahedron_local_id_t tet1{tri2tets_ab2b[tri2tets_a2ab[tri] + 1]};
            if (mesh.isOwned(tet0) and mesh.isOwned(tet1)) {
                tris.container().emplace_back(tri);
            }
        }
    }

    return tris;
}

util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t>
TauLeapingDiffusionOperator::border_tetrahedrons(
    const util::strongid_vector<mesh::triangle_internal_id_t,
                                TauLeapingDiffusionOperator::BorderTriangle>& border_triangles) {
    util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t> tets;
    for (const auto& t: border_triangles) {
        if (t.itet0.get() >= tets.size()) {
            tets.container().resize(t.itet0.get() + 1);
        }
        tets[t.itet0] = t.tet0;
        if (t.itet1.get() >= tets.size()) {
            tets.container().resize(t.itet1.get() + 1);
        }
        tets[t.itet1] = t.tet1;
    }
    return tets;
}

osh::LOs TauLeapingDiffusionOperator::diff_species_per_border_tetrahedron(
    const util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t>&
        border_tets,
    const osh::LOs spec_per_elems) {
    osh::Write<osh::LO> spec_per_tet(border_tets.size());
    const auto fill_specs = [&](osh::LO i) {
        mesh::tetrahedron_internal_id_t itet(i);
        spec_per_tet[itet.get()] = spec_per_elems[border_tets[itet].get()];
    };
    osh::parallel_for(border_tets.size(), fill_specs);

    return spec_per_tet;
}

osh::LOs TauLeapingDiffusionOperator::diff_species_per_border_triangle(
    const util::strongid_vector<mesh::triangle_internal_id_t,
                                TauLeapingDiffusionOperator::BorderTriangle>& border_triangles,
    const osh::LOs spec_per_elems) {
    osh::Write<osh::LO> spec_per_tri(border_triangles.size());
    for (const auto itri: border_triangles.range()) {
        const auto& info = border_triangles[itri];
        spec_per_tri[itri.get()] = spec_per_elems[info.tet0.get()];
    }
    return spec_per_tri;
}

}  // namespace steps::dist

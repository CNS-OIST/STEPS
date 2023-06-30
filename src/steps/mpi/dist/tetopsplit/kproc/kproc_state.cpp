
#include "kproc_state.hpp"

#include <stdexcept>

#include <Omega_h_for.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>

#include "diffusions.hpp"
#include "mpi/dist/tetopsplit/definition/compdef.hpp"
#include "mpi/dist/tetopsplit/definition/statedef.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist::kproc {

////------------------------------------------------------------------
KProcState::KProcState(const Statedef& statedef,
                       DistMesh& mesh,
                       MolState& mol_state,
                       bool independent_kprocs)
    : mol_state_(mol_state)
    , reactions_(statedef, mesh, mol_state)
    , surface_reactions_(statedef, mesh, mol_state)
    , vdep_surface_reactions_(statedef, mesh, mol_state)
    , ghk_surface_reactions_(statedef, mesh, mol_state) {
    setupDependencies();
    // extract connected components of the Gibson-Bruck dependency graph
    setupGroups(independent_kprocs);
}

//------------------------------------------------------------------

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
typename propensity_function_traits::value KProcState::propensityFun() const {
    return [this](KProcID k_id, const MolState& mol_state) {
        osh::Real rate;
        switch (k_id.type()) {
        case KProcType::Reac:
            rate = reactions_.computeRate(mol_state, k_id.id());
            break;
        case KProcType::SReac:
            rate = surface_reactions_.computeRate(mol_state, k_id.id());
            break;
        case KProcType::VDepSReac:
            rate = vdep_surface_reactions_.computeRate(mol_state, k_id.id());
            break;
        case KProcType::GHKSReac:
            rate = ghk_surface_reactions_.computeRate(mol_state, k_id.id());
            break;
        case KProcType::Diff:
            throw std::logic_error("Unhandled kinetic process type: Diff");
        }

        if (rate == std::numeric_limits<double>::infinity()) {
            throw std::logic_error("The reaction rate for kproc id: " + std::to_string(k_id.id()) +
                                   " is infinite. Either"
                                   " a constant is 0 or "
                                   "the voltage spiked "
                                   "to very high values");
        }

        return rate;
    };
}
#pragma GCC diagnostic pop

//------------------------------------------------------------------

template <typename KineticProcesses>
void KProcState::countDependencies(const KineticProcesses& processes,
                                   osh::Write<osh::LO>& dep_map_elems_num_data,
                                   osh::Write<osh::LO>& dep_map_bnds_num_data) const {
    for (const auto& process: processes) {
        for (const auto& dependency: process.getPropensityDependency()) {
            auto species = std::get<1>(dependency);
            std::visit(
                [&](auto& entity) {
                    using T = std::decay_t<decltype(entity)>;

                    if constexpr (std::is_same_v<T, mesh::tetrahedron_id_t>) {
                        const auto idx = mol_state_.moleculesOnElements().ab(entity, species);
                        ++dep_map_elems_num_data[idx];
                    } else if constexpr (std::is_same_v<T, mesh::triangle_id_t>) {
                        const auto idx = mol_state_.moleculesOnPatchBoundaries().ab(entity,
                                                                                    species);
                        ++dep_map_bnds_num_data[idx];
                    } else {
                        static_assert(util::always_false_v<T>, "unmanaged entity type");
                    }
                },
                std::get<0>(dependency));
        }
    }
}

//------------------------------------------------------------------

template <typename KineticProcesses>
void KProcState::fillDependencies(const KineticProcesses& processes,
                                  osh::Write<osh::LO>& elems_curr_counters,
                                  osh::Write<osh::LO>& bnds_curr_counters) {
    for (const auto& process: processes) {
        auto kid = KProcID(processes.getKProcType(), process.getIndex());
        for (const auto& dependency: process.getPropensityDependency()) {
            auto species = std::get<1>(dependency);
            std::visit(
                [&](auto& entity) {
                    using T = std::decay_t<decltype(entity)>;

                    if constexpr (std::is_same_v<T, mesh::tetrahedron_id_t>) {
                        const auto idx = mol_state_.moleculesOnElements().ab(entity, species);
                        // elems_curr_counters is always off by one. -- goes first
                        dependency_map_elems_(idx, --elems_curr_counters[idx]) = kid.data();
                    } else if constexpr (std::is_same_v<T, mesh::triangle_id_t>) {
                        const auto idx = mol_state_.moleculesOnPatchBoundaries().ab(entity,
                                                                                    species);
                        // dependency_map_bnds_ is always off by one. -- goes first
                        dependency_map_bnds_(idx, --bnds_curr_counters[idx]) = kid.data();
                    } else {
                        static_assert(util::always_false_v<T>, "unmanaged entity type");
                    }
                },
                std::get<0>(dependency));
        }
    }
}

//------------------------------------------------------------------

template <typename KineticProcesses>
void KProcState::cacheDependencies(const KineticProcesses& processes,
                                   dependencies_t& dependencies) {
    std::vector<std::set<KProcID>> unique_deps(processes.size());
    osh::Write<osh::LO> sizes(static_cast<osh::LO>(processes.size()));

    // compute dependencies in a temporary datastructure that avoid taking the
    // same dependency more than once
    for (const auto& process: processes) {
        const auto& mol_state_elements_to_update = process.getMolStateElementsUpdates();
        for (const auto& mol_state_element: mol_state_elements_to_update) {
            auto species = std::get<1>(mol_state_element);
            std::visit(
                [&](auto& entity) {
                    using T = std::decay_t<decltype(entity)>;

                    if constexpr (std::is_same_v<T, mesh::tetrahedron_id_t>) {
                        const auto idx = mol_state_.moleculesOnElements().ab(entity, species);
                        auto& deps = unique_deps[process.getIndex()];
                        for (auto elem: dependency_map_elems_[idx]) {
                            deps.emplace(elem);
                        }
                    } else if constexpr (std::is_same_v<T, mesh::triangle_id_t>) {
                        const auto idx = mol_state_.moleculesOnPatchBoundaries().ab(entity,
                                                                                    species);
                        auto& deps = unique_deps[process.getIndex()];
                        for (auto bound: dependency_map_bnds_[idx]) {
                            deps.emplace(bound);
                        }
                    } else {
                        static_assert(util::always_false_v<T>, "unmanaged entity type");
                    }
                },
                std::get<0>(mol_state_element));
            sizes[static_cast<osh::LO>(process.getIndex())] = static_cast<osh::LO>(
                unique_deps[process.getIndex()].size());
        }
    }

    dependencies.reshape(sizes);
    for (auto pid = 0u; pid < processes.size(); ++pid) {
        const auto& udeps = unique_deps[pid];
        std::transform(udeps.begin(),
                       udeps.end(),
                       dependencies[static_cast<osh::LO>(pid)].begin(),
                       [](const KProcID& kid) { return static_cast<osh::LO>(kid.data()); });
    }
}

//------------------------------------------------------------------

#pragma GCC diagnostic push
#if !defined(__clang__) && (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

std::ostream& KProcState::write_dependency_graph(std::ostream& ostr) const {
    // print dependency graph only if we have a filename
    const auto& [grd, labels] = getDependenciesGraphAndLabels();
    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::vecS> gd;
    boost::copy_graph(grd, gd);
    std::vector<unsigned> color(boost::num_vertices(grd));
    auto it = boost::make_iterator_property_map(color.begin(), get(boost::vertex_index, grd));
    boost::connected_components(grd, it);
    const std::vector<std::string> col_vals{"red", "blue", "green", "yellow", "magenta", "grey"};
    auto node_fmt = [labels = &labels, &color, &col_vals](std::ostream& out, unsigned v) {
        out << " [label=\"" << (*labels)[v] << "\"]" << std::endl;
        out << " [color=" << col_vals[color[v] % col_vals.size()] << "]" << std::endl;
    };
    boost::write_graphviz(ostr, gd, node_fmt);
    return ostr;
}

#pragma GCC diagnostic pop

void KProcState::setupDependencies() {
    osh::Write<osh::LO> elems_a2ab(mol_state_.moleculesOnElements().num_data(), 0);
    osh::Write<osh::LO> bnds_a2ab(mol_state_.moleculesOnPatchBoundaries().num_data(), 0);

    // counting all the dependencies to create/reshape the flat multimaps for the dependencies
    countDependencies(reactions(), elems_a2ab, bnds_a2ab);
    countDependencies(surfaceReactions(), elems_a2ab, bnds_a2ab);
    countDependencies(vDepSurfaceReactions(), elems_a2ab, bnds_a2ab);
    countDependencies(ghkSurfaceReactions(), elems_a2ab, bnds_a2ab);

    // reshaping/creation
    dependency_map_elems_.reshape(elems_a2ab);
    dependency_map_bnds_.reshape(bnds_a2ab);

    // filling of the dep maps. Reverse order (compared to counting) for max efficiency.
    // elems_a2ab and bnds_a2ab get gradually emptied. All their elements must be 0 at the end.
    // No checks for improved efficiency, this bug is improbable once it is set correctly
    fillDependencies(reactions(), elems_a2ab, bnds_a2ab);
    fillDependencies(surfaceReactions(), elems_a2ab, bnds_a2ab);
    fillDependencies(vDepSurfaceReactions(), elems_a2ab, bnds_a2ab);
    fillDependencies(ghkSurfaceReactions(), elems_a2ab, bnds_a2ab);

    // dependency. It says from which kprocids a certain kproc depends
    // dep_map.first depends on all the seconds
    // caching the kinetic processes dependencies
    // necessary for the ssa operator
    cacheDependencies(reactions(), reactions_dependencies_);
    cacheDependencies(surfaceReactions(), surface_reactions_dependencies_);
    cacheDependencies(vDepSurfaceReactions(), vdep_surface_reactions_dependencies_);
    cacheDependencies(ghkSurfaceReactions(), ghk_surface_reactions_dependencies_);
}


//------------------------------------------------------------------

typename KProcState::Graph KProcState::getDependenciesGraph(
    const Propensities<>& propensities) const {
    // propensities needed to extract a kproc index
    auto num_edges = [](const dependencies_t& deps) {
        return std::accumulate(deps.begin(), deps.end(), 0u, [](size_t acc, const auto& x) {
            return acc + static_cast<size_t>(x.size());
        });
    };
    assert(reactions().size() + surfaceReactions().size() + vDepSurfaceReactions().size() +
               ghkSurfaceReactions().size() <=
           static_cast<size_t>(std::numeric_limits<osh::LO>::max()));


    const size_t tot_num_edges = num_edges(reactions_dependencies_) +
                                 num_edges(surface_reactions_dependencies_) +
                                 num_edges(vdep_surface_reactions_dependencies_) +
                                 num_edges(ghk_surface_reactions_dependencies_);
    using Edge = std::pair<unsigned, unsigned>;
    std::vector<Edge> edges;
    edges.reserve(tot_num_edges);
    auto populate_edges = [&propensities, &edges](KProcType ty, const dependencies_t& deps) {
        for (osh::LO k = 0; k < deps.size(); ++k) {
            auto klo = propensities.ab(KProcID(ty, static_cast<unsigned>(k)));
            std::transform(deps[k].begin(),
                           deps[k].end(),
                           std::back_inserter(edges),
                           [&](const auto dep) -> Edge {
                               return {klo, propensities.ab(KProcID(static_cast<unsigned>(dep)))};
                           });
        }
    };
    populate_edges(KProcType::Reac, reactions_dependencies_);
    populate_edges(KProcType::SReac, surface_reactions_dependencies_);
    populate_edges(KProcType::VDepSReac, vdep_surface_reactions_dependencies_);
    populate_edges(KProcType::GHKSReac, ghk_surface_reactions_dependencies_);

    return {edges.begin(),
            edges.end(),
            reactions().size() + surfaceReactions().size() + vDepSurfaceReactions().size() +
                ghkSurfaceReactions().size()};
}

//------------------------------------------------------------------

void KProcState::setupGroups(bool independent_kprocs) {
    // propensities needed to extract a kproc. index
    Propensities propensities;
    propensities.init(handledKProcsClassesAndSizes());
    // extract the number of kproc dependencies
    if (independent_kprocs) {
        auto graph = getDependenciesGraph(propensities);
        std::vector<unsigned> color(boost::num_vertices(graph));
        auto it = boost::make_iterator_property_map(color.begin(), get(boost::vertex_index, graph));
        size_t num_components = boost::connected_components(graph, it);

        std::vector<std::vector<osh::LO>> disjoint_k_procs_w(num_components);
        osh::Write<osh::LO> sizes(static_cast<osh::LO>(num_components), 0);
        for (size_t k = 0; k < color.size(); ++k) {
            sizes[static_cast<osh::LO>(color[k])]++;
            disjoint_k_procs_w[color[k]].push_back(
                static_cast<osh::LO>(propensities.kProcId(k).data()));
        }
        disjoint_kprocs_.reshape(sizes);
        for (size_t k = 0; k < num_components; ++k) {
            std::copy(disjoint_k_procs_w[k].begin(),
                      disjoint_k_procs_w[k].end(),
                      disjoint_kprocs_[static_cast<osh::LO>(k)].begin());
        }
    } else {
        disjoint_kprocs_.reshape({static_cast<osh::LO>(propensities.size())});
        size_t k = 0;
        for (auto it = disjoint_kprocs_[0].begin(); it != disjoint_kprocs_[0].end(); it++, k++) {
            *it = static_cast<osh::LO>(propensities.kProcId(k).data());
        }
    }
}

//------------------------------------------------------------------

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
const std::vector<MolStateElementID>& KProcState::updateMolStateAndOccupancy(
    MolState& mol_state,
    const osh::Real event_time,
    const KProcID& event) const {
    switch (event.type()) {
    case KProcType ::Reac:
        return reactions().updateMolStateAndOccupancy(mol_state, event.id(), event_time);
    case KProcType::SReac:
        return surfaceReactions().updateMolStateAndOccupancy(mol_state, event.id(), event_time);
    case KProcType::VDepSReac:
        return vDepSurfaceReactions().updateMolStateAndOccupancy(mol_state, event.id(), event_time);
    case KProcType::GHKSReac:
        return ghkSurfaceReactions().updateMolStateAndOccupancy(mol_state, event.id(), event_time);
    case KProcType::Diff:
        throw std::logic_error("Unhandled kinetic process: Diffusion");
    }
}

#pragma GCC diagnostic ignored "-Wreturn-type"
void KProcState::report(std::ostream& report_stream, KProcID kid) const {
    switch (kid.type()) {
    case KProcType ::Reac:
        return reactions().report(report_stream, kid.id());
    case KProcType::SReac:
        return surfaceReactions().report(report_stream, kid.id());
    case KProcType::VDepSReac:
        return vDepSurfaceReactions().report(report_stream, kid.id());
    case KProcType::GHKSReac:
        return ghkSurfaceReactions().report(report_stream, kid.id());
    case KProcType::Diff:
        throw std::logic_error("Unhandled kinetic process: Diffusion");
    }
}
#pragma GCC diagnostic pop

//------------------------------------------------------------------

std::pair<typename KProcState::Graph, std::vector<std::string>>
KProcState::getDependenciesGraphAndLabels() const {
    Propensities<> propensities;
    propensities.init(handledKProcsClassesAndSizes());
    std::vector<std::string> labels(propensities.size());
    std::transform(boost::make_counting_iterator(size_t{}),
                   boost::make_counting_iterator(propensities.size()),
                   labels.begin(),
                   [&](auto x) {
                       std::ostringstream s;
                       report(s, propensities.kProcId(x));
                       return s.str();
                   });
    return {getDependenciesGraph(propensities), labels};
}

//--------------------------------------------------------

// explicit template instantiation definitions
template void KProcState::countDependencies(const Reactions& processes,
                                            osh::Write<osh::LO>& dep_map_elems_num_data,
                                            osh::Write<osh::LO>& dep_map_bnds_num_data) const;
template void KProcState::countDependencies(const SurfaceReactions& processes,
                                            osh::Write<osh::LO>& dep_map_elems_num_data,
                                            osh::Write<osh::LO>& dep_map_bnds_num_data) const;
template void KProcState::countDependencies(const VDepSurfaceReactions& processes,
                                            osh::Write<osh::LO>& dep_map_elems_num_data,
                                            osh::Write<osh::LO>& dep_map_bnds_num_data) const;
template void KProcState::countDependencies(const GHKSurfaceReactions& processes,
                                            osh::Write<osh::LO>& dep_map_elems_num_data,
                                            osh::Write<osh::LO>& dep_map_bnds_num_data) const;

template void KProcState::fillDependencies(const Reactions& processes,
                                           osh::Write<osh::LO>& elems_curr_counters,
                                           osh::Write<osh::LO>& bnds_curr_counters);
template void KProcState::fillDependencies(const SurfaceReactions& processes,
                                           osh::Write<osh::LO>& elems_curr_counters,
                                           osh::Write<osh::LO>& bnds_curr_counters);
template void KProcState::fillDependencies(const VDepSurfaceReactions& processes,
                                           osh::Write<osh::LO>& elems_curr_counters,
                                           osh::Write<osh::LO>& bnds_curr_counters);
template void KProcState::fillDependencies(const GHKSurfaceReactions& processes,
                                           osh::Write<osh::LO>& elems_curr_counters,
                                           osh::Write<osh::LO>& bnds_curr_counters);

}  // namespace steps::dist::kproc

#include "surface_reactions.hpp"

#include <cassert>
#include <functional>

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_shape.hpp>
#include <stdexcept>
#include <string>

#include "geom/dist/distmesh.hpp"
#include "math/constants.hpp"
#include "model/complex.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist::kproc {


//-------------------------------------------------------

template <typename RDefT>
SurfaceReactionsBase<RDefT>::SurfaceReactionsBase(const Statedef& statedef,
                                                  DistMesh& mesh,
                                                  MolState& mol_state)
    : state_def_(statedef) {
    if (!statedef.patchdefs().container().empty()) {
        /// tet_owned informs whether mesh element is owned
        const auto& tet_owned = mesh.owned_elems_mask();
        /// tri2tets mapping a2ab and ab2b from a triangle id to its two neighbouring tets
        const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
        const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());
        /// and work through all patches
        for (const auto& patch: statedef.patchdefs()) {
            patch2ind_.container().emplace_back(boundary_id_.size(),
                                                patch->reacdefs<RDefT>().size());

            model::patch_id patch_id = patch->getID();
            model::compartment_id inner_compartment_id = patch->getInnerCompId();
            const std::optional<model::compartment_id>& outer_compartment_id =
                patch->getOuterCompId();

            osh::Write<osh::LO> inner_comp_mask(tet_owned.size(), 0);
            for (auto elem: mesh.getEntities(inner_compartment_id)) {
                inner_comp_mask[elem.get()] = 1;
            }

            osh::Write<osh::LO> outer_comp_mask(outer_compartment_id ? tet_owned.size() : 0, 0);
            if (outer_compartment_id) {
                for (auto elem: mesh.getEntities(*outer_compartment_id)) {
                    outer_comp_mask[elem.get()] = 1;
                }
            }

            /// patch_owned_boundaries is a vector of patch boundary ids owned by the
            /// process.
            const auto& patch_owned_boundaries = mesh.getOwnedEntities(patch_id);

            for (const auto boundary: patch_owned_boundaries) {
                /// adjacent tet ids of patch element 'element' are stored in
                /// tritets.a2ab
                osh::LO noffsets = tri2tets_a2ab[boundary.get() + 1] -
                                   tri2tets_a2ab[boundary.get()];
                if (noffsets > 2) {
                    std::ostringstream oss;
                    oss << "Geometry: wrong number of adjacent elements of dimension " << mesh.dim()
                        << " for boundary " << boundary << " of dimension " << mesh.dim() - 1
                        << ": " << noffsets;

                    throw std::logic_error(oss.str());
                }
                std::optional<mesh::tetrahedron_id_t> outer_compartment_elmt,
                    inner_compartment_elmt;
                /// work out the compartment associated with adjacent tets
                for (osh::LO k = 0; k < noffsets; ++k) {
                    auto tet_id = tri2tets_ab2b[tri2tets_a2ab[boundary.get()] + k];
                    /// check whether it is in the inner compartment
                    if (inner_comp_mask[tet_id] == 1) {
                        inner_compartment_elmt = mesh::tetrahedron_id_t(tet_id);
                    } else {
                        /// it should be associated with the outer compartment
                        if (outer_compartment_id) {
                            if (outer_comp_mask[tet_id] == 0) {
                                std::ostringstream oss;
                                oss << "Geometry : the outer compartment cannot be resolved "
                                       "for patch "
                                    << patch_id;
                                throw std::logic_error(oss.str());
                            }
                        }
                        outer_compartment_elmt = mesh::tetrahedron_id_t(tet_id);
                    }
                }
                /// inner comp needs to be resolved by now.
                if (!inner_compartment_elmt) {
                    std::ostringstream oss;
                    oss << "Inner compartment of element " << boundary << " cannot be resolved.";
                    std::logic_error(oss.str());
                }
                /// are these tets owned by the process?
                if (tet_owned[inner_compartment_elmt->get()] == 0) {
                    std::ostringstream oss;
                    oss << "SReac : the inner compartment element " << *inner_compartment_elmt
                        << " is not in the same process as the patch " << patch_id;
                    throw std::logic_error(oss.str());
                }
                if (outer_compartment_id && (tet_owned[outer_compartment_elmt->get()] == 0)) {
                    std::ostringstream oss;
                    oss << "SReac : the outer compartment element " << *outer_compartment_elmt
                        << " is not in the same process as the patch " << patch_id;
                    throw std::logic_error(oss.str());
                }
                for (const auto& reacdef: patch->reacdefs<RDefT>()) {
                    reacdefs_.push_back(*reacdef);
                    boundary_id_.push_back(boundary);
                    inner_compartment_element_id_.push_back(*inner_compartment_elmt);
                    outer_compartment_element_id_.push_back(outer_compartment_elmt);
                    /// resolve the stoichiometry of the reaction in respect of the mesh.
                    {
                        const auto& [elmts, stoichiometry, region_id] =
                            reactionMolStateDependencyAndStoichiometry<PoolChangeType::LHS>(
                                *reacdef,
                                boundary,
                                *inner_compartment_elmt,
                                outer_compartment_elmt);

                        reaction_lhs_.push_back(elmts);
                        stoichiometry_lhs_.push_back(stoichiometry);
                    }

                    /// region_id is necessary to identify the region of species involved
                    /// in the reaction, to check whether the species is diffused.
                    const auto& [elmts, stoichiometry, region_id] =
                        reactionMolStateDependencyAndStoichiometry<PoolChangeType::UPD>(
                            *reacdef, boundary, *inner_compartment_elmt, outer_compartment_elmt);

                    for (size_t k = 0; k < elmts.size(); ++k) {
                        const auto spec = std::get<1>(elmts[k]);
                        std::visit(
                            [&spec, &statedef, &mol_state](const auto& v, const auto& e) {
                                if (statedef.getDefinition(v).isDiffused(spec)) {
                                    mol_state.track_occupancy_rd(e, spec);
                                }
                            },
                            region_id[k],
                            std::get<0>(elmts[k]));
                    }
                    reaction_upd_.push_back(elmts);
                    stoichiometry_upd_.push_back(stoichiometry);
                    kinConstantGeomFactor_.push_back(
                        kinConstantGeomFactor(mesh, kinConstantGeomFactor_.size()));
                    // compute the propensity rate constant.
                    ccsts_.push_back(std::numeric_limits<osh::Real>::quiet_NaN());
                }
            }
        }
    }
    currents_ = osh::Write<osh::Real>(this->size(), 0.0);
}

//-------------------------------------------------------

template <typename RDefT>
const std::vector<MolStateElementID>& SurfaceReactionsBase<RDefT>::updateMolStateAndOccupancy(
    MolState& mol_state,
    rng::RNG& /*rng*/,
    size_t reaction_index,
    const osh::Real event_time,
    int nb) const {
    const auto& upd = stoichiometry_upd_[reaction_index];
    const auto& mol_state_elements = reaction_upd_[reaction_index];
    reacdefs_[reaction_index].get().add_extent(nb);

    for (size_t k = 0; k < mol_state_elements.size(); ++k) {
        const auto& elmt_id = mol_state_elements[k];
        const auto s = upd[k] * nb;
        mol_state.add_and_update_occupancy(elmt_id, static_cast<molecules_t>(s), event_time);
    }
    return mol_state_elements;
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real SurfaceReactionsBase<RDefT>::kinConstantGeomFactor(const DistMesh& mesh,
                                                             size_t reaction_index) const {
    const RDefT& reacdef = reacdefs_[reaction_index];
    const mesh::triangle_id_t boundary{boundary_id_[reaction_index]};
    const mesh::tetrahedron_id_t inner_compartment_element{
        inner_compartment_element_id_[reaction_index]};
    const std::optional<mesh::tetrahedron_id_t>& outer_compartment_element =
        outer_compartment_element_id_[reaction_index];
    if (reacdef.isSurfaceSurfaceReaction()) {
        /// the reaction happens on the surface of the patch.
        const auto area = mesh.getTri(boundary).area;
        const double ascale = area * math::AVOGADRO;
        const osh::I64 o1 = reacdef.getOrder() - 1;
        return std::pow(ascale, static_cast<osh::Real>(-o1));
    } else {
        osh::Real vol;
        if (reacdef.isInnerCompartmentReaction()) {
            vol = mesh.getTet(inner_compartment_element).vol;
        } else {
            vol = mesh.getTet(*outer_compartment_element).vol;
        }
        const double vscale = 1.0e3 * vol * math::AVOGADRO;
        const osh::I64 o1 = reacdef.getOrder() - 1;
        return std::pow(vscale, static_cast<osh::Real>(-o1));
    }
}

//-------------------------------------------------------

template <typename RDefT>
size_t SurfaceReactionsBase<RDefT>::getInd(container::patch_id patch,
                                           container::triangle_id tri,
                                           ContReacID reac) const {
    auto [ind, nbreacs] = patch2ind_.at(patch);
    assert(tri.valid());
    assert(reac.valid());
    return ind + nbreacs * tri.get() + reac.get();
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real SurfaceReactionsBase<RDefT>::computeRate(const MolState& mol_state,
                                                   size_t reaction_index) const {
    const auto& lhs = stoichiometry_lhs_[reaction_index];
    const auto& mol_state_elements = reaction_lhs_[reaction_index];
    osh::Real h_mu = 1.0;
    for (size_t k = 0; k < mol_state_elements.size(); ++k) {
        osh::I64 lhs_s = -lhs[k];
        auto pool_s = mol_state(mol_state_elements[k]);
        if (lhs_s > pool_s) {
            return 0.0;
        }
        switch (lhs_s) {
        case 4: {
            h_mu *= static_cast<osh::Real>(pool_s - 3);
            OMEGA_H_FALLTHROUGH;
        }
        case 3: {
            h_mu *= static_cast<osh::Real>(pool_s - 2);
            OMEGA_H_FALLTHROUGH;
        }
        case 2: {
            h_mu *= static_cast<osh::Real>(pool_s - 1);
            OMEGA_H_FALLTHROUGH;
        }
        case 1: {
            h_mu *= static_cast<osh::Real>(pool_s);
            break;
        }
        default: {
            throw std::runtime_error("Reaction rate computation error");
        }
        }
    }
    return h_mu * ccsts_[reaction_index];
}

//-------------------------------------------------------

template <typename RDefT>
template <PoolChangeType PoolChangeTy>
[[nodiscard]] std::tuple<std::vector<MolStateElementID>,
                         typename SurfaceReactionsBase<RDefT>::Stoichiometry,
                         std::vector<model::region_id>>
SurfaceReactionsBase<RDefT>::reactionMolStateDependencyAndStoichiometry(
    const RDefT& reacdef,
    mesh::triangle_id_t patch_element_id,
    mesh::tetrahedron_id_t inner_compartment_element_id,
    const std::optional<mesh::tetrahedron_id_t>& outer_compartment_element) const {
    using SL = steps::model::Location;
    const size_t size{reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_SURF).size() +
                      reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_IN).size() +
                      reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_OUT).size()};
    std::vector<MolStateElementID> elements_and_species_ids;
    std::vector<osh::I64> stoichiometry;
    /// region_ids contain information about the compartment of the related
    /// species. Needed to compute occupancy dependency.
    std::vector<model::region_id> region_ids;
    elements_and_species_ids.reserve(size);
    stoichiometry.reserve(size);
    region_ids.reserve(size);
    const auto& mol_changes_patch = reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_SURF);
    for (const auto& c: mol_changes_patch) {
        elements_and_species_ids.emplace_back(patch_element_id, c.first);
        stoichiometry.push_back(c.second);
        region_ids.push_back(reacdef.patchdef().getID());
    }
    const auto& mol_changes_inner = reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_IN);
    for (const auto& c: mol_changes_inner) {
        elements_and_species_ids.emplace_back(inner_compartment_element_id, c.first);
        stoichiometry.push_back(c.second);
        region_ids.push_back(reacdef.patchdef().getInnerCompId());
    }
    const auto& mol_changes_outer = reacdef.template getStoichiometry<PoolChangeTy>(SL::PATCH_OUT);
    for (const auto& c: mol_changes_outer) {
        elements_and_species_ids.emplace_back(*outer_compartment_element, c.first);
        stoichiometry.push_back(c.second);
        region_ids.push_back(*reacdef.patchdef().getOuterCompId());
    }
    return std::make_tuple(elements_and_species_ids, stoichiometry, region_ids);
}

//-------------------------------------------------------

template <typename RDefT>
void SurfaceReactionsBase<RDefT>::resetCurrents() {
    osh::fill(currents_, 0.0);
}

//-------------------------------------------------------

template <typename RDefT>
void SurfaceReactionsBase<RDefT>::updateChargeFlow(size_t reaction_index, int nb) {
    const auto& info = this->reacdefs_[reaction_index].get().getInfo();
    currents_[reaction_index] += static_cast<double>(info.charge) * nb * math::E_CHARGE;
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real SurfaceReactionsBase<RDefT>::getCurrent(container::patch_id patch,
                                                  container::triangle_id tri) const {
    osh::Real res = 0;
    size_t start = getInd(patch, tri, ContReacID(0));
    size_t end = getInd(patch, container::triangle_id(tri.get() + 1), ContReacID(0));
    for (size_t i = start; i < end; ++i) {
        res += currents_[i];
    }
    return res;
}

//-------------------------------------------------------

template <typename RDefT>
std::string SurfaceReactionsBase<RDefT>::name() {
    if constexpr (std::is_same_v<RDefT, SReacdef>) {
        return "SurfaceReaction";
    } else if constexpr (std::is_same_v<RDefT, ComplexSReacdef>) {
        return "ComplexSurfaceReaction";
    } else if constexpr (std::is_same_v<RDefT, VDepComplexSReacdef>) {
        return "VDepComplexSurfaceReaction";
    } else if constexpr (std::is_same_v<RDefT, VDepSReacdef>) {
        return "VDepSurfaceReaction";
    } else if constexpr (std::is_same_v<RDefT, GHKSReacdef>) {
        return "GHKSurfaceReaction";
    } else if constexpr (std::is_same_v<RDefT, ComplexGHKSReacdef>) {
        return "ComplexGHKSurfaceReaction";
    } else {
        static_assert(util::always_false_v<RDefT>, "Unmanaged reaction type");
    }
}

//-------------------------------------------------------

template <typename RDefT>
void SurfaceReactionsBase<RDefT>::report(std::ostream& report_stream, size_t reaction_index) const {
    using SR = PoolChangeType;
    const auto& [elmts, stoichiometry, region_id] =
        reactionMolStateDependencyAndStoichiometry<SR::UPD>(
            reacdefs_[reaction_index],
            boundary_id_[reaction_index],
            inner_compartment_element_id_[reaction_index],
            outer_compartment_element_id_[reaction_index]);

    auto f = [&](std::ostream& ostr, const MolStateElementID& el, const model::region_id& region) {
        std::string name = std::visit(
            [&](const auto id) -> std::string {
                return state_def_.getSpecID(
                    state_def_.getDefinition(id).getSpecModelIdx(std::get<1>(el)));
            },
            region);
        switch (std::get<0>(el).index()) {
        case 0:
            ostr << name << "[Tet_" << std::get<mesh::tetrahedron_id_t>(std::get<0>(el)) << ']';
            break;
        case 1:
            ostr << name << "[Tri_" << std::get<mesh::triangle_id_t>(std::get<0>(el)) << ']';
            break;
        default:
            throw std::logic_error("Unknown region");
        }
    };

    report_stream << "Type: " << name() << '\n';
    std::ostringstream stream;
    size_t l{}, m{};
    for (size_t k = 0; k < elmts.size(); ++k) {
        if (stoichiometry[k] <= 0) {
            report_stream << ((m > 0) ? " + " : "")
                          << ((stoichiometry[k] != -1) ? (std::to_string(-stoichiometry[k]) + " * ")
                                                       : "");
            f(report_stream, elmts[k], region_id[k]);
            m++;
        } else {
            stream << ((l == 0) ? "" : " + ")
                   << ((stoichiometry[k] != 1) ? (std::to_string(stoichiometry[k]) + " * ") : "");
            f(stream, elmts[k], region_id[k]);
            l++;
        }
    }
    report_stream << " -> " << stream.str() << '\n';
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real ConstantSurfaceReactions<RDefT>::get_Kcst(container::patch_id patch,
                                                    container::triangle_id tri,
                                                    ContReacID reac) const {
    return get_Kcst(this->getInd(patch, tri, reac));
}

//-------------------------------------------------------

template <typename RDefT>
void ConstantSurfaceReactions<RDefT>::set_Kcst(container::patch_id patch,
                                               container::triangle_id tri,
                                               ContReacID reac,
                                               osh::Real kcst) {
    size_t ind = this->getInd(patch, tri, reac);
    kcsts_[ind] = kcst;
    this->ccsts_[ind] = this->kinConstantGeomFactor_[ind] * kcst;
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real ConstantSurfaceReactions<RDefT>::get_Kcst(size_t ind) const {
    auto it = kcsts_.find(ind);
    if (it != kcsts_.end()) {
        return it->second;
    }
    return this->reacdefs_[ind].get().getInfo().kCst;
}

//-------------------------------------------------------

template <typename RDefT>
void ConstantSurfaceReactions<RDefT>::clear_Kcst(container::patch_id patch) {
    std::vector<size_t> inds;
    for (auto& [ind, _]: kcsts_) {
        if (this->reacdefs_[ind].get().patchdef().getIdx() == patch) {
            inds.emplace_back(ind);
        }
    }
    for (auto ind: inds) {
        kcsts_.erase(ind);
    }
}

//-------------------------------------------------------

template <typename RDefT>
void ConstantSurfaceReactions<RDefT>::updateCcst() {
    for (size_t k = 0; k < this->size(); k++) {
        this->ccsts_[k] = this->kinConstantGeomFactor_[k] * get_Kcst(k);
        // If you reach these asserts probably your mesh is too big/too small or
        // kcst is an invalid value
        assert(this->ccsts_[k] != std::numeric_limits<osh::Real>::infinity());
        assert(this->ccsts_[k] >= 0);
    }
}

//-------------------------------------------------------

template <typename RDefT>
ConstantSurfaceReactions<RDefT>::ConstantSurfaceReactions(const Statedef& statedef,
                                                          DistMesh& mesh,
                                                          MolState& mol_state)
    : SurfaceReactionsBase<RDefT>(statedef, mesh, mol_state) {
    updateCcst();
}

//-------------------------------------------------------

template <typename RDefT>
VDepSurfaceReactionsBase<RDefT>::VDepSurfaceReactionsBase(const Statedef& statedef,
                                                          DistMesh& mesh,
                                                          MolState& mol_state)
    : SurfaceReactionsBase<RDefT>(statedef, mesh, mol_state)
    , tri2verts_(mesh.ask_verts_of(DistMesh::dim() - 1)) {}

//-------------------------------------------------------

template <typename RDefT>
void VDepSurfaceReactionsBase<RDefT>::updatePotential(osh::Reals& potential_on_verts) {
    for (size_t reaction_idx = 0; reaction_idx < this->size(); reaction_idx++) {
        const vdep_propensity_fun_t& f = this->reacdefs_[reaction_idx].get().getInfo().kCstFun;
        const auto& verts =
            Omega_h::gather_verts<DistMesh::dim()>(tri2verts_,
                                                   this->boundary_id_[reaction_idx].get());
        auto avg_potential = std::accumulate(verts.begin(),
                                             verts.end(),
                                             0.0,
                                             [&](const auto acc, const auto h) {
                                                 return acc + potential_on_verts[h];
                                             }) /
                             DistMesh::dim();

        this->ccsts_[reaction_idx] = this->kinConstantGeomFactor_[reaction_idx] * f(avg_potential);

        // If you fail one of these asserts it means that your voltage
        // grew/decreased too much. Probably the mesh is gigantic or minuscule or
        // there is an order of magnitude missmatch and constants do not really make
        // sense
        assert(this->ccsts_[reaction_idx] != std::numeric_limits<double>::infinity());
        assert(this->ccsts_[reaction_idx] >= 0);
    }
}

//-------------------------------------------------------

template <typename RDefT>
void GHKSurfaceReactionsBase<RDefT>::updateChargeFlow(size_t reaction_index, int nb) {
    const GHKInfo& info = this->reacdefs_[reaction_index].get().getInfo();
    SurfaceReactionsBase<RDefT>::updateChargeFlow(reaction_index,
                                                  (2 * static_cast<int>(info.in2out) - 1) * nb);
}

//-------------------------------------------------------

template <typename RDefT>
void GHKSurfaceReactionsBase<RDefT>::updatePotential(osh::Reals& potential_on_verts) {
    for (size_t reaction_idx = 0; reaction_idx < this->size(); reaction_idx++) {
        const auto& verts =
            Omega_h::gather_verts<DistMesh::dim()>(tri2verts_,
                                                   this->boundary_id_[reaction_idx].get());
        auto avg_potential = std::accumulate(verts.begin(),
                                             verts.end(),
                                             0.0,
                                             [&](const auto acc, const auto h) {
                                                 return acc + potential_on_verts[h];
                                             }) /
                             DistMesh::dim();
        potential_on_boundary_[reaction_idx] = avg_potential;
    }
}

//-------------------------------------------------------

template <typename RDefT>
GHKSurfaceReactionsBase<RDefT>::GHKSurfaceReactionsBase(const Statedef& statedef,
                                                        DistMesh& mesh,
                                                        MolState& mol_state)
    : SurfaceReactionsBase<RDefT>(statedef, mesh, mol_state)
    , tri2verts_(mesh.ask_verts_of(DistMesh::dim() - 1))
    , potential_on_boundary_(this->size(), std::numeric_limits<osh::Real>::quiet_NaN()) {
    osh::Write<osh::Real> inner_element_vol_w(this->size()), outer_element_vol_w(this->size());
    std::transform(this->inner_compartment_element_id_.begin(),
                   this->inner_compartment_element_id_.end(),
                   inner_element_vol_w.begin(),
                   [&](const auto el) { return mesh.getTet(el).vol; });
    std::transform(this->outer_compartment_element_id_.begin(),
                   this->outer_compartment_element_id_.end(),
                   outer_element_vol_w.begin(),
                   [&](const auto el) {
                       return el ? mesh.getTet(*el).vol
                                 : std::numeric_limits<osh::Real>::quiet_NaN();
                   });
    inner_element_vol_ = inner_element_vol_w;
    outer_element_vol_ = outer_element_vol_w;
}

//-------------------------------------------------------

template <typename RDefT>
osh::Real GHKSurfaceReactionsBase<RDefT>::ghkRate(const GHKInfo& info,
                                                  size_t index,
                                                  double conc_i,
                                                  double conc_o,
                                                  double nbOpenChan) const {
    // nuFoRT: valence [1] *V [V] * FARADAY [J/(V*mol)]/(R [J/(K*mol)] * T [K]) =
    // [1]
    double nuFoRT = static_cast<double>(info.charge) * potential_on_boundary_[index] *
                    steps::math::FARADAY / (steps::math::GAS_CONSTANT * this->state_def_.getTemp());
    if (nuFoRT >= std::numeric_limits<double>::max_exponent10 * 2.30 ||
        nuFoRT <= std::numeric_limits<double>::min_exponent10 * 2.30) {
        throw std::runtime_error("Overflow encountered, nuFoRT: " + std::to_string(nuFoRT));
    }
    // eNuFoRT: [1]
    double eNuFoRT = std::exp(-nuFoRT);

    // rate: permeability [m/s] * nuFoRT [1] * (conc_i - conc_o * eNuFoRT)
    // [n_molecules/m^3] * n_channels [1] = J [A/m^2]/(valence [1] * Q_charge [C])
    //
    // This comes from a comparison with the GHK current formula on wikipedia. In
    // fact, confronting the 2 equations we have: rate = J [A/m^2] * N_avogadro
    // [mol] /(valence [1] * FARADAY [A/mol]). Since Q_charge [C] * N_avogadro
    // [mol] = FARADAY we have the aforementioned dimensional formula for the
    // rate: J [A/m^2]/(valence [1] * Q_charge [C])
    double rate = (std::abs(nuFoRT) > std::numeric_limits<double>::epsilon())
                      ? info.permeability * nuFoRT * (conc_i - conc_o * eNuFoRT) / (1.0 - eNuFoRT) *
                            nbOpenChan  // number of channels: [1]
                      : info.permeability * (conc_i - conc_o) *
                            nbOpenChan;  // number of channels: [1]

    //  Split result for the 2 twin reactions
    if (info.in2out) {
        return std::max(rate, {});
    } else {
        return -std::min(rate, {});
    }
}


//-------------------------------------------------------

osh::Real GHKSurfaceReactions::computeRate(const MolState& mol_state, size_t index) const {
    const GHKInfo& info = this->reacdefs_[index].get().getInfo();
    size_t idx_in{1}, idx_out{!info.inner_conc ? 2u : 1u};
    const auto& reaction_lhs = this->reaction_lhs_[index];
    // conc_i: [n_molecules/m^3]
    double conc_i = info.inner_conc ? *info.inner_conc
                                    : static_cast<double>(mol_state(reaction_lhs[idx_in])) /
                                          inner_element_vol_[index];
    // conc_o: [n_molecules/m^3]
    double conc_o = info.outer_conc ? *info.outer_conc
                                    : static_cast<double>(mol_state(reaction_lhs[idx_out])) /
                                          outer_element_vol_[index];
    return ghkRate(info, index, conc_i, conc_o, static_cast<double>(mol_state(reaction_lhs[0])));
}

//-------------------------------------------------------

ComplexGHKSurfaceReactions::ComplexGHKSurfaceReactions(const Statedef& statedef,
                                                       DistMesh& mesh,
                                                       MolState& mol_state)
    : GHKSurfaceReactionsBase<ComplexGHKSReacdef>::GHKSurfaceReactionsBase(statedef,
                                                                           mesh,
                                                                           mol_state) {
    for (uint ridx = 0; ridx < this->size(); ++ridx) {
        const auto& reacdef = reacdefs_[ridx].get();
        const auto& tri = this->boundary_id_[ridx];
        auto cId = reacdef.complexID();
        std::vector<MolStateComplexElementID> comp_deps;
        for (const auto& sus: reacdef.complexDEPSET()) {
            comp_deps.emplace_back(tri, cId, sus);
        }
        complex_reactions_deps.push_back(comp_deps);

        surface_candidates.emplace_back(cId, tri);
        surface_candidates.back().addEvent(reacdef.updEvent(),
                                           mol_state.moleculesOnPatchBoundaries());
    }
}

//-------------------------------------------------------

osh::Real ComplexGHKSurfaceReactions::computeRate(const MolState& mol_state, size_t index) const {
    const GHKInfo& info = this->reacdefs_[index].get().getInfo();
    size_t idx_in{0}, idx_out{!info.inner_conc ? 1u : 0u};
    const auto& reaction_lhs = this->reaction_lhs_[index];
    double conc_i = info.inner_conc ? *info.inner_conc
                                    : static_cast<double>(mol_state(reaction_lhs[idx_in])) /
                                          inner_element_vol_[index];
    double conc_o = info.outer_conc ? *info.outer_conc
                                    : static_cast<double>(mol_state(reaction_lhs[idx_out])) /
                                          outer_element_vol_[index];
    return ghkRate(info,
                   index,
                   conc_i,
                   conc_o,
                   surface_candidates[index].rateMult(mol_state.moleculesOnPatchBoundaries()));
}

//-------------------------------------------------------

template <typename surfReacT>
void ComplexSurfaceReactionsBase<surfReacT>::addComplexSurfaceReaction(
    const typename surfReacT::_RDefT& reacdef,
    MolState& mol_state,
    const mesh::tetrahedron_id_t& tet_in,
    const std::optional<mesh::tetrahedron_id_t>& tet_out,
    const mesh::triangle_id_t& tri) {
    std::vector<MolStateComplexElementID> comp_deps;
    const auto& deps = reacdef.complexDEPMAP();

    for (const auto& [cId, substates]: deps.at(steps::model::Location::PATCH_SURF)) {
        for (const auto& sus: substates) {
            comp_deps.emplace_back(tri, cId, sus);
        }
    }
    for (const auto& [cId, substates]: deps.at(steps::model::Location::PATCH_IN)) {
        for (const auto& sus: substates) {
            comp_deps.emplace_back(tet_in, cId, sus);
        }
    }
    if (tet_out.has_value()) {
        for (const auto& [cId, substates]: deps.at(steps::model::Location::PATCH_OUT)) {
            for (const auto& sus: substates) {
                comp_deps.emplace_back(tet_out.value(), cId, sus);
            }
        }
    }
    complex_reactions_deps.push_back(comp_deps);

    std::vector<MolStateComplexElementID> comp_upds;
    const auto& upds = reacdef.complexUPDMAP();

    for (const auto& [cId, substates]: upds.at(steps::model::Location::PATCH_SURF)) {
        for (const auto& sus: substates) {
            comp_upds.emplace_back(tri, cId, sus);
        }
    }
    for (const auto& [cId, substates]: upds.at(steps::model::Location::PATCH_IN)) {
        for (const auto& sus: substates) {
            comp_upds.emplace_back(tet_in, cId, sus);
        }
    }
    if (tet_out.has_value()) {
        for (const auto& [cId, substates]: upds.at(steps::model::Location::PATCH_OUT)) {
            for (const auto& sus: substates) {
                comp_upds.emplace_back(tet_out.value(), cId, sus);
            }
        }
    }
    complex_reactions_upds.push_back(comp_upds);

    {
        surface_candidates.emplace_back();
        auto& lhscands = surface_candidates.back();
        std::map<model::complex_id, uint> added;
        for (auto ev: reacdef.lhsEvents(steps::model::Location::PATCH_SURF)) {
            auto it = added.find(ev->complexIdx());
            if (it == added.end()) {
                it = added.insert({ev->complexIdx(), lhscands.size()}).first;
                lhscands.emplace_back(ev->complexIdx(), tri);
            }
            lhscands[it->second].addEvent(ev, mol_state.moleculesOnPatchBoundaries());
        }
    }
    {
        inner_candidates.emplace_back();
        auto& lhscands = inner_candidates.back();
        std::map<model::complex_id, uint> added;
        for (auto ev: reacdef.lhsEvents(steps::model::Location::PATCH_IN)) {
            auto it = added.find(ev->complexIdx());
            if (it == added.end()) {
                it = added.insert({ev->complexIdx(), lhscands.size()}).first;
                lhscands.emplace_back(ev->complexIdx(), tet_in);
            }
            lhscands[it->second].addEvent(ev, mol_state.moleculesOnElements());
        }
    }
    outer_candidates.emplace_back();
    if (tet_out.has_value()) {
        auto& lhscands = outer_candidates.back();
        std::map<model::complex_id, uint> added;
        for (auto ev: reacdef.lhsEvents(steps::model::Location::PATCH_OUT)) {
            auto it = added.find(ev->complexIdx());
            if (it == added.end()) {
                it = added.insert({ev->complexIdx(), lhscands.size()}).first;
                lhscands.emplace_back(ev->complexIdx(), tet_out.value());
            }
            lhscands[it->second].addEvent(ev, mol_state.moleculesOnElements());
        }
    }
}

template <typename surfReacT>
osh::Real ComplexSurfaceReactionsBase<surfReacT>::computeRate(const MolState& mol_state,
                                                              size_t index) const {
    auto specRate = surfReacT::computeRate(mol_state, index);
    if (specRate == 0.0) {
        return 0;
    }
    // Get the rates for the complex reaction part
    double cmult = 1.0;
    for (auto& cand: surface_candidates[index]) {
        cmult *= cand.rateMult(mol_state.moleculesOnPatchBoundaries());
    }

    if (cmult == 0.0) {
        return 0;
    }
    for (auto& cand: inner_candidates[index]) {
        cmult *= cand.rateMult(mol_state.moleculesOnElements());
    }

    if (cmult == 0.0) {
        return 0;
    }
    for (auto& cand: outer_candidates[index]) {
        cmult *= cand.rateMult(mol_state.moleculesOnElements());
    }

    return cmult * specRate;
}

template <typename surfReacT>
const std::vector<MolStateElementID>&
ComplexSurfaceReactionsBase<surfReacT>::updateMolStateAndOccupancy(MolState& mol_state,
                                                                   rng::RNG& rng,
                                                                   size_t index,
                                                                   osh::Real event_time) const {
    // Species part of the complex reaction
    auto& upd = surfReacT::updateMolStateAndOccupancy(mol_state, rng, index, event_time, 1);

    const auto& tri = this->boundary_id_[index];
    const auto& tet_in = this->inner_compartment_element_id_[index];
    const auto& tet_out = this->outer_compartment_element_id_[index];
    auto& entmols_tets = mol_state.moleculesOnElements();
    auto& entmols_tris = mol_state.moleculesOnPatchBoundaries();
    const auto& creEvents = this->reacdefs_[index].get().creEvents();
    // Updates on surface
    for (auto& cands: surface_candidates[index]) {
        for (auto& event: cands.selectEvents(entmols_tris, rng)) {
            if (event.first->type() == UPDEvent) {
                const auto ev = std::dynamic_pointer_cast<const ComplexUpdateEventdef>(event.first);
                const auto& state = entmols_tris.complexStates(cands.complexIdx()).at(event.second);
                entmols_tris.updateComplexUpdateOccupancy(cands.complexIdx(),
                                                          event.second,
                                                          ev->getUpdate(state, rng),
                                                          event_time);

                switch (ev->destLoc()) {
                case steps::model::Location::PATCH_IN: {
                    entmols_tets.addComplexUpdateOccupancy(
                        tet_in, cands.complexIdx(), event.second, state, event_time);
                    entmols_tris.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              event_time);
                    break;
                }
                case steps::model::Location::PATCH_OUT: {
                    AssertLog(tet_out.has_value());
                    entmols_tets.addComplexUpdateOccupancy(
                        tet_out.value(), cands.complexIdx(), event.second, state, event_time);
                    entmols_tris.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              event_time);
                    break;
                }
                default:
                    break;
                }
            } else {
                // Deletions
                entmols_tris.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                          event.second,
                                                          event_time);
            }
        }
    }
    // Creations on surface
    for (auto& ce: creEvents.at(steps::model::Location::PATCH_SURF)) {
        entmols_tris.createComplexUpdateOccupancy(tri, ce->complexIdx(), ce->init(), event_time);
    }

    // Updates in inner compartment
    for (auto& cands: inner_candidates[index]) {
        for (auto& event: cands.selectEvents(entmols_tets, rng)) {
            if (event.first->type() == UPDEvent) {
                const auto ev = std::dynamic_pointer_cast<const ComplexUpdateEventdef>(event.first);
                const auto& state = entmols_tets.complexStates(cands.complexIdx()).at(event.second);
                entmols_tets.updateComplexUpdateOccupancy(cands.complexIdx(),
                                                          event.second,
                                                          ev->getUpdate(state, rng),
                                                          event_time);

                switch (ev->destLoc()) {
                case steps::model::Location::PATCH_SURF: {
                    entmols_tris.addComplexUpdateOccupancy(
                        tri, cands.complexIdx(), event.second, state, event_time);
                    entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              event_time);
                    break;
                }
                case steps::model::Location::PATCH_OUT: {
                    AssertLog(tet_out.has_value());
                    entmols_tets.addComplexUpdateOccupancy(
                        tet_out.value(), cands.complexIdx(), event.second, state, event_time);
                    entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              event_time);
                    break;
                }
                default:
                    break;
                }
            } else {
                // Deletions
                entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                          event.second,
                                                          event_time);
            }
        }
    }
    // Creations in inner compartment
    for (auto& ce: creEvents.at(steps::model::Location::PATCH_IN)) {
        entmols_tets.createComplexUpdateOccupancy(tet_in, ce->complexIdx(), ce->init(), event_time);
    }

    if (tet_out.has_value()) {
        // Updates in outer compartment
        for (auto& cands: outer_candidates[index]) {
            for (auto& event: cands.selectEvents(entmols_tets, rng)) {
                if (event.first->type() == UPDEvent) {
                    const auto ev = std::dynamic_pointer_cast<const ComplexUpdateEventdef>(
                        event.first);
                    const auto& state =
                        entmols_tets.complexStates(cands.complexIdx()).at(event.second);
                    entmols_tets.updateComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              ev->getUpdate(state, rng),
                                                              event_time);

                    switch (ev->destLoc()) {
                    case steps::model::Location::PATCH_SURF: {
                        entmols_tris.addComplexUpdateOccupancy(
                            tri, cands.complexIdx(), event.second, state, event_time);
                        entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                                  event.second,
                                                                  event_time);
                        break;
                    }
                    case steps::model::Location::PATCH_IN: {
                        entmols_tets.addComplexUpdateOccupancy(
                            tet_in, cands.complexIdx(), event.second, state, event_time);
                        entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                                  event.second,
                                                                  event_time);
                        break;
                    }
                    default:
                        break;
                    }
                } else {
                    // Deletions
                    entmols_tets.removeComplexUpdateOccupancy(cands.complexIdx(),
                                                              event.second,
                                                              event_time);
                }
            }
        }
        // Creations in outer compartment
        for (auto& ce: creEvents.at(steps::model::Location::PATCH_OUT)) {
            entmols_tets.createComplexUpdateOccupancy(tet_out.value(),
                                                      ce->complexIdx(),
                                                      ce->init(),
                                                      event_time);
        }
    }

    // TODO Take complex changes into account for the return value, so RSSA can be used
    return upd;
}

template class SurfaceReactionsBase<SReacdef>;
template class SurfaceReactionsBase<ComplexSReacdef>;
template class SurfaceReactionsBase<VDepComplexSReacdef>;
template class SurfaceReactionsBase<VDepSReacdef>;
template class SurfaceReactionsBase<GHKSReacdef>;
template class SurfaceReactionsBase<ComplexGHKSReacdef>;

template class ConstantSurfaceReactions<SReacdef>;
template class ConstantSurfaceReactions<ComplexSReacdef>;

template class VDepSurfaceReactionsBase<VDepSReacdef>;
template class VDepSurfaceReactionsBase<VDepComplexSReacdef>;

template class GHKSurfaceReactionsBase<GHKSReacdef>;
template class GHKSurfaceReactionsBase<ComplexGHKSReacdef>;

template class ComplexSurfaceReactionsBase<ConstantSurfaceReactions<ComplexSReacdef>>;
template class ComplexSurfaceReactionsBase<VDepSurfaceReactionsBase<VDepComplexSReacdef>>;

}  // namespace steps::dist::kproc

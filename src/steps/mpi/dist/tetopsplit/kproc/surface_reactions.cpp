#include "surface_reactions.hpp"

#include <cassert>
#include <functional>

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_shape.hpp>
#include <boost/optional/optional.hpp>

#include "diffusions.hpp"
#include "kproc_state.hpp"
#include "geom/dist/distmesh.hpp"
#include "math/constants.hpp"
#include "mpi/dist/tetopsplit/definition/sreacdef.hpp"
#include "mpi/dist/tetopsplit/kproc/surface_reactions.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

namespace kproc {

//-------------------------------------------------------

template <typename PropensityType>
template <typename NumMolecules>
SurfaceReactionsBase<PropensityType>::SurfaceReactionsBase(const Statedef& statedef,
                                                           DistMesh& mesh,
                                                           MolState<NumMolecules>& mol_state)
    : state_def_(statedef) {
    if (!statedef.patchdefs().empty()) {
        /// tet_owned informs whether mesh element is owned
        const auto& tet_owned = mesh.owned_elems_mask();
        /// tri2tets mapping a2ab and ab2b from a triangle id to its two neighbouring tets
        const auto tri2tets_a2ab = mesh.bounds2elems_a2ab(mesh.dim() - 1, mesh.dim());
        const auto tri2tets_ab2b = mesh.bounds2elems_ab2b(mesh.dim() - 1, mesh.dim());
        /// and work through all patches
        for (const auto& patch: statedef.patchdefs()) {
            model::patch_id patch_id = patch->getID();
            model::compartment_id inner_compartment_id = patch->getInnerCompId();
            boost::optional<model::compartment_id> outer_compartment_id = patch->getOuterCompId();

            osh::Write<osh::LO> inner_comp_mask(tet_owned.size(), 0);
            for (auto elem: mesh.getEntities(inner_compartment_id)) {
                inner_comp_mask[elem.get()] = 1;
            }

            osh::Write<osh::LO> outer_comp_mask(outer_compartment_id ? tet_owned.size() : 0, 0);
            if (outer_compartment_id) {
                for (auto elem: mesh.getEntities(outer_compartment_id.get())) {
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
                boost::optional<mesh::tetrahedron_id_t> outer_compartment_elmt = boost::none,
                                                        inner_compartment_elmt = boost::none;
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
                if (!tet_owned[inner_compartment_elmt->get()]) {
                    std::ostringstream oss;
                    oss << "SReac : the inner compartment element " << *inner_compartment_elmt
                        << " is not in the same process as the patch " << patch_id;
                    throw std::logic_error(oss.str());
                }
                if (outer_compartment_id && !tet_owned[outer_compartment_elmt->get()]) {
                    std::ostringstream oss;
                    oss << "SReac : the outer compartment element " << *outer_compartment_elmt
                        << " is not in the same process as the patch " << patch_id;
                    throw std::logic_error(oss.str());
                }
                for (const auto& reacdef: patch->getContainer<PropensityType>()) {
                    reacdefs_.push_back(*reacdef);
                    boundary_id_.push_back(boundary);
                    inner_compartment_element_id_.push_back(*inner_compartment_elmt);
                    outer_compartment_element_id_.push_back(outer_compartment_elmt);
                    /// resolve the stoichiometry of the reaction in respect of the mesh.
                    std::vector<MolStateElementID> elmts;
                    Stoichiometry stoichiometry;
                    std::vector<model::region_id> region_id;
                    std::tie(elmts, stoichiometry, region_id) =
                        reactionMolStateDependencyAndStoichiometry<
                            SReacdefBase<PropensityType>::PoolChangeType::LHS>(
                            *reacdef, boundary, *inner_compartment_elmt, outer_compartment_elmt);

                    reaction_lhs_.push_back(elmts);
                    stoichiometry_lhs_.push_back(stoichiometry);

                    /// region_id is necessary to identify the region of species involved
                    /// in the reaction, to check whether the species is diffused.
                    std::tie(elmts, stoichiometry, region_id) =
                        reactionMolStateDependencyAndStoichiometry<
                            SReacdefBase<PropensityType>::PoolChangeType::UPD>(
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
                    // compute the propensity rate constant.
                    ccsts_.push_back(std::numeric_limits<osh::Real>::quiet_NaN());
                }
            }
        }
    }
}

//-------------------------------------------------------

template <typename PropensityType>
template <typename NumMolecules>
void SurfaceReactionsBase<PropensityType>::apply(
    MolState<NumMolecules> &mol_state, size_t reaction_index) const {
  const auto &upd = stoichiometry_upd_[reaction_index];
  const auto &mol_state_elements = reaction_upd_[reaction_index];
  for (size_t k = 0; k < mol_state_elements.size(); ++k) {
    assert(mol_state(mol_state_elements[k]) <=
           std::numeric_limits<NumMolecules>::max() -
               std::max(0, static_cast<osh::LO>(upd[k])));
    assert(mol_state(mol_state_elements[k]) >= -static_cast<osh::LO>(upd[k]));
    mol_state(mol_state_elements[k]) += static_cast<osh::LO>(upd[k]);
  }
}

//-------------------------------------------------------

template <typename PropensityType>
template <typename NumMolecules>
const std::vector<MolStateElementID>&
SurfaceReactionsBase<PropensityType>::updateMolStateAndOccupancy(MolState<NumMolecules>& mol_state,
                                                                 size_t reaction_index,
                                                                 const osh::Real event_time) const {
    const auto& upd = stoichiometry_upd_[reaction_index];
    const auto& mol_state_elements = reaction_upd_[reaction_index];

    for (size_t k = 0; k < mol_state_elements.size(); ++k) {
        const auto& elmt_id = mol_state_elements[k];
        assert(mol_state(elmt_id) <= std::numeric_limits<NumMolecules>::max() -
                                         std::max(static_cast<NumMolecules>(upd[k]), {}));
        assert(mol_state(elmt_id) >= -static_cast<NumMolecules>(upd[k]));
        mol_state.add_and_update_occupancy(elmt_id, static_cast<NumMolecules>(upd[k]), event_time);
    }
    return mol_state_elements;
}

//-------------------------------------------------------

template <typename PropensityType>
osh::Real SurfaceReactionsBase<PropensityType>::kinConstantGeomFactor(
    const DistMesh &mesh, size_t reaction_index) const {
  const SReacdefBase<PropensityType> &reacdef = reacdefs_[reaction_index];
  const mesh::triangle_id_t boundary{boundary_id_[reaction_index]};
  const mesh::tetrahedron_id_t inner_compartment_element{
      inner_compartment_element_id_[reaction_index]};
  const boost::optional<mesh::tetrahedron_id_t> &outer_compartment_element =
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
template <typename PropensityType>
template <typename NumMolecules>
osh::Real SurfaceReactionsBase<PropensityType>::computeRate(
    const MolState<NumMolecules> &mol_state, size_t reaction_index) const {
  const auto &lhs = stoichiometry_lhs_[reaction_index];
  const auto &mol_state_elements = reaction_lhs_[reaction_index];
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

template <typename PropensityType>
template <typename SReacdefBase<PropensityType>::PoolChangeType PoolChangeTy>
std::tuple<std::vector<MolStateElementID>,
           typename SurfaceReactionsBase<PropensityType>::Stoichiometry,
           std::vector<model::region_id>>
SurfaceReactionsBase<PropensityType>::
    reactionMolStateDependencyAndStoichiometry(
        const SReacdefBase<PropensityType> &reacdef,
        mesh::triangle_id_t patch_element_id,
        mesh::tetrahedron_id_t inner_compartment_element_id,
        const boost::optional<mesh::tetrahedron_id_t>
            &outer_compartment_element) const {
  using SL = typename SReacdefBase<PropensityType>::SpecieLocation;
  const size_t size{
      reacdef.template getStoichiometry<PoolChangeTy, SL::Patch>().size() +
      reacdef.template getStoichiometry<PoolChangeTy, SL::InnerCompartment>()
          .size() +
      reacdef.template getStoichiometry<PoolChangeTy, SL::OuterCompartment>()
          .size()};
  std::vector<MolStateElementID> elements_and_species_ids;
  std::vector<osh::I64> stoichiometry;
  /// region_ids contain information about the compartment of the related
  /// species. Needed to compute occupancy dependency.
  std::vector<model::region_id> region_ids;
  elements_and_species_ids.reserve(size);
  stoichiometry.reserve(size);
  region_ids.reserve(size);
  const auto &mol_changes_patch =
      reacdef.template getStoichiometry<PoolChangeTy, SL::Patch>();
  for (const auto &c : mol_changes_patch) {
    elements_and_species_ids.emplace_back(patch_element_id, c.first);
    stoichiometry.push_back(c.second);
    region_ids.push_back(reacdef.patchdef().getID());
  }
  const auto &mol_changes_inner =
      reacdef.template getStoichiometry<PoolChangeTy, SL::InnerCompartment>();
  for (const auto &c : mol_changes_inner) {
    elements_and_species_ids.emplace_back(inner_compartment_element_id,
                                          c.first);
    stoichiometry.push_back(c.second);
    region_ids.push_back(reacdef.patchdef().getInnerCompId());
  }
  const auto &mol_changes_outer =
      reacdef.template getStoichiometry<PoolChangeTy, SL::OuterCompartment>();
  for (const auto &c : mol_changes_outer) {
    elements_and_species_ids.emplace_back(*outer_compartment_element, c.first);
    stoichiometry.push_back(c.second);
    region_ids.push_back(*reacdef.patchdef().getOuterCompId());
  }
  return std::make_tuple(elements_and_species_ids, stoichiometry, region_ids);
}

//-------------------------------------------------------

template <> std::string SurfaceReactionsBase<SReacInfo>::name() {
  return "SurfaceReaction";
}

template <> std::string SurfaceReactionsBase<VDepInfo>::name() {
  return "VDepSurfaceReaction";
}

template <> std::string SurfaceReactionsBase<GHKInfo>::name() {
  return "GHKSurfaceReaction";
}

//-------------------------------------------------------

template <typename PropensityType>
void SurfaceReactionsBase<PropensityType>::report(std::ostream &report_stream,
                                                  size_t reaction_index) const {
  std::vector<MolStateElementID> elmts;
  Stoichiometry stoichiometry;
  std::vector<model::region_id> region_id;
  using SR = typename SReacdefBase<PropensityType>::PoolChangeType;
  std::tie(elmts, stoichiometry, region_id) =
      reactionMolStateDependencyAndStoichiometry<SR::UPD>(
          reacdefs_[reaction_index], boundary_id_[reaction_index],
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
                    << ((stoichiometry[k] != -1)
                            ? (std::to_string(-stoichiometry[k]) + " * ")
                            : "");
      f(report_stream, elmts[k], region_id[k]);
      m++;
    } else {
      stream << ((l == 0) ? "" : " + ")
             << ((stoichiometry[k] != 1)
                     ? (std::to_string(stoichiometry[k]) + " * ")
                     : "");
      f(stream, elmts[k], region_id[k]);
      l++;
    }
  }
  report_stream << " -> " << stream.str() << '\n';
}

//-------------------------------------------------------

void SurfaceReactions::updateCcst(DistMesh &mesh) {
  for (size_t k = 0; k < size(); k++) {
    ccsts_[k] =
        kinConstantGeomFactor(mesh, k) * reacdefs_[k].get().getInfo().kCst;
  }
}

//-------------------------------------------------------
template <typename NumMolecules>
SurfaceReactions::SurfaceReactions(const Statedef& statedef,
                                   DistMesh& mesh,
                                   MolState<NumMolecules>& mol_state)
    : SurfaceReactionsBase<SReacInfo>(statedef, mesh, mol_state) {
    updateCcst(mesh);
}

//-------------------------------------------------------
template <typename NumMolecules>
VDepSurfaceReactions::VDepSurfaceReactions(const Statedef& statedef,
                                           DistMesh& mesh,
                                           MolState<NumMolecules>& mol_state)
    : SurfaceReactionsBase<VDepInfo>(statedef, mesh, mol_state)
    , kinConstantGeomFactor_(reacdefs_.size())
    , tri2verts_(mesh.ask_verts_of(DistMesh::dim() - 1)) {
    for (size_t reaction_idx = 0; reaction_idx < size(); reaction_idx++) {
        kinConstantGeomFactor_[reaction_idx] = kinConstantGeomFactor(mesh, reaction_idx);
    }
}

//-------------------------------------------------------

void VDepSurfaceReactions::kCstUpdate(osh::Reals &potential_on_verts) {
  for (size_t reaction_idx = 0; reaction_idx < size(); reaction_idx++) {
    const vdep_propensity_fun_t &f =
        reacdefs_[reaction_idx].get().getInfo().kCstFun;
    const auto &verts = Omega_h::gather_verts<DistMesh::dim()>(
        tri2verts_, boundary_id_[reaction_idx].get());
    auto avg_potential = std::accumulate(verts.begin(), verts.end(), 0.0,
                                         [&](const auto acc, const auto h) {
                                           return acc + potential_on_verts[h];
                                         }) /
                         DistMesh::dim();
    ccsts_[reaction_idx] =
        kinConstantGeomFactor_[reaction_idx] * f(avg_potential);
  }
}

//-------------------------------------------------------

template <typename NumMolecules>
osh::Real
GHKSurfaceReactions::computeRate(const MolState<NumMolecules> &mol_state,
                                 size_t index) const {
  const GHKInfo &info = reacdefs_[index].get().getInfo();
  size_t idx_in{1}, idx_out{!info.inner_conc ? 2u : 1u};
  const auto &reaction_lhs = reaction_lhs_[index];
  // conc_i: [n_molecules/m^3]
  double conc_i = info.inner_conc
                      ? *info.inner_conc
                      : static_cast<double>(mol_state(reaction_lhs[idx_in])) /
                            inner_element_vol_[index];
  // conc_o: [n_molecules/m^3]
  double conc_o = info.outer_conc
                      ? *info.outer_conc
                      : static_cast<double>(mol_state(reaction_lhs[idx_out])) /
                            outer_element_vol_[index];
  // nuFoRT: valence [1] *V [V] * FARADAY [J/(V*mol)]/(R [J/(K*mol)] * T [K]) =
  // [1]
  double nuFoRT = static_cast<double>(info.valence) * potential_on_boundary_[index] *
                  steps::math::FARADAY / (steps::math::GAS_CONSTANT * state_def_.getTemp());
  if (nuFoRT >= std::numeric_limits<double>::max_exponent10 * 2.30 ||
      nuFoRT <= std::numeric_limits<double>::min_exponent10 * 2.30) {
    throw std::runtime_error("Overflow encountered, nuFoRT: " +
                             std::to_string(nuFoRT));
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
  double rate =
      (std::abs(nuFoRT) > std::numeric_limits<double>::epsilon())
          ? info.permeability * nuFoRT * (conc_i - conc_o * eNuFoRT) /
                (1.0 - eNuFoRT) *
                static_cast<double>(
                    mol_state(reaction_lhs[0])) // number of channels: [1]
          : info.permeability * (conc_i - conc_o) *
                static_cast<double>(
                    mol_state(reaction_lhs[0])); // number of channels: [1]

  //  Split result for the 2 twin reactions
  if (info.in2out) {
    return std::max(rate, {});
  } else {
    return -std::min(rate, {});
  }
}

//-------------------------------------------------------

void GHKSurfaceReactions::resetCurrents() {
  std::fill(currents_.begin(), currents_.end(), 0.0);
}

//-------------------------------------------------------

void GHKSurfaceReactions::updateChargeFlow(size_t reaction_index) {
    const GHKInfo& info = reacdefs_[reaction_index].get().getInfo();
    currents_[reaction_index] += (2 * static_cast<int>(info.in2out) - 1) *
                                 static_cast<double>(info.valence) * math::E_CHARGE;
}

//-------------------------------------------------------

void GHKSurfaceReactions::updatePotential(osh::Reals &potential_on_verts) {
  for (size_t reaction_idx = 0; reaction_idx < size(); reaction_idx++) {
    const auto &verts = Omega_h::gather_verts<DistMesh::dim()>(
        tri2verts_, boundary_id_[reaction_idx].get());
    auto avg_potential = std::accumulate(verts.begin(), verts.end(), 0.0,
                                         [&](const auto acc, const auto h) {
                                           return acc + potential_on_verts[h];
                                         }) /
                         DistMesh::dim();
    potential_on_boundary_[reaction_idx] = avg_potential;
  }
}

//-------------------------------------------------------
template <typename NumMolecules>
GHKSurfaceReactions::GHKSurfaceReactions(const Statedef& statedef,
                                         DistMesh& mesh,
                                         MolState<NumMolecules>& mol_state)
    : SurfaceReactionsBase<GHKInfo>(statedef, mesh, mol_state)
    , tri2verts_(mesh.ask_verts_of(DistMesh::dim() - 1))
    , currents_(size(), std::numeric_limits<osh::Real>::quiet_NaN())
    , potential_on_boundary_(size(), std::numeric_limits<osh::Real>::quiet_NaN()) {
    osh::Write<osh::Real> inner_element_vol_w(size()), outer_element_vol_w(size());
    std::transform(inner_compartment_element_id_.begin(),
                   inner_compartment_element_id_.end(),
                   inner_element_vol_w.begin(),
                   [&](const auto el) { return mesh.getTet(el.get()).vol; });
    std::transform(outer_compartment_element_id_.begin(),
                   outer_compartment_element_id_.end(),
                   outer_element_vol_w.begin(),
                   [&](const auto el) {
                       return el ? mesh.getTet(el->get()).vol
                                 : std::numeric_limits<osh::Real>::quiet_NaN();
                   });
    inner_element_vol_ = inner_element_vol_w;
    outer_element_vol_ = outer_element_vol_w;

    // Fill the curr2tri2reac_ mapping
    for (uint ridx = 0; ridx < size(); ++ridx) {
        const auto &curr_id = reacdefs_[ridx].get().getInfo().curr_id;
        const auto &tri_id = boundary_id_[ridx];
        auto tri2reac_it = curr2tri2reac_.find(curr_id);
        if (tri2reac_it == curr2tri2reac_.end()) {
            tri2reac_it =
                curr2tri2reac_
                    .emplace(curr_id,
                             osh::Write<osh::GO>(
                                 mesh.owned_bounds_mask().size() * rpt(), -1))
                    .first;
        }
        bool added{false};
        for (uint i = 0; i < rpt(); ++i) {
            auto &val = tri2reac_it->second[tri_id.get() * rpt() + i];
            if (val == -1) {
                val = ridx;
                added = true;
                break;
            }
        }
        if (not added) {
            std::ostringstream msg;
            msg << "Expected " << rpt()
                << " but got a higher number of GHK reactions associated with "
                   "GHK current "
                << curr_id << " in triangle " << tri_id.get();
            throw std::logic_error(msg.str());
        }
    }
}

//-------------------------------------------------------

const osh::Write<osh::GO> &
GHKSurfaceReactions::getTri2Curr(const model::ghk_current_id &curr_id) const {
    const auto curr_it = curr2tri2reac_.find(curr_id);
    if (curr_it == curr2tri2reac_.end()) {
        std::ostringstream msg;
        msg << "GHK current " << curr_id << " is not defined.";
        throw std::logic_error(msg.str());
    }
    return curr_it->second;
}

//-------------------------------------------------------

// explicit template instantiation declarations
template osh::Real
SurfaceReactionsBase<SReacInfo>::computeRate(const MolState<osh::LO> &mol_state,
                                             size_t index) const;
template osh::Real
SurfaceReactionsBase<VDepInfo>::computeRate(const MolState<osh::LO> &mol_state,
                                            size_t index) const;
template osh::Real
SurfaceReactionsBase<SReacInfo>::computeRate(const MolState<osh::GO> &mol_state,
                                             size_t index) const;
template osh::Real
SurfaceReactionsBase<VDepInfo>::computeRate(const MolState<osh::GO> &mol_state,
                                            size_t index) const;
template osh::Real
GHKSurfaceReactions::computeRate(const MolState<osh::LO> &mol_state,
                                 size_t index) const;
template osh::Real
GHKSurfaceReactions::computeRate(const MolState<osh::GO> &mol_state,
                                 size_t index) const;
template osh::Real
SurfaceReactionsBase<SReacInfo>::kinConstantGeomFactor(const DistMesh &mesh,
                                                       size_t index) const;

template const std::vector<MolStateElementID>&
SurfaceReactionsBase<SReacInfo>::updateMolStateAndOccupancy(MolState<osh::I32>& mol_state,
                                                            size_t index,
                                                            const osh::Real event_time) const;
template const std::vector<MolStateElementID>&
SurfaceReactionsBase<VDepInfo>::updateMolStateAndOccupancy(MolState<osh::I32>& mol_state,
                                                           size_t index,
                                                           const osh::Real event_time) const;
template const std::vector<MolStateElementID>&
SurfaceReactionsBase<SReacInfo>::updateMolStateAndOccupancy(MolState<osh::I64>& mol_state,
                                                            size_t index,
                                                            const osh::Real event_time) const;
template const std::vector<MolStateElementID>&
SurfaceReactionsBase<VDepInfo>::updateMolStateAndOccupancy(MolState<osh::I64>& mol_state,
                                                           size_t index,
                                                           const osh::Real event_time) const;
template const std::vector<MolStateElementID>&
SurfaceReactionsBase<GHKInfo>::updateMolStateAndOccupancy(MolState<osh::I32>& mol_state,
                                                          size_t index,
                                                          const osh::Real event_time) const;
template const std::vector<MolStateElementID>&
SurfaceReactionsBase<GHKInfo>::updateMolStateAndOccupancy(MolState<osh::I64>& mol_state,
                                                          size_t index,
                                                          const osh::Real event_time) const;

template void
SurfaceReactionsBase<VDepInfo>::report(std::ostream &report_stream,
                                       size_t index) const;
template void
SurfaceReactionsBase<SReacInfo>::report(std::ostream &report_stream,
                                        size_t index) const;
template void SurfaceReactionsBase<GHKInfo>::report(std::ostream &report_stream,
                                                    size_t index) const;

template SurfaceReactionsBase<SReacInfo>::SurfaceReactionsBase(const Statedef& statedef,
                                                               DistMesh& mesh,
                                                               MolState<osh::LO>& mol_state);
template SurfaceReactionsBase<SReacInfo>::SurfaceReactionsBase(const Statedef& statedef,
                                                               DistMesh& mesh,
                                                               MolState<osh::GO>& mol_state);

template SurfaceReactions::SurfaceReactions(const Statedef& statedef,
                                            DistMesh& mesh,
                                            MolState<osh::LO>& mol_state);
template SurfaceReactions::SurfaceReactions(const Statedef& statedef,
                                            DistMesh& mesh,
                                            MolState<osh::GO>& mol_state);

template VDepSurfaceReactions::VDepSurfaceReactions(const Statedef& statedef,
                                                    DistMesh& mesh,
                                                    MolState<osh::LO>& mol_state);
template VDepSurfaceReactions::VDepSurfaceReactions(const Statedef& statedef,
                                                    DistMesh& mesh,
                                                    MolState<osh::GO>& mol_state);

template GHKSurfaceReactions::GHKSurfaceReactions(const Statedef& statedef,
                                                  DistMesh& mesh,
                                                  MolState<osh::LO>& mol_state);
template GHKSurfaceReactions::GHKSurfaceReactions(const Statedef& statedef,
                                                  DistMesh& mesh,
                                                  MolState<osh::GO>& mol_state);


} // namespace kproc
} // namespace dist
} // namespace steps

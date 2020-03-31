#include <boost/optional/optional.hpp>
#include <cassert>

#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_shape.hpp>
#include <hadoken/format/format.hpp>

#include "opsplit/sreacdef.hpp"
#include "opsplit/vocabulary.hpp"
#include "surface_reactions.hpp"


namespace zee {

namespace kproc {

//-------------------------------------------------------

template <osh::Int Dim>
SurfaceReactions::SurfaceReactions(const Statedef& statedef, OmegaHMesh<Dim>& mesh, bool discovery)
    : num_species_per_boundary_({}) {
    // Retrieve the class ids of all tets
    if (statedef.patchdefs().size() > 0) {
        // initialize a vector to record the number of species owned by a patch element and owned by
        // the process
        osh::Write<osh::LO> num_species_per_boundary_element(mesh.getMesh().nents(Dim - 1), 0);
        // tet_owned informs whether mesh element is owned
        const auto& tet_owned = mesh.getOwnedElemsMask();
        // tri2tets is a mapping from a triangle id to its two neighbouring tets
        auto tri2tets = mesh.getMesh().ask_up(Dim - 1, Dim);
        // and work through all patches
        for (const auto& patch: statedef.patchdefs()) {
            // fetch all elements on the patch and look up
            // neighbouring tets
            model::patch_id patch_id = patch->getID();
            model::compartment_id inner_compartment_id = patch->getInnerCompId();
            boost::optional<model::compartment_id> outer_compartment_id = patch->getOuterCompId();

            // class pairs is a vector containing class ids associated with mesh elements.
            // This vector identifies elements of a given physical entity.
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

            // patch_elements_owned is a vector of patch boundary ids owned by the process.
            const auto& patch_owned_boundaries = mesh.getOwnedEntities(patch_id);

            for (const auto boundary: patch_owned_boundaries) {
                // record the number of species on that element
                num_species_per_boundary_element[boundary.get()] = patch->getNSpecs();
                if (discovery) {
                    continue;
                }
                // adjacent tet ids of patch element 'element' are stored in tritets.a2ab
                osh::LO noffsets = tri2tets.a2ab[boundary.get() + 1] -
                                   tri2tets.a2ab[boundary.get()];
                if (noffsets > 2) {
                    throw std::logic_error(
                        hadoken::scat("Geometry: wrong number of adjacent elements of dimension ",
                                      Dim,
                                      " for boundary ",
                                      boundary,
                                      " of dimension ",
                                      Dim - 1,
                                      ": ",
                                      noffsets));
                }
                boost::optional<mesh::element_id> outer_compartment_elmt = boost::none,
                                                  inner_compartment_elmt = boost::none;
                // work out the compartment associated with adjacent tets
                for (osh::LO k = 0; k < noffsets; ++k) {
                    auto tet_id = tri2tets.ab2b[tri2tets.a2ab[boundary.get()] + k];
                    // check whether it is in the inner compartment by working through its class ids
                    if (inner_comp_mask[tet_id] == 1) {
                        inner_compartment_elmt = mesh::element_id(tet_id);
                    } else {
                        // it should be associated with the outer compartment
                        if (outer_compartment_id) {
                            if (outer_comp_mask[tet_id] == 0) {
                                throw std::logic_error(hadoken::scat(
                                    "Geometry : the outer compartment cannot be resolved "
                                    "for patch ",
                                    patch_id));
                            }
                        }
                        outer_compartment_elmt = mesh::element_id(tet_id);
                    }
                }
                // inner comp needs to be resolved by now.
                if (!inner_compartment_elmt) {
                    throw std::logic_error(hadoken::scat("Inner compartment of element ",
                                                         boundary,
                                                         " cannot be resolved."));
                }
                // are these tets owned by the process?
                if (!tet_owned[inner_compartment_elmt->get()]) {
                    throw std::logic_error(
                        hadoken::scat("SReac : the inner compartment element ",
                                      *inner_compartment_elmt,
                                      " is not in the same process as the patch ",
                                      patch_id));
                }
                if (outer_compartment_id && !tet_owned[outer_compartment_elmt->get()]) {
                    throw std::logic_error(
                        hadoken::scat("SReac : the outer compartment element ",
                                      *outer_compartment_elmt,
                                      " is not in the same process as the patch ",
                                      patch_id));
                }

                for (const auto& reacdef: patch->reacdefs()) {
                    reacdefs_.push_back(*reacdef);
                    inner_compartment_element_id_.push_back(*inner_compartment_elmt);
                    outer_compartment_element_id_.push_back(outer_compartment_elmt);
                    // resolve the stoichiometry of the reaction in respect of the mesh.
                    std::vector<MolState::ElementID> elmts;
                    std::vector<PetscInt> stoichiometry;
                    std::vector<model::region_id> region_id;
                    std::tie(elmts, stoichiometry, region_id) =
                        reactionMolStateDependencyAndStoichiometry<SReacdef::PoolChangeType::LHS>(
                            *reacdef, boundary, *inner_compartment_elmt, outer_compartment_elmt);

                    reaction_lhs_.push_back(elmts);
                    stoichiometry_lhs_.push_back(stoichiometry);

                    // region_id is necessary to identify the region of species involved in the
                    // reaction, to check whether the specie is diffused.
                    std::tie(elmts, stoichiometry, region_id) =
                        reactionMolStateDependencyAndStoichiometry<SReacdef::PoolChangeType::UPD>(
                            *reacdef, boundary, *inner_compartment_elmt, outer_compartment_elmt);
                    reaction_upd_.push_back(elmts);
                    stoichiometry_upd_.push_back(stoichiometry);
                    region_ids_.push_back(region_id);

                    // compute the propensity rate constant.
                    ccsts_.push_back(compute_ccst(
                        mesh, *reacdef, boundary, *inner_compartment_elmt, outer_compartment_elmt));
                    rates_.push_back(0.);
                }
            }
        }
        num_species_per_boundary_ = num_species_per_boundary_element;
    }
}

//-------------------------------------------------------

void SurfaceReactions::apply(MolState& mol_state, size_t index) const {
    const auto& upd = stoichiometry_upd_[index];
    const auto& mol_state_elements = reaction_upd_[index];
    for (size_t k = 0; k < mol_state_elements.size(); ++k) {
        assert(mol_state(mol_state_elements[k]) <= std::numeric_limits<MolState::elem_type>::max() -
                                                       std::max(0, static_cast<osh::LO>(upd[k])) &&
               mol_state(mol_state_elements[k]) >= -static_cast<osh::LO>(upd[k]));
        mol_state(mol_state_elements[k]) += static_cast<osh::LO>(upd[k]);
    }
}

//-------------------------------------------------------

PetscScalar SurfaceReactions::updateRate(const MolState& mol_state, size_t index) {
    rates_[index] = computeRate(mol_state, index);
    return rates_[index];
}

//-------------------------------------------------------

PetscScalar SurfaceReactions::computeRate(const MolState& mol_state, size_t index) const {
    const auto& lhs = stoichiometry_lhs_[index];
    const auto& mol_state_elements = reaction_lhs_[index];
    PetscScalar h_mu = 1.0;
    for (size_t k = 0; k < mol_state_elements.size(); ++k) {
        PetscInt lhs_s = -lhs[k];
        auto pool_s = mol_state(mol_state_elements[k]);
        if (lhs_s > pool_s) {
            return 0.0;
        }
        switch (lhs_s) {
        case 4: {
            h_mu *= (pool_s - 3);
            OMEGA_H_FALLTHROUGH;
        }
        case 3: {
            h_mu *= (pool_s - 2);
            OMEGA_H_FALLTHROUGH;
        }
        case 2: {
            h_mu *= (pool_s - 1);
            OMEGA_H_FALLTHROUGH;
        }
        case 1: {
            h_mu *= pool_s;
            break;
        }
        default: {
            throw std::runtime_error("Reaction rate computation error");
        }
        }
    }
    return h_mu * ccsts_[index];
}

//-------------------------------------------------------

template <osh::Int Dim>
PetscScalar SurfaceReactions::compute_ccst(
    OmegaHMesh<Dim>& mesh,
    const SReacdef& reacdef,
    mesh::boundary_id boundary,
    const mesh::element_id inner_compartment_element,
    const boost::optional<mesh::element_id>& outer_compartment_element) {
    if (reacdef.isSurfaceSurfaceReaction()) {
        // the reaction happens on the surface of the patch.
        osh::Reals area = osh::measure_ents_real(&mesh.getMesh(),
                                                 Dim - 1,
                                                 {boundary.get()},
                                                 mesh.getMesh().coords());
        double ascale = area[0] * AVOGADRO;
        PetscInt o1 = reacdef.getOrder() - 1;
        return reacdef.getKcst() * pow(ascale, static_cast<PetscScalar>(-o1));
    } else {
        osh::Reals vol;
        if (reacdef.isInnerCompartmentReaction()) {
            vol = osh::measure_ents_real(&mesh.getMesh(),
                                         Dim,
                                         {inner_compartment_element.get()},
                                         mesh.getMesh().coords());
        } else {
            vol = osh::measure_ents_real(&mesh.getMesh(),
                                         Dim,
                                         {outer_compartment_element->get()},
                                         mesh.getMesh().coords());
        }
        double vscale = 1.0e3 * vol[0] * AVOGADRO;
        PetscInt o1 = reacdef.getOrder() - 1;
        return reacdef.getKcst() * std::pow(vscale, static_cast<PetscScalar>(-o1));
    }
}

//-------------------------------------------------------

template <SReacdef::PoolChangeType PoolChangeTy>
std::tuple<std::vector<MolState::ElementID>,
           SurfaceReactions::Stoichiometry,
           std::vector<model::region_id>>
SurfaceReactions::reactionMolStateDependencyAndStoichiometry(
    const SReacdef& reacdef,
    mesh::boundary_id patch_element_id,
    mesh::element_id inner_compartment_element_id,
    const boost::optional<mesh::element_id>& outer_compartment_element) const {
    size_t size =
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation ::Patch>().size() +
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation ::InnerCompartment>()
            .size() +
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation ::OuterCompartment>()
            .size();
    std::vector<MolState::ElementID> elements_and_species_ids;
    std::vector<PetscInt> stoichiometry;
    // region_ids contain information about the compartment of the related specie.
    // Needed to compute occupancy dependency.
    std::vector<model::region_id> region_ids;
    elements_and_species_ids.reserve(size);
    stoichiometry.reserve(size);
    region_ids.reserve(size);
    const auto& mol_changes_patch =
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation::Patch>();
    for (const auto& c: mol_changes_patch) {
        elements_and_species_ids.push_back(MolState::mkBoundaryElement(patch_element_id, c.first));
        stoichiometry.push_back(c.second);
        region_ids.push_back(reacdef.patchdef().getID());
    }
    const auto& mol_changes_inner =
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation::InnerCompartment>();
    for (const auto& c: mol_changes_inner) {
        elements_and_species_ids.push_back(
            MolState::mkVolumeElement(inner_compartment_element_id, c.first));
        stoichiometry.push_back(c.second);
        region_ids.push_back(reacdef.patchdef().getInnerCompId());
    }
    const auto& mol_changes_outer =
        reacdef.getStoichiometry<PoolChangeTy, SReacdef::SpecieLocation::OuterCompartment>();
    for (const auto& c: mol_changes_outer) {
        elements_and_species_ids.push_back(
            MolState::mkVolumeElement(*outer_compartment_element, c.first));
        stoichiometry.push_back(c.second);
        region_ids.push_back(*reacdef.patchdef().getOuterCompId());
    }
    return std::make_tuple(elements_and_species_ids, stoichiometry, region_ids);
}


// explicit template instantiation declarations
template SurfaceReactions::SurfaceReactions(const zee::Statedef& statedef,
                                            zee::OmegaHMesh<2>& mesh,
                                            bool discovery);
template SurfaceReactions::SurfaceReactions(const zee::Statedef& statedef,
                                            zee::OmegaHMesh<3>& mesh,
                                            bool discovery);

};  // namespace kproc

};  // namespace zee

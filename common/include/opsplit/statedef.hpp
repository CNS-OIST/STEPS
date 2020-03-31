#pragma once

#include <boost/optional.hpp>
#include <cassert>
#include <map>
#include <memory>
#include <vector>

#include <petscsys.h>

#include "fwd.hpp"
#include "vocabulary.hpp"


namespace zee {

/**
 * \brief State definition of the biochemical container.
 *
 * The Statedef class provide the state definition and indexing system
 * of a biochemical container, such as reaction and diffusion for each
 * geometry point.
 * Note that the compartment/patch indexing is different from the
 * label id of the compartment/patch in the DistMesh class.
 * As in STEPS, biochemical container definition is isolated from
 * geometry, and the association is resolved when the simulation solver
 * is created.
 * This class corresponds to the solver::Statedef class in STEPS.
 */


class Statedef {
  public:
    /**
     * Add the species to the biochemical state definition and return its container index.
     * If the species has been added before, return its container index in record,
     * otherwise add the species to the record and return its new container index.
     */
    model::specie_id addSpec(const model::specie_name& name);
    model::specie_id getSpecModelIdx(const model::specie_name& name) const;
    container::compartment_id addComp(const model::compartment_id& compartment);
    container::patch_id addPatch(
        const model::patch_id& patchId,
        const model::compartment_id& inner_compartment_id,
        const boost::optional<model::compartment_id>& outer_compartment_id = boost::none);
    container::compartment_id getCompModelIdx(const model::compartment_id& compartment) const
        noexcept;

    /**
     * Add the species to the compartment definition and return its comparmental index.
     * If the species has been added before, return its lidx in record,
     * otherwise add the species to the record and return its new lidx.
     */
    model::specie_id addCompSpec(const model::compartment_id& compartment,
                                 const model::specie_name& specie);

    std::vector<model::specie_id> addCompSpecs(const model::compartment_id& compartment,
                                               const std::vector<model::specie_name>& species);
    std::vector<model::specie_id> addPatchSpecs(const model::patch_id& patchId,
                                                const std::vector<model::specie_name>& species);
    model::specie_id addPatchSpec(const model::patch_id& patch_id,
                                  const model::specie_name& specie);
    container::specie_id getCompSpecContainerIdx(const model::compartment_id& compartment,
                                                 const model::specie_name& specie) const;

    /**
     * Register a diffusion in a given compartment
     * \param compartment the compartment identifier where to register the diffusion
     * \param specie_name the diffusing chemical specie
     * \param dcst diffusion constant
     * \return the diffusion identifier
     */
    container::diffusion_id addCompDiff(const model::compartment_id& compartment,
                                        const model::specie_name& specie_name,
                                        PetscScalar dcst);

    /**
     * Register a chemical reaction to a compartment
     * \param compartment the compartment identifier
     * \param reactants substances initially involved in the chemical reaction
     * \param products substances yield by the chemical reaction
     * \param kcst the reaction constant
     * \return the reaction identifier
     */
    container::reaction_id addCompReac(const model::compartment_id& compartment,
                                       const std::vector<model::specie_name>& reactants,
                                       const std::vector<model::specie_name>& products,
                                       PetscScalar kcst);

    container::surface_reaction_id addSurfReac(const model::patch_id& patchId,
                                               const std::vector<model::specie_name>& reactants_i,
                                               const std::vector<model::specie_name>& reactants_s,
                                               const std::vector<model::specie_name>& reactants_o,
                                               const std::vector<model::specie_name>& products_i,
                                               const std::vector<model::specie_name>& products_s,
                                               const std::vector<model::specie_name>& products_o,
                                               PetscScalar kcst);

    Compdef& getCompdef(container::compartment_id compartment) const;

    Compdef& getCompdef(const model::compartment_id& compartment) const;

    inline Compdef& getDefinition(const model::compartment_id& compartmentId) const {
        return getCompdef(compartmentId);
    }

    Patchdef& getPatchdef(const container::patch_id& patchId) const;

    Patchdef& getPatchdef(const model::patch_id& patchId) const;

    inline Patchdef& getDefinition(const model::patch_id& patchId) const {
        return getPatchdef(patchId);
    }

    inline PetscInt getNComps() const noexcept {
        return static_cast<PetscInt>(compdefPtrs.size());
    }

    inline PetscInt getNumberOfSpecies() const noexcept {
        return static_cast<PetscInt>(specIDs.size());
    }

    inline const std::vector<std::unique_ptr<Compdef>>& compdefs() const noexcept {
        return compdefPtrs;
    }

    inline const std::vector<std::unique_ptr<Patchdef>>& patchdefs() const noexcept {
        return patchdefPtrs;
    }

    inline const std::map<model::specie_name, model::specie_id>& getSpecModelIdxs() const noexcept {
        return specModelIdxs;
    }

    inline const model::specie_name& getSpecID(model::specie_id spec_model_idx) const {
        assert(static_cast<size_t>(spec_model_idx.get()) < specIDs.size());
        return specIDs[static_cast<size_t>(spec_model_idx.get())];
    }

    std::string createReport() const;

  private:
    std::map<model::specie_name, model::specie_id> specModelIdxs;
    std::vector<model::specie_name> specIDs;
    std::map<model::compartment_id, container::compartment_id> compModelIdxs;
    std::map<model::patch_id, container::patch_id> patchModelIdxs;
    std::vector<std::unique_ptr<Compdef>> compdefPtrs;
    std::vector<std::unique_ptr<Patchdef>> patchdefPtrs;
};

}  // namespace zee

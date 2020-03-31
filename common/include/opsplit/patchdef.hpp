#pragma once


#include <cassert>
#include <memory>
#include <numeric>
#include <set>
#include <unordered_map>
#include <vector>


#include <boost/optional.hpp>
#include <petscsys.h>

#include "fwd.hpp"
#include "vocabulary.hpp"


namespace zee {

/**
 * \brief State definition of a patch.
 *
 * The Patchdef class defines the sub biochemical container of a patch.
 * It provides the global and local indexing for species, reactions
 * and diffusions in the compartment.
 *
 */
class Patchdef {
  public:
    Patchdef(const Statedef& statedef,
             model::patch_id t_model_patch,
             model::compartment_id inner_compartment,
             const boost::optional<model::compartment_id>& outer_compartment);


    inline const model::patch_id& getID() const noexcept {
        return model_patch_;
    }

    inline container::patch_id getModelIdx() const noexcept {
        return container_patch_;
    }

    inline model::compartment_id getInnerCompId() const noexcept {
        return inner_compartment_id_;
    }

    inline const Compdef& getInnerComp() const noexcept;

    inline const boost::optional<model::compartment_id>& getOuterCompId() const noexcept {
        return outer_compartment_id_;
    }

    container::specie_id addSpec(model::specie_id specie);

    model::specie_id getSpecModelIdx(container::specie_id specie) const;

    inline osh::LO getNSpecs() const noexcept {
        return static_cast<osh::LO>(specC2M_.size());
    }

    container::surface_reaction_id addSurfaceReac(
        const std::vector<container::specie_id>& reactants_i,
        const std::vector<container::specie_id>& reactants_s,
        const std::vector<container::specie_id>& reactants_o,
        const std::vector<container::specie_id>& product_i,
        const std::vector<container::specie_id>& product_s,
        const std::vector<container::specie_id>& product_o,
        PetscScalar kcst);

    const SReacdef& getReac(container::surface_reaction_id reaction_id) const;

    /**
     * \return the reaction definitions
     */
    inline const std::vector<std::unique_ptr<SReacdef>>& reacdefs() const noexcept {
        return reacdefPtrs_;
    }

    inline const Statedef& statedef() const noexcept {
        return pStatedef_;
    }

    inline const std::set<container::specie_id>& getAllSpeciesDiffused() const noexcept {
        return species_diffused_;
    }

    inline bool isDiffused(const container::specie_id& specie) const {
        return std::find(species_diffused_.begin(), species_diffused_.end(), specie) !=
               species_diffused_.end();
    }

    inline PetscInt getNReacs() const;

    inline PetscInt getNKProcs() const;

    inline KProcType getKProcType(container::kproc_id kproc) const;

    container::specie_id getSpecPatchIdx(model::specie_id specie) {
        auto it = specM2C_.find(specie);
        if (it != specM2C_.end()) {
            return (*it).second;
        }
        throw std::logic_error("Unregistered specie id " + std::to_string(specie.get()));
    }

  private:
    // compartment KProc order: Reac then Diff
    const Statedef& pStatedef_;
    model::patch_id model_patch_;
    model::compartment_id inner_compartment_id_;
    boost::optional<model::compartment_id> outer_compartment_id_;
    container::patch_id container_patch_;
    std::unordered_map<model::specie_id, container::specie_id> specM2C_;
    std::vector<model::specie_id> specC2M_;
    PetscInt nKProcs_;
    std::vector<std::unique_ptr<SReacdef>> reacdefPtrs_;

    // This is in preparation for SReac
    std::set<container::specie_id> species_diffused_;
};

}  // namespace zee

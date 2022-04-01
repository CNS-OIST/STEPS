#pragma once

#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#include "fwd.hpp"

#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

/**
 * \brief State definition of a compartment.
 *
 * The Compdef class defines the sub biochemical container of a compartment.
 * It provides the global and local indexing for species, reactions
 * and diffusions in the compartment.
 *
 * This class corresponds to the solver::Compdef class in STEPS.
 */
class Compdef {
  public:
    Compdef(const Statedef& statedef,
            model::compartment_id t_model_compartment,
            container::compartment_id t_container_compartment);

    inline const model::compartment_id& getID() const noexcept {
        return model_compartment;
    }

    inline container::compartment_id getModelIdx() const noexcept {
        return container_compartment;
    }
    /**
     * Add the species to the compartment definition and return its local index.
     * If the species has been added before, return its lidx in record,
     * otherwise add the species to the record and return its new lidx.
     */
    container::species_id addSpec(model::species_id species);
    container::species_id getSpecContainerIdx(model::species_id species) const;
    model::species_id getSpecModelIdx(container::species_id species) const;

    /**
     * \return number of chemical species in the compartment
     */
    inline osh::I32 getNSpecs() const noexcept {
      return static_cast<osh::I32>(specC2M.size());
    }

    container::reaction_id
    addReac(const std::vector<container::species_id> &reactants,
            const std::vector<container::species_id> &products, osh::Real kcst);
    Reacdef& getReac(container::reaction_id reaction) const;

    /**
     * Register a diffusion process in the compartment
     * \param species the diffusion chemical specie
     * \param dcst the diffusion constant
     * \return the diffusion identifier
     */
    container::diffusion_id addDiff(container::species_id species,
                                    osh::Real dcst);

    /**
     * Get internal species identifier of the given diffusion
     * \param diffusion the diffusion identifier
     * \return the internal chemical species identifier
     */
    container::species_id
    getDiffSpecContainerIdx(container::diffusion_id diffusion);

    /**
     * Get global species identifier of the given diffusion
     * \param diffusion the diffusion identifier
     * \return the global chemical species identifier
     */
    model::species_id getDiffSpecModelIdx(container::diffusion_id diffusion);

    /**
     * \param diffusion internal diffusion identifier
     * \return definition of the specified diffusion
     */
    Diffdef& getDiff(container::diffusion_id diffusion);

    /**
     * Get diffusion definition
     * \param kproc internal kproc identifier of the diffusion
     * \return the diffusion definition
     */
    Diffdef& getDiffByKProcContainerIdx(container::kproc_id kproc);

    /**
     * Check if a KProc with local index kproc_lidx depends on
     * the species with local index spec_lidx.
     * \param kproc Local index of the KProc.
     * \param species Local index of the Species.
     * return True if there is dependency, false if not.
     */
    bool KProcDepSpec(container::kproc_id kproc,
                      container::species_id species) const;

    /**
     * \return number of kinetic processes defined in the compartment
     */
    inline osh::I64 getNKProcs() const noexcept { return nKProcs; }

    /**
     * \return number of reactions defined in the compartment
     */
    inline osh::I64 getNReacs() const noexcept {
      return static_cast<osh::I64>(reacdefPtrs.size());
    }

    /**
     * \return number of diffusions defined in the compartment
     */
    inline osh::I64 getNDiffs() const noexcept {
      return static_cast<osh::I64>(diffdefPtrs.size());
    }

    /**
     * \return the diffusion definitions
     */
    inline const std::vector<std::unique_ptr<Diffdef>>& diffdefs() const noexcept {
        return diffdefPtrs;
    }

    inline const std::set<container::species_id> &getAllSpeciesDiffused() const
        noexcept {
      return species_diffused_;
    }

    inline bool isDiffused(const container::species_id &species) const {
      return std::find(species_diffused_.begin(), species_diffused_.end(),
                       species) != species_diffused_.end();
    }

    /**
     * \return the reaction definitions
     */
    inline const std::vector<std::unique_ptr<Reacdef>>& reacdefs() const noexcept {
        return reacdefPtrs;
    }

    inline const Statedef& statedef() const noexcept {
        return pStatedef;
    }

    inline kproc::KProcType getKProcType(container::kproc_id kproc) const {
      if (kproc.get() >= 0 && kproc < getNReacs()) {
        return kproc::KProcType::Reac;
      } else if (kproc >= getNReacs() && kproc < getNKProcs()) {
        return kproc::KProcType::Diff;
      } else {
        throw std::out_of_range("KProc local index error.");
      }
    }

    void report(std::ostream& ostr) const;

  private:
    // compartment KProc order: Reac then Diff
    const Statedef& pStatedef;
    model::compartment_id model_compartment;
    container::compartment_id container_compartment;
    std::set<container::species_id> species_diffused_;
    std::unordered_map<model::species_id, container::species_id> specM2C;
    std::vector<model::species_id> specC2M;

    osh::I64 nKProcs{};
    std::vector<std::unique_ptr<Reacdef>> reacdefPtrs;
    std::vector<std::unique_ptr<Diffdef>> diffdefPtrs;
};

}  // namespace dist
}  // namespace steps

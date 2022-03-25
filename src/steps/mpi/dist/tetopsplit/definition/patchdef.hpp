#pragma once


#include <cassert>
#include <memory>
#include <numeric>
#include <set>
#include <unordered_map>
#include <vector>


#include <boost/optional.hpp>

#include "fwd.hpp"
#include "sreacdef.hpp"
#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

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
    Patchdef(const Statedef &statedef, model::patch_id t_model_patch,
             container::patch_id container_patch_id,
             model::compartment_id inner_compartment,
             const boost::optional<model::compartment_id> &outer_compartment);

    inline const model::patch_id& getID() const noexcept {
        return model_patch_;
    }

    inline container::patch_id getModelIdx() const noexcept {
      return container_patch_id_;
    }

    inline model::compartment_id getInnerCompId() const noexcept {
        return inner_compartment_id_;
    }

    inline const Compdef& getInnerComp() const noexcept;

    inline const boost::optional<model::compartment_id>& getOuterCompId() const noexcept {
        return outer_compartment_id_;
    }

    container::species_id addSpec(model::species_id species);

    model::species_id getSpecModelIdx(container::species_id species) const;

    inline osh::LO getNSpecs() const noexcept {
        return static_cast<osh::LO>(specC2M_.size());
    }


    const SReacdef& getReac(container::surface_reaction_id reaction_id) const;

    /**
     * \return the reaction definitions
     */
    inline const std::vector<std::unique_ptr<SReacdef>>& reacdefs() const noexcept {
        return reacdefPtrs_;
    }

    inline const std::vector<std::unique_ptr<VDepSReacdef>> &
    vDepReacdefs() const noexcept {
      return vdepSReacPtrs_;
    }
    
    inline const std::vector<std::unique_ptr<GHKSReacdef>> &
    ghkSReacdefs() const noexcept {
      return ghkSReacPtrs_;
    }    

    inline const Statedef& statedef() const noexcept {
        return pStatedef_;
    }

    inline const std::set<container::species_id> &getAllSpeciesDiffused() const
        noexcept {
      return species_diffused_;
    }

    inline bool isDiffused(const container::species_id &species) const {
      return species_diffused_.find(species) != species_diffused_.end();
    }

    inline osh::I64 getNReacs() const;

    inline osh::I64 getNKProcs() const;

    inline kproc::KProcType getKProcType(container::kproc_id kproc) const;

    container::species_id getSpecPatchIdx(model::species_id species) const;

    template <typename PropensityType>
    container::surface_reaction_id
    addSurfaceReacImpl(const std::vector<container::species_id> &reactants_i,
                       const std::vector<container::species_id> &reactants_s,
                       const std::vector<container::species_id> &reactants_o,
                       const std::vector<container::species_id> &products_i,
                       const std::vector<container::species_id> &products_s,
                       const std::vector<container::species_id> &products_o,
                       PropensityType kcst);

    template <typename T>
    std::vector<std::unique_ptr<SReacdefBase<T>>> &getContainer();
    //-------------------------------------------------------

  private:
    // compartment KProc order: Reac then Diff
    const Statedef& pStatedef_;
    model::patch_id model_patch_;
    model::compartment_id inner_compartment_id_;
    boost::optional<model::compartment_id> outer_compartment_id_;
    container::patch_id container_patch_id_;
    std::unordered_map<model::species_id, container::species_id> specM2C_;
    std::vector<model::species_id> specC2M_;
    osh::I64 nKProcs_;
    std::vector<std::unique_ptr<SReacdef>> reacdefPtrs_;
    std::vector<std::unique_ptr<VDepSReacdef>> vdepSReacPtrs_;
    std::vector<std::unique_ptr<GHKSReacdef>> ghkSReacPtrs_;
    // This is in preparation for SReac
    std::set<container::species_id> species_diffused_;
};

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<SReacdef>> &Patchdef::getContainer<SReacInfo>();

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<VDepSReacdef>> &Patchdef::getContainer<VDepInfo>();

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<GHKSReacdef>> &Patchdef::getContainer<GHKInfo>();

//-------------------------------------------------------

extern template std::vector<std::unique_ptr<SReacdefBase<SReacInfo>>> &
Patchdef::getContainer();

//-------------------------------------------------------

extern template std::vector<std::unique_ptr<SReacdefBase<VDepInfo>>> &
Patchdef::getContainer();

//-------------------------------------------------------

extern template std::vector<std::unique_ptr<SReacdefBase<GHKInfo>>> &
Patchdef::getContainer();

//-------------------------------------------------------

extern template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id> &reactants_i,
    const std::vector<container::species_id> &reactants_s,
    const std::vector<container::species_id> &reactants_o,
    const std::vector<container::species_id> &products_i,
    const std::vector<container::species_id> &products_s,
    const std::vector<container::species_id> &products_o, SReacInfo kcst);

//-------------------------------------------------------

extern template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id> &reactants_i,
    const std::vector<container::species_id> &reactants_s,
    const std::vector<container::species_id> &reactants_o,
    const std::vector<container::species_id> &products_i,
    const std::vector<container::species_id> &products_s,
    const std::vector<container::species_id> &products_o, VDepInfo kcst);

//-------------------------------------------------------

extern template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id> &reactants_i,
    const std::vector<container::species_id> &reactants_s,
    const std::vector<container::species_id> &reactants_o,
    const std::vector<container::species_id> &products_i,
    const std::vector<container::species_id> &products_s,
    const std::vector<container::species_id> &products_o, GHKInfo kcst);

//-------------------------------------------------------

}  // namespace dist
}  // namespace steps

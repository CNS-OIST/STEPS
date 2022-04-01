#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <vector>

#include <boost/optional.hpp>

#include "efield.hpp"
#include "fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "model/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

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
     * \brief Construct a Statedef object from steps::model::Model and steps::dist::DistMesh.
     *
     * \attention Parallelism: Distribute
     *
     * \param model Reference to the steps::model::Model object.
     * \param mesh Reference to the steps::dist::DistMesh object.
     */
    Statedef(const steps::model::Model &model,
             const steps::dist::DistMesh &mesh);

    /**
     * Add the species to the biochemical state definition and return its container index.
     * If the species has been added before, return its container index in record,
     * otherwise add the species to the record and return its new container index.
     */
    model::species_id addSpec(const model::species_name &name);
    model::species_id getSpecModelIdx(const model::species_name &name) const;
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
    model::species_id addCompSpec(const model::compartment_id &compartment,
                                  const model::species_name &species);

    std::vector<model::species_id>
    addCompSpecs(const model::compartment_id &compartment,
                 const std::vector<model::species_name> &species);
    std::vector<model::species_id>
    addPatchSpecs(const model::patch_id &patchId,
                  const std::vector<model::species_name> &species);
    model::species_id addPatchSpec(const model::patch_id &patch_id,
                                   const model::species_name &species);
    container::species_id
    getCompSpecContainerIdx(const model::compartment_id &compartment,
                            const model::species_name &species) const;

    /**
     * Register a diffusion in a given compartment
     * \param compartment the compartment identifier where to register the
     * diffusion \param species_name the diffusing chemical specie \param dcst
     * diffusion constant \return the diffusion identifier
     */
    container::diffusion_id
    addCompDiff(const model::compartment_id &compartment,
                const model::species_name &species_name, osh::Real dcst);

    /**
     * Register a chemical reaction to a compartment
     * \param compartment the compartment identifier
     * \param reactants substances initially involved in the chemical reaction
     * \param products substances yield by the chemical reaction
     * \param kcst the reaction constant
     * \return the reaction identifier
     */
    container::reaction_id
    addCompReac(const model::compartment_id &compartment,
                const std::vector<model::species_name> &reactants,
                const std::vector<model::species_name> &products,
                osh::Real kcst);

    container::surface_reaction_id
    addSurfReac(const model::patch_id &patchId,
                const std::vector<model::species_name> &reactants_i,
                const std::vector<model::species_name> &reactants_s,
                const std::vector<model::species_name> &reactants_o,
                const std::vector<model::species_name> &products_i,
                const std::vector<model::species_name> &products_s,
                const std::vector<model::species_name> &products_o,
                osh::Real kcst);

    container::surface_reaction_id
    addVDepSurfReac(const model::patch_id &patchId,
                    const std::vector<model::species_name> &reactants_i,
                    const std::vector<model::species_name> &reactants_s,
                    const std::vector<model::species_name> &reactants_o,
                    const std::vector<model::species_name> &products_i,
                    const std::vector<model::species_name> &products_s,
                    const std::vector<model::species_name> &products_o,
                    const std::function<osh::Real(osh::Real)> &kcst);

    Compdef &getCompdef(container::compartment_id compartment) const noexcept;

    Compdef &getCompdef(const model::compartment_id &compartment) const
        noexcept;

    inline Compdef &
    getDefinition(const model::compartment_id &compartmentId) const noexcept {
      return getCompdef(compartmentId);
    }

    Patchdef &getPatchdef(const container::patch_id &patchId) const noexcept;

    Patchdef &getPatchdef(const model::patch_id &patchId) const noexcept;

    inline Patchdef &getDefinition(const model::patch_id &patchId) const
        noexcept {
      return getPatchdef(patchId);
    }

    inline osh::I64 getNComps() const noexcept {
      return static_cast<osh::I64>(compdefPtrs.size());
    }

    inline osh::I64 getNumberOfSpecies() const noexcept {
      return static_cast<osh::I64>(specIDs.size());
    }

    inline const std::vector<std::unique_ptr<Compdef>>& compdefs() const noexcept {
        return compdefPtrs;
    }

    inline const std::vector<std::unique_ptr<Patchdef>>& patchdefs() const noexcept {
        return patchdefPtrs;
    }

    inline const std::map<model::ohmic_current_id,
                          std::unique_ptr<OhmicCurrent>> &
    ohmicCurrents() const noexcept {
        return ohmicCurrPtrs;
    }

    inline const std::map<model::ghk_current_id, std::unique_ptr<GHKCurrent>> &
    ghkCurrents() const noexcept {
        return ghkCurrPtrs;
    }

    inline const std::map<model::membrane_id, std::unique_ptr<Membrane>> &
    membranes() const noexcept {
      return membranePtrs;
    }

    inline void addMembrane(const model::membrane_id &membrane,
                            const model::patch_id &patch, double capacitance) {
      membranePtrs.emplace(membrane,
                           std::make_unique<Membrane>(patch, capacitance));
    }

    void addChannel(const model::membrane_id &membrane,
                    const model::channel_id &channel,
                    const std::vector<model::species_name> &channel_states);

    void
    addOhmicCurrent(const model::ohmic_current_id &curr,
                    const model::membrane_id &membrane,
                    const model::channel_id &channel,
                    const boost::optional<model::species_name> &species_name,
                    double conductance, double reversal_potential);

    inline void setStimulus(const model::membrane_id& membrane, osh::Real current) {
        membranePtrs[membrane]->setStimulus([current](auto) { return current; });
    }

    void addCompartmentConductivity(const model::compartment_id &comp_id,
                                    osh::Real conductivity);

    osh::Real
    getCompartmentConductivity(const model::compartment_id &compartment) const {
      const auto it = compartment_conductivity_.find(compartment);
      if (it != compartment_conductivity_.end()) {
        return it->second;
      } else {
        throw std::invalid_argument("No conductivity defined for " + compartment);
      }
    }

    /**
     * Add a GHK current surface reaction
     *
     * Check the wiki or this online calculator for more info:
     * https://www.physiologyweb.com/calculators/ghk_equation_calculator.html
     *
     * \param membrane: on which the reaction is
     * \param ion_channel: the channel that follows the GHK dynamics
     * \param ion_channel_state: the open state of the ion channel
     * \param ion_name: the ion interested by this reaction
     * \param permeability: usually computed by math::permeability in real life scenarios
     * \param valence: of the ion
     *
     * \return the reaction identifier
     */
    void addGHKCurrentSurfReac(const model::ghk_current_id& curr_id,
                               const model::membrane_id& membrane,
                               const model::channel_id& ion_channel,
                               const model::species_name& ion_channel_state,
                               const model::species_name& ion_name,
                               osh::Real permeability,
                               osh::I64 valence,
                               boost::optional<osh::Real> outer_conc = boost::none,
                               boost::optional<osh::Real> inner_conc = boost::none);

    inline const std::map<model::species_name, model::species_id> &
    getSpecModelIdxs() const noexcept {
      return specModelIdxs;
    }

    const container::surface_reaction_id &
    getSReacIdx(const model::surface_reaction_id &reac) const;

    inline const model::species_name &
    getSpecID(model::species_id spec_model_idx) const noexcept {
      assert(static_cast<size_t>(spec_model_idx.get()) < specIDs.size());
      return specIDs[static_cast<size_t>(spec_model_idx.get())];
    }

    std::string createReport() const;

    /**
     * \return true if at least one membrane is defined, false otherwise.
     */
    inline bool is_efield_enabled() const noexcept {
      return !membranes().empty() && efield_enabled_;
    }

    void disableEField() { efield_enabled_ = false; }

    /// Get temperature
    inline osh::Real getTemp() const noexcept {
        return temperature;
    }
    /// Set temperature
    inline void setTemp(const osh::Real temp) noexcept {
        temperature = temp;
    }

  private:
    template <typename PropensityType>
    container::surface_reaction_id
    addSurfReacImpl(const model::patch_id &patchId,
                    const std::vector<model::species_name> &reactants_i,
                    const std::vector<model::species_name> &reactants_s,
                    const std::vector<model::species_name> &reactants_o,
                    const std::vector<model::species_name> &products_i,
                    const std::vector<model::species_name> &products_s,
                    const std::vector<model::species_name> &products_o,
                    PropensityType kcst);

    std::map<model::species_name, model::species_id> specModelIdxs;
    std::vector<model::species_name> specIDs;
    std::map<model::compartment_id, container::compartment_id> compModelIdxs;
    std::map<model::patch_id, container::patch_id> patchModelIdxs;
    std::map<model::surface_reaction_id, container::surface_reaction_id> surfReacIdxs;
    std::vector<std::unique_ptr<Compdef>> compdefPtrs;
    std::vector<std::unique_ptr<Patchdef>> patchdefPtrs;
    std::map<model::membrane_id, std::unique_ptr<Membrane>> membranePtrs;
    std::map<model::compartment_id, osh::Real> compartment_conductivity_;
    std::map<model::ohmic_current_id, std::unique_ptr<OhmicCurrent>> ohmicCurrPtrs;
    std::map<model::ghk_current_id, std::unique_ptr<GHKCurrent>> ghkCurrPtrs;
    bool efield_enabled_ = true;

    /// Global temperature. Default = 20c as in src/steps/mpi/tetopsplit/tetopsplit.cpp
    osh::Real temperature{293.15};
};

}  // namespace dist
}  // namespace steps

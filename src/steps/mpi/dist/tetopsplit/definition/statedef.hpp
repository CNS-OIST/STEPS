#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "complexdef.hpp"
#include "efield.hpp"
#include "fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "model/fwd.hpp"
#include "model/model.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

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
    Statedef(const steps::model::Model& model, const steps::dist::DistMesh& mesh);

    /// Reset the def-object values to model defaults
    void reset();

    /**
     * Add the species to the biochemical state definition and return its container index.
     * If the species has been added before, return its container index in record,
     * otherwise add the species to the record and return its new container index.
     */
    model::species_id getSpecModelIdx(const model::species_name& name) const;
    model::species_id getSpecModelIdx(const steps::model::Spec& spec) const;

    model::complex_id getComplexModelIdx(const model::complex_name& name) const;
    container::compartment_id getCompModelIdx(
        const model::compartment_id& compartment) const noexcept;

    model::complex_id addComplex(const steps::model::Complex& cmplx);

    container::species_id getCompSpecContainerIdx(const model::compartment_id& compartment,
                                                  const model::species_name& species) const;

    container::surface_reaction_id addVDepSurfReac(
        const model::patch_id& patchId,
        const std::vector<model::species_name>& reactants_i,
        const std::vector<model::species_name>& reactants_s,
        const std::vector<model::species_name>& reactants_o,
        const std::vector<model::species_name>& products_i,
        const std::vector<model::species_name>& products_s,
        const std::vector<model::species_name>& products_o,
        const std::function<osh::Real(osh::Real)>& kcst);

    Compdef& getCompdef(container::compartment_id compartment) const noexcept;

    Compdef& getCompdef(const model::compartment_id& compartment) const noexcept;

    inline Compdef& getDefinition(const model::compartment_id& compartmentId) const noexcept {
        return getCompdef(compartmentId);
    }

    inline Compdef& getDefinition(const container::compartment_id& compartmentId) const noexcept {
        return getCompdef(compartmentId);
    }

    Patchdef& getPatchdef(const container::patch_id& patchId) const noexcept;

    Patchdef& getPatchdef(const model::patch_id& patchId) const;

    Membrane& getMembrane(const model::membrane_id& membraneId) const;

    inline Patchdef& getDefinition(const model::patch_id& patchId) const noexcept {
        return getPatchdef(patchId);
    }

    inline Patchdef& getDefinition(const container::patch_id& patchId) const noexcept {
        return getPatchdef(patchId);
    }

    inline osh::I64 getNComps() const noexcept {
        return static_cast<osh::I64>(compdefPtrs.size());
    }

    inline osh::I64 getNumberOfSpecies() const noexcept {
        return static_cast<osh::I64>(specIDs.size());
    }

    inline const util::strongid_vector<container::compartment_id, std::unique_ptr<Compdef>>&
    compdefs() const noexcept {
        return compdefPtrs;
    }

    inline const util::strongid_vector<container::patch_id, std::unique_ptr<Patchdef>>& patchdefs()
        const noexcept {
        return patchdefPtrs;
    }

    inline const util::strongid_vector<container::membrane_id, std::unique_ptr<Membrane>>&
    membranes() const noexcept {
        return membranePtrs;
    }

    void setStimulus(const model::membrane_id& membrane, osh::Real current) const;

    void setResistivity(const model::membrane_id& membrane, osh::Real resistivity) const;

    osh::Real getResistivity(const model::membrane_id& membrane) const;

    void setReversalPotential(const model::membrane_id& membrane,
                              osh::Real reversal_potential) const;

    osh::Real getReversalPotential(const model::membrane_id& membrane) const;

    inline const std::map<model::species_name, model::species_id>& getSpecModelIdxs()
        const noexcept {
        return specModelIdxs;
    }

    container::surface_reaction_id getSReacIdx(const model::patch_id& patchId,
                                               const model::surface_reaction_id& reac) const;

    inline const model::species_name& getSpecID(model::species_id spec_model_idx) const noexcept {
        assert(spec_model_idx < specIDs.size());
        return specIDs[spec_model_idx];
    }

    std::string createReport() const;

    /**
     * \return true if at least one membrane is defined, false otherwise.
     */
    inline bool is_efield_enabled() const noexcept {
        return !membranes().container().empty() && efield_enabled_;
    }

    void disableEField() {
        efield_enabled_ = false;
    }

    /// Get temperature
    inline osh::Real getTemp() const noexcept {
        return temperature;
    }
    /// Set temperature
    inline void setTemp(const osh::Real temp) noexcept {
        temperature = temp;
    }

    [[nodiscard]] const steps::model::Model& model() const noexcept {
        return pModel;
    }

    [[nodiscard]] const DistMesh& mesh() const noexcept {
        return pMesh;
    }

    osh::LO num_complexes() const noexcept {
        return complexdefPtrs.size();
    }

    osh::LOs substates_per_complexes() const noexcept {
        osh::Write<osh::LO> ret(complexdefPtrs.size());
        std::transform(complexdefPtrs.begin(),
                       complexdefPtrs.end(),
                       ret.begin(),
                       [](const auto& complexdef) { return complexdef->nbSubStates(); });
        return ret;
    }

  private:
    model::species_id addSpec(const steps::model::Spec& spec);
    container::compartment_id addComp(const DistComp& compartment);
    container::patch_id addPatch(const DistPatch& patch);
    container::membrane_id addMembrane(const DistMemb& membrane);

    template <typename PropensityType>
    container::surface_reaction_id addSurfReacImpl(
        const model::patch_id& patchId,
        std::optional<model::surface_reaction_id> sreacId,
        const std::vector<model::species_name>& reactants_i,
        const std::vector<model::species_name>& reactants_s,
        const std::vector<model::species_name>& reactants_o,
        const std::vector<model::species_name>& products_i,
        const std::vector<model::species_name>& products_s,
        const std::vector<model::species_name>& products_o,
        PropensityType kcst);

    const steps::model::Model& pModel;
    const DistMesh& pMesh;

    std::map<model::species_name, model::species_id> specModelIdxs;
    std::map<model::complex_name, model::complex_id> complexModelIdxs;
    std::map<model::compartment_id, container::compartment_id> compModelIdxs;
    std::map<model::patch_id, container::patch_id> patchModelIdxs;
    std::map<model::membrane_id, container::membrane_id> membraneModelIdxs;

    util::strongid_vector<model::species_id, model::species_name> specIDs;
    util::strongid_vector<model::complex_id, model::complex_name> complexIDs;
    util::strongid_vector<container::compartment_id, std::unique_ptr<Compdef>> compdefPtrs;
    util::strongid_vector<container::patch_id, std::unique_ptr<Patchdef>> patchdefPtrs;
    util::strongid_vector<container::membrane_id, std::unique_ptr<Membrane>> membranePtrs;
    util::strongid_vector<model::complex_id, std::unique_ptr<Complexdef>> complexdefPtrs;

    bool efield_enabled_ = true;

    /// Global temperature. Default = 20c as in src/steps/mpi/tetopsplit/tetopsplit.cpp
    osh::Real temperature{293.15};
};

}  // namespace steps::dist

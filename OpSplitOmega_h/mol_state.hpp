#pragma once

#include <boost/optional.hpp>

#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>

#include "common.hpp"
#include "flat_multimap.hpp"

namespace zee {
/**
 * Keep track of molecules per specie per element
 */

struct NumMolecules {
    friend class MolState;
    using elem_type = osh::LO;

    explicit NumMolecules(const osh::LOs& t_species_per_elements)
        : pools_(t_species_per_elements)
        , species_per_elements_(t_species_per_elements) {}

    inline elem_type& operator()(osh::LO element, container::specie_id specie) noexcept {
        assert(specie.get() < numSpecies(element));
        return pools_(element, specie.get());
    }
    inline elem_type operator()(osh::LO element, container::specie_id specie) const noexcept {
        assert(specie.get() < numSpecies(element));
        return pools_(element, specie.get());
    }
    inline bool empty(osh::LO element, container::specie_id specie) const noexcept {
        assert(specie.get() < numSpecies(element));
        return pools_(element, specie.get()) == 0;
    }
    inline void reset() {
        pools_.assign(0);
    }

    inline osh::LO numElements() const noexcept {
        return pools_.size();
    }

    inline osh::LO numSpecies(osh::LO element) const noexcept {
        return pools_.size(element);
    }

    inline elem_type sumNumMolecules(container::specie_id specie) const {
        elem_type num_molecules{};
        for (osh::LO elem = 0; elem < numElements(); ++elem) {
            if (specie.get() < numSpecies(elem)) {
                num_molecules += this->operator()(elem, specie);
            }
        }
        return num_molecules;
    }

    inline const osh::LOs& speciesPerElements() const noexcept {
        return species_per_elements_;
    }

  private:
    zee::flat_multimap<elem_type, 1, zee::OSH> pools_;
    osh::LOs species_per_elements_;
};


class MolState {
  public:
    using elem_type = NumMolecules::elem_type;
    enum class Location : char { volume, boundary };
    /// field 1: kind of entity
    /// field 2: entity identifier
    /// field 3: specie identifier
    using ElementID =
        std::tuple<boost::variant<mesh::element_id, mesh::boundary_id>, container::specie_id>;
    static inline ElementID mkBoundaryElement(mesh::boundary_id element,
                                              container::specie_id specie) {
        return std::make_pair(element, specie);
    }
    static inline ElementID mkVolumeElement(mesh::element_id element, container::specie_id specie) {
        return std::make_pair(element, specie);
    }

    MolState(const osh::LOs& t_species_per_elements,
             const boost::optional<osh::LOs>& t_species_per_boundary_element = boost::none)
        : molecules_on_elements_(t_species_per_elements)
        , molecules_on_patch_boundaries_(t_species_per_boundary_element
                                             ? (*t_species_per_boundary_element)
                                             : osh::Write<osh::LO>(1, 0)) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"

    inline osh::LO& operator()(const ElementID& elementId) noexcept {
        auto element = std::get<0>(elementId);
        auto specie = std::get<1>(elementId);
        return boost::apply_visitor([this, specie](auto entity)
                                        -> osh::LO& { return this->operator()(entity, specie); },
                                    std::get<0>(elementId));
    }

    inline osh::LO operator()(const ElementID& elementId) const noexcept {
        auto element = std::get<0>(elementId);
        auto specie = std::get<1>(elementId);
        return boost::apply_visitor([this, specie](auto entity)
                                        -> osh::LO { return this->operator()(entity, specie); },
                                    std::get<0>(elementId));
    }
#pragma GCC diagnostic pop

    inline osh::LO& operator()(mesh::element_id element, container::specie_id specie) noexcept {
        return molecules_on_elements_(element.get(), specie);
    }
    inline osh::LO operator()(mesh::element_id element, container::specie_id specie) const
        noexcept {
        return molecules_on_elements_(element.get(), specie);
    }

    inline osh::LO& operator()(mesh::boundary_id element, container::specie_id specie) noexcept {
        return molecules_on_patch_boundaries_(element.get(), specie);
    }
    inline osh::LO operator()(mesh::boundary_id element, container::specie_id specie) const
        noexcept {
        return molecules_on_patch_boundaries_(element.get(), specie);
    }

    inline bool empty(mesh::element_id element, container::specie_id specie) const noexcept {
        return molecules_on_elements_(element.get(), specie) == 0;
    }
    inline void reset() {
        molecules_on_elements_.pools_.assign(0);
    }

    inline osh::LO numElements() const noexcept {
        return molecules_on_elements_.pools_.size();
    }

    inline osh::LO numSpecies(mesh::element_id element) const noexcept {
        return molecules_on_elements_.pools_.size(element.get());
    }

    /**
     * \copybrief molecules_on_elements_
     */
    inline const NumMolecules& moleculesOnElements() const noexcept {
        return molecules_on_elements_;
    }

    /**
     * \copybrief molecules_on_elements_
     */
    inline NumMolecules& moleculesOnElements() noexcept {
        return molecules_on_elements_;
    }

    /**
     * \copybrief NumMolecules::molecules_on_patch_boundaries_
     */
    inline const NumMolecules& moleculesOnPatchBoundaries() const noexcept {
        return molecules_on_patch_boundaries_;
    }

    /**
     * \copybrief NumMolecules::molecules_on_patch_boundaries_
     */
    inline NumMolecules& moleculesOnPatchBoundaries() noexcept {
        return molecules_on_patch_boundaries_;
    }

    /**
     * \copybrief NumMolecules::molecules_on_patch_boundaries_
     */
    NumMolecules& moleculesCountOnPatchBoundaries() noexcept {
        return molecules_on_patch_boundaries_;
    }

    inline const osh::LOs& species_per_elements() const noexcept {
        return molecules_on_elements_.species_per_elements_;
    }

  private:
    /**
     * \brief Container providing the number of molecules of every specie
     * within the elements of the local mesh.
     */
    NumMolecules molecules_on_elements_;

    /**
     * \brief Container providing the number of molecules of every specie
     * within the boundaries of the local mesh that belong to a patch.
     */
    NumMolecules molecules_on_patch_boundaries_;
};

}  // namespace zee

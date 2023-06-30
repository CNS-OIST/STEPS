#pragma once

#include <set>
#include <string>

#include "geom/dist/distmesh.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

class DistMemb {
  public:
    /**
     * \brief Add a membrane to a distributed mesh.
     *
     * \attention Parallelism: Collective
     *
     * \param membrane Name id of the membrane.
     * \param mesh The distributed mesh instance
     * \param patches The list of patches included in the membrane
     * \param capac Capacitance per unit area
     */
    DistMemb(const model::membrane_id& membrane,
             DistMesh& mesh,
             std::set<model::patch_id> patches,
             double capac);

    /**
     * \brief Return a pointer to the distmesh container object.
     *
     * \attention Parallelism: Collective
     */
    inline DistMesh* getContainer() const noexcept {
        return &meshRef;
    }

    /**
     * \brief Return the membrane id.
     *
     * \attention Parallelism: Collective
     */
    inline std::string const& getID() const noexcept {
        return pID;
    }

    /**
     * \brief Return the patches.
     *
     * \attention Parallelism: Collective
     */
    inline std::set<model::patch_id> const& patches() const noexcept {
        return pPatches;
    }

    /**
     * \brief Return the membrane capacitance.
     *
     * \attention Parallelism: Collective
     */
    inline double getCapacitance() const noexcept {
        return pCapacitance;
    }

    /**
     * \brief Set the membrane capacitance.
     *
     * \attention Parallelism: Collective
     */
    inline void setCapacitance(double capacitance) noexcept {
        pCapacitance = capacitance;
    }

  private:
    std::string pID;
    DistMesh& meshRef;
    std::set<model::patch_id> pPatches;
    double pCapacitance;
};

}  // namespace steps::dist

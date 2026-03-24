#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <variant>

#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>

#include "definition/compdef.hpp"
#include "definition/complexeventsdef.hpp"
#include "definition/patchdef.hpp"
#include "definition/statedef.hpp"
#include "geom/dist/distmesh.hpp"
#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "solver/fwd.hpp"
#include "util/collections.hpp"
#include "util/flat_multimap.hpp"
#include "util/strong_ra.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

/** Small object that implements occupancy
 *
 * Let us define occupancy as the integral average, over an certain period of time, of a quantity.
 * At the moment, there are two kinds of occupancies that need to be tracked:
 * - particle counts
 * - channel state counts
 *
 * Particle occupancy
 *
 * Diffusions depend on the number of particle present in certain element (tet/tri). Since the
 * reaction-diffusion loop operates advancing first the reactions and then the diffusions (over the
 * time step \f rd_{dt}\f), the diffusion operator accounts for the fact that the final number of
 * particles was not present for the whole rd_dt using the occupancy and not the final particle
 * count. Particle occupancy is only needed when there are reactions with particles that can
 * diffuse.
 *
 * Channel state occupancy
 *
 * The efield operator requires the number of channels in open state to compute their relative ohmic
 * currents. Similar to what we said before, channel states open and close during the various
 * reaction-diffusion loops. In order to take into account that the final number of channel states
 * was not open for the whole duration of the \f ef_{dt}\f we consider instead the channel state
 * occupancy (this time over the \f ef_{dt} \f).
 *
 * Implementation
 *
 * The occupancy \f o \f
 *
 * \f o = n-c \f
 * \f \sum_e v_e (t_e - t_s) \f
 *
 * where:
 * - \f n\f is the final (at the end of the appropriate dt) number of particles/channel states
 * - \f c\f is a corrective term. It is reset to 0 at the beginning of the appropriate time step
 * (rd_dt for a particle count, ef_dt for a channel state)
 * - \f v_e\f is the total valence of particles that are added/removed in the event e (usually it
 * is: \f n_{particlesvalenceion}\f) \f t_e\f is the event time \f t_s\f is the start time of the
 * recording. It is updated at the beginning of the appropriate dt (\f rd_{dt} \f for a particle
 * count, \f ef_{dt} \f for a channel state) The advantages for computing it in this way is that we
 * do not need to forecast when the time step will end and in all the functions we provide the
 * current time or the event time.
 *
 * There are a few ways to store the corrective term and reset the occupancy values every time step.
 * Since the main bottleneck for STEPS 4 is time complexity, we chose the fastest option (that
 * requires marginally more or equal space). We use:
 *
 * - one correction vector of doubles to track occupancy
 * - a vector of indexes of the particles/channel states to loop only over them to speed-up the
 * reset phase
 *
 * \Note: as convention, a nan in the correction vector means that we are not interested in tracking
 * that occupancy. To activate tracking we can use the track function which sets the occupancy to 0
 * and starts tracking. Thus, having nans in the correction vector is totally normal: it means that
 * we are not interested in that occupancy and if we still ask for the occupancy in that case we
 * just return the final value of the pools.
 *
 * FAQ:
 * - When can we safely disregard occupancy?
 *
 * For particles that do not diffuse or particle that do not react we are not interested in the
 * reaction-diffusion occupancy. In the same way, if there are no ohmic currents, we are not
 * interested in channel state tracking for the efield loop.
 *
 * - Why are not we throwing an error if we ask for the occupancy of a particle that is not marked
 * for tracking?
 *
 * Let us focus on particles and the reaction-diffusion loop. Imagine a particle that diffuses but
 * does not react. In that case we are sure that its number does not change during the reaction part
 * of the loop. The occupancy is always equal to the final number in pools.
 *
 */
class Occupancy {
  public:
    /// Ctor
    Occupancy(const size_t correction_size = 0, const osh::Real start_time = 0.0)
        : corrections_(correction_size, std::numeric_limits<osh::Real>::quiet_NaN())
        , start_time_(start_time) {}

    /// Size. To throw out_of_range errors
    inline size_t size() const noexcept {
        return static_cast<size_t>(corrections_.size());
    }

    /// Empty. To verify if we are tracking or not
    inline bool empty() const noexcept {
        return size() == 0;
    }

    /// Track occupancy of the molecule that corresponds to index
    void track(const size_t index) {
        assert(index < size());

        // If nan it means that we were not tracking before. Start now. This prevents doubles in
        // ids_
        auto& corr = corrections_[index];
        if (ignore_correction(corr)) {
            corr = 0.0;
            ids_.emplace_back(index);
        }
    }

    /** Full reset for
     *
     * \Note: it does not reset occupancy tracking
     *
     * @param start_time: every integral type starts at start_time
     */
    inline void reset(const osh::Real start_time) {
        start_time_ = start_time;
        osh::parallel_for(
            ids_.size(), OMEGA_H_LAMBDA(osh::LO index) { corrections_[ids_[index]] = 0.0; });
    }

    /** Some new molecules are added to the pools
     *
     * We record the amount of time these molecules were NOT present in the pool (times the
     * number of molecules).
     *
     * @param index: hash of elem id and species
     * @param val: molecule change
     * @param event_time: simulation time at the event
     */
    inline void add_correction(size_t index, const molecules_t val, const osh::Real event_time) {
        assert(index < size());
        auto& corr = corrections_[index];

        if (ignore_correction(corr)) {
            return;
        }
        // the - is there so that to get the correct occupancy at the end we sum the values
        corr -= val * (event_time - start_time_);
    }

    /** Get the occupancy
     *
     * The occupancy is just the current pool value + the correction value divided by the time
     * that passed since last reset. Notice that we always subtract every positive correction
     * (and add every negative) so that this is an addition
     *
     * @param index: hash of elem id and species
     * @param pool: molecule quantity (taken from pools)
     * @param end_time: end_time of the integral so that the timestep is end_time - start_time
     * @return occupancy
     */
    inline osh::Real get_occupancy(const size_t index,
                                   const molecules_t pool,
                                   const osh::Real end_time) const {
        assert(index < size());
        const auto corr = corrections_[index];

        if (ignore_correction(corr)) {
            return pool;
        }

        const auto dt = end_time - start_time_;

        // in case the dt ~ 0 we could end up with the case 0/0. We avoid it here
        if (util::almost_equal(dt, 0.0)) {
            return pool;
        }

        assert(dt > 0);

        return pool + corr / dt;
    }

    /// Pretty print
    friend std::ostream& operator<<(std::ostream& os, const Occupancy& o);

  private:
    /// \return true if the given correction should be ignore, false otherwise
    inline bool ignore_correction(const osh::Real correction) const noexcept {
        return std::isnan(correction);
    }

    /// Corrective terms that added to the pools and divided by the integration time gives the
    /// occupancy
    osh::Write<osh::Real> corrections_;

    /// Indexes of correction_ that need to be reset every new step
    std::vector<size_t> ids_;

    /// Integration starting time
    osh::Real start_time_;
};


/**
 * \brief Rank and comm size
 */
struct RankInfo {
    RankInfo(int r, int s)
        : rank(r)
        , size(s) {}
    int rank;
    int size;
};

/**
 * \brief Class for generating unique identifiers for complexes
 */
class ComplexIndexer {
  public:
    ComplexIndexer(RankInfo info)
        : nextComplexInd(std::numeric_limits<model::complex_individual_id::value_type>::max() /
                         info.size * info.rank)
        , maxComplexInd(std::numeric_limits<model::complex_individual_id::value_type>::max() /
                        info.size * (info.rank + 1)) {}

    model::complex_individual_id getNextComplexInd() {
        assert(nextComplexInd.get() < maxComplexInd.get());
        return nextComplexInd++;
    }

  private:
    model::complex_individual_id nextComplexInd;
    model::complex_individual_id maxComplexInd;
};

/** Keeps track of molecules per species per entity (element or boundary)
 *
 * \tparam Entity a strong_id type, \c mesh::element_id or \c
 * mesh::triangle_id_t for instance
 */
template <typename Entity>
struct EntityMolecules {
    explicit EntityMolecules(const DistMesh& mesh_,
                             const Statedef& statedef_,
                             const osh::LOs& t_species_per_elements,
                             const osh::LOs& t_substates_per_complexes,
                             ComplexIndexer& indexer,
                             const bool with_occupancy = true)
        : mesh(mesh_)
        , statedef(statedef_)
        , pools_(t_species_per_elements)
        , clamped_(pools_.num_data(), false)
        , species_per_elements_(t_species_per_elements)
        , substates_per_complexes(t_substates_per_complexes)
        , substates_shifts(t_substates_per_complexes.size() + 1)
        , pComplexStates(substates_per_complexes.size())
        , pComplexLocations(substates_per_complexes.size())
        , pFilters(substates_per_complexes.size())
        , pFiltersMap(substates_per_complexes.size())
        , occupancy_rd_(with_occupancy ? pools_.num_data() : 0)
        , occupancy_ef_(with_occupancy ? pools_.num_data() : 0)
        , complex_occupancy_filters_per_elements(t_species_per_elements.size())
        , complex_occupancy_shifts(t_species_per_elements.size() + 1)
        , complex_indexer(indexer) {
        std::partial_sum(substates_per_complexes.begin(),
                         substates_per_complexes.end(),
                         substates_shifts.begin() + 1);
    }

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_rd(const Entity entity, const container::species_id species) {
        if (!occupancy_rd_.empty()) {
            const size_t index = pools_.ab(entity.get(), species.get());
            occupancy_rd_.track(index);
        } else {
            throw std::logic_error(
                "You asked for occupancy tracking in an EntityMolecules built without "
                "occupancy "
                "enabled");
        }
    }

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_ef(const Entity entity, const container::species_id species) {
        if (!occupancy_ef_.empty()) {
            const size_t index = pools_.ab(entity.get(), species.get());
            occupancy_ef_.track(index);
        } else {
            throw std::logic_error(
                "You asked for occupancy tracking in an EntityMolecules built without "
                "occupancy "
                "enabled");
        }
    }

    /// Activate occupancy tracking for the pair entity/complex filter
    inline model::complex_filter_occupancy_id track_complex_occupancy_ef(
        const Entity entity,
        const ComplexFilterDescr& filt) {
        auto filter = GetComplexFilter(filt.complexId, filt.filters);
        if (not filter->occupancy_id()) {
            filter->setOccupancyId(
                model::complex_filter_occupancy_id(complex_occupancy_filters.size()));
            complex_occupancy_filters.container().emplace_back(filt.complexId, filter->id());
        }
        auto fid = *filter->occupancy_id();
        complex_occupancy_filters_per_elements[entity].push_back(fid);
        return fid;
    }

    /** Add molecule quantity
     *
     * Necessary for diffusions where we do not care for occupancy
     *
     * @param entity
     * @param species
     * @param val
     */
    template <bool CheckRegionClamped = true>
    inline void add(const Entity entity,
                    const container::species_id species,
                    const molecules_t val) noexcept {
        assert(species.get() < numSpecies(entity));
        if (clamped_[pools_.ab(entity.get(), species.get())]) {
            return;
        }
        if constexpr (CheckRegionClamped) {
            if (check_region_clamped(entity, species)) {
                return;
            }
        }
        pools_(entity.get(), species.get()) += val;
    }

    /** Add molecule quantity and update occupancy
     *
     * Check Occupancy for more information
     *
     * This should be called only if occupancy exists
     *
     * @param entity
     * @param species
     * @param val
     * @param event_time
     */
    template <bool CheckRegionClamped = true>
    inline void add_and_update_occupancy(const Entity entity,
                                         const container::species_id species,
                                         const molecules_t val,
                                         const osh::Real event_time) noexcept {
        const size_t index = pools_.ab(entity.get(), species.get());
        if (clamped_[index]) {
            return;
        }
        if constexpr (CheckRegionClamped) {
            if (check_region_clamped(entity, species)) {
                return;
            }
        }
        pools_(entity.get(), species.get()) += val;
        assert(!occupancy_rd_.empty());
        assert(!occupancy_ef_.empty());

        occupancy_rd_.add_correction(index, val, event_time);
        occupancy_ef_.add_correction(index, val, event_time);
    }

    /// Get the clamped flag for a species in an entity
    bool get_clamped(const Entity entity, const container::species_id species) const noexcept {
        assert(species.get() < numSpecies(entity));
        return clamped_[pools_.ab(entity.get(), species.get())];
    }

    /// Set the clamped flag for a species in an entity
    void set_clamped(const Entity entity,
                     const container::species_id species,
                     bool clamped) noexcept {
        assert(species.get() < numSpecies(entity));
        clamped_[pools_.ab(entity.get(), species.get())] = clamped;
    }

    /// Return whether a species in clamped in the region associated to the entity
    inline bool check_region_clamped(const Entity entity,
                                     const container::species_id species) const noexcept {
        auto regionId = mesh.getRegionContID(entity);
        if (regionId.valid()) {
            const auto& regiondef = statedef.getDefinition(regionId);
            return regiondef.getSpecClamped(species);
        }
        return false;
    }

    /** Return whether there are negative counts
     */
    inline bool has_negative_counts() const noexcept {
        for (const auto& elem: entities()) {
            for (const auto& spec: species(elem)) {
                if (pools_(elem.get(), spec.get()) < 0) {
                    return true;
                }
            }
        }
        return false;
    }

    /** Get occupancy based on the reaction-diffusion time step
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_rd(const Entity entity,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        assert(!occupancy_rd_.empty());
        const size_t index = pools_.ab(entity.get(), species.get());
        return occupancy_rd_.get_occupancy(index, pools_(entity.get(), species.get()), end_time);
    }

    /** Get occupancy based on the efield time step
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_ef(const Entity entity,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        assert(!occupancy_ef_.empty());
        const size_t index = pools_.ab(entity.get(), species.get());
        return occupancy_ef_.get_occupancy(index, pools_(entity.get(), species.get()), end_time);
    }

    inline osh::Real get_occupancy_ef(const Entity entity,
                                      const model::complex_filter_occupancy_id& fid,
                                      const osh::Real end_time) const {
        assert(complex_occupancy_ef_ and !complex_occupancy_ef_->empty());
        const size_t index = complex_occupancy_shifts[entity] + fid.get();
        const auto fullId = complex_occupancy_filters[fid];
        auto filter = pFilters[fullId.complexId][fullId.filterId];
        return complex_occupancy_ef_->get_occupancy(index,
                                                    (*this)(entity, fullId.complexId, *filter),
                                                    end_time);
    }

    /** Assign value to the pools
     *
     * We do not need to update/invalidate occupancy because this function cannot/should not be
     * used in the middle of a time step
     *
     * Warning: do not use in the middle of a time step
     *
     * @param entity
     * @param species
     * @param val
     */
    inline void assign(const Entity entity,
                       const container::species_id species,
                       const molecules_t val) noexcept {
        assert(species.get() < numSpecies(entity));
        pools_(entity.get(), species.get()) = val;
    }

    void assign(const Entity entity,
                const model::complex_id complex,
                const util::strongid_vector<model::complex_substate_id, uint>& i,
                const molecules_t val) {
        std::vector<model::complex_individual_id> existingInds;
        auto [bgn, end] = pComplexLocations[complex].equal_range(entity);
        for (auto it = bgn; it != end; ++it) {
            auto it2 = pComplexStates[complex].find(it->second);
            if (it2 != pComplexStates[complex].end() and it2->second == i) {
                existingInds.push_back(it->second);
            }
        }
        // If we need to remove some of the complexes
        if (static_cast<molecules_t>(existingInds.size()) > val) {
            uint nbToDel = existingInds.size() - val;
            for (auto it = existingInds.begin(); it != existingInds.begin() + nbToDel; ++it) {
                removeComplex(complex, *it);
            }
        } else {
            for (uint nb = 0; nb < val - existingInds.size(); ++nb) {
                addComplex(entity, complex, complex_indexer.getNextComplexInd(), i);
            }
        }
    }

    void removeComplex(model::complex_id cmplIdx, model::complex_individual_id stIdx) {
        removeComplexUpdateOccupancy(cmplIdx, stIdx, -1);
    }

    void removeComplexUpdateOccupancy(model::complex_id cmplIdx,
                                      model::complex_individual_id stIdx,
                                      const osh::Real event_time) {
        assert(cmplIdx.get() < pComplexStates.container().size());
        assert(pComplexStates[cmplIdx].count(stIdx) > 0);

        // Remove the complex state
        auto entity = pComplexStates[cmplIdx][stIdx].location();
        auto [bgn, end] = pComplexLocations[cmplIdx].equal_range(entity);
        for (auto it = bgn; it != end; ++it) {
            if (it->second == stIdx) {
                pComplexLocations[cmplIdx].erase(it);
                break;
            }
        }
        pComplexStates[cmplIdx].erase(stIdx);

        // Compute occupancy changes
        std::set<model::complex_filter_id> processedFilters;
        if (event_time >= 0 and complex_occupancy_ef_) {
            for (const auto& fid: complex_occupancy_filters_per_elements[entity]) {
                const auto fullId = complex_occupancy_filters[fid];
                if (fullId.complexId == cmplIdx) {
                    const size_t index = complex_occupancy_shifts[entity] + fid.get();
                    auto filter = pFilters[fullId.complexId][fullId.filterId];
                    if (filter->matches(stIdx)) {
                        // Matched before removal
                        filter->removeMatch(stIdx);

                        complex_occupancy_ef_->add_correction(index, -1, event_time);
                    }
                    processedFilters.insert(fullId.filterId);
                }
            }
        }

        // Only add updates to filters that were not processed during the occupancy update
        for (auto filt: pFilters[cmplIdx]) {
            if (processedFilters.count(filt->id()) == 0) {
                filt->toUpdate(stIdx);
            }
        }
    }

    void addComplex(const Entity entity,
                    model::complex_id cmplIdx,
                    model::complex_individual_id stIdx,
                    const util::strongid_vector<model::complex_substate_id, uint>& init) {
        addComplexUpdateOccupancy(entity, cmplIdx, stIdx, init, -1);
    }

    void addComplexUpdateOccupancy(
        const Entity entity,
        model::complex_id cmplIdx,
        model::complex_individual_id stIdx,
        const util::strongid_vector<model::complex_substate_id, uint>& init,
        const osh::Real event_time) {
        assert(cmplIdx.get() < pComplexStates.container().size());

        // Add the complex state
        auto it = pComplexStates[cmplIdx].emplace(stIdx, ComplexState(init, stIdx, entity)).first;
        pComplexLocations[cmplIdx].emplace(entity, stIdx);

        // Compute occupancy changes
        std::set<model::complex_filter_id> processedFilters;
        if (event_time >= 0 and complex_occupancy_ef_) {
            for (const auto& fid: complex_occupancy_filters_per_elements[entity]) {
                const auto fullId = complex_occupancy_filters[fid];
                if (fullId.complexId == cmplIdx) {
                    const size_t index = complex_occupancy_shifts[entity] + fid.get();
                    auto filter = pFilters[fullId.complexId][fullId.filterId];
                    if (filter->computeMatch(it->second)) {
                        filter->addMatch(stIdx);

                        complex_occupancy_ef_->add_correction(index, 1, event_time);
                    }
                    processedFilters.insert(fullId.filterId);
                }
            }
        }

        // Only add updates to filters that were not processed during the occupancy update
        for (auto filt: pFilters[cmplIdx]) {
            if (processedFilters.count(filt->id()) == 0) {
                filt->toUpdate(stIdx);
            }
        }
    }

    void createComplex(const Entity entity,
                       model::complex_id cmplIdx,
                       const util::strongid_vector<model::complex_substate_id, uint>& init) {
        addComplexUpdateOccupancy(entity, cmplIdx, complex_indexer.getNextComplexInd(), init, -1);
    }

    void createComplexUpdateOccupancy(
        const Entity entity,
        model::complex_id cmplIdx,
        const util::strongid_vector<model::complex_substate_id, uint>& init,
        const osh::Real event_time) {
        addComplexUpdateOccupancy(
            entity, cmplIdx, complex_indexer.getNextComplexInd(), init, event_time);
    }

    void updateComplex(model::complex_id cmplIdx,
                       model::complex_individual_id stIdx,
                       const util::strongid_vector<model::complex_substate_id, int>& upd) {
        updateComplexUpdateOccupancy(cmplIdx, stIdx, upd, -1);
    }

    void updateComplexUpdateOccupancy(
        model::complex_id cmplIdx,
        model::complex_individual_id stIdx,
        const util::strongid_vector<model::complex_substate_id, int>& upd,
        const osh::Real event_time) {
        assert(cmplIdx.get() < pComplexStates.container().size());
        auto it = pComplexStates[cmplIdx].find(stIdx);
        assert(it != pComplexStates[cmplIdx].end());

        // Apply the update to the state
        auto& state = it->second;
        std::vector<model::complex_substate_id> modifInds;
        for (auto sus: upd.range()) {
            if (upd[sus] != 0) {
                modifInds.push_back(sus);
                state[sus] += upd[sus];
            }
        }

        // Compute occupancy changes
        std::set<model::complex_filter_id> processedFilters;
        Entity entity = state.location();
        if (event_time >= 0 and complex_occupancy_ef_) {
            for (const auto& fid: complex_occupancy_filters_per_elements[entity]) {
                const auto fullId = complex_occupancy_filters[fid];
                if (fullId.complexId == cmplIdx) {
                    const size_t index = complex_occupancy_shifts[entity] + fid.get();
                    auto filter = pFilters[fullId.complexId][fullId.filterId];
                    bool couldChange =
                        filter->matchAll() or
                        std::any_of(modifInds.begin(), modifInds.end(), [&](auto sus) {
                            return filter->dependsOnSus(sus);
                        });
                    if (couldChange) {
                        if (filter->matches(stIdx)) {
                            // Matched before the changes...
                            if (not filter->computeMatch(state)) {
                                // ... but not after
                                filter->removeMatch(stIdx);

                                complex_occupancy_ef_->add_correction(index, -1, event_time);
                            }
                        } else if (filter->computeMatch(state)) {
                            // Matched after the changes but not before
                            filter->addMatch(stIdx);

                            complex_occupancy_ef_->add_correction(index, 1, event_time);
                        }
                    }
                    processedFilters.insert(fullId.filterId);
                }
            }
        }

        // Only add updates to filters that were not processed during the occupancy update
        if (not modifInds.empty()) {
            for (auto filt: pFilters[cmplIdx]) {
                if (processedFilters.count(filt->id()) == 0 and
                    std::any_of(modifInds.begin(), modifInds.end(), [&](auto sus) {
                        return filt->dependsOnSus(sus);
                    })) {
                    filt->toUpdate(stIdx);
                }
            }
        }
    }

    std::shared_ptr<ComplexFilter<Entity>> GetComplexFilter(
        const model::complex_id& cmplIdx,
        const std::vector<util::strongid_vector<model::complex_substate_id,
                                                steps::model::SubunitStateFilter>>& filts) const {
        auto it = pFiltersMap[cmplIdx].find(filts);
        if (it == pFiltersMap[cmplIdx].end()) {
            model::complex_filter_id fid(pFilters[cmplIdx].size());
            pFilters[cmplIdx].container().emplace_back(
                new ComplexFilter(filts, fid, pComplexStates[cmplIdx]));
            it = pFiltersMap[cmplIdx].emplace(filts, pFilters[cmplIdx].back()).first;
        }
        return it->second;
    }

    const ComplexFilter<Entity> updatedComplexFilter(
        const model::complex_id& cmplIdx,
        const std::vector<util::strongid_vector<model::complex_substate_id,
                                                steps::model::SubunitStateFilter>>& filts) const {
        const auto filt = GetComplexFilter(cmplIdx, filts);
        filt->processUpdates(pComplexStates[cmplIdx]);
        return *filt;
    }

    const std::unordered_map<model::complex_individual_id, ComplexState<Entity>>& complexStates(
        model::complex_id cidx) const {
        return pComplexStates[cidx];
    }

    const std::multimap<Entity, model::complex_individual_id>& complexLocations(
        model::complex_id cidx) const {
        return pComplexLocations[cidx];
    }

    /// Get a copy of the molecule quantity
    inline molecules_t operator()(Entity entity, container::species_id species) const noexcept {
        assert(species.get() < numSpecies(entity));
        return pools_(entity.get(), species.get());
    }

    molecules_t operator()(Entity entity,
                           model::complex_id complex,
                           const ComplexFilter<Entity>& filt) const noexcept {
        assert(filt.isUpdated());
        auto [bgn, end] = pComplexLocations[complex].equal_range(entity);
        return std::count_if(bgn, end, [&filt](auto& loc) { return filt.matches(loc.second); });
    }

    molecules_t operator()(Entity entity,
                           model::complex_id complex,
                           const ComplexFilter<Entity>& filt,
                           model::complex_substate_id sus) const noexcept {
        assert(filt.isUpdated());
        auto [bgn, end] = pComplexLocations[complex].equal_range(entity);
        const auto& states = pComplexStates[complex];
        return std::accumulate(bgn, end, 0, [&](auto acc, auto& loc) -> molecules_t {
            if (filt.matches(loc.second)) {
                return acc + states.at(loc.second)[sus];
            } else {
                return acc;
            }
        });
    }

    /// Check if a pool is empty
    inline bool empty(Entity entity, container::species_id species) const noexcept {
        assert(species.get() < numSpecies(entity));
        return pools_(entity.get(), species.get()) == 0;
    }

    /** Full Reset
     *
     * All integrals start from state_time
     */
    inline void reset(const osh::Real state_time) {
        pools_.assign(0);
        std::fill(clamped_.begin(), clamped_.end(), false);
        reset_occupancy_rd(state_time);
        reset_occupancy_ef(state_time);
        for (const auto& cmplIdx: pComplexStates.range()) {
            pComplexStates[cmplIdx].clear();
            pComplexLocations[cmplIdx].clear();
            for (auto& filt: pFilters[cmplIdx]) {
                filt->reset();
            }
        }
    }

    /** Reset reaction-diffusion based occupancy
     *
     * @param current_time: necessary for occupancy. The time integrals start from here
     */
    inline void reset_occupancy_rd(const osh::Real current_time) {
        occupancy_rd_.reset(current_time);
    }

    /** Reset efield based occupancy
     *
     * @param current_time: necessary for occupancy. The time integrals start from here
     */
    inline void reset_occupancy_ef(const osh::Real current_time) {
        occupancy_ef_.reset(current_time);
        if (complex_occupancy_ef_) {
            complex_occupancy_ef_->reset(current_time);
            // Update the filters that are used with occupancy
            for (const auto& fullId: complex_occupancy_filters) {
                auto filter = pFilters[fullId.complexId][fullId.filterId];
                filter->processUpdates(pComplexStates[fullId.complexId]);
            }
        }
    }

    inline osh::LO numEntities() const noexcept {
        return pools_.size();
    }

    inline osh::LO numSpecies(Entity entity) const noexcept {
        return pools_.size(entity.get());
    }

    osh::LO num_complex_data() const {
        return numEntities() * substates_shifts.back();
    }

    inline molecules_t sumNumMolecules(container::species_id species) const {
        molecules_t num_molecules{};
        for (osh::LO elem = 0; elem < numEntities(); ++elem) {
            if (species.get() < numSpecies(elem)) {
                num_molecules += this->operator()(elem, species);
            }
        }
        return num_molecules;
    }

    inline const osh::LOs& species() const noexcept {
        return species_per_elements_;
    }

    inline const osh::LOs& substates() const noexcept {
        return substates_per_complexes;
    }

    inline auto entities() const noexcept {
        return util::EntityIterator<Entity, typename Entity::value_type>(numEntities());
    }

    inline auto species(Entity entity) const noexcept {
        const auto num_species = numSpecies(entity);
        return util::EntityIterator<container::species_id, container::species_id::value_type>(
            num_species);
    }

    inline auto ab(Entity entity, container::species_id species) const noexcept {
        return pools_.ab(entity.get(), species.get());
    }

    inline auto ab(Entity entity,
                   model::complex_id complex,
                   model::complex_substate_id sus) const noexcept {
        return substates_shifts.back() * entity.get() + substates_shifts[complex.get()] + sus.get();
    }

    inline auto num_data() const noexcept {
        return pools_.num_data();
    }

    void finalize_complex_occupancy() {
        std::transform(complex_occupancy_filters_per_elements.begin(),
                       complex_occupancy_filters_per_elements.end(),
                       complex_occupancy_shifts.begin() + 1,
                       [](auto& v) { return v.size(); });
        std::partial_sum(complex_occupancy_shifts.begin(),
                         complex_occupancy_shifts.end(),
                         complex_occupancy_shifts.begin());
        complex_occupancy_ef_ = Occupancy(complex_occupancy_shifts.back());
        for (const auto& element_id: complex_occupancy_filters_per_elements.range()) {
            for (const auto& fid: complex_occupancy_filters_per_elements[element_id]) {
                complex_occupancy_ef_->track(complex_occupancy_shifts[element_id] + fid.get());
            }
        }
    }

  private:
    const DistMesh& mesh;
    const Statedef& statedef;

    /** Number of molecules/channels per element and species **/
    util::flat_multimap<molecules_t, 1> pools_;
    std::vector<bool> clamped_;

    osh::LOs species_per_elements_;

    osh::LOs substates_per_complexes;
    std::vector<osh::LO> substates_shifts;

    // Complexes
    util::strongid_vector<model::complex_id,
                          std::unordered_map<model::complex_individual_id, ComplexState<Entity>>>
        pComplexStates;

    util::strongid_vector<model::complex_id, std::multimap<Entity, model::complex_individual_id>>
        pComplexLocations;

    mutable util::strongid_vector<
        model::complex_id,
        util::strongid_vector<model::complex_filter_id, std::shared_ptr<ComplexFilter<Entity>>>>
        pFilters;

    mutable util::strongid_vector<
        model::complex_id,
        std::unordered_map<std::vector<util::strongid_vector<model::complex_substate_id,
                                                             steps::model::SubunitStateFilter>>,
                           std::shared_ptr<ComplexFilter<Entity>>,
                           FilterHash>>
        pFiltersMap;

    /** Occupancy based on molecules (based on rd dt)
     *
     * The dt is given as end_time (given when you ask for the occupancy) - start_time (given
     * when occupancy is reset)
     */
    Occupancy occupancy_rd_;

    /** Occupancy based on channel states (based on ef dt)
     *
     * The dt is given as end_time (given when you ask for the occupancy) - start_time (given
     * when occupancy is reset)
     */
    Occupancy occupancy_ef_;

    std::optional<Occupancy> complex_occupancy_ef_;
    util::strongid_vector<model::complex_filter_occupancy_id, model::complex_filter_full_id>
        complex_occupancy_filters;
    util::strongid_vector<Entity, std::vector<model::complex_filter_occupancy_id>>
        complex_occupancy_filters_per_elements;
    util::strongid_vector<Entity, osh::LO> complex_occupancy_shifts;

    ComplexIndexer& complex_indexer;

    static constexpr bool is_osh_integral = std::is_same<molecules_t, osh::LO>::value ||
                                            std::is_same<molecules_t, osh::GO>::value;
    static_assert(is_osh_integral, "Expected type Omega_h::LO or Omega_h::GO");
};


/////////////////////////////////////////

using ElementsMolecules = EntityMolecules<mesh::tetrahedron_id_t>;
using BoundariesMolecules = EntityMolecules<mesh::triangle_id_t>;

using MolStateElementID =
    std::tuple<std::variant<mesh::tetrahedron_id_t, mesh::triangle_id_t>, container::species_id>;
using MolStateComplexElementID =
    std::tuple<std::variant<mesh::tetrahedron_id_t, mesh::triangle_id_t>,
               model::complex_id,
               model::complex_substate_id>;

using MolStateElementQuantity = std::
    tuple<std::variant<mesh::tetrahedron_id_t, mesh::triangle_id_t>, container::species_id, int>;

class MolState {
  public:
    enum class Location : char { volume, boundary };
    /// field 1: kind of entity
    /// field 2: entity identifier
    /// field 3: species identifier
    using ElementID = MolStateElementID;

    explicit MolState(const DistMesh& _mesh,
                      const Statedef& _statedef,
                      const osh::LOs& t_species_per_elements,
                      const osh::LOs& substates_per_complexes,
                      RankInfo info,
                      const bool with_occupancy = true,
                      const std::optional<osh::LOs>& t_species_per_boundary_element = std::nullopt)
        : mesh_(_mesh)
        , statedef_(_statedef)
        , complex_indexer(info)
        , molecules_on_elements_(_mesh,
                                 _statedef,
                                 t_species_per_elements,
                                 substates_per_complexes,
                                 complex_indexer,
                                 with_occupancy)
        , molecules_on_patch_boundaries_(_mesh,
                                         _statedef,
                                         t_species_per_boundary_element
                                             ? (*t_species_per_boundary_element)
                                             : osh::Write<osh::LO>(1, 0),
                                         substates_per_complexes,
                                         complex_indexer,
                                         with_occupancy)
        , elements(molecules_on_elements_.entities())
        , boundaries(molecules_on_patch_boundaries_.entities()) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"

    /// Add molecule quantity and update occupancy
    template <bool CheckRegionClamped = true>
    inline void add(const ElementID& elementId, const molecules_t val) noexcept {
        const auto species = std::get<1>(elementId);
        std::visit([this, species, val](auto entity)
                       -> void { this->add<CheckRegionClamped>(entity, species, val); },
                   std::get<0>(elementId));
    }

    /** Add molecule quantity and update occupancy
     *
     * Check Occupancy for more information
     *
     * @param entity
     * @param species
     * @param val
     * @param event_time
     */
    template <bool CheckRegionClamped = true>
    inline void add_and_update_occupancy(const ElementID& elementId,
                                         const molecules_t val,
                                         const osh::Real event_time) noexcept {
        const auto species = std::get<1>(elementId);
        std::visit(
            [this, species, val, event_time](auto entity) -> void {
                this->add_and_update_occupancy<CheckRegionClamped>(entity,
                                                                   species,
                                                                   val,
                                                                   event_time);
            },
            std::get<0>(elementId));
    }

    /// Get the clamped flag for a species in an entity
    bool get_clamped(const ElementID& elementId) const noexcept {
        const auto species = std::get<1>(elementId);
        return std::visit([this, species](
                              auto entity) -> bool { return this->get_clamped(entity, species); },
                          std::get<0>(elementId));
    }

    /// Set the clamped flag for a species in an entity
    void set_clamped(const ElementID& elementId, bool clamped) noexcept {
        const auto species = std::get<1>(elementId);
        std::visit([this, species, clamped](
                       auto entity) -> void { this->set_clamped(entity, species, clamped); },
                   std::get<0>(elementId));
    }

    /** Return whether there are negative counts
     */
    inline bool has_negative_counts() const noexcept {
        return molecules_on_elements_.has_negative_counts() or
               molecules_on_patch_boundaries_.has_negative_counts();
    }

    /** Assign molecule species on entity
     *
     * Warning: do not use in the middle of a time step
     *
     * @param element
     * @param species
     * @param val
     */
    inline void assign(const ElementID& elementId, const molecules_t val) noexcept {
        const auto species = std::get<1>(elementId);
        std::visit([this, species, val](
                       auto entity) -> void { this->assign(entity, species, val); },
                   std::get<0>(elementId));
    }

    /** Get occupancy reaction-diffusion based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_rd(const ElementID& elementId, const osh::Real end_time) const {
        const auto species = std::get<1>(elementId);
        return std::visit(
            [this, species, end_time](auto entity) -> osh::Real {
                return this->get_occupancy_rd(entity, species, end_time);
            },
            std::get<0>(elementId));
    }

    /** Get occupancy efield based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_ef(const ElementID& elementId, const osh::Real end_time) const {
        const auto species = std::get<1>(elementId);
        return std::visit(
            [this, species, end_time](auto entity) -> osh::Real {
                return this->get_occupancy_ef(entity, species, end_time);
            },
            std::get<0>(elementId));
    }

    /// Returns a copy of the pool
    inline molecules_t operator()(const ElementID& elementId) const noexcept {
        auto species = std::get<1>(elementId);
        return std::visit([this, species](auto entity)
                              -> molecules_t { return this->operator()(entity, species); },
                          std::get<0>(elementId));
    }
#pragma GCC diagnostic pop

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_rd(const mesh::tetrahedron_id_t element,
                                   const container::species_id species) {
        molecules_on_elements_.track_occupancy_rd(element, species);
    }

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_ef(const mesh::tetrahedron_id_t element,
                                   const container::species_id species) {
        molecules_on_elements_.track_occupancy_ef(element, species);
    }

    /// Add to the pools without updating occupancy
    template <bool CheckRegionClamped = true>
    inline void add(const mesh::tetrahedron_id_t element,
                    const container::species_id species,
                    const molecules_t val) noexcept {
        molecules_on_elements_.add<CheckRegionClamped>(element, species, val);
    }

    /** Add molecule quantity and update occupancy
     *
     * Check Occupancy for more information
     *
     * @param element
     * @param species
     * @param val
     * @param event_time
     */
    template <bool CheckRegionClamped = true>
    inline void add_and_update_occupancy(const mesh::tetrahedron_id_t element,
                                         const container::species_id species,
                                         const molecules_t val,
                                         const osh::Real event_time) noexcept {
        molecules_on_elements_.add_and_update_occupancy<CheckRegionClamped>(element,
                                                                            species,
                                                                            val,
                                                                            event_time);
    }

    /// Get the clamped flag for a species in an entity
    bool get_clamped(const mesh::tetrahedron_id_t element,
                     const container::species_id species) const noexcept {
        return molecules_on_elements_.get_clamped(element, species);
    }

    /// Set the clamped flag for a species in an entity
    void set_clamped(const mesh::tetrahedron_id_t element,
                     const container::species_id species,
                     bool clamped) noexcept {
        molecules_on_elements_.set_clamped(element, species, clamped);
    }

    /** Assign molecule species on tetrahedron
     *
     * Warning: do not use in the middle of a time step
     *
     * @param element
     * @param species
     * @param val
     */
    inline void assign(const mesh::tetrahedron_id_t element,
                       const container::species_id species,
                       const molecules_t val) noexcept {
        molecules_on_elements_.assign(element, species, val);
    }

    inline void assign(const mesh::tetrahedron_id_t element,
                       const model::complex_id complex,
                       const util::strongid_vector<model::complex_substate_id, uint>& i,
                       const molecules_t val) noexcept {
        molecules_on_elements_.assign(element, complex, i, val);
    }

    /** Get occupancy reaction-diffusion based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_rd(const mesh::tetrahedron_id_t element,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        return molecules_on_elements_.get_occupancy_rd(element, species, end_time);
    }

    /** Get occupancy efield based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_ef(const mesh::tetrahedron_id_t element,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        return molecules_on_elements_.get_occupancy_ef(element, species, end_time);
    }

    inline molecules_t operator()(mesh::tetrahedron_id_t element,
                                  container::species_id species) const noexcept {
        return molecules_on_elements_(element, species);
    }

    inline molecules_t operator()(mesh::tetrahedron_id_t element,
                                  model::complex_id complex,
                                  const ComplexFilter<mesh::tetrahedron_id_t>& f) const noexcept {
        return molecules_on_elements_(element, complex, f);
    }

    molecules_t operator()(mesh::tetrahedron_id_t element,
                           model::complex_id complex,
                           const ComplexFilter<mesh::tetrahedron_id_t>& f,
                           model::complex_substate_id sus) const noexcept {
        return molecules_on_elements_(element, complex, f, sus);
    }

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_rd(const mesh::triangle_id_t element,
                                   const container::species_id species) {
        molecules_on_patch_boundaries_.track_occupancy_rd(element, species);
    }

    /// Activate occupancy tracking for the pair entity/species
    inline void track_occupancy_ef(const mesh::triangle_id_t element,
                                   const container::species_id species) {
        molecules_on_patch_boundaries_.track_occupancy_ef(element, species);
    }

    /// Activate occupancy tracking for complex filters
    inline model::complex_filter_occupancy_id track_complex_occupancy_ef(
        mesh::triangle_id_t element,
        const ComplexFilterDescr& filt) {
        return molecules_on_patch_boundaries_.track_complex_occupancy_ef(element, filt);
    }

    /// Add to the pools without updating occupancy
    template <bool CheckRegionClamped = true>
    inline void add(const mesh::triangle_id_t element,
                    const container::species_id species,
                    const molecules_t val) noexcept {
        molecules_on_patch_boundaries_.add<CheckRegionClamped>(element, species, val);
    }
    /** Add molecule quantity and update occupancy
     *
     * Check Occupancy for more information
     *
     * @param element
     * @param species
     * @param val
     * @param event_time
     */
    template <bool CheckRegionClamped = true>
    inline void add_and_update_occupancy(const mesh::triangle_id_t element,
                                         const container::species_id species,
                                         const molecules_t val,
                                         const osh::Real event_time) noexcept {
        molecules_on_patch_boundaries_.add_and_update_occupancy<CheckRegionClamped>(element,
                                                                                    species,
                                                                                    val,
                                                                                    event_time);
    }

    /// Get the clamped flag for a species in an entity
    bool get_clamped(const mesh::triangle_id_t element,
                     const container::species_id species) const noexcept {
        return molecules_on_patch_boundaries_.get_clamped(element, species);
    }

    /// Set the clamped flag for a species in an entity
    void set_clamped(const mesh::triangle_id_t element,
                     const container::species_id species,
                     bool clamped) noexcept {
        molecules_on_patch_boundaries_.set_clamped(element, species, clamped);
    }

    /** Assign molecule species on triangle
     *
     * Warning: do not use in the middle of a time step
     *
     * @param element
     * @param species
     * @param val
     */
    inline void assign(const mesh::triangle_id_t element,
                       const container::species_id species,
                       const molecules_t val) noexcept {
        molecules_on_patch_boundaries_.assign(element, species, val);
    }

    inline void assign(const mesh::triangle_id_t element,
                       const model::complex_id complex,
                       const util::strongid_vector<model::complex_substate_id, uint>& i,
                       const molecules_t val) noexcept {
        molecules_on_patch_boundaries_.assign(element, complex, i, val);
    }


    /** Get occupancy reaction-diffusion based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_rd(const mesh::triangle_id_t element,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        return molecules_on_patch_boundaries_.get_occupancy_rd(element, species, end_time);
    }


    /** Get occupancy efield based
     *
     * The occupancy is the average integral of the molecule count (per entity and species) over
     * the time step.
     *
     * An entity is typically a tet or tri
     *
     * @param entity
     * @param species
     * @param end_time: time stamp at the end of the time step
     * @return occupancy
     */
    inline osh::Real get_occupancy_ef(const mesh::triangle_id_t element,
                                      const container::species_id species,
                                      const osh::Real end_time) const {
        return molecules_on_patch_boundaries_.get_occupancy_ef(element, species, end_time);
    }
    inline osh::Real get_occupancy_ef(const mesh::triangle_id_t element,
                                      const model::complex_filter_occupancy_id& filter,
                                      const osh::Real end_time) const {
        return molecules_on_patch_boundaries_.get_occupancy_ef(element, filter, end_time);
    }

    /// Returns a copy of the pool
    inline molecules_t operator()(mesh::triangle_id_t element,
                                  container::species_id species) const noexcept {
        return molecules_on_patch_boundaries_(element, species);
    }

    inline molecules_t operator()(mesh::triangle_id_t element,
                                  model::complex_id complex,
                                  const ComplexFilter<mesh::triangle_id_t>& f) const noexcept {
        return molecules_on_patch_boundaries_(element, complex, f);
    }

    molecules_t operator()(mesh::triangle_id_t element,
                           model::complex_id complex,
                           const ComplexFilter<mesh::triangle_id_t>& f,
                           model::complex_substate_id sus) const noexcept {
        return molecules_on_patch_boundaries_(element, complex, f, sus);
    }


    inline bool empty(mesh::tetrahedron_id_t element,
                      container::species_id species) const noexcept {
        return molecules_on_elements_(element, species) == 0;
    }

    inline void reset(const osh::Real state_time) {
        molecules_on_elements_.reset(state_time);
        molecules_on_patch_boundaries_.reset(state_time);
    }

    inline void reset_occupancy_rd(const osh::Real state_time) {
        molecules_on_elements_.reset_occupancy_rd(state_time);
        molecules_on_patch_boundaries_.reset_occupancy_rd(state_time);
    }

    inline void reset_occupancy_ef(const osh::Real state_time) {
        molecules_on_elements_.reset_occupancy_ef(state_time);
        molecules_on_patch_boundaries_.reset_occupancy_ef(state_time);
    }

    inline osh::LO numElements() const noexcept {
        return molecules_on_elements_.numEntities();
    }

    inline osh::LO numBoundaries() const noexcept {
        return molecules_on_patch_boundaries_.numEntities();
    }

    inline osh::LO numSpecies(mesh::tetrahedron_id_t element) const noexcept {
        return molecules_on_elements_.numSpecies(element);
    }

    /**
     * \copybrief molecules_on_elements_
     */
    inline const auto& moleculesOnElements() const noexcept {
        return molecules_on_elements_;
    }

    inline auto& moleculesOnElements() noexcept {
        return molecules_on_elements_;
    }

    /**
     * \copybrief NumMolecules::molecules_on_patch_boundaries_
     */
    inline const auto& moleculesOnPatchBoundaries() const noexcept {
        return molecules_on_patch_boundaries_;
    }

    inline auto& moleculesOnPatchBoundaries() noexcept {
        return molecules_on_patch_boundaries_;
    }

    inline const osh::LOs& species_per_elements() const noexcept {
        return molecules_on_elements_.species();
    }

    inline const osh::LOs& species_per_boundaries() const noexcept {
        return molecules_on_patch_boundaries_.species();
    }

    auto species(mesh::tetrahedron_id_t element) const noexcept {
        return moleculesOnElements().species(element);
    }

    auto species(mesh::triangle_id_t boundary) const noexcept {
        return moleculesOnPatchBoundaries().species(boundary);
    }

    /**
     * \brief Setup occupancy tracking for complexes
     */
    void finalize_complex_occupancy() {
        molecules_on_elements_.finalize_complex_occupancy();
        molecules_on_patch_boundaries_.finalize_complex_occupancy();
    }

    const DistMesh& mesh() const noexcept {
        return mesh_;
    }

    const Statedef& statedef() const noexcept {
        return statedef_;
    }

  private:
    const DistMesh& mesh_;
    const Statedef& statedef_;

    ComplexIndexer complex_indexer;
    /**
     * \brief Container providing the number of molecules of every species
     * within the elements of the local mesh.
     */
    ElementsMolecules molecules_on_elements_;

    /**
     * \brief Container providing the number of molecules of every specie
     * within the boundaries of the local mesh that belong to a patch.
     */
    BoundariesMolecules molecules_on_patch_boundaries_;

  public:
    const util::EntityIterator<mesh::tetrahedron_id_t, mesh::tetrahedron_id_t::value_type> elements;
    const util::EntityIterator<mesh::triangle_id_t, mesh::triangle_id_t::value_type> boundaries;
};

}  // namespace steps::dist

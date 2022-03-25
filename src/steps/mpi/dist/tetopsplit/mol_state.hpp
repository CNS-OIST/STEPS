#pragma once

#include <boost/optional.hpp>
#include <variant>

#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>

#include "mpi/dist/tetopsplit/definition/fwd.hpp"
#include "util/collections.hpp"
#include "util/flat_multimap.hpp"
#include "util/vocabulary.hpp"


namespace steps {
namespace dist {

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
    template <class NumMolecules>
    inline void add_correction(size_t index, const NumMolecules val, const osh::Real event_time) {
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
    template <class NumMolecules>
    inline osh::Real get_occupancy(const size_t index,
                                   const NumMolecules pool,
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

    /** Keeps track of molecules per species per entity (element or boundary)
     *
     * \tparam Entity a strong_id type, \c mesh::element_id or \c
     * mesh::triangle_id_t for instance
     *
     * \tparam NumMolecules integral type used to
     * store a number of molecules. May be \c Omega_h::LO or \c Omega_h::GO
     */
    template <typename Entity, typename NumMolecules>
    struct EntityMolecules {
        using entity_t = Entity;
        using molecules_t = NumMolecules;

        explicit EntityMolecules(const osh::LOs& t_species_per_elements,
                                 const bool with_occupancy = true)
            : pools_(t_species_per_elements)
            , species_per_elements_(t_species_per_elements)
            , occupancy_rd_(with_occupancy ? pools_.num_data() : 0)
            , occupancy_ef_(with_occupancy ? pools_.num_data() : 0) {}

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


        /** Add molecule quantity
         *
         * Necessary for diffusions where we do not care for occupancy
         *
         * @param entity
         * @param species
         * @param val
         */
        inline void add(const Entity entity,
                        const container::species_id species,
                        const molecules_t val) noexcept {
            assert(species.get() < numSpecies(entity));
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
        inline void add_and_update_occupancy(const Entity entity,
                                             const container::species_id species,
                                             const molecules_t val,
                                             const osh::Real event_time) noexcept {
            add(entity.get(), species.get(), val);
            assert(!occupancy_rd_.empty());
            assert(!occupancy_ef_.empty());
            const size_t index = pools_.ab(entity.get(), species.get());

            occupancy_rd_.add_correction(index, val, event_time);
            occupancy_ef_.add_correction(index, val, event_time);
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
            return occupancy_rd_.get_occupancy(index,
                                               pools_(entity.get(), species.get()),
                                               end_time);
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
            return occupancy_ef_.get_occupancy(index,
                                               pools_(entity.get(), species.get()),
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

        /// Get a copy of the molecule quantity
        inline molecules_t operator()(Entity entity, container::species_id species) const noexcept {
            assert(species.get() < numSpecies(entity));
            return pools_(entity.get(), species.get());
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
            reset_occupancy_rd(state_time);
            reset_occupancy_ef(state_time);
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
        }


        inline osh::LO numEntities() const noexcept {
            return pools_.size();
        }

        inline osh::LO numSpecies(Entity entity) const noexcept {
            return pools_.size(entity.get());
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

        inline auto entities() const noexcept {
            return util::EntityIterator<Entity, typename Entity::value_type>(numEntities());
        }

        inline auto species(Entity entity) const noexcept {
            const auto num_species = numSpecies(entity);
            return util::EntityIterator<container::species_id, container::species_id::value_type>(
                num_species);
        }

      private:
        /** Number of molecules/channels per element and species **/
        util::flat_multimap<molecules_t, 1> pools_;

        osh::LOs species_per_elements_;

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

        static constexpr bool is_osh_integral = std::is_same<molecules_t, osh::LO>::value ||
                                                std::is_same<molecules_t, osh::GO>::value;
        static_assert(is_osh_integral, "Expected type Omega_h::LO or Omega_h::GO");
    };


/////////////////////////////////////////

template <typename NumMolecules>
using ElementsMolecules = EntityMolecules<mesh::tetrahedron_id_t, NumMolecules>;
template <typename NumMolecules>
using BoundariesMolecules = EntityMolecules<mesh::triangle_id_t, NumMolecules>;

using MolStateElementID =
    std::tuple<std::variant<mesh::tetrahedron_id_t, mesh::triangle_id_t>, container::species_id>;

static inline MolStateElementID mkVolumeElement(mesh::tetrahedron_id_t element,
                                                container::species_id species) {
  return std::make_pair(element, species);
}

template <typename NumMolecules> class MolState {
public:
  enum class Location : char { volume, boundary };
  /// field 1: kind of entity
  /// field 2: entity identifier
  /// field 3: species identifier
  using ElementID = MolStateElementID;
  using molecules_t = NumMolecules;
  static inline ElementID mkBoundaryElement(mesh::triangle_id_t element,
                                            container::species_id species) {
    return std::make_pair(element, species);
  }

  explicit MolState(const osh::LOs& t_species_per_elements,
                    const bool with_occupancy = true,
                    const boost::optional<osh::LOs>& t_species_per_boundary_element = boost::none)
      : molecules_on_elements_(t_species_per_elements, with_occupancy)
      , molecules_on_patch_boundaries_(t_species_per_boundary_element
                                           ? (*t_species_per_boundary_element)
                                           : osh::Write<osh::LO>(1, 0),
                                       with_occupancy)
      , elements(molecules_on_elements_.entities())
      , boundaries(molecules_on_patch_boundaries_.entities()) {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"

  /// Activate occupancy tracking for the pair entity/species
  inline bool track_occupancy(const ElementID& elementId) {
      const auto species = std::get<1>(elementId);
      return std::visit([this,
                         species](auto entity) -> bool { this->track_occupancy(entity, species); },
                        std::get<0>(elementId));
  }

  /// Add molecule quantity and update occupancy
  inline void add(const ElementID& elementId, const molecules_t val) noexcept {
      const auto species = std::get<1>(elementId);
      std::visit([this, species, val](auto entity) -> void { this->add(entity, species, val); },
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
  inline void add_and_update_occupancy(const ElementID& elementId,
                                       const molecules_t val,
                                       const osh::Real event_time) noexcept {
      const auto species = std::get<1>(elementId);
      std::visit([this, species, val, event_time](auto entity)
                     -> void { this->add_and_update_occupancy(entity, species, val, event_time); },
                 std::get<0>(elementId));
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
      std::visit([this, species, val](auto entity) -> void { this->assign(entity, species, val); },
                 std::get<0>(elementId));
  }

  /** Get occupancy reaction-diffusion based
   *
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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
      return std::visit([this, species, end_time](auto entity)
                            -> osh::Real { this->get_occupancy_rd(entity, species, end_time); },
                        std::get<0>(elementId));
  }

  /** Get occupancy efield based
   *
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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
      return std::visit([this, species, end_time](auto entity)
                            -> osh::Real { this->get_occupancy_ef(entity, species, end_time); },
                        std::get<0>(elementId));
  }

  /// Returns a copy of the pool
  inline molecules_t operator()(const ElementID &elementId) const noexcept {
    auto species = std::get<1>(elementId);
    return std::visit([this, species](
                          auto entity) -> molecules_t { return this->operator()(entity, species); },
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
  inline void add(const mesh::tetrahedron_id_t element,
                  const container::species_id species,
                  const molecules_t val) noexcept {
      molecules_on_elements_.add(element, species, val);
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
  inline void add_and_update_occupancy(const mesh::tetrahedron_id_t element,
                                       const container::species_id species,
                                       const molecules_t val,
                                       const osh::Real event_time) noexcept {
      molecules_on_elements_.add_and_update_occupancy(element, species, val, event_time);
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

  /** Get occupancy reaction-diffusion based
   *
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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

  /// Add to the pools without updating occupancy
  inline void add(const mesh::triangle_id_t element,
                  const container::species_id species,
                  const molecules_t val) noexcept {
      molecules_on_patch_boundaries_.add(element, species, val);
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
  inline void add_and_update_occupancy(const mesh::triangle_id_t element,
                                       const container::species_id species,
                                       const molecules_t val,
                                       const osh::Real event_time) noexcept {
      molecules_on_patch_boundaries_.add_and_update_occupancy(element, species, val, event_time);
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

  /** Get occupancy reaction-diffusion based
   *
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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
   * The occupancy is the average integral of the molecule count (per entity and species) over the
   * time step.
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

  /// Returns a copy of the pool
  inline molecules_t operator()(mesh::triangle_id_t element,
                                container::species_id species) const noexcept {
    return molecules_on_patch_boundaries_(element, species);
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

  inline osh::LO numSpecies(mesh::tetrahedron_id_t element) const noexcept {
    return molecules_on_elements_.numSpecies(element);
  }

  /**
   * \copybrief molecules_on_elements_
   */
  inline const auto &moleculesOnElements() const noexcept {
    return molecules_on_elements_;
  }

  /**
   * \copybrief NumMolecules::molecules_on_patch_boundaries_
   */
  inline const auto &moleculesOnPatchBoundaries() const noexcept {
    return molecules_on_patch_boundaries_;
  }

  inline const osh::LOs &species_per_elements() const noexcept {
    return molecules_on_elements_.species();
  }

  inline const osh::LOs &species_per_boundaries() const noexcept {
    return molecules_on_patch_boundaries_.species();
  }

  auto species(mesh::tetrahedron_id_t element) const noexcept {
    return moleculesOnElements().species(element);
  }

  auto species(mesh::triangle_id_t boundary) const noexcept {
    return moleculesOnPatchBoundaries().species(boundary);
  }

private:
  /**
   * \brief Container providing the number of molecules of every specie
   * within the elements of the local mesh.
   */
  ElementsMolecules<NumMolecules> molecules_on_elements_;

  /**
   * \brief Container providing the number of molecules of every specie
   * within the boundaries of the local mesh that belong to a patch.
   */
  BoundariesMolecules<NumMolecules> molecules_on_patch_boundaries_;

public:
  const util::EntityIterator<mesh::tetrahedron_id_t, osh::LO> elements;
  const util::EntityIterator<mesh::triangle_id_t, osh::LO> boundaries;
};

} // namespace dist
} // namespace steps

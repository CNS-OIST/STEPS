#pragma once

#include <random>

#include <Omega_h_array.hpp>
#include <Omega_h_shape.hpp>

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_state.hpp"
#include "../mol_state.hpp"
#include "geom/dist/fwd.hpp"
#include "rng/rng.hpp"
#include "util/flat_multimap.hpp"
#include "util/strong_ra.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

class DiffusionOperator {
  public:
    DiffusionOperator(DistMesh& mesh,
                      rng::RNG& t_rng,
                      MolState& t_pools,
                      kproc::Diffusions& t_diffusions,
                      kproc::KProcState& t_kproc_state);
    ~DiffusionOperator() = default;

    void initialize();

    inline bool has_active_diffusions() const {
        return active_diffusions;
    }

    inline osh::Real getDt() const {
        return time_delta;
    }

    void operator()(osh::Real opsplit_period, osh::Real state_time);

    /**
     * \brief Get debugging information
     *
     * Return the following values:
     *   - total_diff_steps: The total number of diffusion events
     *
     * \return a map from value name to value
     */
    std::map<std::string, double> getDebugInfo() const;

    inline osh::I64 getExtent() const noexcept {
        return num_diffusions_;
    }

    inline void setBinomialThreshold(osh::I64 threshold) noexcept {
        diffusion_threshold_ = threshold;
    }

    inline auto getBinomialThreshold() const noexcept {
        return diffusion_threshold_;
    }

    void reset();

  protected:
    /** Compute leaving species on all owned elements
     *
     * @param opsplit_period
     * @param state_time
     */
    void species_leaving_elements(osh::Real opsplit_period, osh::Real state_time);

    /**
     * Update state of owned elements to take into account entering species
     *
     * no need for opsplit_period apparently
     */
    void species_entering_elements();

    /** Compute species leaving a given element (triangle/tetrahedron)
     *
     * @param element
     * @param species
     * @param num_molecules
     * @param opsplit_period
     * @param state_time
     */
    void species_leaving_element(mesh::tetrahedron_id_t element,
                                 container::species_id species,
                                 molecules_t num_molecules,
                                 osh::Real opsplit_period,
                                 osh::Real state_time);

    void species_leaving_element_standard(mesh::tetrahedron_id_t element,
                                          container::species_id species,
                                          molecules_t delta_pool_total,
                                          osh::Real scaled_dcst);

    void species_leaving_element_binomial(mesh::tetrahedron_id_t element,
                                          container::species_id species,
                                          molecules_t delta_pool_total,
                                          osh::Real scaled_dcst);

    /** Compute number of molecules leaving an element (triangle/tetrahedron)
     *
     * \param elem an element identifier
     * \param species the molecule specie
     * @param num_molecules
     * @param sum_rates
     * @param opsplit_period
     * @param state_time
     * \return number of molecules leaving
     */
    molecules_t get_leaving_molecules(mesh::tetrahedron_id_t elem,
                                      container::species_id species,
                                      molecules_t num_molecules,
                                      osh::Real sum_rates,
                                      osh::Real opsplit_period,
                                      osh::Real state_time);

    const DistMesh& mesh;
    rng::RNG& rng;
    MolState& pools;
    kproc::Diffusions& diffusions_;
    // Needed to access the dependency_map and relate (elem_id,species) with KProcIDs (optimization
    // on updating the propensities).
    kproc::KProcState& kproc_state_;
    osh::I64 num_diffusions_{};

    bool active_diffusions{};
    osh::Real time_delta{};
    osh::I64 diffusion_threshold_{10};
    std::uniform_real_distribution<double> ur_distribution;
};


/**
 * \brief Diffusion operator inspired by tau-leaping
 *
 * This class implements a diffusion operator that has a variable diffusion dt.
 * The diffusion dt is computed in the same way as tau in tau-leaping methods but taking
 * into account diffusion fluxes only. In practice, it selects a diffusion dt such that the
 * mean and variance of the change in species population after the diffusion event is kept
 * below some threshold. If that resulting diffusion dt is lower than the default diffusion dt from
 * the normal diffusion operator, it selects the normal diffusion dt instead.
 *
 * The diffusion dt is computed according to:
 * Cao, Y., Gillespie, D. T., & Petzold, L. R. (2006). Efficient step size selection for the
 * tau-leaping simulation method. The Journal of chemical physics, 124(4).
 *
 * Since the resulting diffusion dt can be substantially higher than the normal one, we cannot use
 * the same method for computing the actual diffusion (the probability of going out of each
 * tetrahedron would just plateau to 1 for all tetrahedrons). The diffusion is instead computed
 * by sampling the net flux across all triangles.
 *
 * A dedicated method is used to prevent the apparition of negative species after the sampling of
 * net fluxes.
 */
class TauLeapingDiffusionOperator {
  public:
    TauLeapingDiffusionOperator(DistMesh& mesh,
                                rng::RNG& t_rng,
                                MolState& t_pools,
                                kproc::Diffusions& t_diffusions,
                                kproc::KProcState& t_kproc_state);
    ~TauLeapingDiffusionOperator() = default;

    void initialize();

    void reset();

    void operator()(osh::Real opsplit_period, osh::Real state_time);

    osh::Real getDt();

    inline bool has_active_diffusions() const {
        return active_diffusions;
    }

    inline osh::I64 getExtent() const noexcept {
        return num_diffusions_;
    }

    inline osh::Real getDefaultDt() const {
        return const_dt;
    }

    /**
     * \brief Get debugging information
     *
     * Return the following values:
     *   - total_diff_steps: The total number of net diffusion events
     *   - total_diffleaping_steps: The total number of diffusion steps
     *   - diffleaping_cum_rel_step_sizes: The sum of relative step sizes (actual diffusion dt
     *                                     divided by the standard diffusion dt) over all diffusion
     *                                     steps
     *   - diffleaping_cum_squared_rel_step_sizes: The sum of squared relative step sizes (see
     *                                             above)
     *
     * \return a map from value name to value
     */
    std::map<std::string, double> getDebugInfo() const;

    inline osh::Real getTolerance() const noexcept {
        return tolerance_;
    }
    inline void setTolerance(osh::Real tolerance) {
        ArgErrLogIf(tolerance <= 0, "Tolerance should be strictly positive");
        tolerance_ = tolerance;
    }

    inline osh::Real getNormalApproximationThreshold() const noexcept {
        return tolerance_;
    }
    inline void setNormalApproximationThreshold(osh::Real threshold) {
        ArgErrLogIf(threshold < 0, "The normal approximation threshold should be positive");
        normal_approx_thresh = threshold;
    }

    inline osh::Real getCrankNicolsonThreshold() const noexcept {
        return crank_nicolson_thresh;
    }
    inline void setCrankNicolsonThreshold(osh::Real threshold) {
        ArgErrLogIf(threshold <= 0, "The Crank-Nicolson threshold should be strictly positive");
        crank_nicolson_thresh = threshold;
    }

    inline uint getLeapThreshold() const noexcept {
        return leap_threshold;
    }

    inline void setLeapThreshold(uint leap_thresh) noexcept {
        leap_threshold = leap_thresh;
    }

    inline uint getMaxDtSkips() const noexcept {
        return max_dt_skips;
    }

    inline void setMaxDtSkips(uint max_skips) noexcept {
        max_dt_skips = max_skips;
    }

    inline uint getMinDtFactor() const noexcept {
        return min_dt_factor;
    }

    inline void setMinDtFactor(osh::Real factor) noexcept {
        min_dt_factor = factor;
    }

  protected:
    // Information about a border triangle (i.e. a triangle which sits in between tetrahedrons owned
    // by different ranks)
    struct BorderTriangle {
        BorderTriangle(mesh::triangle_local_id_t tri_,
                       bool owned_,
                       mesh::tetrahedron_local_id_t tet0_,
                       mesh::tetrahedron_local_id_t tet1_,
                       mesh::tetrahedron_internal_id_t itet0_,
                       mesh::tetrahedron_internal_id_t itet1_,
                       unsigned short face0_,
                       unsigned short face1_)
            : tri(tri_)
            , owned(owned_)
            , tet0(tet0_)
            , tet1(tet1_)
            , itet0(itet0_)
            , itet1(itet1_)
            , face0(face0_)
            , face1(face1_) {}

        mesh::triangle_local_id_t tri;  // Local triangle index
        bool owned;                     // Whether the triangle is owned by the current rank

        // Neihboring tetrahedron indexes
        mesh::tetrahedron_local_id_t tet0;
        mesh::tetrahedron_local_id_t tet1;

        // Indexes of neighboring tetrahedrons in the list of border tetrahedrons
        mesh::tetrahedron_internal_id_t itet0;
        mesh::tetrahedron_internal_id_t itet1;

        // Face index of the triangle in each neighboring tetrahedron
        unsigned short face0;
        unsigned short face1;
    };

    // Information about a border tetrahedron (i.e. a tetrahedron that has at least one face that is
    // a border triangle)
    struct BorderTetrahedron {
        int owned_border_faces;
        int total_border_faces;
        short border_faces;
    };

    void sync_out_fluxes(osh::Real end_time);
    void compute_available_border_pop();
    template <bool Border>
    void sample_triangle_fluxes(osh::Real dt);
    void apply_border_diffusion();

    // Utility static methods for setting up data structures
    static osh::LOs diff_species_per_triangle(DistMesh& mesh, const osh::LOs spec_per_elems);
    static osh::LOs neighbs_per_triangle(DistMesh& mesh);
    static util::strongid_vector<mesh::triangle_internal_id_t, BorderTriangle>
    border_triangles_info(DistMesh& mesh);
    static util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t>
    border_triangles_ids(const util::strongid_vector<mesh::triangle_internal_id_t, BorderTriangle>&
                             border_tris_info);
    static util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t>
    internal_triangles(DistMesh& mesh);
    static util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t>
    border_tetrahedrons(
        const util::strongid_vector<mesh::triangle_internal_id_t,
                                    TauLeapingDiffusionOperator::BorderTriangle>& border_triangles);
    static osh::LOs diff_species_per_border_tetrahedron(
        const util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t>&
            border_tets,
        const osh::LOs spec_per_elems);
    static osh::LOs diff_species_per_border_triangle(
        const util::strongid_vector<mesh::triangle_internal_id_t,
                                    TauLeapingDiffusionOperator::BorderTriangle>& border_triangles,
        const osh::LOs spec_per_elems);

    DistMesh& mesh;
    rng::RNG& rng;
    MolState& pools;
    kproc::Diffusions& diffusions_;
    // Needed to access the dependency_map and relate (elem_id,species) with KProcIDs (optimization
    // on updating the propensities).
    kproc::KProcState& kproc_state_;
    // Total number of diffusion events that happened
    osh::I64 num_diffusions_{};
    // Number of diffusion steps that happened
    osh::I64 num_steps{};
    // Jump sizes relative to const_dt for each diffusion step
    std::vector<osh::Real> relative_step_sizes{};

    // Diffusive fluxes going out of the faces of each tetrahedron
    kproc::PoolsOutFluxes out_fluxes;
    // Structure for keeping track of outbound species movement and diffusion propensities
    util::flat_multimap<osh::Real, 2> tet_outbound;

    // Map that associates (triangle, tet neigbor index) to (face index in the tet, compartment id
    // of the tet)
    util::flat_multimap<osh::LO, 2> tri2face;
    // Border triangle information
    util::strongid_vector<mesh::triangle_internal_id_t, BorderTriangle> border_tris_info;
    // Border triangles
    util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t> border_tris;

    // Border tetrahedron indexes
    util::strongid_vector<mesh::tetrahedron_internal_id_t, mesh::tetrahedron_local_id_t>
        border_tets;
    // Available species populations in border tetrahedrons
    kproc::BorderTetrahedronPopulation border_tet_pop;
    // Border tetrahedron information
    util::strongid_vector<mesh::tetrahedron_internal_id_t, BorderTetrahedron> border_tet_info;

    // Sampled net fluxes going across border triangles
    kproc::BorderTriangleDiff border_tri_diff;

    // Internal triangles (tetrahedrons on both sides are owned)
    util::strongid_vector<mesh::triangle_internal_id_t, mesh::triangle_local_id_t> internal_tris;

    // Maximum fraction change in species population for computing diffusion dt
    osh::Real tolerance_{0.1};
    // Minimum number of average transfer per diffusion period for using the normal approximation to
    // the skellam distribution
    double normal_approx_thresh{1.0};
    // Fraction of the default diffusion dt above which Crank-Nicolson should be used
    osh::Real crank_nicolson_thresh{2};

    // Whether there are active diffusion rules
    bool active_diffusions{};
    // The default diffusion dt
    osh::Real const_dt{};
    // Species threshold below which the diffusion dt is set to the default value
    uint leap_threshold{10};
    // How many times we skip the computation dt after a default diffusion dt has been selected
    uint max_dt_skips{5};
    // This value determines how the minimum diffusion dt is computed. 0 means that the lowest
    // diffusion propensity out of all tetrahedron is used to compute the diffusion dt. 1 means that
    // the highest value is used, and 0.5 means that the average value is used.
    // Values between 0 and 0.5 are linearly interpolated, same thing for values between 0.5 and 1.
    osh::Real min_dt_factor{1.0};
};

}  // namespace steps::dist

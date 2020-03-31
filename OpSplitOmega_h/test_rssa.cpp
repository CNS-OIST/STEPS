#include "test_rssa.hpp"

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include <Omega_h_array.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

// clang-format off
#include "opsplit/test/scenario.hpp"
// clang-format on

#include "mesh.hpp"
#include "mesh_utils.hpp"
#include "opsplit/diffdef.hpp"
#include "simulation.hpp"

#undef MPI_Allreduce

namespace zee {

namespace {
/**
 * Main diffusion simulation wrapper
 */
template <osh::Int Dim>
class InternalSimulation {
  public:
    using rng_type = std::mt19937;
    template <SSAMethod SSA>
    using sim_data_type = SimulationData<Dim, SSA, rng_type>;
    using sim_input_type = SimulationInput<rng_type>;

    InternalSimulation(OmegaHMesh<Dim>& t_mesh, rng_type& gen)
        : mesh(t_mesh)
        , rank(mesh.getMesh().comm()->rank())
        , measures(osh::measure_elements_real(&mesh.getMesh()))
        , elems2verts(mesh.getMesh().ask_elem_verts())
        , coords(mesh.getMesh().coords())
        , rng_(gen)
        , elems2boundary_(
              mesh.getMesh().ask_down(mesh.getMesh().dim(), mesh.getMesh().dim() - 1).ab2b)
        , boundary2elems_(mesh.getMesh().ask_up(mesh.getMesh().dim() - 1, mesh.getMesh().dim())) {}

    InternalSimulation(const Simulation<std::mt19937>&) = delete;

    template <SSAMethod SSA>
    void run(sim_input_type& input, const Statedef& statedef) {
        sim_data_type<SSA> data(mesh, statedef, input, false, rng_);
        extent = 0;
        for (auto i = 0; i < input.num_iterations; ++i) {
            data.time_delta = 2.0 / 17.0;
            this->run_step(data);
        }
    }

    void report(const std::string& filename, const MolState& pools) {
        // FIXME(TCL): assume they all have same number of species per element
        for (container::specie_id specie(0); specie < pools.numSpecies(mesh::element_id(0));
             ++specie) {
            osh::Write<osh::Real> concentration(mesh.getMesh().nelems());
            for (auto elem = 0; elem < mesh.getMesh().nelems(); ++elem) {
                concentration[elem] = pools(mesh::element_id(elem), specie) / measures[elem];
            }
            std::ostringstream oss;
            oss << "conc" << specie;
            const std::string name = oss.str();
            mesh.getMesh().add_tag(mesh.getMesh().dim(), name, 1, osh::Reals(concentration));
        }
        osh::vtk::write_parallel(filename, &mesh.getMesh());
    }

    OmegaHMesh<Dim>& mesh;
    const osh::I32 rank;
    const osh::Reals measures;
    const osh::LOs elems2verts;
    const osh::Reals coords;

  private:
    template <SSAMethod SSA>
    void run_step(sim_data_type<SSA>& data) {
        // SSA
        data.ssaOp.run(data.time_delta);  // comment this for pure diffusion
        // no diffusion
    }

    /// random number generator used by simulation
    std::mt19937& rng_;

    const osh::LOs elems2boundary_;
    const osh::Adj boundary2elems_;
    int extent{};
};


#if 0
// zee
/**
 * Provide coordinates of the point where molecules are inserted at the beginning
 * of the simulation
 * \tparam Dim mesh dimension
 */
template <osh::Int Dim>
osh::Vector<Dim> seed_point() {
    return {initializer_list<Dim>(0.5)};
}
#endif  // 0

/**
 * Two empirical distributions v1 and v2.
 * \param v1 a probability distribution
 * \param v2 another probability distribution
 * \return The Kolmogorov-Smirnov stat.
 *   https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
 */
PetscScalar kolmogorovSmirnovStat(const std::vector<osh::LO>& v1, const std::vector<osh::LO>& v2) {
    // Since the distributions of mols are discrete,
    // the stat will not be entirely accurate.
    assert(v1.size() == v2.size() && !v1.empty());
    std::map<osh::LO, PetscScalar> m;
    // here we remove duplicates and aggregate prob. masses
    for (auto v: v1) {
        if (m.count(v)) {
            m[v] += 1.0 / static_cast<PetscScalar>(v1.size());
        } else {
            m[v] = 1.0 / static_cast<PetscScalar>(v1.size());
        }
    }
    for (auto v: v2) {
        if (m.count(v)) {
            m[v] -= 1.0 / static_cast<PetscScalar>(v2.size());
        } else {
            m[v] = -1.0 / static_cast<PetscScalar>(v2.size());
        }
    }
    std::vector<PetscScalar> v1_v2_diff(m.size() + 1);
    size_t count = 0;
    for (auto el: m) {
        // keys are sorted in increasing order
        v1_v2_diff[count + 1] = v1_v2_diff[count] + el.second;
        count++;
    }
    std::transform(v1_v2_diff.begin(), v1_v2_diff.end(), v1_v2_diff.begin(), [](const auto l) {
        return std::abs(l);
    });
    return *std::max_element(v1_v2_diff.begin(), v1_v2_diff.end());
}

}  // namespace

template <osh::Int Dim, SSAMethod SSA>
auto opsplit_simulation(InternalSimulation<Dim>& simulation,
                        const Statedef& statedef,
                        const osh::LO num_time_steps,
                        osh::LO num_simulations) {
    const auto num_elements = simulation.mesh.getMesh().nelems();
    std::vector<std::vector<std::vector<osh::LO>>> sims(static_cast<size_t>(num_elements));
    std::mt19937 rng;
    SimulationInput<decltype(rng)> input(osh::LOs(simulation.mesh.getMesh().nelems(),
                                                  static_cast<osh::LO>(
                                                      statedef.getNumberOfSpecies())),
                                         num_time_steps,
                                         rng);

    // Allocate nested containers of `sims` to prevent reallocation in block below
    for (auto elem = 0u; elem < sims.size(); ++elem) {
        sims[elem].resize(static_cast<size_t>(
            input.pools.numSpecies(mesh::element_id(static_cast<osh::LO>(elem)))));
        for (auto& specie_data: sims[elem]) {
            specie_data.reserve(static_cast<size_t>(num_simulations));
        }
    }

    while (num_simulations-- > 0) {
        for (mesh::element_id elem{}; elem < static_cast<osh::LO>(sims.size()); elem++) {
            // scatter 'random' initial populations
            input.pools(elem, container::specie_id(0)) =
                std::max<osh::LO>(1000000 + 100 * (elem.get() - 500), 0);
            input.pools(elem, container::specie_id(1)) =
                std::max<osh::LO>(2000000 - 100 * (elem.get() - 500), 0);
            input.pools(elem, container::specie_id(2)) =
                std::max<osh::LO>(1700000 + 100 * (elem.get() - 500), 0);
            input.pools(elem, container::specie_id(3)) =
                std::max<osh::LO>(1700000 - 200 * (elem.get() - 500), 0);
        }
        simulation.template run<SSA>(input, statedef);
        for (mesh::element_id elem{}; elem < input.pools.numElements(); elem++) {
            for (container::specie_id specie(0); specie.get() < statedef.getNumberOfSpecies();
                 ++specie) {
                sims[static_cast<size_t>(elem.get())][static_cast<size_t>(specie.get())].push_back(
                    input.pools(elem, specie));
            }
        }
    }
    return sims;
}

template <osh::Int Dim>
int rssa_opsplit_test(OmegaHMesh<Dim>& mesh, const osh::LO num_simulation_paths) {
    // const osh::LO center_index = point_to_elem(mesh.getMesh(), seed_point<Dim>());
    thread_local std::mt19937 gen{std::random_device()()};  // NOLINT(cert-msc51-cpp, cert-msc32-c)
    InternalSimulation<Dim> simulation(mesh, gen);

    Statedef statedef;
    statedef.addComp("comp1");
    mesh.addComp("comp1", model::compartment_label(0));
    statedef.addCompSpecs("comp1", {"A", "B", "C", "D"});
    statedef.addCompReac("comp1", {"A", "B"}, {"C"}, 1e12);
    statedef.addCompReac("comp1", {"A", "C"}, {"D"}, 2e12);
    statedef.addCompDiff("comp1", "A", 1e-6);
    statedef.addCompDiff("comp1", "B", 1e-6);
    statedef.addCompDiff("comp1", "C", 1e-6);
    statedef.addCompDiff("comp1", "D", 1e-6);

    const size_t num_time_steps = 3;
    // run with and without RSSA
    const auto& sims_rssa = opsplit_simulation<Dim, SSAMethod::RSSA>(simulation,
                                                                     statedef,
                                                                     num_time_steps,
                                                                     num_simulation_paths);
    const auto& sims_no_rssa = opsplit_simulation<Dim, SSAMethod::SSA>(simulation,
                                                                       statedef,
                                                                       num_time_steps,
                                                                       num_simulation_paths);

    assert(sims_rssa.size() == sims_no_rssa.size());
    // compare the distributions on every element,
    // return failure if the stat indicates rejection of identical distributions
    for (auto k = 0u; k < sims_rssa.size(); k++) {
        for (container::specie_id specie(0); specie.get() < statedef.getNumberOfSpecies();
             ++specie) {
            auto specie_idx = static_cast<size_t>(specie.get());
            double ks_stat = kolmogorovSmirnovStat(sims_no_rssa[k][specie_idx],
                                                   sims_rssa[k][specie_idx]);
            double rejection_level = 3.0 *
                                     std::sqrt(2.0 / num_simulation_paths);  // it is supposed to be
                                                                             // pretty conservative
            if (ks_stat > rejection_level) {
                std::clog << "\nRSSA TEST FAILURE\n";
                std::clog << "With RSSA:\n";
                for (auto l: sims_rssa[k][specie_idx]) {
                    std::clog << l << " ";
                }
                std::clog << "\nWithout RSSA:\n";
                for (auto l: sims_no_rssa[k][specie_idx]) {
                    std::clog << l << " ";
                }
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}

// explicit template instantiation definitions
template int rssa_opsplit_test(OmegaHMesh<2>& mesh, const osh::LO num_iterations);
template int rssa_opsplit_test(OmegaHMesh<3>& mesh, const osh::LO num_iterations);

};  // namespace zee

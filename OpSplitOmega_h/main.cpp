#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

// clang-format off
#include "opsplit/test/scenario.hpp"
// clang-format on

#include "debug.hpp"
#include "mesh.hpp"
#include "mesh_utils.hpp"
#include "mpitools.hpp"
#include "opsplit/diffdef.hpp"
#include "simulation.hpp"
#include "test_rssa.hpp"

namespace zee {
// see README.md for description

// Execution example on 8 cores. 4 physical CPU, 1 socket with OpenMP 3.0
// export OMP_NUM_THREADS=2
// export OMP_PROC_BIND=true
// mpirun --bind-to core -np 4 ./OpSplitOmega_h

#undef MPI_Allreduce


template <osh::Int Dim>
struct MeshTraits {};

template <>
struct MeshTraits<2> {
    static const int boundary_type = osh::EDGE;
};

template <>
struct MeshTraits<3> {
    static const int boundary_type = osh::FACE;
};


/**
 * Main diffusion simulation wrapper
 */
template <osh::Int Dim, typename RNG>
class InternalSimulation {
  public:
    InternalSimulation(OmegaHMesh<Dim>& t_mesh, RNG& t_rng)
        : mesh(t_mesh)
        , rank(mesh.getMesh().comm()->rank())
        , measures(osh::measure_elements_real(&mesh.getMesh()))
        , elems2verts(mesh.getMesh().ask_elem_verts())
        , coords(mesh.getMesh().coords())
        , rng_(t_rng)
        , elems2boundary_(
              mesh.getMesh().ask_down(mesh.getMesh().dim(), mesh.getMesh().dim() - 1).ab2b)
        , boundary2elems_(mesh.getMesh().ask_up(mesh.getMesh().dim() - 1, mesh.getMesh().dim())) {}

    InternalSimulation(const InternalSimulation<Dim, RNG>&) = delete;

    template <SSAMethod SSA>
    void setup(SimulationData<Dim, SSA, RNG>& data) {
        const auto& n_boundary = mesh.getMesh().nents(MeshTraits<Dim>::boundary_type);
        const auto& boundary_measures = measure_ents_real(&mesh.getMesh(),
                                                          MeshTraits<Dim>::boundary_type,
                                                          osh::LOs(n_boundary, 0, 1),
                                                          this->coords);  // use reference instead
        osh::parallel_for(
            this->mesh.getMesh().nelems(),
            [&data, &boundary_measures,
             this](osh::LO elem_id) __attribute__((always_inline, flatten)) {
                mesh::element_id element(elem_id);
                auto elem2verts = osh::gather_verts<Dim + 1>(elems2verts, element.get());
                auto elem2x = osh::gather_vectors<Dim + 1, Dim>(coords, elem2verts);
                const auto center = barycenter(elem2x);
                const auto& elem2neighbours =
                    gather_neighbours<Dim>(boundary2elems_, elems2boundary_, element.get());
                // SM: the function below is also called above in gather_neighbours
                const auto& elem2boundary = osh::gather_down<3>(elems2boundary_, element.get());
                const auto& elem2boundary_measures = osh::gather_scalars(boundary_measures,
                                                                         elem2boundary);
                for (container::specie_id specie(0); specie < data.pools.numSpecies(element);
                     ++specie) {
                    data.diffusions.rates_sum(element, specie) = 0;
                    for (auto e = 0; e < elem2neighbours.size(); ++e) {
                        const auto neighbour_elem = elem2neighbours[e];
                        if (neighbour_elem >= 0) {
                            elem2verts = osh::gather_verts<Dim + 1>(elems2verts, neighbour_elem);
                            elem2x = osh::gather_vectors<Dim + 1, Dim>(coords, elem2verts);
                            data.diffusions.ith_rate(element, specie, e) =
                                0.000001 /* TODO(TCL) read in diffdef */ *
                                elem2boundary_measures[e] / measures[element.get()] /
                                osh::norm(center - barycenter(elem2x));
                            data.diffusions.rates_sum(element, specie) +=
                                data.diffusions.ith_rate(element, specie, e);
                        }
                    }
                }
            });
    }

    template <SSAMethod SSA>
    void run(SimulationInput<RNG>& input, const Statedef& statedef) {
        SimulationData<Dim, SSA, RNG> data(mesh, statedef, input, false, rng_);
        extent = 0;
        this->setup(data);
        for (auto i = 0; i < input.num_iterations; ++i) {
            if (i % 100 == 0 && rank == 0) {
                std::clog << i << std::endl;
            }
            this->run_step(data);
        }
    }

    void report(const std::string& filename, const MolState& pools, bool volume = true) {
        const NumMolecules& molecules = volume ? pools.moleculesOnElements()
                                               : pools.moleculesOnPatchBoundaries();
        osh::Int dim = volume ? mesh.getMesh().dim() : mesh.getMesh().dim() - 1;
        std::vector<osh::Write<osh::Real>> concentrations;
        for (auto elem = 0; elem < molecules.numElements(); ++elem) {
            size_t sz = static_cast<size_t>(molecules.numSpecies(elem));
            if (concentrations.size() < sz) {
                concentrations.resize(sz, osh::Write<osh::Real>(molecules.numElements(), 0.0));
            }
            for (size_t specie = 0; specie < sz; ++specie) {
                concentrations[specie][elem] =
                    molecules(elem, container::specie_id(static_cast<osh::LO>(specie))) /
                    (volume ? measures[elem] : 1.0);
            }
        }

        for (size_t k = 0; k < concentrations.size(); ++k) {
            std::ostringstream oss;
            oss << (volume ? "conc" : "nummols") << k;
            const std::string name = oss.str();
            std::cout << name << "\n";
            mesh.getMesh().add_tag(dim, name, 1, osh::Reals(concentrations[k]));
        }
        if (volume) {
            osh::vtk::write_parallel(filename, &mesh.getMesh());
        }
    }

    OmegaHMesh<Dim>& mesh;
    const osh::I32 rank;
    const osh::Reals measures;
    const osh::LOs elems2verts;
    const osh::Reals coords;

  private:
    /**
     * Compute number of molecules leaving an element (triangle/tetrahedron)
     * \param data container holding all simulation data
     * \param elem an element identifier
     * \param specie the molecule specie
     * \return number of molecules leaving
     */
    template <SSAMethod SSA>
    int get_leaving_molecules(SimulationData<Dim, SSA, RNG>& data,
                              mesh::element_id elem,
                              container::specie_id specie) {
        PetscScalar mean_population = data.diffusions.dv().occupancy(elem, specie) /
                                      data.time_delta;
        if (mean_population > data.pools(elem, specie)) {
            mean_population = data.pools(elem, specie);  // force this for pure diffusion
        }
        return data.diffusions.total_leaving()(static_cast<int>(mean_population),
                                               data.diffusions.rates_sum(elem, specie),
                                               data.time_delta);
    }

    /**
     * Compute species leaving a given element (triangle/tetrahedron)
     * \param data the simulation state
     * \param elem an element identifier
     */
    template <SSAMethod SSA>
    void species_leaving_element(SimulationData<Dim, SSA, RNG>& data, mesh::element_id elem) {
        if (!this->mesh.isOwned(elem)) {
            return;
        }
        for (container::specie_id specie(0); specie < data.pools.numSpecies(elem); ++specie) {
            if (!data.pools.empty(elem, specie)) {
                auto delta_pool_total = this->get_leaving_molecules(data, elem, specie);
                data.pools(elem, specie) -= delta_pool_total;
                extent += delta_pool_total;
                if (delta_pool_total <= 10) {  // threshold) {
                    std::uniform_real_distribution<double> distribution(0.0, 1.0);
                    for (; delta_pool_total > 0; --delta_pool_total) {
                        const auto selector = distribution(this->rng_) *
                                              data.diffusions.rates_sum(elem, specie);
                        auto partial_sum_scaled_dcst = 0.0;
                        for (auto e = 0; e < Dim + 1; ++e) {  // loop over boundary/faces
                            // for (auto direction = 0; direction < n_tet_neighbors; direction++) {
                            const auto current_scaled_dcst =
                                data.diffusions.ith_rate(elem, specie, e);
                            if (current_scaled_dcst == 0.0) {
                                continue;
                            }
                            partial_sum_scaled_dcst += current_scaled_dcst;
                            if (selector < partial_sum_scaled_dcst) {
                                data.diffusions.leaving_molecules().increment_ith_delta_pool(elem,
                                                                                             specie,
                                                                                             e,
                                                                                             1);
                                break;
                            }
                        }
                    }
                } else {
                    auto partial_diffusion_sum = data.diffusions.rates_sum(elem, specie);
                    for (auto e = 0; e < Dim + 1; ++e) {  // loop over boundary/faces
                        const auto& probability_e = data.diffusions.ith_rate(elem, specie, e) /
                                                    partial_diffusion_sum;
                        if (probability_e > 0) {
                            int delta_pool_e = delta_pool_total;
                            if (probability_e < 1) {
                                std::binomial_distribution<> leaving_through_e(delta_pool_total,
                                                                               probability_e);
                                delta_pool_e = leaving_through_e(this->rng_);
                                delta_pool_total -= delta_pool_e;
                                partial_diffusion_sum -= data.diffusions.ith_rate(elem, specie, e);
                            }
                            data.diffusions.leaving_molecules().increment_ith_delta_pool(
                                elem, specie, e, delta_pool_e);
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute species entering a given element (triangle/tetrahedron)
     * \param data the simulation data
     * \param elem an element identifier
     */
    template <SSAMethod SSA>
    void species_entering_element(SimulationData<Dim, SSA, RNG>& data, mesh::element_id elem) {
        const auto elem2neighbours =
            gather_neighbours<Dim>(boundary2elems_, elems2boundary_, elem.get());
        for (const auto& neighbour_elem: elem2neighbours) {
            if (neighbour_elem >= 0) {
                const auto neighbour2neighbours =
                    gather_neighbours<Dim>(boundary2elems_, elems2boundary_, neighbour_elem);
                for (auto e2 = 0; e2 < neighbour2neighbours.size(); ++e2) {
                    if (elem == neighbour2neighbours[e2]) {
                        for (container::specie_id specie(0); specie < data.pools.numSpecies(elem);
                             ++specie) {
                            const auto mols = data.diffusions.leaving_molecules().ith_delta_pool(
                                neighbour_elem, specie, e2);
                            if ((data.pools(elem, specie) >
                                 std::numeric_limits<MolState::elem_type>::max() -
                                     std::max(mols, 0)) ||
                                data.pools(elem, specie) < -mols) {
                                throw std::logic_error(
                                    hadoken::scat("Bounds on molecules count violated."));
                            }
                            data.pools(elem, specie) += mols;
                        }
                        break;
                    }
                }
            }
        }
    }

    template <SSAMethod SSA>
    void run_step(SimulationData<Dim, SSA, RNG>& data) {
        // SSA
        data.updateIterationTimeStep();
        data.ssaOp.run(data.time_delta);  // comment this for pure diffusion

        run_diffusion(data);
    }

    template <SSAMethod SSA>
    void run_diffusion(SimulationData<Dim, SSA, RNG>& data) {
        osh::parallel_for(
            this->mesh.getMesh().nelems(),
            [&data, this](osh::LO elem) __attribute__((always_inline, flatten)) {
                this->species_leaving_element(data, mesh::element_id(elem));
            });

        data.diffusions.leaving_molecules().sync_delta_pools();

        osh::parallel_for(
            this->mesh.getMesh().nelems(),
            [&data, this](osh::LO elem) __attribute__((always_inline, flatten)) {
                this->species_entering_element(data, mesh::element_id(elem));
            });
        data.reset();
        data.kproc_state.updateAllPropensities(data.pools);
    }

    /// random number generator used by simulation
    RNG& rng_;

    const osh::LOs elems2boundary_;
    const osh::Adj boundary2elems_;
    int extent{};
};

template <osh::Int Dim, typename RNG>
static void report_num_molecules_in_concentric_rings(const InternalSimulation<Dim, RNG>& simulation,
                                                     const MolState& pools,
                                                     mesh::element_id center_index) {
    std::vector<osh::LO> num_molecules;
    // FIXME(TCL): assume same number of species per element
    const osh::LO num_species = pools.numSpecies(mesh::element_id(0));
    const osh::LO loop_count{10};
    if (simulation.rank == 0) {
        num_molecules.reserve(static_cast<unsigned>(loop_count * num_species));
    }
    const auto& center_elem2verts = osh::gather_verts<Dim + 1>(simulation.elems2verts,
                                                               center_index.get());
    const auto& center_elem2x = osh::gather_vectors<Dim + 1, Dim>(simulation.coords,
                                                                  center_elem2verts);
    const auto& center = barycenter(center_elem2x);
    for (auto loop = 1; loop <= loop_count; ++loop) {
        osh::Write<osh::LO> molecule_counts(num_species, 0);
        const auto radius = loop * 0.01;

        for (auto elem: simulation.mesh.getOwnedElems()) {
            const auto& elem2verts = osh::gather_verts<Dim + 1>(simulation.elems2verts, elem.get());
            const auto& elem2x = osh::gather_vectors<Dim + 1, Dim>(simulation.coords, elem2verts);
            const auto& distance = osh::norm(center - barycenter(elem2x));
            if (distance < radius || loop == loop_count) {
                for (auto specie = 0; specie < pools.numSpecies(elem); ++specie) {
                    molecule_counts[specie] += pools(elem, container::specie_id(specie));
                }
            }
        }
        for (auto molecule_count: molecule_counts) {
            int sum_count = 0;
            MPI_Allreduce(&molecule_count,
                          &sum_count,
                          1,
                          MPI_INT,
                          MPI_SUM,
                          simulation.mesh.getMesh()
                              .comm()
                              ->get_impl());  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)*/
            if (simulation.rank == 0) {
                num_molecules.push_back(sum_count);
            }
        }
    }
    if (simulation.rank == 0) {
        for (auto specie = 0; specie < num_species; ++specie) {
            for (auto loop = 0; loop < loop_count; ++loop) {
                const auto id = static_cast<unsigned>(loop * num_species) +
                                static_cast<unsigned>(specie);
                std::clog << num_molecules[id] << ' ';
            }
            std::clog << std::endl;
        }
    }
}


/**
 * Provide coordinates of the point where molecules are inserted at the beginning
 * of the simulation
 * \tparam Dim mesh dimension
 */
template <osh::Int Dim>
osh::Vector<Dim> seed_point() {
    return {initializer_list<Dim>(0.5)};
}


template <osh::Int Dim, SSAMethod SSA>
int opsplit(OmegaHMesh<Dim>& mesh, osh::LO num_iterations) {
    const auto center_index = point_to_elem(mesh.getMesh(), seed_point<Dim>());
    thread_local std::mt19937 rng;  // NOLINT(cert-msc51-cpp, cert-msc32-c)
    InternalSimulation<Dim, decltype(rng)> simulation(mesh, rng);

    Statedef statedef;
    statedef.addComp("comp1");
    mesh.addComp("comp1", model::compartment_label(0));
    statedef.addCompSpecs("comp1", {"A", "B", "C"});
    statedef.addCompReac("comp1", {"A", "B"}, {"C"}, 1e12);
    statedef.addCompDiff("comp1", "A", 1e-6);
    statedef.addCompDiff("comp1", "B", 1e-6);

    SimulationInput<decltype(rng)> input(osh::LOs(mesh.getMesh().nelems(), 3), num_iterations, rng);
    input.pools(center_index, container::specie_id(0)) = 1000000;
    input.pools(center_index, container::specie_id(1)) = 2000000;

    simulation.template run<SSA>(input, statedef);
    std::ostringstream report_name;
    report_name << "OpSplitOmega_h_" << Dim << 'D';
    simulation.report(report_name.str(), input.pools);
    report_num_molecules_in_concentric_rings(simulation, input.pools, center_index);
    return EXIT_SUCCESS;
}


template <osh::Int Dim>
int surface_reactions_test(OmegaHMesh<Dim>& mesh, int num_iterations) {
    thread_local std::mt19937 gen{std::random_device()()};  // NOLINT(cert-msc51-cpp, cert-msc32-c)

    Statedef statedef;
    statedef.addComp("Compartment");
    mesh.addComp("Compartment", model::compartment_label(100));
    statedef.addPatch("Patch", "Compartment");
    statedef.addPatchSpecs("Patch", {"F", "G", "H"});
    statedef.addCompSpecs("Compartment", {"A", "B", "C", "D", "E", "F"});
    statedef.addCompReac("Compartment", {"A", "B"}, {"C"}, 1e4);
    statedef.addCompReac("Compartment", {"A", "C"}, {"D"}, 2e4);
    statedef.addSurfReac("Patch", {"A", "A"}, {"F", "G"}, {}, {"B"}, {"H"}, {}, 1.e40);

    statedef.addSurfReac("Patch", {}, {"H"}, {}, {}, {}, {}, 1.0e-10);

    statedef.addCompDiff("Compartment", "A", 1e-6);
    statedef.addCompDiff("Compartment", "B", 1e-6);
    statedef.addCompDiff("Compartment", "C", 1e-6);
    statedef.addCompDiff("Compartment", "D", 1e-6);

    if (statedef.getNumberOfSpecies() != 8) {
        return EXIT_FAILURE;
    }
    if (statedef.patchdefs().size() != 1) {
        return EXIT_FAILURE;
    }

    for (const auto& p: statedef.patchdefs()) {
        if (p->reacdefs().size() != 2) {
            return EXIT_FAILURE;
        }
        const auto& reac = p->reacdefs()[0];
        container::specie_id specieId =
            statedef.getCompdef("Compartment").getSpecContainerIdx(statedef.getSpecModelIdx("A"));
        std::map<container::specie_id, osh::LO> stoi =
            reac->getStoichiometry<SReacdef::PoolChangeType ::LHS,
                                   SReacdef::SpecieLocation ::InnerCompartment>();
        if (stoi.size() != 1) {
            return EXIT_FAILURE;
        } else {
            auto h = *stoi.begin();
            if (h.first != specieId) {
                return EXIT_FAILURE;
            }
            if (h.second != -2) {
                return EXIT_FAILURE;
            }
        }
        stoi = reac->getStoichiometry<SReacdef::PoolChangeType ::LHS,
                                      SReacdef::SpecieLocation ::Patch>();

        if (stoi.size() != 2 || stoi[container::specie_id(0)] != -1 ||
            stoi[container::specie_id(1)] != -1) {
            return EXIT_FAILURE;
        }

        stoi = reac->getStoichiometry<SReacdef::PoolChangeType ::LHS,
                                      SReacdef::SpecieLocation ::OuterCompartment>();

        if (!stoi.empty()) {
            return EXIT_FAILURE;
        }

        stoi = reac->getStoichiometry<SReacdef::PoolChangeType ::UPD,
                                      SReacdef::SpecieLocation ::Patch>();

        if (stoi.size() != 3 || stoi[container::specie_id(0)] != -1 ||
            stoi[container::specie_id(1)] != -1 || stoi[container::specie_id(2)] != 1) {
            return EXIT_FAILURE;
        }

        stoi = reac->getStoichiometry<SReacdef::PoolChangeType ::UPD,
                                      SReacdef::SpecieLocation ::InnerCompartment>();

        if (stoi.size() != 2 || stoi[container::specie_id(0)] != -2 ||
            stoi[container::specie_id(1)] != 1) {
            return EXIT_FAILURE;
        }

        stoi = reac->getStoichiometry<SReacdef::PoolChangeType ::UPD,
                                      SReacdef::SpecieLocation ::OuterCompartment>();

        if (!stoi.empty()) {
            return EXIT_FAILURE;
        }
    }
    // const osh::LO center_index = point_to_elem(mesh.getMesh(), );
    thread_local std::mt19937 rng;  // NOLINT(cert-msc51-cpp, cert-msc32-c)

    // Inefficient but for the sake of the test
    KProcState k_proc_state(statedef, mesh);
    auto species_per_owned_elements = k_proc_state.getNumberOfSpeciesPerOwnedElement();
    const auto& species_per_elements = k_proc_state.getNumberOfSpeciesPerElement();
    SimulationInput<decltype(rng)> input(species_per_owned_elements.first,
                                         *species_per_owned_elements.second,
                                         species_per_elements,
                                         num_iterations,
                                         rng);

    for (osh::LO k = 0; k < input.pools.moleculesOnElements().numElements(); k++) {
        for (osh::LO l = 0; l < input.pools.moleculesOnElements().numSpecies(k); ++l) {
            input.pools.moleculesOnElements()(k, container::specie_id(l)) = 100000;
        }
    }
    for (osh::LO k = 0; k < input.pools.moleculesCountOnPatchBoundaries().numElements(); k++) {
        for (osh::LO l = 0; l < input.pools.moleculesCountOnPatchBoundaries().numSpecies(k); ++l) {
            input.pools.moleculesCountOnPatchBoundaries()(k, container::specie_id(l)) = 100000;
        }
    }

    InternalSimulation<Dim, decltype(rng)> simulation(mesh, rng);
    simulation.template run<SSAMethod ::SSA>(input, statedef);
    std::ostringstream report_name;
    report_name << "OpSplitOmega_h_" << Dim << 'D';
    simulation.report(report_name.str(), input.pools, false);
    simulation.report(report_name.str(), input.pools);
    return EXIT_SUCCESS;
}


static timemory_fixture& ti_full_timer() {
    static timemory_fixture timer("OpSplitOmega_h::main", false, false);
    return timer;
}

template <osh::Int Dim>
int splitting_operator(osh::Mesh& mesh, Omega_h::CmdLine& cmdline) {
    {
        int debugRank{-1};
        if (cmdline.parsed("--debug-rank")) {
            debugRank = cmdline.get<int>("--debug-rank", "value");
            if (debugRank >= 0) {
                const auto rank = zee::mpi_comm_rank();
                if (rank == debugRank) {
                    wait_for_gdb();
                }
            }
        }
    }

    const auto mesh_file = cmdline.get<std::string>("square.msh");
    int scenario = -1;
    if (cmdline.parsed("--test")) {
        scenario = cmdline.get<int>("--test", "number");
    }
    if (scenario < 0) {
        PetscScalar scale{1.0};
        if (cmdline.parsed("--scale")) {
            scale = cmdline.get<PetscScalar>("--scale", "value");
        }
        int num_iterations{3000};
        if (cmdline.parsed("--num-iterations")) {
            num_iterations = cmdline.get<int>("--num-iterations", "value");
        }
        OmegaHMesh<Dim> mesh_wrapper(mesh, mesh_file, scale, false);

        if (scenario == -1) {
            if (cmdline.parsed("--use-rssa")) {
                return opsplit<Dim, SSAMethod::RSSA>(mesh_wrapper, num_iterations);
            } else {
                return opsplit<Dim, SSAMethod::SSA>(mesh_wrapper, num_iterations);
            }
        } else if (scenario == -3) {
            return surface_reactions_test<Dim>(mesh_wrapper, num_iterations);
        } else if (scenario == -2) {
            return rssa_opsplit_test<Dim>(mesh_wrapper, num_iterations);
        } else {
            throw std::invalid_argument(hadoken::scat("Unrecognized scenario number: ", scenario));
        }
    }
    PetscScalar threshold = 10.0;
    if (cmdline.parsed("--threshold")) {
        threshold = cmdline.get<int>("--threshold", "value");
    }

    PetscScalar scale{1.0e-6};
    if (cmdline.parsed("--scale")) {
        scale = cmdline.get<PetscScalar>("--scale", "value");
    }

    int comp_label{100};
    if (cmdline.parsed("--comp-label")) {
        comp_label = cmdline.get<int>("--comp-label", "value");
    }

    PetscScalar num_mols_factor{1};
    if (cmdline.parsed("--num-mols-factor")) {
        num_mols_factor = cmdline.get<PetscScalar>("--num-mols-factor", "value");
    }

    bool log_mesh_report{false};
    if (cmdline.parsed("--log-mesh-report")) {
        log_mesh_report = true;
    }

    bool log_statedef_report{false};
    if (cmdline.parsed("--log-statedef-report")) {
        log_statedef_report = true;
    }

    bool log_state_report{false};
    if (cmdline.parsed("--log-state-report")) {
        log_state_report = true;
    }

    if (cmdline.parsed("--log-all-reports")) {
        log_mesh_report = true;
        log_statedef_report = true;
        log_state_report = true;
    }

    bool molecules_pools_force_dist_for_variable_sized = false;
    if (cmdline.parsed("--molecules-pools-force-dist-for-variable-sized")) {
        molecules_pools_force_dist_for_variable_sized = true;
    }

    unsigned rng_seed = std::random_device()();
    if (cmdline.parsed("--rng-seed")) {
        rng_seed = static_cast<unsigned>(cmdline.get<int>("--rng-seed", "value"));
    }

    PetscScalar do_interval{1e-7};
    if (cmdline.parsed("--do-interval")) {
        do_interval = cmdline.get<PetscScalar>("--do-interval", "value");
    }

    PetscScalar end_time{20};
    if (cmdline.parsed("--end-time")) {
        end_time = cmdline.get<PetscScalar>("--end-time", "value");
    }

    int do_inject_in_elt = -1;
    if (cmdline.parsed("--do-inject-in-elt")) {
        do_inject_in_elt = cmdline.get<int>("--do-inject-in-elt", "value");
    }

    int expected_num_reactions{-1};
    if (cmdline.parsed("--expected-reactions")) {
        expected_num_reactions = cmdline.get<int>("--expected-reactions", "value");
    }

    PetscScalar expected_num_reactions_tolerance{0.01};
    if (cmdline.parsed("--expected-reactions-tolerance")) {
        expected_num_reactions_tolerance =
            cmdline.get<PetscScalar>("--expected-reactions-tolerance", "value");
    }

    int expected_num_diffusions{-1};
    if (cmdline.parsed("--expected-diffusions")) {
        expected_num_diffusions = cmdline.get<int>("--expected-diffusions", "value");
    }

    PetscScalar expected_num_diffusions_tolerance{0.02};
    if (cmdline.parsed("--expected-diffusions-tolerance")) {
        expected_num_diffusions_tolerance =
            cmdline.get<PetscScalar>("--expected-diffusions-tolerance", "value");
    }

    OmegaHMesh<Dim> mesh_wrapper(mesh, mesh_file, scale, false);

    const ScenarioInput input(mesh_file,
                              threshold,
                              num_mols_factor,
                              scale,
                              model::compartment_label(comp_label),
                              do_interval,
                              end_time,
                              expected_num_reactions,
                              expected_num_reactions_tolerance,
                              expected_num_diffusions,
                              expected_num_diffusions_tolerance,
                              do_inject_in_elt,
                              log_mesh_report,
                              log_statedef_report,
                              log_state_report,
                              molecules_pools_force_dist_for_variable_sized);

    const auto run_sim = [&](auto& sim) {
        timemory_fixture ti_sim_timer("test_splitting_operator");
        const auto& result = test_splitting_operator(sim, scenario, input);
        result.ti_init_exec_timer = ti_sim_timer.stop();
        result.ti_full_timer = ti_full_timer().stop();
        result.log_summary(sim, input);
        return result.status;
    };

    std::mt19937 rng(rng_seed);

    if (cmdline.parsed("--use-rssa")) {
        OmegaHSimulation<Dim, SSAMethod::RSSA, decltype(rng)> simulation(input, mesh_wrapper, rng);
        return run_sim(simulation);
    } else {
        OmegaHSimulation<Dim, SSAMethod::SSA, decltype(rng)> simulation(input, mesh_wrapper, rng);
        return run_sim(simulation);
    }
}

}  // namespace zee

static inline bool is_power_of_2(int value) {
    return !(value == 0) && !(value & (value - 1));
}

int main(int argc, char** argv) {
    zee::timemory_cout_output(false);
    zee::ti_full_timer().start();

    auto lib = zee::osh::Library(&argc, &argv);
    if (lib.world()->size() > 1 && !is_power_of_2(lib.world()->size())) {
        if (lib.world()->rank() == 0) {
            std::cerr << "Error: Omega_h requires the number of "
                         "processes to be a power of 2.\nAbort.\n";
        }
        return EXIT_FAILURE;
    }
    //    auto mesh = osh::gmsh::read("square.msh", lib.world());
    Omega_h::CmdLine cmdline;
    cmdline.add_arg<std::string>("square.msh");

    auto& test_flag = cmdline.add_flag("--test", "scenario");
    test_flag.add_arg<int>("number");

    auto& threshold_flag = cmdline.add_flag("--threshold", "threshold");
    threshold_flag.add_arg<int>("value");

    auto& scale_flag = cmdline.add_flag("--scale", "scale");
    scale_flag.add_arg<PetscScalar>("value");

    auto& comp_label_flag = cmdline.add_flag("--comp-label", "compLabel");
    comp_label_flag.add_arg<int>("value");

    auto& num_mols_factor_flag = cmdline.add_flag("--num-mols-factor", "numMolsFactor");
    num_mols_factor_flag.add_arg<PetscScalar>("value");

    cmdline.add_flag("--log-mesh-report", "logging");
    cmdline.add_flag("--log-statedef-report", "logging");
    cmdline.add_flag("--log-state-report", "logging");
    cmdline.add_flag("--log-all-reports", "logging");
    cmdline.add_flag("--molecules-pools-force-dist-for-variable-sized", "mpi");

    auto& num_iterations = cmdline.add_flag("--num-iterations", "numIterations");
    num_iterations.add_arg<int>("value");

    auto& rng_seed_flag = cmdline.add_flag("--rng-seed", "RngSeed");
    rng_seed_flag.add_arg<int>("value");

    auto& expected_num_diffusions_flag = cmdline.add_flag("--expected-diffusions", "numDiffusions");
    expected_num_diffusions_flag.add_arg<int>("value");

    auto& expected_num_diffusions_tolerance_flag =
        cmdline.add_flag("--expected-diffusions-tolerance", "tolerance");
    expected_num_diffusions_tolerance_flag.add_arg<PetscScalar>("value");

    auto& expected_num_reactions_flag = cmdline.add_flag("--expected-reactions", "numReactions");
    expected_num_reactions_flag.add_arg<int>("value");

    auto& expected_num_reactions_tolerance_flag = cmdline.add_flag("--expected-reactions-tolerance",
                                                                   "tolerance");
    expected_num_reactions_tolerance_flag.add_arg<PetscScalar>("value");

    // diffusionOnly parameters
    auto& do_interval_flag = cmdline.add_flag("--do-interval", "diffusionOnlyInterval");
    do_interval_flag.add_arg<PetscScalar>("value");

    auto& endtime_flag = cmdline.add_flag("--end-time", "simulationEndTime");
    endtime_flag.add_arg<PetscScalar>("value");

    auto& do_inject_in_elt_flag = cmdline.add_flag("--do-inject-in-elt", "InjectInElement");
    do_inject_in_elt_flag.add_arg<int>("value");

    cmdline.add_flag("--use-rssa", "UseRSSA");

    auto& debug_rank_flag = cmdline.add_flag("--debug-rank", "mpiRank");
    debug_rank_flag.add_arg<int>("value");

    if (!cmdline.parse_final(lib.world(), &argc, argv)) {
        return EXIT_FAILURE;
    }
    // end diffusionOnly parameters

    const auto mesh_file = cmdline.get<std::string>("square.msh");
    auto mesh = zee::osh::gmsh::read(mesh_file, lib.world());
    switch (mesh.dim()) {
    case 2:
        return zee::splitting_operator<2>(mesh, cmdline);
    case 3:
        return zee::splitting_operator<3>(mesh, cmdline);
    default:
        std::cerr << "Unexpected mesh dimension: " << mesh.dim() << '\n';
        return EXIT_FAILURE;
    }
}

#include "steps/mpi/dist/tetopsplit/mol_state.hpp"

#include <iostream>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace osh = Omega_h;


TEST_CASE("MolState_Occupancy") {
    size_t id = 3;
    osh::Real n_molecules = 1;
    int corr = 0;
    osh::Real t = 1.0;
    osh::Real val;

    steps::dist::Occupancy occ(id + 2);
    REQUIRE_NOTHROW(occ.track(id));

    // no correction, result same as the pool
    t += 4.0;
    val = occ.get_occupancy(id, n_molecules, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(n_molecules, 4));

    occ.reset(t);

    // the correction happens at the time we ask for the integral. The last value is not taken into
    // account
    t += 2.0;
    corr = 2;
    n_molecules += corr;
    occ.add_correction(id, corr, t);

    val = occ.get_occupancy(id, n_molecules, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(n_molecules - 2, 4));

    // the previous addition kicks in and has effects
    t += 2.0;
    val = occ.get_occupancy(id, n_molecules, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(2.0, 4));

    // subtraction
    t += 1.0;
    corr = -2;
    n_molecules += corr;
    occ.add_correction(id, corr, t);

    // time to let the integral accumulate and check double result
    t += 3.0;
    val = occ.get_occupancy(id, n_molecules, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(1.75, 4));
    // this id was not marked for tracking, we return the pool. No error on purpose.
    val = occ.get_occupancy(id + 1, n_molecules, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(n_molecules, 4));
}

TEST_CASE("MolState_EntityMolecules") {
    osh::LOs structure = {1, 3, 11, 5};
    steps::dist::EntityMolecules<steps::dist::mesh::tetrahedron_id_t> en_mol(structure);
    osh::Real val;

    steps::dist::mesh::tetrahedron_id_t elem_rd(3);
    steps::dist::container::species_id species_rd(2);
    en_mol.track_occupancy_rd(elem_rd, species_rd);
    steps::dist::mesh::tetrahedron_id_t elem_ef(2);
    steps::dist::container::species_id species_ef(1);
    en_mol.track_occupancy_ef(elem_ef, species_ef);

    // start from 0
    REQUIRE(en_mol(elem_rd, species_rd) == 0);

    osh::Real t(4.0);
    // reset only rd
    en_mol.reset_occupancy_rd(t);

    t += 2.0;
    // add rd event
    en_mol.add_and_update_occupancy(elem_rd, species_rd, 2, t);
    // add ef event
    en_mol.add_and_update_occupancy(elem_ef, species_ef, 2, t);

    t += 2.0;
    // get integral after 2 more seconds
    val = en_mol.get_occupancy_rd(elem_rd, species_rd, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(1.f, 4));
    // get integral after 2 more seconds
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(0.5f, 4));

    // reset ef and check that the integrals are at 0
    en_mol.reset_occupancy_ef(t);
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(2.f, 4));

    // check full reset
    en_mol.add_and_update_occupancy(elem_ef, species_ef, 5, t);
    en_mol.reset(t);
    val = en_mol(elem_rd, species_rd);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(0.f, 4));
    val = en_mol.get_occupancy_rd(elem_rd, species_rd, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(0.f, 4));
    val = en_mol(elem_ef, species_ef);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(0.f, 4));
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    REQUIRE_THAT(val, Catch::Matchers::WithinULP(0.f, 4));
}

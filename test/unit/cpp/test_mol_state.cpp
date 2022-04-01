#include "steps/mpi/dist/tetopsplit/mol_state.hpp"

#include "gtest/gtest.h"

#include <iostream>

namespace osh = Omega_h;


TEST(MolState, Occupancy) {
    size_t id = 3;
    size_t n_molecules = 1;
    int corr = 0;
    osh::Real t = 1.0;
    osh::Real val;

    steps::dist::Occupancy occ(id + 2);
    EXPECT_NO_THROW(occ.track(id));

    // no correction, result same as the pool
    t += 4.0;
    val = occ.get_occupancy(id, n_molecules, t);
    EXPECT_DOUBLE_EQ(val, n_molecules);

    occ.reset(t);

    // the correction happens at the time we ask for the integral. The last value is not taken into
    // account
    t += 2.0;
    corr = 2;
    n_molecules += corr;
    occ.add_correction(id, corr, t);

    val = occ.get_occupancy(id, n_molecules, t);
    EXPECT_DOUBLE_EQ(val, n_molecules - 2);

    // the previous addition kicks in and has effects
    t += 2.0;
    val = occ.get_occupancy(id, n_molecules, t);
    EXPECT_DOUBLE_EQ(val, 2.0);

    // subtraction
    t += 1.0;
    corr = -2;
    n_molecules += corr;
    occ.add_correction(id, corr, t);

    // time to let the integral accumulate and check double result
    t += 3.0;
    val = occ.get_occupancy(id, n_molecules, t);
    EXPECT_DOUBLE_EQ(val, 1.75);
    // this id was not marked for tracking, we return the pool. No error on purpose.
    val = occ.get_occupancy(id + 1, n_molecules, t);
    EXPECT_DOUBLE_EQ(val, n_molecules);
}

TEST(MolState, EntityMolecules) {
    osh::LOs structure = {1, 3, 11, 5};
    steps::dist::EntityMolecules<steps::dist::mesh::tetrahedron_id_t, osh::LO> en_mol(structure);
    osh::Real val;

    steps::dist::mesh::tetrahedron_id_t elem_rd(3);
    osh::LO species_rd(2);
    en_mol.track_occupancy_rd(elem_rd, species_rd);
    steps::dist::mesh::tetrahedron_id_t elem_ef(2);
    osh::LO species_ef(1);
    en_mol.track_occupancy_ef(elem_ef, species_ef);

    // start from 0
    EXPECT_EQ(en_mol(elem_rd, species_rd), 0);

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
    EXPECT_DOUBLE_EQ(val, 1);
    // get integral after 2 more seconds
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    EXPECT_DOUBLE_EQ(val, 0.5);

    // reset ef and check that the integrals are at 0
    en_mol.reset_occupancy_ef(t);
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    EXPECT_DOUBLE_EQ(val, 2);

    // check full reset
    en_mol.add_and_update_occupancy(elem_ef, species_ef, 5, t);
    en_mol.reset(t);
    val = en_mol(elem_rd, species_rd);
    EXPECT_DOUBLE_EQ(val, 0);
    val = en_mol.get_occupancy_rd(elem_rd, species_rd, t);
    EXPECT_DOUBLE_EQ(val, 0);
    val = en_mol(elem_ef, species_ef);
    EXPECT_DOUBLE_EQ(val, 0);
    val = en_mol.get_occupancy_ef(elem_ef, species_ef, t);
    EXPECT_DOUBLE_EQ(val, 0);
}


int main(int argc, char* argv[]) {
    int r = 0;
    ::testing::InitGoogleTest(&argc, argv);
    r = RUN_ALL_TESTS();
    return r;
}

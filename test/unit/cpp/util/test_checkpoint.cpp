#include <fstream>
#include <vector>

#include "solver/fwd.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
#include "util/init.hpp"
#include "util/vocabulary.hpp"

#include "gtest/gtest.h"

class CheckpointTestEnvironment: public ::testing::Environment {
  public:
    virtual void SetUp() {
        steps::init();
    }
};

TEST(CheckpointTest, ScalarRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    int cp_a{10};
    double cp_b{0.0001};
    steps::vertex_id_t cp_v{20};

    steps::util::checkpoint(cp_file, cp_a);
    steps::util::checkpoint(cp_file, cp_b);
    steps::util::checkpoint(cp_file, cp_v);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    int rs_a;
    double rs_b;
    steps::vertex_id_t rs_v;
    steps::util::restore(cp_file, rs_a);
    ASSERT_EQ(cp_a, rs_a);
    steps::util::restore(cp_file, rs_b);
    ASSERT_EQ(cp_b, rs_b);
    steps::util::restore(cp_file, rs_v);
    ASSERT_EQ(cp_v, rs_v);
    cp_file.close();
}

TEST(CheckpointTest, ScalarCompare) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    int cp_a{10};
    double cp_b{0.0001};
    steps::vertex_id_t cp_v{20};

    steps::util::checkpoint(cp_file, cp_a);
    steps::util::checkpoint(cp_file, cp_b);
    steps::util::checkpoint(cp_file, cp_v);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    int rs_a{20};
    double rs_b{1.0};
    steps::vertex_id_t rs_v{4};
    EXPECT_THROW(steps::util::compare(cp_file, rs_a), steps::CheckpointErr);
    EXPECT_THROW(steps::util::compare(cp_file, rs_b), steps::CheckpointErr);
    EXPECT_THROW(steps::util::compare(cp_file, rs_v), steps::CheckpointErr);
    cp_file.close();
}

TEST(CheckpointTest, VectorRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::vector<double> double_vec = {0.1, 0.5, 0.4, 1.4};
    std::vector<steps::triangle_id_t> strong_id_vec{steps::triangle_id_t(1),
                                                    steps::triangle_id_t(5),
                                                    steps::triangle_id_t(7),
                                                    steps::triangle_id_t(8),
                                                    steps::triangle_id_t(100)};

    steps::util::checkpoint(cp_file, double_vec);
    steps::util::checkpoint(cp_file, strong_id_vec);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::vector<double> restore_double_vec;
    std::vector<steps::triangle_id_t> restore_strong_id_vec;

    steps::util::restore(cp_file, restore_double_vec);

    ASSERT_EQ(double_vec.size(), restore_double_vec.size());
    for (size_t i = 0; i < double_vec.size(); ++i) {
        EXPECT_EQ(double_vec[i], restore_double_vec[i]);
    }

    steps::util::restore(cp_file, restore_strong_id_vec);
    ASSERT_EQ(strong_id_vec.size(), restore_strong_id_vec.size());
    for (size_t i = 0; i < restore_strong_id_vec.size(); ++i) {
        EXPECT_EQ(strong_id_vec[i], restore_strong_id_vec[i]);
    }
    cp_file.close();
}

TEST(CheckpointTest, MapRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::map<steps::tetrahedron_id_t, long> cp_map = {{steps::tetrahedron_id_t(2), 0.44},
                                                      {steps::tetrahedron_id_t(5), 2.43},
                                                      {steps::tetrahedron_id_t(3), 6.9}};

    steps::util::checkpoint(cp_file, cp_map);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::map<steps::tetrahedron_id_t, long> restore_map;

    steps::util::restore(cp_file, restore_map);

    ASSERT_EQ(cp_map.size(), restore_map.size());
    for (auto& e: cp_map) {
        auto result = restore_map.find(e.first);
        EXPECT_NE(result, restore_map.end());
        EXPECT_EQ(result->second, e.second);
    }
    cp_file.close();
}

TEST(CheckpointTest, CRDataRestore) {
    struct CRKProcData {
        bool recorded{false};
        int pow{0};
        unsigned pos{0};
        double rate{0.0};
    };

    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    CRKProcData data{true, 10, 5, 0.441};

    steps::util::checkpoint(cp_file, data);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    CRKProcData restored_data;

    steps::util::restore(cp_file, restored_data);

    ASSERT_EQ(data.recorded, restored_data.recorded);
    ASSERT_EQ(data.pow, restored_data.pow);
    ASSERT_EQ(data.pos, restored_data.pos);
    ASSERT_EQ(data.rate, restored_data.rate);
    cp_file.close();
}

TEST(CheckpointTest, CRDataVecRestore) {
    struct CRKProcData {
        bool recorded{false};
        int pow{0};
        unsigned pos{0};
        double rate{0.0};
    };

    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::vector<CRKProcData> data_vec;
    data_vec.push_back({true, 4, 3, 0.241});
    data_vec.push_back({false, 2, 7, 1.3441});
    data_vec.push_back({true, 977, 45345, 987.441});

    steps::util::checkpoint(cp_file, data_vec);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::vector<CRKProcData> rs_data_vec;

    steps::util::restore(cp_file, rs_data_vec);

    for (size_t e = 0; e < 3; e++) {
        ASSERT_EQ(data_vec[e].recorded, rs_data_vec[e].recorded);
        ASSERT_EQ(data_vec[e].pow, rs_data_vec[e].pow);
        ASSERT_EQ(data_vec[e].pos, rs_data_vec[e].pos);
        ASSERT_EQ(data_vec[e].rate, rs_data_vec[e].rate);
    }
    cp_file.close();
}

TEST(CheckpointTest, CArrayRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    int c_array[4] = {0, 1, 2, 3};

    steps::util::checkpoint(cp_file, c_array, 4);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    int restored_array[4];

    steps::util::restore(cp_file, restored_array, 4);

    for (size_t i = 0; i < 4; ++i) {
        ASSERT_EQ(c_array[i], restored_array[i]);
    }
    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);
    int wrong_size_array[6];
    EXPECT_THROW(steps::util::restore(cp_file, wrong_size_array, 6), steps::CheckpointErr);
    cp_file.close();
}

TEST(CheckpointTest, StrongVecRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    steps::util::strongid_vector<steps::solver::spec_local_id, uint> counts;
    counts.container().resize(4);
    counts.at(steps::solver::spec_local_id(0)) = 4;
    counts.at(steps::solver::spec_local_id(1)) = 45;
    counts.at(steps::solver::spec_local_id(2)) = 2;
    counts.at(steps::solver::spec_local_id(3)) = 999;
    steps::util::checkpoint(cp_file, counts);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    steps::util::strongid_vector<steps::solver::spec_local_id, uint> restored_counts;

    steps::util::restore(cp_file, restored_counts);

    for (uint i = 0; i < 4; ++i) {
        ASSERT_EQ(counts.at(steps::solver::spec_local_id(i)),
                  restored_counts.at(steps::solver::spec_local_id(i)));
    }
    cp_file.close();
}

TEST(CheckpointTest, StrongMapOfStrongVecsRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::map<steps::solver::vesicle_global_id,
             steps::util::strongid_vector<steps::solver::spec_global_id, double>>
        data;

    for (uint i = 0; i < 4; ++i) {
        data[steps::solver::vesicle_global_id(i)].container().resize(4);
        data[steps::solver::vesicle_global_id(i)].at(steps::solver::spec_global_id(0)) = 1.0 * i +
                                                                                         0.1;
        data[steps::solver::vesicle_global_id(i)].at(steps::solver::spec_global_id(1)) = 1.0 * i +
                                                                                         0.2;
        data[steps::solver::vesicle_global_id(i)].at(steps::solver::spec_global_id(2)) = 1.0 * i +
                                                                                         0.3;
        data[steps::solver::vesicle_global_id(i)].at(steps::solver::spec_global_id(3)) = 1.0 * i +
                                                                                         0.4;
    }

    steps::util::checkpoint(cp_file, data);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::map<steps::solver::vesicle_global_id,
             steps::util::strongid_vector<steps::solver::spec_global_id, double>>
        restored_data;

    steps::util::restore(cp_file, restored_data);

    ASSERT_EQ(data.size(), restored_data.size());

    for (auto const& e: data) {
        auto key = e.first;
        auto restore_result = restored_data.find(key);
        ASSERT_NE(restore_result, restored_data.end());
        ASSERT_EQ(e.second.size(), restore_result->second.size());
        for (uint i = 0; i < 4; ++i) {
            ASSERT_EQ(e.second.at(steps::solver::spec_global_id(i)),
                      restore_result->second.at(steps::solver::spec_global_id(i)));
        }
    }
    cp_file.close();
}

TEST(CheckpointTest, STLArrayRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::array<double, 4> stl_array{0.343, 232.232, 245.22, 5343.33};

    steps::util::checkpoint(cp_file, stl_array);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::array<double, 4> restored_array;

    steps::util::restore(cp_file, restored_array);

    for (size_t i = 0; i < 4; ++i) {
        ASSERT_EQ(stl_array[i], restored_array[i]);
    }
    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);
    std::array<double, 10> wrong_size_array;
    EXPECT_THROW(steps::util::restore(cp_file, wrong_size_array), steps::CheckpointErr);
    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);
    std::array<int, 4> wrong_type_array;
    EXPECT_THROW(steps::util::restore(cp_file, wrong_type_array), steps::CheckpointErr);
    cp_file.close();
}

TEST(CheckpointTest, MapOfVecsRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::map<uint, std::vector<double>> map_vecs;
    map_vecs[1] = {4.43, 0.444, 23.343};
    map_vecs[99] = {5., 343.};
    map_vecs[12234] = {564.};

    steps::util::checkpoint(cp_file, map_vecs);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::map<uint, std::vector<double>> restore_map_vecs;

    steps::util::restore(cp_file, restore_map_vecs);

    ASSERT_EQ(map_vecs.size(), restore_map_vecs.size());
    for (auto const& e: map_vecs) {
        auto key = e.first;
        auto original_vec = e.second;
        auto restore_result = restore_map_vecs.find(key);
        ASSERT_NE(restore_result, restore_map_vecs.end());
        auto restore_vec = restore_result->second;
        ASSERT_EQ(original_vec.size(), restore_vec.size());
        for (size_t i = 0; i < original_vec.size(); i++) {
            ASSERT_EQ(original_vec[i], restore_vec[i]);
        }
    }
    cp_file.close();
}

TEST(CheckpointTest, SetRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::set<double> cp_set = {10.33, 2.886, 0.33, 42424., 76565.};

    steps::util::checkpoint(cp_file, cp_set);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::set<double> restore_set;

    steps::util::restore(cp_file, restore_set);

    ASSERT_EQ(cp_set.size(), restore_set.size());

    for (auto const& e: cp_set) {
        ASSERT_NE(restore_set.find(e), restore_set.end());
    }
    cp_file.close();
}

TEST(CheckpointTest, MapOfSetsRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::map<double, std::set<double>> map_sets;

    map_sets[0.4] = {10.33, 2.886, 0.33, 42424., 76565.};
    map_sets[53.4554] = {4534.67, 56.65, 7.554};
    map_sets[23] = {344.3};

    steps::util::checkpoint(cp_file, map_sets);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::map<double, std::set<double>> restore_map_sets;

    steps::util::restore(cp_file, restore_map_sets);

    ASSERT_EQ(map_sets.size(), restore_map_sets.size());

    for (auto const& e: map_sets) {
        auto key = e.first;
        auto original_set = e.second;
        auto restore_result = restore_map_sets.find(key);
        ASSERT_NE(restore_result, restore_map_sets.end());
        auto restore_set = restore_result->second;
        ASSERT_EQ(original_set.size(), restore_set.size());
        for (auto const& set_elem: original_set) {
            ASSERT_NE(restore_set.find(set_elem), restore_set.end());
        }
    }

    cp_file.close();
}

TEST(CheckpointTest, MapOfMapsRestore) {
    std::fstream cp_file;
    cp_file.open("serial_cp_test.bin",
                 std::fstream::out | std::fstream::binary | std::fstream::trunc);

    std::map<double, std::map<double, double>> map_maps;

    map_maps[0.4] = {{1.3, 23232.2}, {10.33, 2.886}, {42424., 76565.}};
    map_maps[53.4554] = {{4534.67, 56.65}, {7.554, 434}};
    map_maps[23] = {{344.3, 75676.43}};

    steps::util::checkpoint(cp_file, map_maps);

    cp_file.close();

    cp_file.open("serial_cp_test.bin", std::fstream::in | std::fstream::binary);
    cp_file.seekg(0);

    std::map<double, std::map<double, double>> restore_map_maps;

    steps::util::restore(cp_file, restore_map_maps);

    ASSERT_EQ(map_maps.size(), restore_map_maps.size());

    for (auto const& sub_map: map_maps) {
        auto key = sub_map.first;
        auto original_map = sub_map.second;
        auto restore_result = restore_map_maps.find(key);
        ASSERT_NE(restore_result, restore_map_maps.end());
        auto restore_map = restore_result->second;
        ASSERT_EQ(original_map.size(), restore_map.size());
        for (auto const& elem: original_map) {
            auto restore_elem_result = restore_map.find(elem.first);
            ASSERT_NE(restore_elem_result, restore_map.end());
            ASSERT_EQ(elem.second, restore_elem_result->second);
        }
    }

    cp_file.close();
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new CheckpointTestEnvironment);
    return RUN_ALL_TESTS();
}

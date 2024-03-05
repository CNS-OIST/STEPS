#include "mpi/dist/tetopsplit/kproc/propensities.hpp"

#include "boost/iostreams/device/null.hpp"
#include "boost/iostreams/stream.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("propensities", "[steps4]") {
    boost::iostreams::stream<boost::iostreams::null_sink> null_ostr(
        (boost::iostreams::null_sink()));
    steps::dist::kproc::Propensities propensities;

    for (const auto& group: propensities.groups()) {
        null_ostr << group << '\n';
    }
}

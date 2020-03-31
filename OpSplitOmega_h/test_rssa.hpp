#pragma once

#include "common.hpp"

namespace zee {

template <osh::Int Dim>
int rssa_opsplit_test(OmegaHMesh<Dim>& mesh, const osh::LO num_simulation_paths);

// explicit template instantiation declarations
extern template int rssa_opsplit_test(OmegaHMesh<2>& mesh, const osh::LO num_simulation_paths);
extern template int rssa_opsplit_test(OmegaHMesh<3>& mesh, const osh::LO num_simulation_paths);

}  // namespace zee

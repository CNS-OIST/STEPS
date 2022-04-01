#pragma once

#include "mpi/dist/tetopsplit/fwd.hpp"

namespace steps {
namespace dist {

class EFieldOperator;

template <typename RNG, typename NumMolecules> class DiffusionOperator;

template <typename RNG, typename NumMolecules, NextEventSearchMethod SearchMethod>
class SSAOperator;

template <typename RNG, typename NumMolecules>
class RSSAOperator;

} // namespace dist
} // namespace steps

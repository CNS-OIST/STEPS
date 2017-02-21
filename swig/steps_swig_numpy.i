/*
 ___license_placeholder___
 */


#define WITH_NUMPY
%{
#define WITH_NUMPY
%}

%module steps_swig_numpy

%{
// SWIG bug workaround: http://sourceforge.net/p/swig/bugs/1187/
#include <cstddef>
%} 

%include "model.i"
%include "geom.i"
%include "rng.i"
%include "numpy.i"
%include "solver.i"
#ifdef WITH_MPI
%include "mpi_solver.i"
#endif

%{
#include "steps/init.hpp"
#include "steps/finish.hpp"
%}

#ifdef WITH_MPI
%{
#include "steps/mpi/mpi_init.hpp"
#include "steps/mpi/mpi_finish.hpp"
%}
#endif
////////////////////////////////////////////////////////////////////////////////

namespace steps
{
    void init(void);
    void finish(void);
    
} // end steps namesapce

#ifdef WITH_MPI
namespace steps
{
namespace mpi
{
    void mpiInit(void);
    int getRank(void);
    int getNHosts(void);
    void mpiFinish(void);
}
} // end steps namesapce
#endif

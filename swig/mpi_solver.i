/*
 ___license_placeholder___
 */


%module mpi_solver_swig

%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

%template(map_unsigned_unsigned) std::map<unsigned int, unsigned int>;
%template(vector_unsigned) std::vector<unsigned int>;

#ifdef WITH_NUMPY
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"

%init %{
    import_array();
%}

%apply (unsigned int* IN_ARRAY1, int DIM1) {
    (unsigned int* indices, int input_size)
}
%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* output, int output_size)
}

%import "unchecked_stl_seq.i"
UNCHECKED_STL_SEQ_CONVERT(std::vector<unsigned int>,push_back,PyInt_AsUnsignedLongMask)
UNCHECKED_STL_DICT_CONVERT(%arg(std::map<unsigned int,unsigned int>),insert,PyInt_AsUnsignedLongMask,PyInt_AsUnsignedLongMask)

#endif

////////////////////////////////////////////////////////////////////////////////

%include "error.i"
%import "steps/common.h"
%{
#include "steps/solver/api.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
%}

%feature("autodoc", "1");

%exception 
{
    try {
        $action
    } catch (steps::ArgErr & ae) {
        PyErr_SetString(PyExc_NameError, ae.getMsg());
        return NULL;
    } catch (steps::NotImplErr & nie) {
        PyErr_SetString(PyExc_NotImplementedError, nie.getMsg());
        return NULL;
    } catch (steps::ProgErr & pe){
        PyErr_SetString(PyExc_RuntimeError, pe.getMsg());
        return NULL;
    }
}

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace mpi
{
namespace tetopsplit
{
class TetOpSplitP : public steps::solver::API
{    

public:
    %feature("autodoc", "1");
/* note: interface with default argyments is not used with the python wrapper, and SWIG
   (as of version 3.0.7) does not correctly generate the typechecks for overload resolution
   with typecheck typemaps on constructors, which prevents the workaround for numpy 1.4.1's
   improper integral type registration from being triggered.

    TetOpSplitP(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r, int calcMembPot = EF_NONE, std::vector<uint> const & tet_hosts = std::vector<uint>(), std::map<uint, uint> const & tri_hosts = std::map<uint, uint>(), std::vector<uint> const & wm_hosts = std::vector<uint>());
*/
    TetOpSplitP(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r, int calcMembPot, std::vector<uint> const & tet_hosts, std::map<uint, uint> const & tri_hosts, std::vector<uint> const & wm_hosts);
    %feature("autodoc", "1");
    ~TetOpSplitP(void);

/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Returns a string of the solver's name.

Syntax::
    
    getSolverName()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverName(void) const;
    
/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Returns a string giving a short description of the solver.

Syntax::
    
    getSolverDesc()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverDesc(void) const;
    
/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Returns a string of the solver authors names.

Syntax::
    
    getSolverAuthors()
    
Arguments:
    None

Return:
    string
");
    virtual std::string getSolverAuthors(void) const;
    
/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Returns a string giving the author's email address.

Syntax::
    
    getSolverEmail()
    
Arguments:
    None

Return:
    string
");

    virtual std::string getSolverEmail(void) const;
    
/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Checkpoint data to a file. Not yet implemented.
    
Syntax::
    
    checkpoint(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");

    virtual void checkpoint(std::string const & file_name);
    
/////////////------------------------------------------------------////////////////
    
    %feature("autodoc", 
"
Restore data from a file. Not yet implemented.
    
Syntax::
    
    restore(file_name)
    
Arguments:
    string file_name
    
Return:
    None
");
    virtual void restore(std::string const & file_name);
    
/////////////------------------------------------------------------////////////////
    
    %feature("autodoc", 
"
Reset the simulation to the state the solver was initialised to.

Syntax::
    
    reset()
    
Arguments:
    None

Return:
    None
");
    virtual void reset(void);
    
/////////////------------------------------------------------------////////////////

    %feature("autodoc", 
"
Advance the simulation until endtime (given in seconds) is reached. 
The endtime must be larger or equal to the current simulation time.

Syntax::
    
    run(endtime)
    
Arguments:
    float endtime

Return:
    None
");
    virtual void run(double endtime);
    
/////////////------------------------------------------------------////////////////
    
    %feature("autodoc", 
"
Returns the current simulation time in seconds.

Syntax::
    
    getTime()
    
Arguments:
    None

Return:
    float
");
    virtual double getTime(void) const;
    
/////////////------------------------------------------------------////////////////
    
    %feature("autodoc", 
"
Return the update period tau of the Operator-Splitting solution.
See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.

Syntax::
    
    getUpdPeriod()
    
Arguments:
    None

Return:
    float
");
    double getUpdPeriod(void);
    
/////////////------------------------------------------------------////////////////
    
    %feature("autodoc", 
"
Set the threshold for using binomial distribution for molecule diffusion instead of
single molecule diffusion.

If the number of molecules in a tetrahedron await for diffusion is higher than this
threshold, the solver will use binomial function to distribute these molecules to
each neighboring tetrahedron. Otherwise the molecules will diffuse one by one.

The default threshold is 10.

Syntax::
    
    setDiffApplyThreshold(threshold)
    
Arguments:
    int threshold

Return:
    None
");
    void setDiffApplyThreshold(int threshold);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the number of reaction events happened in the simulation.

if all processes call this function, it will return the accumulated
result accross all processes. It can also be called in individual process with
the local argument set to true, in which case it returns the local result of this process.

By default it is called globally and return the accumlated result.

Syntax::
    
    getReacExtent(local)
    
Arguments:
    bool local (default: false)

Return:
    float
");
    double getReacExtent(bool local = false);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the number of diffusion events happened in the simulation.

if all processes call this function, it will return the accumulated
result accross all processes. It can also be called in individual process with
the local argument set to true, in which case it returns the local result of this process.

By default it is called globally and return the accumlated result.

Syntax::
    
    getDiffExtent(local)
    
Arguments:
    bool local (default: false)

Return:
    float
");
    double getDiffExtent(bool local = false);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the number of Operator-Splitting iterations happened in the simulation.

See (Hepburn et al, 2016) and (Chen et al, 2017) for more detail.

This function can be called locally.

Syntax::
    
    getNIteration()
    
Arguments:
    None

Return:
    float
");
    double getNIteration(void);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated sum of species s in a batch of tetrahedrons.

This function requires NumPy array as input, and called globally in all processes.

Syntax::
    
    sumBatchTetCountsNP(tet_list, s)
    
Arguments:
    numpy.array tet_list
    string s

Return:
    float
");
    double sumBatchTetCountsNP(unsigned int* indices, int input_size, std::string const & s);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated sum of species s in a batch of triangles.

This function requires NumPy array as input, and called globally in all processes.

Syntax::
    
    sumBatchTriCountsNP(tri_list, s)
    
Arguments:
    numpy.array tri_list
    string s

Return:
    float
");
    double sumBatchTriCountsNP(unsigned int* indices, int input_size, std::string const & s);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated sum of GHK currents in a batch of triangles.

This function requires NumPy array as input, and called globally in all processes.

Syntax::
    
    sumBatchTriGHKIsNP(tri_list, ghk)
    
Arguments:
    numpy.array tri_list
    string ghk

Return:
    float
");
    double sumBatchTriGHKIsNP(unsigned int* indices, int input_size, std::string const & ghk);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated sum of Ohmic currents in a batch of triangles.

This function requires NumPy array as input, and called globally in all processes.

Syntax::
    
    sumBatchTriOhmicIsNP(tri_list, oc)
    
Arguments:
    numpy.array tri_list
    string oc

Return:
    float
");
    double sumBatchTriOhmicIsNP(unsigned int* indices, int input_size, std::string const & oc);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Repartition and reset the simulation.

Note: It is possible to repartitioning the mesh using any subset of the processes,
in which case processes with no assigned subvolumes will mostly idle for the working processes.

Syntax::
    
    repartitionAndReset(tet_hosts, tri_hosts, wm_hosts)
    
Arguments:
    list tet_hosts
    dict tri_hosts (default: {})
    dict wm_hosts (default: {})

Return:
    None
");
    void repartitionAndReset(std::vector<uint> const &tet_hosts, std::map<uint, uint> const &tri_hosts  = std::map<uint, uint>(), std::vector<uint> const &wm_hosts = std::vector<uint>());

/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated computation time of the process.

To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
line in the src/CmakeLists.txt file. This function is always called and return result locally.

See (Chen et al, 2017) for more detail.

Syntax::
    
    getCompTime()
    
Arguments:
    None
    
Return:
    float
");
    double getCompTime(void);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated synchronization time of the process.

To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
line in the src/CmakeLists.txt file. This function is always called and return result locally.

See (Chen et al, 2017) for more detail.

Syntax::
    
    getSyncTime()
    
Arguments:
    None
    
Return:
    float
");
    double getSyncTime(void);
    
/////////////------------------------------------------------------////////////////
    %feature("autodoc", 
"
Return the accumulated idle time of the process.

To use the funtion, it is necessary to enable the add_definitions(-DMPI_PROFILING=1)
line in the src/CmakeLists.txt file. This function is always called and return result locally.

See (Chen et al, 2017) for more detail.

Syntax::
    
    getSyncTime()
    
Arguments:
    None
    
Return:
    float
");
    double getIdleTime(void);
    
/////////////------------------------------------------------------////////////////
};
    
////////////////////////////////////////////////////////////////////////////////
} // end namespace tetopsplit
} // end namespace mpi
} // end namespace steps


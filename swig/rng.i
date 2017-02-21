/*
 ___license_placeholder___
 */


%module rng

%import "steps/common.h"



%include "python/std_string.i"
%include "error.i"
%{
#include "steps/rng/rng.hpp"
#include "steps/rng/create.hpp"
%}


///////////////////////////////////////////////////////////////////////////////    

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

///////////////////////////////////////////////////////////////////////////////    

namespace steps
{
    namespace rng
    {
        
        ////////////////////////////////////////////////////////////////////////////////
        
        %feature("autodoc", 
"
Proxy of C++ random number generator.

");    
        class RNG
        {
            
        public:
            
            RNG(unsigned int bufsize);
            virtual ~RNG(void);
            
            %feature("autodoc", 
"
Initialize the random number generator with given seed value.


Syntax::
    
    initialize(seed)
    
Arguments:
    uint seed

Return:
    None
");
            void initialize(unsigned long const & seed);
            
            unsigned int get(void);
            
            double getUnfII(void);
            double getUnfIE(void);
            double getUnfEE(void);
            double getUnfIE53(void);
            
            float getStdExp(void);
            double getExp(double lambda);
            long getPsn(float lambda);
            float getStdNrm(void);
            
        protected:
            // Mark this class as abstract.
            virtual void concreteInitialize(ulong seed) = 0;
            virtual void concreteFillBuffer(void) = 0;
            
        };
        
        ////////////////////////////////////////////////////////////////////////////////
        %feature("autodoc", 
"
Equivalent to: create('mt19937', buffer_size)

Syntax::

    create_mt19937(buffer_size)
    
Arguments:
    uint buffer_size

Return:
    steps.rng.RNG
");        
        RNG * create_mt19937(unsigned int bufsize);
        
        %feature("autodoc", 
"
Creates and returns a reference to a steps.rng.RNG random number generator object, 
which is specified by type and pre-allocates a buffer list with size of buffer_size.

Syntax::
    
    create(type, buffer_size)

Arguments:
    * string type
    * uint buffer_size

Return:
    steps.rng.RNG
");    
        RNG * create(std::string rng_name, unsigned int bufsize);
        
        ////////////////////////////////////////////////////////////////////////////////
        
    } // end namespace rng
} // end namespace steps

// END


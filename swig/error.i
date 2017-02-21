/*
 ___license_placeholder___
 */


%module error_swig

%{
// Autotools definitions.
#include "steps/error.hpp"
%}

////////////////////////////////////////////////////////////////////////////////

%feature("autodoc", "1");

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
	
struct Err
{
	Err(std::string const & msg = "");
	const char * getMsg(void);
	
};

struct ArgErr
: public Err
{
	ArgErr(std::string const & msg = "")
	: Err(msg) { }
};


struct NotImplErr
: public Err
{
    NotImplErr(std::string const & msg = "")
    : Err(msg) { }
};
    
struct ProgErr
: public Err
{
    ProgErr(std::string const & msg = "")
    : Err(msg) { }
};

} // end steps namesapce

////////////////////////////////////////////////////////////////////////////////


// END

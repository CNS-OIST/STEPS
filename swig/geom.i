/*
 ___license_placeholder___
 */


%module geom_swig

%include "std_map.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_set.i"

%template(vector_unsigned) std::vector<unsigned int>;
%template(set_unsigned) std::set<unsigned int>;

#ifdef WITH_NUMPY
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (unsigned int* IN_ARRAY1, int DIM1) {
    (unsigned int* indices, int input_size),
    (unsigned int* t_indices, int input_size),
    (unsigned int* indices, int index_size)
}

%apply (double* INPLACE_ARRAY1, int DIM1) {
    (double* centres, int output_size),
    (double* cords, int cord_size),
    (double* coordinates, int output_size),
    (double* volumes, int volume_size),
    (double* areas, int area_size)
}

%apply (unsigned int* INPLACE_ARRAY1, int DIM1) {
    (unsigned int* v_indices, int output_size),
    (unsigned int* t_vertices, int t_vertices_size),
    (unsigned int* v_set, int v_set_size),
    (unsigned int* point_counts, int count_size)
}

%import "unchecked_stl_seq.i"
UNCHECKED_STL_SEQ_CONVERT(std::vector<unsigned int>,push_back,PyInt_AsUnsignedLongMask)
UNCHECKED_STL_SEQ_CONVERT(std::set<unsigned int>,insert,PyInt_AsUnsignedLongMask)

#endif

///////////////////////////////////////////////////////////////////////////////

%include "error.i"
%import "steps/common.h"

%{
#include "steps/geom/geom.hpp"
#include "steps/geom/comp.hpp"
#include "steps/geom/patch.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/geom/memb.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/geom/tetmesh_rw.hpp"
#include "steps/error.hpp"
#include "steps/model/model.hpp"
#include "steps/model/spec.hpp"
#include "steps/model/reac.hpp"
#include "steps/model/diff.hpp"
#include "steps/geom/diffboundary.hpp"
%}

namespace std
{
    %template(set_str)     set<std::string>;
    %template(set_ptc)     set<steps::wm::Patch *>;
    %template(vector_ptc)  vector<steps::wm::Patch *>;
    %template(vector_cmp)  vector<steps::wm::Comp *>;
    %template(vector_tmp)  vector<steps::tetmesh::TmPatch *>;
    %template(vector_db)   vector<steps::tetmesh::DiffBoundary *>;
    %template(vector_str)  vector<std::string>;
    %template(vector_int)  vector<int>;
    %template(vector_dbl)  vector<double>;
    %template(vector_bool) vector<bool>;
}

namespace steps
{
namespace wm
{

class Comp;
class Patch;
class Geom;

}
}

namespace steps
{
namespace model
{	

class Surfsys;
class Volsys;

}
}

namespace steps
{
namespace tetmesh
{
		
class Tetmesh;
class Tri;
class Tet;
class TmPatch;
class TmComp;
class Memb;
class DiffBoundary;

enum ElementType {ELEM_VERTEX, ELEM_TRI, ELEM_TET, ELEM_UNDEFINED = 99};

struct ROISet {
    ElementType                                         type;
    std::vector<uint>                                   indices;
};

}
}



////////////////////////////////////////////////
// add for mesh object casting from wm objects
// Weiliang 20130521

%inline %{
    steps::tetmesh::TmComp* castToTmComp(steps::wm::Comp* base) {
        return dynamic_cast<steps::tetmesh::TmComp*>(base);
    }
    
    steps::tetmesh::TmPatch* castToTmPatch(steps::wm::Patch* base) {
        return dynamic_cast<steps::tetmesh::TmPatch*>(base);
    }
%}


%feature("autodoc", "1");

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
namespace wm
{

////////////////////////////////////////////////////////////////////////////////
/*
bool isValidID(std::string const & id);
void checkID(std::string const & id);
*/
////////////////////////////////////////////////////////////////////////////////

class Geom
{
	
public:
	
	Geom(void);
	virtual ~Geom(void);
	
    %feature("autodoc", 
"
Returns a reference to the steps.model.Comp compartment object with 
identifier string comp_id (if defined).

Syntax::

    getComp(comp_id)
    
Arguments:
    string comp_id
             
Return:
    steps.model.Comp
");
	steps::wm::Comp * getComp(std::string const & id) const;

    %feature("autodoc", 
"
Removes the steps.geom.Comp object with identifier string comp_id (if defined) 
from the geometry container.

Syntax::

    delComp(comp_id)
    
Arguments:
    string comp_id
             
Return:
    None
");
	void delComp(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of references to all steps.geom.Comp compartment objects in the 
geometry container.

Syntax::

    getAllComps()
    
Arguments:
    None
             
Return:
    list<steps.geom.Comp>
");
	std::vector<steps::wm::Comp *> getAllComps(void) const;
	
    %feature("autodoc", 
"
Removes the steps.geom.Patch object with identifier string patch_id (if defined) 
from the geometry container.

Syntax::

    getPatch(patch_id)
    
Arguments:
    string patch_id
             
Return:
    steps.geom.Patch
");
	steps::wm::Patch * getPatch(std::string const & id) const;
    
    %feature("autodoc", 
"
Removes the steps.geom.Patch object with identifier string patch_id (if defined) 
from the geometry container.

Syntax::

    delPatch(patch_id)
    
Arguments:
    string patch_id
             
Return:
    None
");
	void delPatch(std::string const & id);
    
    
    %feature("autodoc", 
"
Returns a list of references to all steps.geom.Patch patch objects in the 
geometry container.

Syntax::

    getAllPatches()
    
Arguments:
    None
             
Return:
    list<steps.geom.Patch>
");
	std::vector<steps::wm::Patch *> getAllPatches(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") Patch::Patch;

class Patch
{
	
public:
	
	Patch(std::string const & id, steps::wm::Geom * container, 
		  steps::wm::Comp * icomp, steps::wm::Comp * ocomp = 0, double area = 0.0);
	virtual ~Patch(void);
	
    %feature("autodoc", 
"
Get the identifier string of the patch.

Syntax::

    getID()
    
Arguments:
    None
             
Return:
    string
");
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the patch.

Syntax::

    setID(name)
    
Arguments:
    string name
             
Return:
    None
");
	void setID(std::string const & id);
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.geom.Geom container object.

Syntax::

    getContainer()
    
Arguments:
    None
             
Return:
    steps.geom.Geom
");
	steps::wm::Geom * getContainer(void) const;
    
    %feature("autodoc", 
"
Get the area of the patch (in m^2).

Syntax::

    getArea()

Arguments:
    None
             
Return:
    float
");
	double getArea(void) const;
    
    %feature("autodoc", 
"
Set the area of the patch (in m^2).

Syntax::

    setArea(area)
    
Arguments:
    float area
             
Return:
    None
");
	virtual void setArea(double vol);
	
    %feature("autodoc", 
"
Add surface system identifier string surfsys_id to the patch object.

Syntax::

    addSurfsys(surfsys_id)
    
Arguments:
    string surfsys_id
             
Return:
    None
");
	void addSurfsys(std::string const & id);

    %feature("autodoc", 
"
Returns a list of the surface system identifier strings which have 
been added to the patch.

Syntax::

    getSurfsys()
    
Arguments:
    None
             
Return:
    list<string>
");
	std::set<std::string> getSurfsys(void) const;
    
    %feature("autodoc", 
"
Removes surface system identifier string surfsys_id from this patch.

Syntax::

    delSurfsys(surfsys_id)

Arguments:
    string surfsys_id
             
Return:
    None
");
	void delSurfsys(std::string const & id);
	
    %feature("autodoc", 
"
Giving a steps.model.Model, return all species in the compartment.

Syntax::

    getAllSpecs(model)

Arguments:
    steps.model.Model model

Return:
    list<steps.model.Spec>
");
    
    std::vector<steps::model::Spec*> getAllSpecs(steps::model::Model* model);
    
    %feature("autodoc", 
"
Giving a steps.model.Model, return all surface reactions in the compartment.

Syntax::

    getAllSReacs(model)

Arguments:
    steps.model.Model model

Return:
    list<steps.model.SReac>
");
    std::vector<steps::model::SReac*> getAllSReacs(steps::model::Model* model);
    

    %feature("autodoc", 
"
Returns a reference to the steps.geom.Comp compartment object representing
the inner compartment.

Syntax::

    getIComp()

Arguments:
    None
             
Return:
    steps.geom.Comp
");
	steps::wm::Comp * getIComp(void) const;
    
    %feature("autodoc", 
"
Returns a reference to the steps.geom.Comp compartment object representing
the outer compartment.

Syntax::
    
    getOComp()

Arguments:
    None
             
Return:
    steps.geom.Comp
");
	steps::wm::Comp * getOComp(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") Comp::Comp;

class Comp
{
	
public:
	
	Comp(std::string const & id, steps::wm::Geom * container,
		 double vol = 0.0);
	virtual ~Comp(void);
	
    %feature("autodoc", 
"
Get the identifier string of the compartment.

Syntax::

    getID()

Arguments:
    None
             
Return:
    string
");
	std::string getID(void) const;
    
    %feature("autodoc", 
"
Set the identifier string of the compartment.

Syntax::

    setID(name)

Arguments:
    string name
             
Return:
    None
");
	void setID(std::string const & id);
	
    %feature("autodoc", 
"
Returns a reference to the parent steps.geom.Geom container object.

Syntax::

    getContainer()

Arguments:
    None
             
Return:
    steps.geom.Geom
");
	steps::wm::Geom * getContainer(void) const;
	
    %feature("autodoc", 
"
Get the volume of the compartment (in m^3).

Syntax::

    getVol()
    
Arguments:
    None
             
Return:
    float
");
	double getVol(void) const;
    
    %feature("autodoc", 
"
Set the volume of the compartment (in m^3).

Syntax::

    setVol(vol)

Arguments:
    float vol
             
Return:
    None
");

	virtual void setVol(double vol);
	
    %feature("autodoc", 
"
Add volume system identifier string volsys_id to the compartment object.

Syntax::

    addVolsys(volsys_id)

Arguments:
    string volsys_id
             
Return:
    None
");
	void addVolsys(std::string const & id);
    
    %feature("autodoc", 
"
Returns a list of the volume system identifier strings which have been 
added to the compartment.

Syntax::

    getVolsys()

Arguments:
    None
             
Return:
    list<string>
");
	std::set<std::string> getVolsys(void) const;
    
    %feature("autodoc", 
"
Removes volume system identifier string volsys_id from this compartment.

Syntax::

    delVolsys(volsys_id)

Arguments:
    string volsys_id
             
Return:
    None
");
	void delVolsys(std::string const & id);
    
    %feature("autodoc",	
"
Giving a steps.model.Model, return all species in the compartment.
    
Syntax::
    
    getAllSpecs(model)
    
Arguments:
    steps.model.Model model
    
Return:
    list<steps.model.Spec>
");
    
    std::vector<steps::model::Spec*> getAllSpecs(steps::model::Model* model);
    
    %feature("autodoc", 
"
Giving a steps.model.Model, return all reactions in the compartment.

Syntax::

    getAllReacs(model)

Arguments:
    steps.model.Model model

Return:
    list<steps.model.Reac>
");
    std::vector<steps::model::Reac*> getAllReacs(steps::model::Model* model);
    
    %feature("autodoc", 
"
Giving a steps.model.Model, return all diffusions in the compartment.

Syntax::

    getAllDiffs(model)

Arguments:
    steps.model.Model model

Return:
    list<steps.model.Diff>
");
    std::vector<steps::model::Diff*> getAllDiffs(steps::model::Model* model);

    
    %feature("autodoc", 
"
Returns a list of references to steps.geom.Patch patch objects: 
the 'inner' patches.

Syntax::

    getIPatches()

Arguments:
    None
             
Return:
    list<steps.geom.Patch>
");
	std::set<steps::wm::Patch *> getIPatches(void) const;
    
    %feature("autodoc", 
"
Returns a list of references to steps.geom.Patch patch objects: 
the 'outer' patches.

Syntax::

    getOPatches()

Arguments:
    None
             
Return:
    list<steps.geom.Patch>
");
	std::set<steps::wm::Patch *> getOPatches(void) const;
    
};

////////////////////////////////////////////////////////////////////////////////

} // namespace wm
} // namespace steps

////////////////////////////////////////////////////////////////////////////////

namespace steps
{
namespace tetmesh
{	

////////////////////////////////////////////////////////////////////////////////

class Tetmesh : public steps::wm::Geom
{

public:
	
	// Disabling this constructor because of the 'bars' issue.
	// Tetmesh(unsigned int nverts, unsigned int ntets, unsigned int ntris);
	Tetmesh(std::vector<double> const & verts,
			std::vector<unsigned int> const & tets,
			std::vector<unsigned int> const & tris = std::vector<unsigned int>());
	Tetmesh(std::vector<double> const & verts,
			std::vector<unsigned int> const & tris,
			std::vector<double> const & tri_areas,
			std::vector<double> const & tri_norms,
			std::vector<int> const & tri_tet_neighbs,
			std::vector<unsigned int> const & tets,
			std::vector<double> const & tet_vols,
			std::vector<double> const & tet_barycs,
			std::vector<unsigned int> const & tet_tri_neighbs,
			std::vector<int> const & tet_tet_neighbs);
	
	virtual ~Tetmesh(void);
    
    %feature("autodoc", 
"
Returns the coordinates of vertex with index vidx in the container.

Syntax::

    getVertex(vidx)

Arguments:
    uint vidx
             
Return:
    list<float, length = 3>
");
	std::vector<double> getVertex(unsigned int vidx) const;
    
    %feature("autodoc", 
"
Returns the total number of vertices in the mesh.

Syntax::

    countVertices()

Arguments:
    None
             
Return:
    uint
");
	unsigned int countVertices(void) const;
	
	//steps::tetmesh::TmPatch * getTmPatch(std::string const & id) const;
    %feature("autodoc", 
"Returns the vertices of bar with index bidx in the container.

Syntax::
			 
	getBar(bidx)

Arguments:
	uint bidx

Return:
	list<uint, length = 2>
");
	std::vector<unsigned int> getBar(unsigned int bidx) const;
	
	%feature("autodoc",
"
Returns the total nubmer of bars in the mesh.

Syntax::
	
	countBars()
			 
Arguments:
	None

Return:
	uint
");
    unsigned int countBars(void) const;
			 
    %feature("autodoc", 
"
Returns the triangle with index tidx in the container by its three vertex indices.

Syntax::

    getTri(tidx)

Arguments:
    uint tidx
             
Return:
    list<uint, length = 3>
");    
	std::vector<unsigned int> getTri(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the total number of triangles in the mesh.

Syntax::

    countTris()

Arguments:
    None
             
Return:
    uint
");  
	unsigned int countTris(void) const;
    
    %feature("autodoc", 
"
Returns the area of the triangle with index tidx.

Syntax::

    getTriArea(tidx)

Arguments:
    uint tidx
             
Return:
    float
");  
	double getTriArea(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the Cartesian coordinates of the barycenter of triangle with index tidx.

Syntax::

    getTriBarycenter(tidx)

Arguments:
    uint tidx
             
Return:
    list<float, length = 3>
");  
	std::vector<double> getTriBarycenter(unsigned int tidx) const;
    
	%feature("autodoc",
"
Returns the index of the bars that comprise the triangle.
	
Syntax::
			 
	getTriBars(tidx)
			 
Arguments:
	uint tidx

Return:
	list<uint, length = 3>
");	
	std::vector<unsigned int> getTriBars(unsigned int tidx) const;
	
    %feature("autodoc", 
"
Returns the radius-edge-ratio (a quality measurement) of tetrahedron with index tidx.

Syntax::
    
    getTetQualityRER(tidx)

Arguments:
    uint tidx
             
Return:
    float
");
	double getTetQualityRER(unsigned int tidx) const;

    %feature("autodoc", 
"
Returns the normal vector of the triangle with index tidx.

Syntax::

    getTriNorm(tidx)

Arguments:
    uint tidx
             
Return:
    list<float, length = 3>
");
	std::vector<double> getTriNorm(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns a reference to a step.geom.TmPatch object: the patch which triangle 
with index tidx belongs to. Returns None if triangle not assigned to a patch.

Syntax::

    getTriPatch(tidx)

Arguments:
    uint tidx
             
Return:
    steps.geom.TmPatch
");
	steps::tetmesh::TmPatch * getTriPatch(unsigned int tidx) const;

    %feature("autodoc", 
"
Returns a reference to a step.geom.Diffboundary object: the diffusion boundary triangle 
with index tidx belongs to. Returns None if triangle not assigned to a diffusion boundary.

Syntax::
             
    getTriDiffBoundary(tidx)
             
Arguments:
    uint tidx
             
Return:
    steps.geom.DiffBoundary
");
	steps::tetmesh::DiffBoundary * getTriDiffBoundary(unsigned int tidx) const;

    %feature("autodoc", 
"
Returns the indices of the two neighbouring tetrahedrons of triangle with 
index tidx. An index of -1 indicates no neighbour (triangle is on the mesh border). 

Syntax::

    getTriTetNeighb(tidx)

Arguments:
    uint tidx
             
Return:
    list<int, length = 2>
");
	std::vector<int> getTriTetNeighb(unsigned int tidx) const;
    
	%feature("autodoc", 
"
Returns the indices of the neighbouring triangles (that is all triangles that
share a 'bar') of triangle with index tidx. 

Syntax::
	
	getTriTriNeighbs(uint)
	
Arguments:
	uint tidx
	
Returns:
	list<uint>
");
	std::set<unsigned int> getTriTriNeighbs(unsigned int tidx) const;
    
	%feature("autodoc", 
"
Returns the indices of the neighbouring triangles (that is all triangles that
share a 'bar') of triangle with index tidx within the same patch. 
             
Syntax::
             
    getTriTriNeighb(uint)
             
Arguments:
    uint tidx
             
Returns:
    list<uint>
");
	std::vector<int> getTriTriNeighb(unsigned int tidx, TmPatch * tmpatch) const;
	
    // added by Weiliang
    %feature("autodoc", 
"
Returns a list of triangles that form the mesh boundary.
Support function for steps.utilities.visual.

Syntax::

    getTriBoundary()

Arguments:
    None
             
Return:
    list<int>
");
	std::vector<int> getSurfTris(void) const;
	//steps::tetmesh::TmComp * getTmComp(std::string const & id) const;
    
    %feature("autodoc", 
"
Returns the tetrahedron with index tidx in the container by its four vertex indices.

Syntax::
    getTet(tidx)
    
Arguments:
    uint tidx
             
Return:
    list<uint, length = 4>
");
	std::vector<unsigned int> getTet(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the total number of tetrahedrons in the mesh.

Syntax::

    countTets()

Arguments:
    None
             
Return:
    uint
");
	unsigned int countTets(void) const;
    
    
    %feature("autodoc", 
"
Returns the volume of the tetrahedron with index tidx.

Syntax::

    getTetVol(tidx)

Arguments:
    uint tidx
             
Return:
    float
");
	double getTetVol(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the barycenter of the tetrahedron with index tidx.

Syntax::

    getTetBarycenter(tidx)

Arguments:
    uint tidx
             
Return:
    list<float, length = 3>
");
	std::vector<double> getTetBarycenter(unsigned int tidx) const;

    %feature("autodoc", 
"
Returns a reference to a steps.geom.Comp object: the compartment which 
tetrahedron with index tidx belongs to. Returns None if tetrahedron not 
assigned to a compartment.

Syntax::

    getTetComp(tidx)

Arguments:
    uint tidx
             
Return:
    steps.geom.TmComp
");
	steps::tetmesh::TmComp * getTetComp(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the indices of the four neighbouring triangles of tetrahedron with index tidx.

Syntax::

    getTetTriNeighb(tidx)

Arguments:
    uint tidx
             
Return:
    list<uint, length = 4>
");
	std::vector<unsigned int> getTetTriNeighb(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the indices of the four neighbouring tetrahedrons of tetrahedron with index tidx. 
An index of -1 indicates no neighbour (tetrahedron is on the mesh border).

Syntax::

    getTetTetNeighb(tidx)

Arguments:
    uint tidx
             
Return:
    list<int, length = 4>
");
	std::vector<int> getTetTetNeighb(unsigned int tidx) const;
    
    %feature("autodoc", 
"
Returns the index of the tetrahedron which encompasses a given point 
p (given in Cartesian coordinates x,y,z). Returns -1 if p is a position 
outside the mesh.

Syntax::

    findTetByPoint(p)

Arguments:
    list<float, length = 3> p
             
Return:
    int
");
	int findTetByPoint(std::vector<double> p) const;
	
    %feature("autodoc", 
"
Returns the minimal Cartesian coordinate of the rectangular bounding box of the mesh. 

Syntax::

    getBoundMin()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMin(void) const;
    
    %feature("autodoc", 
"
Returns the maximal Cartesian coordinate of the rectangular bounding box of the mesh. 

Syntax::

    getBoundMax()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMax(void) const;
    
    %feature("autodoc", 
"
Returns the total volume of the mesh. 

Syntax::

    getMeshVolume()

Arguments:
    None
             
Return:
    float
");
	double getMeshVolume(void) const;
    
	
    %feature("autodoc", 
"
Returns the barycentres of a list of tetrahedrons.

Syntax::

    getBatchTetBarycentres(tets)

Arguments:
    * list<uint> tets
             
Return:
    list<float, length = len(tets) * 3>
");    
    std::vector<double> getBatchTetBarycentres(std::vector<uint> const & tets) const;
    
    %feature("autodoc", 
"
Returns the barycentres of a list of triangles.

Syntax::

    getBatchTetBarycentres(tris)

Arguments:
    * list<uint> tris
             
Return:
    list<float, length = len(tris) * 3>
");        
    std::vector<double> getBatchTriBarycentres(std::vector<uint> const & tris) const;
    
    %feature("autodoc", 
"
Returns the barycentres of a list of tetrahedrons in numpy arrays.

Syntax::

    import numpy as np
    indices = np.array([0, 1, 2], dtype= np.uint32)
    centres = np.zeros(len(indices) * 3)
    getBatchTetBarycentres(indices, centres)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = len(indices) * 3> centres
    
Return:
    None
"); 
    void getBatchTetBarycentresNP(unsigned int* indices, int input_size, double* centres, int output_size) const;
    
    %feature("autodoc", 
"
Returns the barycentres of a list of triangles in numpy arrays.

Syntax::

    import numpy as np
    indices = np.array([0, 1, 2], dtype= np.uint32)
    centres = np.zeros(len(indices) * 3)
    getBatchTriBarycentres(indices, centres)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = len(indices) * 3> centres
    
Return:
    None
"); 
    void getBatchTriBarycentresNP(unsigned int* indices, int input_size, double* centres, int output_size) const;
    
    %feature("autodoc",
"
Get coordinates of a list of vertices.

Syntax::

    getBatchVertices(verts)

Arguments:
    * list<uint> verts
             
Return:
    list<float, length = len(verts) * 3>
");

    std::vector<double> getBatchVertices(std::vector<uint> const & verts) const;
    
    %feature("autodoc", 
"
Get coordinates of a list of vertices.

Syntax::

    import numpy as np
    indices = np.array([0, 1, 2], dtype= np.uint32)
    coordinates = np.zeros(len(indices) * 3)
    getBatchVertices(indices, coordinates)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = len(indices) * 3> coordinates
    
Return:
    None
"); 
    void getBatchVerticesNP(unsigned int* indices, int input_size, double* coordinates, int output_size) const;

    %feature("autodoc",
"
Get vertex indices of a list of triangles.

Syntax::

    getBatchTris(tris)

Arguments:
    * list<uint> tris
             
Return:
    list<uint, length = len(tris) * 3>
");
    std::vector<uint> getBatchTris(std::vector<uint> const & tris) const;
    
    %feature("autodoc", 
"
Get vertex indices of a list of triangles.

Syntax::

    getBatchTrisNP(t_indices, v_indices)

Arguments:
    * numpy.array<uint> t_indices
    * numpy.array<uint, length = len(t_indices) * 3> v_indices
    
Return:
    None
"); 
    void getBatchTrisNP(unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const;
    
    %feature("autodoc",
"
Get vertex indices of a list of tetrahedrons.

Syntax::

    getBatchTets(tets)

Arguments:
    * list<uint> tets
             
Return:
    list<uint, length = len(tets) * 4>
");
    std::vector<uint> getBatchTets(std::vector<uint> const & tets) const;
    
    %feature("autodoc", 
"
Get vertex indices of a list of triangles.

Syntax::

    getBatchTetsNP(t_indices, v_indices)

Arguments:
    * numpy.array<uint> t_indices
    * numpy.array<uint, length = len(t_indices) * 4> v_indices
    
Return:
    None
"); 
    void getBatchTetsNP(unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const;

    %feature("autodoc", 
"
Return the size of a set with unique vertex indices of a list of triangles,
preparation function for furture numpy data access.

Syntax::

    getTriVerticesSetSizeNP(t_indices)

Arguments:
    * numpy.array<uint> t_indices
    
Return:
    uint
");     

    uint getTriVerticesSetSizeNP(unsigned int* t_indices, int input_size);
    
    %feature("autodoc", 
"
Return the size of a set with unique vertex indices of a list of tetrahedrons,
preparation function for furture numpy data access.

Syntax::

    getTetVerticesSetSizeNP(t_indices)

Arguments:
    * numpy.array<uint> t_indices
    
Return:
    uint
"); 

    uint getTetVerticesSetSizeNP(unsigned int* t_indices, int input_size);
    
    %feature("autodoc", 
"
Get the vertex indices of a list of triangles.
The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
whose size is provided by getTriVerticesSetSizeNP().

Syntax::

    getTriVerticesMappingSetNP(t_indices, t_vertices, v_set)

Arguments:
    * numpy.array<uint> t_indices
    * numpy.array<uint, length = length(t_indices) * 3> t_vertices
    * numpy.array<uint, length = getTriVerticesSetSizeNP(t_indices)> v_set
    
Return:
    None
"); 
    
    void getTriVerticesMappingSetNP(unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size);
    
    %feature("autodoc", 
"
Get the vertex indices of a list of tetrahedrons.
The vertex indices are reindexed, with their oringinal STEPS indices stored in a given array,
whose size is provided by getTriVerticesSetSizeNP().

Syntax::

    getTetVerticesMappingSetNP(t_indices, t_vertices, v_set)

Arguments:
    * numpy.array<uint> t_indices
    * numpy.array<uint, length = length(t_indices) * 4> t_vertices
    * numpy.array<uint, length = getTriVerticesSetSizeNP(t_indices)> v_set
    
Return:
    None
"); 
    void getTetVerticesMappingSetNP(unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size);
    
    %feature("autodoc", 
"
Generate npnts random point coordinates x,y,z within a tetraedron with index tidx, export it to NumPy array cords

Syntax::

    genPointsInTet(t_idx, npnts, coords)

Arguments:
    * unsigned tidx
    * unsigned npnts
    * numpy.array<double, length = npnts * 3> coords
    
Return:
    None
");
    void genPointsInTet(unsigned tidx, unsigned npnts, double* cords, int cord_size);
    
    %feature("autodoc", 
"
Generate npnts random point coordinates x,y,z within a triangle with index tidx, export it to NumPy array cords

Syntax::

    genPointsInTri(t_idx, npnts, coords)

Arguments:
    * unsigned tidx
    * unsigned npnts
    * numpy.array<double, length = npnts * 3> coords
    
Return:
    None
");
    void genPointsInTri(unsigned tidx, unsigned npnts, double* cords, int cord_size);
    
    %feature("autodoc", 
"
For each tetrahedron index in indices, randomly generate a set of point coordinates x,y,z within the tetrahedron, where n is
stored in point_counts. The number of points required to be generated for tetrahedron indices[i] is point_counts[i].
All generated points are stored in cords.

Syntax::

    genTetVisualPointsNP(indices, point_counts, coords)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<uint, length = length(indices)> point_counts
    * numpy.array<double, length = sum(point_counts) * 3> coords
    
Return:
    None
");
    void genTetVisualPointsNP(unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double* cords, int cord_size) const;
    
    %feature("autodoc", 
"
For each triangle index in indices, randomly generate a set of point coordinates x,y,z within the triangle, where n is
stored in point_counts. The number of points required to be generated for triangle indices[i] is point_counts[i].
All generated points are stored in cords.

Syntax::

    genTriVisualPointsNP(indices, point_counts, coords)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<uint, length = length(indices)> point_counts
    * numpy.array<double, length = sum(point_counts) * 3> coords
    
Return:
    None
");
    void genTriVisualPointsNP(unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double* cords, int cord_size) const;

    %feature("autodoc", 
"
Get the volumes of a list of tetrahedrons in indices and stored in volumes.

Syntax::

    getBatchTetVolsNP(indices, volumes)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = length(indices)> volumes
    
Return:
    None
");
    void getBatchTetVolsNP(unsigned int* indices, int index_size, double* volumes, int volume_size) const;
    
    %feature("autodoc", 
"
Get the areas of a list of triangles in indices and stored in areas.

Syntax::

    getBatchTriAreasNP(indices, areas)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = length(indices)> areas
    
Return:
    None
");
    void getBatchTriAreasNP(unsigned int* indices, int index_size, double* areas, int area_size) const;

    %feature("autodoc", 
"
Reduce the number of random point coordinates generated for each tetrahedron in indices so that the point density of the tetrahedron is below max_density. If the density is already below max_density for that tetrahedron, the count stored in point_counts is intacted.

Syntax::

    reduceBatchTetPointCountsNP(indices, point_counts, max_density)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = length(indices)> point_counts
    * double max_density
    
Return:
    None
");
    void reduceBatchTetPointCountsNP(unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density);
    
    %feature("autodoc", 
"
Reduce the number of random point coordinates generated for each triangle in indices so that the point density of the triangle is below max_density. If the density is already below max_density for that triangle, the count stored in point_counts is intacted.

Syntax::

    reduceBatchTriPointCountsNP(indices, point_counts, max_density)

Arguments:
    * numpy.array<uint> indices
    * numpy.array<double, length = length(indices)> point_counts
    * double max_density
    
Return:
    None
");
    void reduceBatchTriPointCountsNP(unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density);

    //std::vector<int> getTetsTetNeighbSet(std::vector<uint> const & t_indices) const;
    
    
    ////////////////////////////////////////////////////////////////////////
    // ROI Recording
    ////////////////////////////////////////////////////////////////////////
             
    %feature("autodoc", 
"
Add a Region of Interest data record with name id to the ROI dataset.
The type of elements stored in the ROI data can be one of the follows:
steps.geom.ELEM_VERTEX, steps.geom.ELEM_TET, steps.geom.ELEM_TRI, steps.geom.ELEM_UNDEFINED.

Syntax::

    addROI(id, type, indices)

Arguments:
    * string id
    * ElementType type
    * list<uint> indices
    
Return:
    None
");
    void addROI(std::string id, steps::tetmesh::ElementType type, std::set<uint> const &indices);
     
    %feature("autodoc", 
"
Remove a Region of Interest data record with name id.

Syntax::

    removeROI(id)

Arguments:
    * string id
    
Return:
    None
");
    void removeROI(std::string id);
     
    %feature("autodoc", 
"
Replace a Region of Interest data record with name id with new data.

Syntax::

    replaceROI(id, type, indices)

Arguments:
    * string id
    * ElementType type
    * list<uint> indices
    
Return:
    None
");
    void replaceROI(std::string id, steps::tetmesh::ElementType type, std::set<uint> const &indices);
     
    %feature("autodoc", 
"
Get the element type of a Region of Interest data record with name id.

Syntax::

    getROIType(id)

Arguments:
    * string id
    
Return:
    ElementType
");
    steps::tetmesh::ElementType getROIType(std::string id) const;
     
    %feature("autodoc", 
"
Get the stored data of a Region of Interest data record with name id.

Syntax::

    getROIData(id)

Arguments:
    * string id
    
Return:
    list<uint>
");
    std::vector<uint> getROIData(std::string id) const;
     
    %feature("autodoc", 
"
Get the number of elements stored in a Region of Interest data record with name id.

Syntax::

    getROIDataSize(id)

Arguments:
    * string id
    
Return:
    uint
");
    uint getROIDataSize(std::string id) const;
     
    %feature("autodoc", 
"
Get the number of Region of Interest data stored in the ROI dataset.

Syntax::

    getNROIs()

Arguments:
    None
    
Return:
    uint
");
    uint getNROIs(void);
     
    %feature("autodoc", 
"
Get a Region of Interest data record with name id.

Syntax::

    getROI(id)

Arguments:
    * string id
    
Return:
    ROISet
");
    ROISet getROI(std::string id) const;
     
    %feature("autodoc", 
"
Get a list of the names of all Region of Interest data stored in ROI dataset.

Syntax::

    getAllROINames(id)

Arguments:
    None
    
Return:
    list<string>
");
    std::vector<std::string> getAllROINames(void);
     
    ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////

    %feature("autodoc", 
"
Get barycentres of elements stored in a tetrahedral ROI.

Syntax::

    getROITetBarycentres(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    list<double>
");
    std::vector<double> getROITetBarycentres(std::string ROI_id) const;

    %feature("autodoc", 
"
Get barycentres of elements stored in a tetrahedral ROI and write to a NumPy array centres.
The size of centres should be the same as the number of elements stored in the ROI.

Syntax::

    getROITetBarycentresNP(ROI_id, centres)

Arguments:
    * string ROI_id
    * numpy.array<double> centres
    
Return:
    None
");
    void getROITetBarycentresNP(std::string ROI_id, double* centres, int output_size) const;

    %feature("autodoc", 
"
Get barycentres of elements stored in a triangular ROI.

Syntax::

    getROITriBarycentres(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    list<double>
");
    std::vector<double> getROITriBarycentres(std::string ROI_id) const;

    %feature("autodoc", 
"
Get barycentres of elements stored in a triangular ROI and write to a NumPy array centres.
The size of centres should be the same as the number of elements stored in the ROI.

Syntax::

    getROITriBarycentresNP(ROI_id, centres)

Arguments:
    * string ROI_id
    * numpy.array<double> centres
    
Return:
    None
");
    void getROITriBarycentresNP(std::string ROI_id, double* centres, int output_size) const;

    %feature("autodoc", 
"
Get coordinates of elements stored in a vertices ROI.

Syntax::

    getROIVertices(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    list<double>
");
    std::vector<double> getROIVertices(std::string ROI_id) const;

    %feature("autodoc", 
"
Get coordinates of elements stored in a vertices ROI and write to a NumPy array coordinates.
The size of coordinates should be the same as the number of elements stored in the ROI.

Syntax::

    getROIVerticesNP(ROI_id, coordinates)

Arguments:
    * string ROI_id
    * numpy.array<double> coordinates
    
Return:
    None
");
    void getROIVerticesNP(std::string ROI_id, double* coordinates, int output_size) const;

    %feature("autodoc", 
"
Get vertices of elements stored in a triangular ROI.

Syntax::

    getROITris(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    list<uint>
");
    std::vector<uint> getROITris(std::string ROI_id) const;

    %feature("autodoc", 
"
Get vertices of elements stored in a triangular ROI and write to a NumPy array v_indices.
The size of v_indices should be 3 * the number of elements stored in the ROI.

Syntax::

    getROITrisNP(ROI_id, v_indices)

Arguments:
    * string ROI_id
    * numpy.array<uint> v_indices
    
Return:
    None
");
    void getROITrisNP(std::string ROI_id, unsigned int* v_indices, int output_size) const;

    %feature("autodoc", 
"
Get vertices of elements stored in a tetrahedral ROI.

Syntax::

    getROITets(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    list<uint>
");
    std::vector<uint> getROITets(std::string ROI_id) const;

    %feature("autodoc", 
"
Get vertices of elements stored in a tetrahedral ROI and write to a NumPy array v_indices.
The size of v_indices should be 3 * the number of elements stored in the ROI.

Syntax::

    getROITetsNP(ROI_id, v_indices)

Arguments:
    * string ROI_id
    * numpy.array<uint> v_indices
    
Return:
    None
");
    void getROITetsNP(std::string ROI_id, unsigned int* v_indices, int output_size) const;

    %feature("autodoc", 
"
Add all vertex indices of a list of triangles in a ROI to a set and return its size.

Syntax::

    getROITriVerticesSetSizeNP(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    uint
");
    uint getROITriVerticesSetSizeNP(std::string ROI_id) const;

    %feature("autodoc", 
"
Add all vertex indices of a list of tetrahedrons in a ROI to a set and return its size.

Syntax::

    getROITetVerticesSetSizeNP(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    uint
");
    uint getROITetVerticesSetSizeNP(std::string ROI_id) const;

    %feature("autodoc", 
"
Add all vertex indices of a list of triangles in a ROI to a set and write it to a NumPy array v_set.
For each of the triangle, t_vertices records the positions of its vertices in v_set.
i.e. For the i triangle in the ROI, the STEPS indices of its vertices are
v_set[t_vertices[3*i]], v_set[t_vertices[3*i + 1]], v_set[t_vertices[3*i + 2]]

Syntax::

    getROITriVerticesMappingSetNP(ROI_id, t_vertices, v_set)

Arguments:
    * string ROI_id
    * numpy.array<uint> v_indices
    * numpy.array<uint> v_set
    
Return:
    None
");
    void getROITriVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;

    %feature("autodoc", 
"
Add all vertex indices of a list of tetrahedrons in a ROI to a set and write it to a NumPy array v_set.
For each of the tetrahedron, t_vertices records the positions of its vertices in v_set.
i.e. For the i tetrahedron in the ROI, the STEPS indices of its vertices are
v_set[t_vertices[4*i]], v_set[t_vertices[4*i + 1]], v_set[t_vertices[4*i + 2]], v_set[t_vertices[4*i + 3]]

Syntax::

    getROITetVerticesMappingSetNP(ROI_id, t_vertices, v_set)

Arguments:
    * string ROI_id
    * numpy.array<uint> v_indices
    * numpy.array<uint> v_set
    
Return:
    None
");
    void getROITetVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;

    %feature("autodoc", 
"
For each tetrahedron index in a ROI, randomly generate a set of point coordinates x,y,z within the tetrahedron, where n is
stored in point_counts. The number of points required to be generated for tetrahedron i in the ROI is point_counts[i].
All generated points are stored in cords.

Syntax::

    genTetVisualPointsNP(ROI_id, point_counts, coords)

Arguments:
    * string ROI_id
    * numpy.array<uint> point_counts
    * numpy.array<double> coords
    
Return:
    None
");
    void genROITetVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* cords, int cord_size) const;

    %feature("autodoc", 
"
For each triangle index in a ROI, randomly generate a set of point coordinates x,y,z within the triangle, where n is
stored in point_counts. The number of points required to be generated for triangle i in the ROI is point_counts[i].
All generated points are stored in cords.

Syntax::

    genROITriVisualPointsNP(ROI_id, point_counts, coords)

Arguments:
    * string ROI_id
    * numpy.array<uint> point_counts
    * numpy.array<double> coords
    
Return:
    None
");
    void genROITriVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* cords, int cord_size) const;

    %feature("autodoc", 
"
Get the volumes of a list of tetrahedrons in a ROI and stored in volumes.

Syntax::

    getROITetVolsNP(ROI_id, volumes)

Arguments:
    * string ROI_id
    * numpy.array<double> volumes
    
Return:
    None
");
    void getROITetVolsNP(std::string ROI_id, double* volumes, int volume_size) const;

    %feature("autodoc", 
"
Get the areas of a list of triangles in a ROI and stored in areas.

Syntax::

    getROITriAreasNP(ROI_id, areas)

Arguments:
    * string ROI_id
    * numpy.array<double> areas
    
Return:
    None
");
    void getROITriAreasNP(std::string ROI_id, double* areas, int area_size) const;

    %feature("autodoc", 
"
Get the total volume of tetrahedral ROI.

Syntax::

    getROIVol(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    Total volume of the ROI
");
    double getROIVol(std::string ROI_id) const;

    %feature("autodoc", 
"
Get the total area of a triangular ROI.

Syntax::

    getROIArea(ROI_id)

Arguments:
    * string ROI_id
    
Return:
    Total area of the ROI
");
    double getROIArea(std::string ROI_id) const;
    
    %feature("autodoc", 
"
Reduce the number of random point coordinates generated for each tetrahedron in a ROI so that the point density of the tetrahedron is below max_density. If the density is already below max_density for that tetrahedron, the count stored in point_counts is intacted.

Syntax::

    reduceROITetPointCountsNP(indices, point_counts, max_density)

Arguments:
    * string ROI_id
    * numpy.array<double> point_counts
    * double max_density
    
Return:
    None
");
    void reduceROITetPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density);

    %feature("autodoc", 
"
Reduce the number of random point coordinates generated for each triangle in a ROI so that the point density of the triangle is below max_density. If the density is already below max_density for that triangle, the count stored in point_counts is intacted.

Syntax::

    reduceROITriPointCountsNP(ROI_id, point_counts, max_density)

Arguments:
    * string ROI_id
    * numpy.array<double, length = length(indices)> point_counts
    * double max_density
    
Return:
    None
");
    void reduceROITriPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density);
};

////////////////////////////////////////////////////////////////////////////////
/*
Tetmesh * loadASCII(std::string pathname);

void saveASCII(std::string pathname, Tetmesh * m);
*/
////////////////////////////////////////////////////////////////////////////////

/* /////////////////////////////////////////////////////////////////////////////
//////// OBJECT REMOVED BECAUSE OF MEMORY ISSUES. SEE TODO NOTE IN C++ /////////
///////////////////////// CONSTRUCTOR FOR DETAILS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class Tri	
{
	
public:
	
	Tri(Tetmesh * mesh, unsigned int tidx);
	~Tri(void);
	
	unsigned int getIdx(void);
	double getArea(void) const;
	std::vector<double> getBarycenter(void) const;
	std::vector<double> getNorm(void) const;
	steps::tetmesh::TmPatch * getPatch(void) const;
	
	Tet getTet0(void) const;
    Tet getTet1(void) const;
    
    
    Tet getInnerTet(void) const;
    Tet getOuterTet(void) const;
    
    
    unsigned int getTet0Idx(void) const;
    unsigned int getTet1Idx(void) const;

    
    unsigned int getInnerTetIdx(void) const;
    unsigned int getOuterTetIdx(void) const;
    

    unsigned int getVertex0Idx(void) const;
    unsigned int getVertex1Idx(void) const;
    unsigned int getVertex2Idx(void) const;
	
};
	
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//////// OBJECT REMOVED BECAUSE OF MEMORY ISSUES. SEE TODO NOTE IN C++ /////////
///////////////////////// CONSTRUCTOR FOR DETAILS //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class Tet
{
	
public:
	
	Tet(Tetmesh * mesh, unsigned int tidx);
	~Tet(void);
	
	unsigned int getIdx(void) const;
	double getVol(void) const;
	std::vector<double> getBarycenter(void) const;
	double getQualityRER(void) const;
	steps::tetmesh::TmComp * getComp(void) const;
	
	Tet getTet0(void) const;
    Tet getTet1(void) const;
    Tet getTet2(void) const;
    Tet getTet3(void) const;

	unsigned int getTet0Idx(void) const;
    unsigned int getTet1Idx(void) const;
    unsigned int getTet2Idx(void) const;
    unsigned int getTet3Idx(void) const;
	
    double getTet0Dist(void) const;
    double getTet1Dist(void) const;
    double getTet2Dist(void) const;
    double getTet3Dist(void) const;
	
    Tri getTri0(void) const;
    Tri getTri1(void) const;
    Tri getTri2(void) const;
    Tri getTri3(void) const;
	
    unsigned int getTri0Idx(void) const;
    unsigned int getTri1Idx(void) const;
    unsigned int getTri2Idx(void) const;
    unsigned int getTri3Idx(void) const;

    double getTri0Dist(void) const;
    double getTri1Dist(void) const;
    double getTri2Dist(void) const;
    double getTri3Dist(void) const;
	
    double getTri0Area(void) const;
    double getTri1Area(void) const;
    double getTri2Area(void) const;
    double getTri3Area(void) const;

    unsigned int getVertex0Idx(void) const;
    unsigned int getVertex1Idx(void) const;
    unsigned int getVertex2Idx(void) const;
    unsigned int getVertex3Idx(void) const;
	
	unsigned int getVertex0(void) const;
    unsigned int getVertex1(void) const;
    unsigned int getVertex2(void) const;
    unsigned int getVertex3(void) const;
	
    bool isInside(std::vector<double> p) const;
	std::vector<double> getRanPnt(steps::rng::RNG * r, unsigned int n = 1) const;

};

*/
////////////////////////////////////////////////////////////////////////////////
	
%feature("kwargs") TmComp::TmComp;

class TmComp : public steps::wm::Comp
{
public:
		
	TmComp(std::string const & id, Tetmesh * container, 
		 std::vector<unsigned int> const & tets);
	virtual ~TmComp(void);
	
    %feature("autodoc", "Obsolete");
	virtual void setVol(double vol);
	
    %feature("autodoc", 
"
Returns a list of indices of all tetrahedrons assigned to the compartment. 

Syntax::

    getAllTetIndices()

Arguments:
    None
             
Return:
    list<uint>
");

	std::vector<unsigned int> getAllTetIndices(void) const;
    
    %feature("autodoc", 
"
Returns the number of tetrahedrons assigned to the compartment. 

Syntax::

    countTets()

Arguments:
    None
             
Return:
    uint
");
    unsigned int countTets(void) const;

    %feature("autodoc", 
"
Returns a list of Booleans describing if tetrahedrons tets are 
assigned to the compartment.

Syntax::
    
    isTetInside(tets)

Arguments:
    list<uint> tets
             
Return:
    list<bool, length = length(tets)>
");
	std::vector<bool> isTetInside(std::vector<unsigned int> tet) const;
	
    %feature("autodoc", 
"
Returns the minimal Cartesian coordinate of the rectangular bounding box 
of the compartment. 

Syntax::

    getBoundMin()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMin(void) const;
    
    %feature("autodoc", 
"
Returns the maximal Cartesian coordinate of the rectangular bounding box 
of the compartment. 

Syntax::

    getBoundMax()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMax(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") TmPatch::TmPatch;

class TmPatch : public steps::wm::Patch
{
public:
	/*
	TmPatch(std::string const & id, Tetmesh * container,
		  std::vector<unsigned int> const & tris, steps::tetmesh::TmComp * icomp,
		  steps::tetmesh::TmComp * ocomp = 0, steps::model::Surfsys * surfsys = 0);
	*/
	TmPatch(std::string const & id, Tetmesh * container,
			std::vector<unsigned int> const & tris, steps::wm::Comp* icomp,
			steps::wm::Comp* ocomp = 0);
	virtual ~TmPatch(void);
		
    %feature("autodoc", 
"
Returns a list of indices of all triangles assigned to the patch.

Syntax::

    getAllTriIndices()

Arguments:
    None
             
Return:
    list<uint>
");
	std::vector<unsigned int> getAllTriIndices(void) const;
    
    %feature("autodoc", 
"
Returns a list of Booleans describing if triangles tris are 
assigned to the patch.

Syntax::

    isTriInside(tris)

Arguments:
    list<uint> tris
             
Return:
    list<bool, length = length(tris)>
");
	std::vector<bool> isTriInside(std::vector<unsigned int> tet) const;
    
	
    %feature("autodoc", 
"
Returns the minimal Cartesian coordinate of the rectangular bounding box 
of the compartment. 

Syntax::

    getBoundMin()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMin(void) const;
    
    %feature("autodoc", 
"
Returns the maximal Cartesian coordinate of the rectangular bounding box 
of the compartment. 

Syntax::

    getBoundMax()

Arguments:
    None
             
Return:
    list<float, length = 3>
");
	std::vector<double> getBoundMax(void) const;
	
};
		
////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") Memb::Memb;

class Memb
{		
public:
	
	Memb(std::string const & id, Tetmesh * container,
		 std::vector<TmPatch *> const & patches, 
         bool verify=false, uint opt_method = 1, double search_percent=100.0, std::string const & opt_file_name= "");
	~Memb(void);

    %feature("autodoc", 
"
Get the identifier string of the membrane.
			 
Syntax::
			 
    getID()
			 
Arguments:
    None
             
Return:
    string
");
	std::string getID(void) const;
	
%feature("autodoc", 
"
Returns a reference to the parent steps.geom.Tetmesh container object.
			 
Syntax::
			 
    getTetmesh()
			 
Arguments:
    None
             
Return:
    steps.tetmesh.Tetmesh
");
	steps::tetmesh::Tetmesh * getContainer(void) const;
	
	%feature("autodoc",
"
Returns a list of Booleans describing if triangles tris are 
assigned to the membrane.

Syntax::

    isTriInside(tris)

Arguments:
    list<uint> tris
             
Return:
    list<bool, length = length(tris)>
");

	std::vector<bool> isTriInside(std::vector<unsigned int> tris) const;
	
    %feature("autodoc", 
"
Returns a list of indices of all triangles assigned to the membrane.

Syntax::

    getAllTriIndices()

Arguments:
    None
             
Return:
    list<uint>
");	
	std::vector<unsigned int> getAllTriIndices(void) const;

    %feature("autodoc", 
"
Returns a list of indices of all tetrahedrons assigned to the conduction volume.

Syntax::
			 
    getAllVolTetIndices()
			 
Arguments:
    None

Return:
    list<uint>
");	
	std::vector<unsigned int> getAllVolTetIndices(void) const;	
	
	%feature("autodoc",
"
Returns a list of all vertices in the conduction volume.

Syntax:
			
    getAllVertices()

Arguments:
	None

Return:
	list<uint>
");
	std::vector<unsigned int> getAllVertIndices(void) const;
	
    %feature("autodoc", 
" 
Returns a list of all virtual triangles for the membrane forming a closed surface.

Syntax:

    getAllVirtTris()

Arguments:
    None

Return:
    list<uint>

");
	std::vector<unsigned int> getAllVirtTriIndices(void) const;
	
	%feature("autodoc",
"
Returns the number of tetrahedrons assigned to the conduction volume.
	
Syntax::
	
	countVolTets()
	
Arguments:
	None

Return:
	uint
	
");
	unsigned int countVolTets(void) const;

    %feature("autodoc", 
" 
Returns the number of virtual triangles for the membrane forming a closed surface.

Syntax:
	
	countTris()
	
Arguments:
	None
	
Return:
	uint

");
	unsigned int countVirtTris(void) const;
    
	%feature("autodoc", 
" 
Returns the number of triangles assigned to the membrane.
			 
Syntax:
			 
    countTris()
			 
Arguments:
    None
			 
Return:
    uint
			 
");
	unsigned int countTris(void) const;	
	
	%feature("autodoc",
" 
Returns the number of vertices in the conduction volume and membrane surface.
	
Syntax:
	
	countVertices()

Arguments:
	None
	
Returns:
	uint
	
");
	unsigned int countVerts(void) const;
	
};

////////////////////////////////////////////////////////////////////////////////

%feature("kwargs") DiffBoundary::DiffBoundary;
    
class DiffBoundary
{
        
public:
        
    DiffBoundary(std::string const & id, Tetmesh * container,
                std::vector<unsigned int> const & tris);
        
    virtual ~DiffBoundary(void);
        
%feature("autodoc", 
"
Get the identifier string of the diffusion boundary.

Syntax::

    getID()

Arguments:
    None

Return:
    string
");
    std::string getID(void) const;
        
%feature("autodoc", 
"
Set the identifier string of the diffusion boundary.

Syntax::

    setID(name)

Arguments:
    string name

Return:
    None
");
    void setID(std::string const & id);
        
%feature("autodoc", 
"
Returns a reference to the parent steps.tetmesh.Tetmesh container object.

Syntax::

    getContainer()

Arguments:
    None

Return:
    steps.tetmesh.Tetmesh
");
    steps::tetmesh::Tetmesh * getContainer(void) const;
        
%feature("autodoc", 
"
Returns a list of the two compartments this diffusion boundary connects.

Syntax::

    getComps()

Arguments:
    None

Return:
    list<steps::wm::Comp, length = 2>
");
        std::vector<steps::wm::Comp *>  getComps(void) const;
        
        %feature("autodoc", 
"
Returns a list of Booleans describing if triangles tris are 
assigned to the Diffusion Boundary.

Syntax::

    isTriInside(tris)

Arguments:
    list<uint> tris

Return:
    list<bool, length = length(tris)>
");
    std::vector<bool> isTriInside(std::vector<unsigned int> tri) const;
        
    %feature("autodoc", 
"
Returns a list of indices of all triangles assigned to the patch.

Syntax::

    getAllTriIndices()

Arguments:
    None

Return:
    list<uint>
");
        std::vector<unsigned int> getAllTriIndices(void) const;                
};	

////////////////////////////////////////////////////////////////////////////////

}	// end namespace tetmesh
}	// end namespace steps

// END
	

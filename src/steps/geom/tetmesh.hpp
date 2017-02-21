/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#ifndef STEPS_TETMESH_TETMESH_HPP
#define STEPS_TETMESH_TETMESH_HPP 1


// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"
#include "steps/math/bbox.hpp"
#include "steps/geom/geom.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/geom/memb.hpp"
#include "steps/geom/diffboundary.hpp"

// STL headers
#include <vector>
#include <map>
#include <set>
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class TmPatch;
class TmComp;
class Memb;
class DiffBoundary;

enum ElementType {ELEM_VERTEX, ELEM_TRI, ELEM_TET, ELEM_UNDEFINED = 99};

struct ROISet {
    ROISet(): type(ELEM_UNDEFINED) {}

    ROISet(ElementType t, std::set<uint> const &i):
        type(t), indices(i.begin(),i.end()) {}
    
    ElementType        type;
    std::vector<uint>  indices;
};

////////////////////////////////////////////////////////////////////////////////

/// The main container class for static tetrahedronl meshes.
/*!
This class stores the vertices points, tetrahedron and boundary triangles
that comprise the tetrahedronl mesh. In addition, it also precomputes
some auxiliary data for the mesh as a whole:

    - Rectangular, axis-aligned bounding box.
    - Overall volume

Auxiliary data is also stored for the tetrahedrons:

    - Volume of each tetrahedron.
    - For each tetrahedron, the indices of its 4 neighbouring tets.
      If there is no neighbour (i.e. if the tetrahedron lies on the
      border), this index will be -1. The sequence of neighbours is
      determined by the following common boundary triangles: (0,1,2);
      (0,1,3); (0,2,3); (1,2,3).
    - For each tetrahedron, the indices of its 4 neighbouring
      boundary triangles. The sequence of neighbours is also
      determined by (0,1,2); (0,1,3); (0,2,3); (1,2,3).
    - The compartment (Comp object) that a tetrahedron belongs to.
      Returns zero pointer if the tetrahedron has not been added to a comp
    - The total number of tetrahedron in the mesh
    - A method of finding which tetrahedron a point given in x,y,z
      coordinates belongs to

And for the triangles:

    - Area of each triangle.
    - Normal vector for each triangle, normalized to length 1.0.
    - For each triangle, the indices of its inside and outside
      neighbouring tetrahedron. If this tetrahedron does not exist
      (because the triangle lies on the outer boundary), this index
      will be -1.
    - The patch (Patch object) that a triangle belongs to. Returns
      zero pointer if triangle has not been added to a patch
    - The total number of triangles in the mesh

And, finally, for the vertices:
    - The total number of vertices in the mesh

NOTES:
    - Adding/deleting/moving vertices, triangles and tetrahedron after
      initiation is currently not implemented

    \warning Methods start with an underscore are not exposed to Python.
*/

class Tetmesh : public steps::wm::Geom
{
public:
    typedef steps::math::point3d point3d;

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param verts List of vertices.
    /// \param tets List of tetrahedrons.
    /// \param tris List of triangles.
    Tetmesh(std::vector<double> const & verts, std::vector<uint> const & tets,
            std::vector<uint> const & tris = std::vector<uint>());

    /// Constructor
    ///
    /// \param verts
    Tetmesh(std::vector<double> const & verts, std::vector<uint> const & tris,
            std::vector<double> const & tri_areas,
            std::vector<double> const & tri_norms,
            std::vector<int> const & tri_tet_neighbs,
            std::vector<uint> const & tets,
            std::vector<double> const & tet_vols,
            std::vector<double> const & tet_barycs,
            std::vector<uint> const & tet_tri_neighbs,
            std::vector<int> const & tet_tet_neighbs);

    /// Destructor
    virtual ~Tetmesh(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): VERTICES
    ////////////////////////////////////////////////////////////////////////

    /// Return the coordinates of a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    std::vector<double> getVertex(uint vidx) const;

    /// Count the vertices in the Temesh.
    ///
    /// \return Number of the vertices.
    inline uint countVertices(void) const
    { return pVertsN; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): BARS
    ////////////////////////////////////////////////////////////////////////

    /// Return the bar with index bidx
    ///
    /// \param bidx Index of the bar.
    /// \return Indices of the two vertices that form the bar.
    std::vector<uint> getBar(uint bidx) const;

    /// Count the bars in the Tetmesh.
    ///
    /// \return Number of bars.
    inline uint countBars(void) const
    { return pBarsN; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): TRIANGLES
    ////////////////////////////////////////////////////////////////////////

    /// Return the triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Indices of the vertices that form the triangle.
    std::vector<uint> getTri(uint tidx) const;

    /// Count the triangles in the Temesh.
    ///
    /// \return Number of the triangles.
    inline uint countTris(void) const
    { return pTrisN; }

    /// Return the area of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Area of the triangle.
    double getTriArea(uint tidx) const;

    /// Return the bars of a triangle
    ///
    /// \param tidx Index of the triangle
    /// \return Bars of the triangle
    std::vector<uint> getTriBars(uint tidx) const;

    /// Return the barycentre of triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Barycentre of the triangle.
    std::vector<double> getTriBarycenter(uint tidx) const;

    /// Return the normalised triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Coordinate of the normalised vertices form the triangle.
    std::vector<double> getTriNorm(uint tidx) const;

    //getTriBarycentre
    /// Return the patch which a triangle associated to.
    ///
    /// \param tidx Index of the triangle.
    /// \return Pointer to the patch.
    steps::tetmesh::TmPatch * getTriPatch(uint tidx) const;

    ///Set the patch which a triangle belongs to.
    ///
    /// \param tidx Index of the triangle.
    /// \param patch Pointer to the associated patch.
    void setTriPatch(uint tidx, steps::tetmesh::TmPatch * patch);

    ///Set the diffusion boundary which a triangle belongs to.
    ///
    /// \param tidx Index of the triangle.
    /// \param patch Pointer to the associated diffusion boundary.
    void setTriDiffBoundary(uint tidx, steps::tetmesh::DiffBoundary * diffb);

    /// Return the diffusion boundary which a triangle is associated to.
    ///
    /// \param tidx Index of the triangle.
    /// \return Pointer to the diffusion boundary.
    steps::tetmesh::DiffBoundary * getTriDiffBoundary(uint tidx) const;

    ///Return the tetrahedron neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Vector of the tetrahedron neighbors.
    ///
    std::vector<int> getTriTetNeighb(uint tidx) const;

    ///Return the 3 triangle neighbors within the same patch of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \param tmpatch Pointer to the patch
    /// \return Vector of the triangle neighbors.
    ///
    std::vector<int> getTriTriNeighb(uint tidx, const TmPatch * tmpatch) const;

    ///Return all the triangle neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Set of the triangle neighbors.
    ///
    std::set<uint> getTriTriNeighbs(uint tidx) const;

    /// Flip the triangle's inner and outer tetrahedron.
    ///
    /// \param tidx Index of the triangle.
    void _flipTriTetNeighb(uint tidx);

    /// Flip the triangle's vertices and recalculate the normal.
    ///
    /// \param Index of the triangle.
    void _flipTriVerts(uint tidx);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): TETRAHEDRA
    ////////////////////////////////////////////////////////////////////////
    /// Return a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the indices of triangles which form the tetrahedron.
    std::vector<uint>  getTet(uint tidx) const;

    /// Count the number of tetrahedrons.
    ///
    /// \return Number of tetrahedrons.
    inline uint countTets(void) const
    { return pTetsN; }

    /// Return the volume of a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Volume of the tetrahedron.

    double getTetVol(uint tidx) const;

    /// Computes the quality of the tetrahedron.
    ///
    /// This method uses the radius-edge ratio (RER) metric for tetrahedron
    /// quality, given by dividing the radius of the tetrahedron's
    /// circumsphere with the length of the shortest edge.
    ///
    /// The smaller this value, the more regular the tetrahedron. The
    /// lowest possible value of this metric is given by computing the
    /// RER for a fully regular tetrahedron:
    ///
    ///    Q = sqrt(6)/4 ~ 0.612
    ///
    /// This is a slightly weaker metric than getQualityAR, because
    /// certain slivers (degenerate tetrahedrons) can still have a fairly
    /// small value.
    ///
    /// \return Quality RER of the tetrahedron.
    double getTetQualityRER(uint tidx) const;

    /// Return the barycentre of the tetrahedron in x,y,z coordinates
    /// \param tidx Index of the tetrahedron
    /// \return Barycentre of the tetrahedron
    std::vector<double> getTetBarycenter(uint tidx) const;

    /// Return the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Pointer to the compartment object.

    steps::tetmesh::TmComp * getTetComp(uint tidx) const;
    ///Set the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param comp Pointer to the compartment object.

    void setTetComp(uint tidx, steps::tetmesh::TmComp * comp);
    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the triangle neighbors.

    std::vector<uint> getTetTriNeighb(uint tidx) const;
    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the tetrahedron neighbors.

    std::vector<int> getTetTetNeighb(uint tidx) const;

    /// Find a tetrahedron which encompasses a given point.
    /// Return the index of the tetrahedron that encompasses point;
    ///  return -1 if point is outside mesh;
    /// if point is on boundary between two or more tetrahedron,
    /// returns first tetrahedron found.
    /// \param p A point given by its coordinates.
    /// \return ID of the found tetrahedron.

    int findTetByPoint(std::vector<double> p) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): MESH
    ////////////////////////////////////////////////////////////////////////

    /// Return the minimal coordinate of the rectangular bounding box.
    ///
    /// \return Minimal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMin(void) const;

    /// Return the maximal coordinate of the rectangular bounding box.
    ///
    /// \return Maximal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMax(void) const;

    /// Return the total volume of the mesh.
    ///
    /// \return Volume of the mesh.
    double getMeshVolume(void) const;

    /// Return the triangles which form the surface boundary of the mesh.
    /// \return Vector of the triangle boundary.
    // Weiliang 2010.02.02
    std::vector<int> getSurfTris(void) const;

    ////////////////////////////////////////////////////////////////////////
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////
    
    /// Get barycentres of a list of tetrahedrons
    std::vector<double> getBatchTetBarycentres(std::vector<uint> const & tets) const;
    
    /// Get barycentres of a list of tetrahedrons
    void getBatchTetBarycentresNP(const unsigned int* indices, int input_size, double* centres, int output_size) const;
    
    /// Get barycentres of a list of triangles
    std::vector<double> getBatchTriBarycentres(std::vector<uint> const & tris) const;
    
    /// Get barycentres of a list of triangles
    void getBatchTriBarycentresNP(const unsigned int* indices, int input_size, double* centres, int output_size) const;
    
    /// Get coordinates of a list of vertices
    std::vector<double> getBatchVertices(std::vector<uint> const & verts) const;
    
    /// Get coordinates of a list of vertices
    void getBatchVerticesNP(const unsigned int* indices, int input_size, double* coordinates, int output_size) const;
    
    /// Get vertex indices of a list of triangles
    std::vector<uint> getBatchTris(std::vector<uint> const & tris) const;
    
    /// Get vertex indices of a list of triangles
    void getBatchTrisNP(const unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const;
    
    /// Get vertex indices of a list of tetrahedrons
    std::vector<uint> getBatchTets(std::vector<uint> const & tets) const;
    
    /// Get vertex indices of a list of tetrahedrons
    void getBatchTetsNP(const unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const;
    
    /// Return the size of a set with unique vertex indices of a list of triangles
    /// preparation function for furture numpy data access
    uint getTriVerticesSetSizeNP(const unsigned int* t_indices, int input_size) const;
    
    /// Return the size of a set with unique vertex indices of a list of tetrahedrons
    /// preparation function for furture numpy data access
    uint getTetVerticesSetSizeNP(const unsigned int* t_indices, int input_size) const;
    
    /// Get the set with unique vertex indices of a list of triangles, write into given 1D array
    void getTriVerticesMappingSetNP(const unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;
    
    /// Get the set with unique vertex indices of a list of tetrahedrons, write into given 1D array
    void getTetVerticesMappingSetNP(const unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;
    
    /// Generate npnts random points inside tetrahedron t_idx and write the coordinates to coords
    void genPointsInTet(unsigned tidx, unsigned npnts, double* coords, int coord_size) const;
    
    /// Generate npnts random points on triangle t_idx and write the coordinates to coords
    void genPointsInTri(unsigned tidx, unsigned npnts, double* coords, int coord_size) const;
    
    /// Generate the random points required in point_counts for tets in indices and store them in coords
    void genTetVisualPointsNP(const unsigned int* indices, int index_size, const unsigned int* point_counts, int count_size, double* coords, int coord_size) const;
    
    /// Generate the random points required in point_counts for tris in indices and store them in coords
    void genTriVisualPointsNP(const unsigned int* indices, int index_size, const unsigned int* point_counts, int count_size, double* coords, int coord_size) const;
    
    /// get the volumes of a list of tetrahedrons
    void getBatchTetVolsNP(const unsigned int* indices, int index_size, double* volumes, int volume_size) const;
    
    /// get the areas of a list of triangles
    void getBatchTriAreasNP(const unsigned int* indices, int index_size, double* areas, int area_size) const;
    
    /// reduce the number of points required to be generated in a list of tets based on maximum point density
    void reduceBatchTetPointCountsNP(const unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density);
    
    /// reduce the number of points required to be generated in a list of tris based on maximum point density
    void reduceBatchTriPointCountsNP(const unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density);
    
    /// Get tet neighbors for a list of tets, no duplication
    ///std::vector<int> getTetsTetNeighbSet(std::vector<uint> const & t_indices) const;
    
    ////////////////////////////////////////////////////////////////////////
    // ROI (Region of Interest) Data
    ////////////////////////////////////////////////////////////////////////
    
    /// Add a ROI data
    void addROI(std::string id, ElementType type, std::set<uint> const &indices);
    
    /// Remove a ROI data
    void removeROI(std::string id);
    
    /// Replace a ROI data with a new set with the same name
    void replaceROI(std::string id, ElementType type, std::set<uint> const &indices);
    
    /// Return the type of a ROI data
    ElementType getROIType(std::string id) const;
    
    /// Return the data of a ROI
    std::vector<uint> getROIData(std::string id) const;
    
    /// Return the data size of a ROI
    uint getROIDataSize(std::string id) const;
    
    /// get the total number of ROI recorded
    uint getNROIs(void);
    
    /// Return a ROI
    ROISet getROI(std::string id) const;
    
    /// get all ROI names
    std::vector<std::string> getAllROINames(void);
    
    /// return pointer of a ROI data
    uint* _getROIData(std::string id) const;
    
    /// check if a ROI enquire is valid
    bool checkROI(std::string id, ElementType type, uint count = 0, bool warning = true) const;
    
    ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////
    
    /// get barycentres of a list of tetrahedrons
    std::vector<double> getROITetBarycentres(std::string ROI_id) const;
    
    /// get barycentres of a list of tetrahedrons
    void getROITetBarycentresNP(std::string ROI_id, double* centres, int output_size) const;
    
    /// get barycentres of a list of triangles
    std::vector<double> getROITriBarycentres(std::string ROI_id) const;
    
    /// get barycentres of a list of triangles
    void getROITriBarycentresNP(std::string ROI_id, double* centres, int output_size) const;
    
    /// get coordinates of a list of vertices
    std::vector<double> getROIVertices(std::string ROI_id) const;
    
    /// get coordinates of a list of vertices
    void getROIVerticesNP(std::string ROI_id, double* coordinates, int output_size) const;
    
    /// get vertex indices of a list of triangles
    std::vector<uint> getROITris(std::string ROI_id) const;
    
    /// get vertex indices of a list of triangles
    void getROITrisNP(std::string ROI_id, unsigned int* v_indices, int output_size) const;
    
    /// get vertex indices of a list of tetrahedrons
    std::vector<uint> getROITets(std::string ROI_id) const;
    
    /// get vertex indices of a list of tetrahedrons
    void getROITetsNP(std::string ROI_id, unsigned int* v_indices, int output_size) const;
    
    /// return the size of a set with unique vertex indices of a list of triangles
    /// preparation function for furture numpy data access
    uint getROITriVerticesSetSizeNP(std::string ROI_id) const;
    
    /// return the size of a set with unique vertex indices of a list of tetrahedrons
    /// preparation function for furture numpy data access
    uint getROITetVerticesSetSizeNP(std::string ROI_id) const;
    
    /// Get the set with unique vertex indices of a list of triangles, write into given 1D array
    void getROITriVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;
    
    /// Get the set with unique vertex indices of a list of tetrahedrons, write into given 1D array
    void getROITetVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const;
    
    /// Generate the random points required in point_counts for tets in indices and store them in coords
    void genROITetVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const;
    
    /// Generate the random points required in point_counts for tris in indices and store them in coords
    void genROITriVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const;
    
    /// get the volumes of a list of tetrahedrons
    void getROITetVolsNP(std::string ROI_id, double* volumes, int volume_size) const;
    
    /// get the areas of a list of triangles
    void getROITriAreasNP(std::string ROI_id, double* areas, int area_size) const;
    
    /// get the total volume of a tetrahedral roi
    double getROIVol(std::string ROI_id) const;
    
    /// get the total area of a triangular roi
    double getROIArea(std::string ROI_id) const;
    
    /// reduce the number of points required to be generated in a list of tets based on maximum point density
    void reduceROITetPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density);
    
    /// reduce the number of points required to be generated in a list of tris based on maximum point density
    void reduceROITriPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (C++ INTERNAL)
    ////////////////////////////////////////////////////////////////////////

    /// Return a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    const point3d &_getVertex(uint vidx) const { return pVerts[vidx]; }

    /// Return the bar with index bidx
    ///
    /// \param bidx Index of the bar.
    /// \return Indices of the two vertices that form the bar.
    const uint *_getBar(uint bidx) const { return &pBars[bidx][0]; }

    /// Return a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return List of the vertices form the triangle.
    const uint * _getTri(uint tidx) const { return &pTris[tidx][0]; }

    /// Return a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return List of the vertices form the tetrahedron.
    const uint * _getTet(uint tidx) const { return &pTets[tidx][0]; }

    ///Return the tetrahedron neighbors of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Array of the tetrahedron neighbors.
    const int * _getTriTetNeighb(uint tidx) const { return &pTri_tet_neighbours[tidx][0]; }

    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the triangle neighbors.
    const uint * _getTetTriNeighb(uint tidx) const { return &pTet_tri_neighbours[tidx][0]; }

    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the tetrahedron neighbors.
    const int * _getTetTetNeighb(uint tidx) const { return &pTet_tet_neighbours[tidx][0]; }

    /// Return the bars of a triangle
    ///
    /// \param tidx Index of the triangle
    /// \return Bars of the triangle
    const uint *_getTriBars(uint tidx) const { return &pTri_bars[tidx][0]; }

    /// Return the normal to the triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Point representing unit normal to triangle.
    const point3d &_getTriNorm(uint tidx) const { return pTri_norms[tidx]; }

    /// Return the barycenter of the triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Barycenter of triangle.
    const point3d &_getTriBarycenter(uint tidx) const { return pTri_barycs[tidx]; }

    /// Return the barycenter of the tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Barycenter of tetrahedron.
    const point3d &_getTetBarycenter(uint tidx) const { return pTet_barycentres[tidx]; }

    ////////////////////////////////////////////////////////////////////////

    /// Check if a membrane id is occupied.
    ///
    /// \param id ID of the membrane.
    void _checkMembID(std::string const & id) const;

    /// Change the id of a membrane.
    ///
    /// \param o Old id of the membrane.
    /// \param n New id of the membrane.
    void _handleMembIDChange(std::string const & o, std::string const & n);

    /// Add a membrane.
    ///
    /// \param memb Pointer to the membrane.
    void _handleMembAdd(steps::tetmesh::Memb * memb);

    /// Delete a membrane.
    ///
    /// \param patch Pointer to the membrane.
    void _handleMembDel(steps::tetmesh::Memb * memb);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the membranes in the tetmesh container.
    ///
    /// \return Number of membranes.
    inline uint _countMembs(void) const
    { return pMembs.size(); }

    /// Return a membrane with index gidx.
    ///
    /// \param gidx Index of the membrane.
    /// \return Pointer to the membrane.
    steps::tetmesh::Memb * _getMemb(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////

    /// Check if a diffusion boundary id is occupied.
    ///
    /// \param id ID of the diffusion boundary.
    void _checkDiffBoundaryID(std::string const & id) const;

    /// Change the id of a diffusion boundary.
    ///
    /// \param o Old id of the diffusion boundary.
    /// \param n New id of the diffusion boundary.
    void _handleDiffBoundaryIDChange(std::string const & o, std::string const & n);

    /// Add a diffusion boundary.
    ///
    /// \param patch Pointer to the diffusion boundary.
    void _handleDiffBoundaryAdd(steps::tetmesh::DiffBoundary * diffb);

    /// Delete a diffusion boundary.
    ///
    /// \param patch Pointer to the diffusion boundary.
    void _handleDiffBoundaryDel(steps::tetmesh::DiffBoundary * diffb);

    /// Count the diffusion boundaries in the tetmesh container.
    ///
    /// \return Number of diffusion boundaries.
    inline uint _countDiffBoundaries(void) const
    { return pDiffBoundaries.size(); }

    /// Return a diffusion boundary with index gidx.
    ///
    /// \param gidx Index of the diffusion boundary.
    /// \return Pointer to the diffusion boundary.
    steps::tetmesh::DiffBoundary * _getDiffBoundary(uint gidx) const;

private:
    typedef std::array<uint,2> bar_verts;
    typedef std::array<uint,3> tri_verts;
    typedef std::array<uint,4> tet_verts;

    typedef std::array<uint,3> tri_bars;
    typedef std::array<int,2>  tri_tets;

    typedef std::array<uint,4> tet_tris;
    typedef std::array<int,4>  tet_tets;

    /// Build pBars, pBarsN, pTri_bars from pTris.
    void buildBarData();

    ///////////////////////// DATA: VERTICES ///////////////////////////////
    ///
    /// The total number of vertices in the mesh
    uint                                pVertsN;
    /// The vertices by x,y,z coordinates
    std::vector<point3d>                pVerts;

    /////////////////////////// DATA: BARS /////////////////////////////////
    ///
    /// The total number of 1D 'bars' in the mesh
    uint                                pBarsN;
    /// The bars by the two vertices index
    std::vector<bar_verts>              pBars;

    ///////////////////////// DATA: TRIANGLES //////////////////////////////
    ///
    /// The total number of triangles in the mesh
    uint                                pTrisN;
    /// The triangles by vertices index
    std::vector<tri_verts>              pTris;
    // The bars of the triangle
    std::vector<tri_bars>               pTri_bars;
    /// The areas of the triangles
    std::vector<double>                 pTri_areas;
    /// The triangle barycentres
    std::vector<point3d>                pTri_barycs;
    /// The triangle normals
    std::vector<point3d>                pTri_norms;
    /// The patch a triangle belongs to
    std::vector<steps::tetmesh::TmPatch *> pTri_patches;

    /// The diffusion boundary a triangle belongs to
    std::vector<steps::tetmesh::DiffBoundary *> pTri_diffboundaries;

    /// The tetrahedron neighbours of each triangle (by index)
    std::vector<tri_tets>               pTri_tet_neighbours;

    ///////////////////////// DATA: TETRAHEDRA /////////////////////////////
    ///
    /// The total number of tetrahedron in the mesh
    uint                                pTetsN;
    /// The tetrahedron by vertices index
    std::vector<tet_verts>              pTets;
    /// The volume of the tetrahedron
    std::vector<double>                 pTet_vols;
    /// The barycentres of the tetrahedra
    std::vector<point3d>                pTet_barycentres;
    /// The compartment a tetrahedron belongs to
    std::vector<steps::tetmesh::TmComp  *> pTet_comps;
    /// The triangle neighbours of each tetrahedron (by index)
    std::vector<tet_tris>               pTet_tri_neighbours;
    /// The tetrahedron neighbours of each tetrahedron (by index)
    std::vector<tet_tets>               pTet_tet_neighbours;

    ////////////////////////////////////////////////////////////////////////

    /// Information about the minimal and maximal boundary values
    steps::math::bounding_box           pBBox;

    ////////////////////////////////////////////////////////////////////////

    // List of contained membranes. Members of this class because they
    // do not belong in a well-mixed geometry description
    std::map<std::string, steps::tetmesh::Memb *>       pMembs;
    std::map<std::string, steps::tetmesh::DiffBoundary *> pDiffBoundaries;
    
    ////////////////////////// ROI Dataset /////////////////////////////////
    std::map<std::string, ROISet>                       mROI;
};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif

// STEPS_TETMESH_TETMESH_HPP
// END

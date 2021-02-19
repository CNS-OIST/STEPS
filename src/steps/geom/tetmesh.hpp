/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/geom/RegionOfInterest.hpp"
#include "steps/math/point.hpp"
#include "steps/math/bbox.hpp"
#include "steps/geom/fwd.hpp"
#include "steps/geom/geom.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/geom/memb.hpp"
#include "steps/geom/diffboundary.hpp"
#include "steps/geom/sdiffboundary.hpp"

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
class SDiffBoundary;

// FIXME TCL: use proper storage according to type of element
enum ElementType {ELEM_VERTEX, ELEM_TRI, ELEM_TET, ELEM_UNDEFINED = 99};

 struct ROISet {
    using data_type = index_t;
    using set_data_type = std::set<data_type>;
    using vector_data_type = std::vector<data_type>;
    ROISet(): type(ELEM_UNDEFINED) {}

    ROISet(ElementType t, const set_data_type& i):
        type(t), indices(i.begin(),i.end()) {}

    ElementType        type;
    vector_data_type  indices;
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
    using point3d = steps::math::point3d;

    /// \brief store the 4 possible tetrahedrons that are neighbors of a tetrahedron,
    /// \a UNKNOWN_TET if there is no such tet.
    using tet_tets = std::array<tetrahedron_id_t, 4>;
    /// \brief store the 2 vertices of a bar
    using bar_verts = std::array<vertex_id_t, 2>;
    /// \brief store the 3 vertices of a triangle
    using tri_verts = std::array<vertex_id_t, 3>;
    /// \brief store the 4 vertices of a tetrahedron
    using tet_verts = std::array<vertex_id_t, 4>;

    /// \brief store the 3 bars which form a triangle
    using tri_bars = std::array<bar_id_t, 3>;

    /// \brief store the 2 possible triangles that are neighbors of a bar
    /// (on a surface diffusion boundary), UINT_MAX if there is no such triangle
    using bar_tris = std::array<triangle_id_t, 2>;

    /// \brief store the 2 possible tetrahedrons that are neighbors of a triangle,
    /// \a UNKNOWN_TET if there is no such tet
    using tri_tets = std::array<tetrahedron_id_t, 2>;

    /// \brief store the 4 triangles that form a tetrahedron
    using tet_tris = std::array<triangle_id_t, 4>;

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param verts List of vertices.
    /// \param tets List of tetrahedrons.
    /// \param tris List of triangles.
    Tetmesh(std::vector<double> const & verts,
            std::vector<index_t> const & tets,
            std::vector<index_t> const & tris = {});

    /// Constructor
    ///
    /// \param verts
    Tetmesh(std::vector<double> const & verts,
            std::vector<vertex_id_t> const & tris,
            std::vector<double> const & tri_areas,
            std::vector<double> const & tri_norms,
            std::vector<tetrahedron_id_t> const & tri_tet_neighbs,
            std::vector<vertex_id_t> const & tets,
            std::vector<double> const & tet_vols,
            std::vector<double> const & tet_barycs,
            std::vector<triangle_id_t> const & tet_tri_neighbs,
            std::vector<tetrahedron_id_t> const & tet_tet_neighbs);

    /// Destructor
    virtual ~Tetmesh();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): VERTICES
    ////////////////////////////////////////////////////////////////////////

    /// Return the coordinates of a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    std::vector<double> getVertex(vertex_id_t vidx) const;


    /// Return the id of tets sharing a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return tets ids.
    std::vector<tetrahedron_id_t> getVertexTetNeighbs(vertex_id_t vidx) const;

    /// Count the vertices in the Tetmesh.
    ///
    /// \return Number of the vertices.
    inline index_t countVertices() const noexcept
    { return pVertsN; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): BARS
    ////////////////////////////////////////////////////////////////////////

    /// Return the bar with index bidx
    ///
    /// \param bidx Index of the bar.
    /// \return Indices of the two vertices that form the bar.
    std::vector<index_t> getBar(bar_id_t bidx) const;

    /// Count the bars in the Tetmesh.
    ///
    /// \return Number of bars.
    inline unsigned long countBars() const noexcept
    { return pBarsN; }

    ///Set the surface diffusion boundary which a bar belongs to.
    ///
    /// \param tidx Index of the bar.
    /// \param sdiffb Pointer to the associated surface diffusion boundary.
    void setBarSDiffBoundary(bar_id_t bidx, SDiffBoundary * sdiffb);

    /// Return the surface diffusion boundary which a triangle is associated to.
    ///
    /// \param bidx Index of the bar.
    /// \return Pointer to the surface diffusion boundary.
    SDiffBoundary * getBarSDiffBoundary(bar_id_t bidx) const;

    ///Return all the triangle neighbors of a bar by its index.
    ///NOTE: This differs from _getBarTriNeighb which returns the 2 triangles that a bar
    ///separates if it is part of a surface diffusion boundary.
    /// \param bidx Index of the bar.
    /// \return Set of the triangle neighbors.
    ///
    std::set<triangle_id_t> getBarTriNeighbs(bar_id_t bidx) const;

    ///Set neighbouring tris to bar in context of surface diffusion boundary
    /// NOTE: Making one unique list for each bar disallows any bar from
    /// being part of more than one surface diffusion boundary.
    void setBarTris(bar_id_t bidx, triangle_id_t itriidx, triangle_id_t otriidx);


    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): TRIANGLES
    ////////////////////////////////////////////////////////////////////////

    /// Return the triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Indices of the vertices that form the triangle.
    std::vector<index_t> getTri(triangle_id_t tidx) const;

    /// Count the triangles in the Temesh.
    ///
    /// \return Number of the triangles.
    inline index_t countTris() const noexcept
    { return pTrisN; }

    /// Return the area of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Area of the triangle.
    double getTriArea(triangle_id_t tidx) const;

    /// Return the bars of a triangle
    ///
    /// \param tidx Index of the triangle
    /// \return Bars of the triangle
    std::vector<index_t> getTriBars(triangle_id_t tidx) const;

    /// Return the barycenter of triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Barycenter of the triangle.
    std::vector<double> getTriBarycenter(triangle_id_t tidx) const;

    /// Return the normalised triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Coordinate of the normalised vertices form the triangle.
    std::vector<double> getTriNorm(triangle_id_t tidx) const;

    /// Return the patch which a triangle associated to.
    ///
    /// \param tidx Index of the triangle.
    /// \return Pointer to the patch.
    TmPatch * getTriPatch(triangle_id_t tidx) const;

    ///Set the patch which a triangle belongs to.
    ///
    /// \param tidx Index of the triangle.
    /// \param patch Pointer to the associated patch.
    void setTriPatch(triangle_id_t tidx, TmPatch * patch);

    ///Set the diffusion boundary which a triangle belongs to.
    ///
    /// \param tidx Index of the triangle.
    /// \param diffb Pointer to the associated diffusion boundary.
    void setTriDiffBoundary(triangle_id_t tidx, DiffBoundary * diffb);

    /// Return the diffusion boundary which a triangle is associated to.
    ///
    /// \param tidx Index of the triangle.
    /// \return Pointer to the diffusion boundary.
    DiffBoundary * getTriDiffBoundary(triangle_id_t tidx) const;

    ///Return the tetrahedron neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Vector of the tetrahedron neighbors.
    ///
    std::vector<index_t> getTriTetNeighb(triangle_id_t tidx) const;

    ///Return the 3 triangle neighbors within the same patch of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \param tmpatch Pointer to the patch
    /// \return Vector of the triangle neighbors.
    ///
    std::vector<index_t> getTriTriNeighb(triangle_id_t tidx, const TmPatch * tmpatch) const;

    ///Return the 3 triangle neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Vector of the triangle neighbors.
    ///
    /// This function differs from the one above as it doesn't exclude triangle neighbors
    /// in other patches

    std::vector<index_t> getTriTriNeighb(triangle_id_t tidx) const;

    ///Return all the triangle neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Set of the triangle neighbors.
    ///
    std::set<index_t> getTriTriNeighbs(triangle_id_t tidx) const;

    /// Flip the triangle's inner and outer tetrahedron.
    ///
    /// \param tidx Index of the triangle.
    void _flipTriTetNeighb(triangle_id_t tidx);

    /// Flip the triangle's vertices and recalculate the normal.
    ///
    /// \param Index of the triangle.
    void _flipTriVerts(triangle_id_t tidx);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): TETRAHEDRA
    ////////////////////////////////////////////////////////////////////////
    /// Return a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the indices of triangles which form the tetrahedron.
    std::vector<index_t>  getTet(tetrahedron_id_t tidx) const;

    /// Count the number of tetrahedrons.
    ///
    /// \return Number of tetrahedrons.
    inline std::size_t countTets() const noexcept
    { return pTetsN; }

    /// Return the volume of a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Volume of the tetrahedron.

    double getTetVol(tetrahedron_id_t tidx) const;

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
    double getTetQualityRER(tetrahedron_id_t tidx) const;

    /// Return the barycenter of the tetrahedron in x,y,z coordinates
    /// \param tidx Index of the tetrahedron
    /// \return Barycenter of the tetrahedron
    std::vector<double> getTetBarycenter(tetrahedron_id_t tidx) const;

    /// Return the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Pointer to the compartment object.

    TmComp * getTetComp(tetrahedron_id_t tidx) const;
    ///Set the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param comp Pointer to the compartment object.

    void setTetComp(tetrahedron_id_t tidx, TmComp * comp);
    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the triangle neighbors.

    std::vector<index_t> getTetTriNeighb(tetrahedron_id_t tidx) const;
    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the tetrahedron neighbors.

    std::vector<index_t> getTetTetNeighb(tetrahedron_id_t tidx) const;

    /// Find a tetrahedron which encompasses a given point.
    /// Return the index of the tetrahedron that encompasses point;
    ///  return -1 if point is outside mesh;
    /// if point is on boundary between two or more tetrahedron,
    /// returns first tetrahedron found.
    /// \param p A point given by its coordinates.
    /// \return ID of the found tetrahedron.

    tetrahedron_id_t findTetByPoint(std::vector<double> const &p) const;

    tetrahedron_id_t findTetByPoint(point3d const &p) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): MESH
    ////////////////////////////////////////////////////////////////////////

    /// Return the minimal coordinate of the rectangular bounding box.
    ///
    /// \return Minimal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMin() const;

    /// Return the maximal coordinate of the rectangular bounding box.
    ///
    /// \return Maximal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMax() const;

    /// Return the total volume of the mesh.
    ///
    /// \return Volume of the mesh.
    double getMeshVolume() const;

    /// Return the triangles which form the surface boundary of the mesh.
    /// \return Vector of the triangle boundary.
    // Weiliang 2010.02.02
    std::vector<index_t> getSurfTris() const;

    /// Computes the percentage of intersection of a segment with the mesh tets.
    ///
    /// \param p_start Beginning of segment.
    /// \param p_end Ending of segment.
    /// \param tet_start Index of tetrahedron containing \p p_start.
    /// \param sampling Number of point to test for monte-carlo method.
    /// \return A vector where each position contains pairs <tet, intersection ratio>.
    std::vector<std::pair<tetrahedron_id_t, double>>
    intersectMontecarlo(const point3d &p_start, const point3d &p_end, const tetrahedron_id_t &tet_start= UNKNOWN_TET, unsigned int sampling = 100)
    const;

    /// Computes the percentage of intersection of a segment with the mesh tets
    ///
    /// \return A vector where each position contains pairs <tet, intersection ratio>
    std::vector<std::pair<tetrahedron_id_t, double>> intersectDeterministic(const point3d &p_start, const point3d &p_end, const tetrahedron_id_t &tet_start= UNKNOWN_TET)
    const;

    // public alias type for segment intersections
    using intersection_list_t = std::vector<std::pair<index_t, double>>;

    /// Computes the percentage of intersection of a line of segments with the mesh tets
    ///
    /// \return A vector of vectors (for each segment) containing pairs <tet, intersection ratio>
    std::vector<intersection_list_t>
    intersect(const double *points, int n_points, int sampling = -1) const;


    ////////////////////////////////////////////////////////////////////////
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////

    /// Get barycenters of a list of tetrahedrons
    std::vector<double> getBatchTetBarycenters(std::vector<tetrahedron_id_t> const & tets) const;

    /// Get barycenters of a list of tetrahedrons
    void getBatchTetBarycentersNP(const tetrahedron_id_t* indices, int input_size, double* centers, int output_size) const;

    /// Get barycenters of a list of triangles
    std::vector<double> getBatchTriBarycenters(std::vector<triangle_id_t> const & tris) const;

    /// Get barycenters of a list of triangles
    void getBatchTriBarycentersNP(const triangle_id_t* indices, int input_size, double* centers, int output_size) const;

    /// Get coordinates of a list of vertices
    std::vector<double> getBatchVertices(std::vector<index_t> const & verts) const;

    /// Get coordinates of a list of vertices
    void getBatchVerticesNP(const index_t* indices, int input_size, double* coordinates, int output_size) const;

    /// Get vertex indices of a list of triangles
    std::vector<index_t> getBatchTris(std::vector<index_t> const & tris) const;

    /// Get vertex indices of a list of triangles
    void getBatchTrisNP(const index_t* t_indices, int input_size, index_t* v_indices, int output_size) const;

    /// Get vertex indices of a list of tetrahedrons
    std::vector<index_t> getBatchTets(std::vector<index_t> const & tets) const;

    /// Get vertex indices of a list of tetrahedrons
    void getBatchTetsNP(const index_t* t_indices, int input_size, index_t* v_indices, int output_size) const;

    /// Return the size of a set with unique vertex indices of a list of triangles
    /// preparation function for future numpy data access
    uint getTriVerticesSetSizeNP(const index_t* t_indices, int input_size) const;

    /// Return the size of a set with unique vertex indices of a list of tetrahedrons
    /// preparation function for future numpy data access
    uint getTetVerticesSetSizeNP(const index_t* t_indices, int input_size) const;

    /// Get the set with unique vertex indices of a list of triangles, write into given 1D array
    void getTriVerticesMappingSetNP(const index_t* t_indices, int input_size, index_t* t_vertices, int t_vertices_size, index_t* v_set, int v_set_size) const;

    /// Get the set with unique vertex indices of a list of tetrahedrons, write into given 1D array
    void getTetVerticesMappingSetNP(const index_t* t_indices, int input_size, index_t* t_vertices, int t_vertices_size, index_t* v_set, int v_set_size) const;

    /// Generate npnts random points inside tetrahedron t_idx and write the coordinates to coords
    void genPointsInTet(tetrahedron_id_t tidx, unsigned npnts, double *coords, unsigned int coord_size) const;

    /// Generate npnts random points on triangle t_idx and write the coordinates to coords
    void genPointsInTri(triangle_id_t tidx, unsigned npnts, double *coords, unsigned int coord_size) const;

    /// Generate the random points required in point_counts for tets in indices and store them in coords
    void genTetVisualPointsNP(const index_t *indices,
                              unsigned int index_size,
                              const unsigned int *point_counts,
                              unsigned int count_size,
                              double *coords,
                              unsigned int coord_size) const;

    /// Generate the random points required in point_counts for tris in indices and store them in coords
    void genTriVisualPointsNP(const index_t *indices,
                              unsigned int index_size,
                              const unsigned int *point_counts,
                              unsigned int count_size,
                              double *coords,
                              unsigned int coord_size) const;

    /// get the volumes of a list of tetrahedrons
    void getBatchTetVolsNP(const index_t* indices, int index_size, double* volumes, int volume_size) const;

    /// get the areas of a list of triangles
    void getBatchTriAreasNP(const index_t* indices, int index_size, double* areas, int area_size) const;

    /// reduce the number of points required to be generated in a list of tets based on maximum point density
    void reduceBatchTetPointCountsNP(const index_t* indices,
                                     unsigned int index_size,
                                     unsigned int *point_counts,
                                     unsigned int count_size,
                                     double max_density);

    /// reduce the number of points required to be generated in a list of tris based on maximum point density
    void reduceBatchTriPointCountsNP(const index_t* indices,
                                     unsigned int index_size,
                                     unsigned int *point_counts,
                                     unsigned int count_size,
                                     double max_density);

    ////////////////////////////////////////////////////////////////////////
    // ROI (Region of Interest) Data
    ////////////////////////////////////////////////////////////////////////

    /// Add a ROI data
    void addROI(std::string const &id, ElementType type, const ROISet::set_data_type& indices);

    /// Remove a ROI data
    [[gnu::deprecated("ROI type will be part of its identifier in future version")]]
    void removeROI(std::string const &id);
    void removeROI(std::string const &id, ElementType type);

    /// Replace a ROI data with a new set with the same name
    void replaceROI(std::string const &id, ElementType type, const ROISet::set_data_type& indices);

    /// Return the type of a ROI data
    [[gnu::deprecated("ROI type will be part of its identifier in future version")]]
    ElementType getROIType(std::string const &id) const;

    /// Return the data of a ROI
    [[gnu::deprecated("type should be specified when retrieving a ROI")]]
    ROISet::vector_data_type getROIData(std::string const &id) const;
    ROISet::vector_data_type const & getROIData(std::string const &id, ElementType type) const;

    /// Return the data size of a ROI
    [[gnu::deprecated("usage lead to bad pattern")]]
    uint getROIDataSize(std::string const &id) const;

    /// get the total number of ROI recorded
    /// get the total number of ROI recorded
    uint getNROIs() const;

    /// Return a ROI
    ROISet getROI(std::string const &id) const;

    /// get all ROI names
    std::vector<std::string> getAllROINames() const;

    /// check if a ROI enquire is valid
    bool checkROI(std::string const &id, ElementType type, uint count = 0, bool warning = true) const;

  ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////

    /// get barycenters of a list of tetrahedrons
    std::vector<double> getROITetBarycenters(std::string const &ROI_id) const;

    /// get barycenters of a list of tetrahedrons
    void getROITetBarycentersNP(std::string const &ROI_id, double* centers, int output_size) const;

    /// get barycenters of a list of triangles
    std::vector<double> getROITriBarycenters(std::string const &ROI_id) const;

    /// get barycenters of a list of triangles
    void getROITriBarycentersNP(std::string const &ROI_id, double* centers, int output_size) const;

    /// get coordinates of a list of vertices
    std::vector<double> getROIVertices(std::string const & ROI_id) const;

    /// get coordinates of a list of vertices
    void getROIVerticesNP(std::string const &ROI_id, double* coordinates, int output_size) const;

    /// get vertex indices of a list of triangles
    std::vector<index_t> getROITris(std::string const &ROI_id) const;

    /// get vertex indices of a list of triangles
    void getROITrisNP(std::string const &ROI_id, index_t* v_indices, int output_size) const;

    /// get vertex indices of a list of tetrahedrons
    std::vector<index_t> getROITets(std::string const &ROI_id) const;

    /// get vertex indices of a list of tetrahedrons
    void getROITetsNP(std::string const &ROI_id, index_t* v_indices, int output_size) const;

    /// return the size of a set with unique vertex indices of a list of triangles
    /// preparation function for future numpy data access
    uint getROITriVerticesSetSizeNP(std::string const &ROI_id) const;

    /// return the size of a set with unique vertex indices of a list of tetrahedrons
    /// preparation function for future numpy data access
    uint getROITetVerticesSetSizeNP(std::string const &ROI_id) const;

    /// Get the set with unique vertex indices of a list of triangles, write into given 1D array
    void getROITriVerticesMappingSetNP(std::string const &ROI_id, index_t* t_vertices, int t_vertices_size, index_t* v_set, int v_set_size) const;

    /// Get the set with unique vertex indices of a list of tetrahedrons, write into given 1D array
    void getROITetVerticesMappingSetNP(std::string const &ROI_id, index_t* t_vertices, int t_vertices_size, index_t* v_set, int v_set_size) const;

    /// Generate the random points required in point_counts for tets in indices and store them in coords
    void genROITetVisualPointsNP(std::string const &ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const;

    /// Generate the random points required in point_counts for tris in indices and store them in coords
    void genROITriVisualPointsNP(std::string const &ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const;

    /// get the volumes of a list of tetrahedrons
    void getROITetVolsNP(std::string const &ROI_id, double* volumes, int volume_size) const;

    /// get the areas of a list of triangles
    void getROITriAreasNP(std::string const &ROI_id, double* areas, int area_size) const;

    /// get the total volume of a tetrahedral roi
    double getROIVol(std::string const &ROI_id) const;

    /// get the total area of a triangular roi
    double getROIArea(std::string const &ROI_id) const;

    /// reduce the number of points required to be generated in a list of tets based on maximum point density
    void reduceROITetPointCountsNP(std::string const &ROI_id, unsigned int* point_counts, int count_size, double max_density);

    /// reduce the number of points required to be generated in a list of tris based on maximum point density
    void reduceROITriPointCountsNP(std::string const &ROI_id, unsigned int* point_counts, int count_size, double max_density);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (C++ INTERNAL)
    ////////////////////////////////////////////////////////////////////////

    /// Return a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    inline const point3d &_getVertex(vertex_id_t vidx) const noexcept { return pVerts[vidx.get()]; }

    /// Return the bar with index bidx
    ///
    /// \param bidx Index of the bar.
    /// \return Indices of the two vertices that form the bar.
    inline const vertex_id_t *_getBar(bar_id_t bidx) const noexcept { return pBars[bidx.get()].data(); }

    ///Return the triangle neighbors of a bar with index bidx.
    ///
    /// \param bidx Index of the bar.
    /// \return Array of the triangle neighbors.
    inline const triangle_id_t* _getBarTriNeighb(bar_id_t bidx) const noexcept { return pBar_tri_neighbours[bidx.get()].data(); }

    /// Return a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return List of the vertices form the triangle.
    inline const vertex_id_t* _getTri(triangle_id_t tidx) const noexcept { return pTris[tidx.get()].data(); }

    /// Return a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return List of the vertices form the tetrahedron.
    inline const vertex_id_t* _getTet(tetrahedron_id_t tidx) const noexcept { return pTets[tidx.get()].data(); }

    ///Return the tetrahedron neighbors of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Array of the tetrahedron neighbors.
    /// \TODO TCL return const tri_tets& instead
    inline const tetrahedron_id_t * _getTriTetNeighb(triangle_id_t tidx) const noexcept { return pTri_tet_neighbours[tidx.get()].data(); }

    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the triangle neighbors.
    /// \TODO TCL use const tet_tris& instead
    inline const triangle_id_t * _getTetTriNeighb(tetrahedron_id_t tidx) const noexcept { return pTet_tri_neighbours[tidx.get()].data(); }

    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the tetrahedron neighbors.
    /// \TODO TCL return const tet_tets& instead
    inline const tetrahedron_id_t * _getTetTetNeighb(tetrahedron_id_t tidx) const noexcept { return pTet_tet_neighbours[tidx.get()].data(); }

    /// Return the bars of a triangle
    ///
    /// \param tidx Index of the triangle
    /// \return Bars of the triangle
    inline const tri_bars& _getTriBars(triangle_id_t tidx) const noexcept { return pTri_bars[tidx.get()]; }

    /// Return the normal to the triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Point representing unit normal to triangle.
    inline const point3d &_getTriNorm(triangle_id_t tidx) const noexcept { return pTri_norms[tidx.get()]; }

    /// Return the barycenter of the triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Barycenter of triangle.
    inline const point3d &_getTriBarycenter(triangle_id_t tidx) const noexcept { return pTri_barycs[tidx.get()]; }

    /// Return the barycenter of the tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Barycenter of tetrahedron.
    inline const point3d &_getTetBarycenter(tetrahedron_id_t tidx) const noexcept { return pTet_barycenters[tidx.get()]; }

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
    void _handleMembAdd(Memb * memb);

    /// Delete a membrane.
    ///
    /// \param patch Pointer to the membrane.
    void _handleMembDel(Memb * memb);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the membranes in the tetmesh container.
    ///
    /// \return Number of membranes.
    inline uint _countMembs() const noexcept
    { return pMembs.size(); }

    /// Return a membrane with index gidx.
    ///
    /// \param gidx Index of the membrane.
    /// \return Pointer to the membrane.
    Memb * _getMemb(uint gidx) const;

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
    /// \param diffb Pointer to the diffusion boundary.
    void _handleDiffBoundaryAdd(DiffBoundary * diffb);

    /// Delete a diffusion boundary.
    ///
    /// \param diffb Pointer to the diffusion boundary.
    void _handleDiffBoundaryDel(DiffBoundary * diffb);

    /// Count the diffusion boundaries in the tetmesh container.
    ///
    /// \return Number of diffusion boundaries.
    inline uint _countDiffBoundaries() const noexcept
    { return pDiffBoundaries.size(); }

    /// Return a diffusion boundary with index gidx.
    ///
    /// \param gidx Index of the diffusion boundary.
    /// \return Pointer to the diffusion boundary.
    DiffBoundary * _getDiffBoundary(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////

    /// Check if a surface diffusion boundary id is occupied.
    ///
    /// \param id ID of the surface diffusion boundary.
    void _checkSDiffBoundaryID(std::string const & id) const;

    /// Change the id of a surface diffusion boundary.
    ///
    /// \param o Old id of the surface diffusion boundary.
    /// \param n New id of the surface diffusion boundary.
    void _handleSDiffBoundaryIDChange(std::string const & o, std::string const & n);

    /// Add a surface diffusion boundary.
    ///
    /// \param sdiffb Pointer to the surface diffusion boundary.
    void _handleSDiffBoundaryAdd(SDiffBoundary * sdiffb);

    /// Delete a diffusion boundary.
    ///
    /// \param sdiffb Pointer to the surface diffusion boundary.
    void _handleSDiffBoundaryDel(SDiffBoundary * sdiffb);

    /// Count the surface diffusion boundaries in the tetmesh container.
    ///
    /// \return Number of surface diffusion boundaries.
    inline uint _countSDiffBoundaries() const noexcept
    { return pSDiffBoundaries.size(); }

    /// Return a surface diffusion boundary with index gidx.
    ///
    /// \param gidx Index of the surface diffusion boundary.
    /// \return Pointer to the surface diffusion boundary.
    SDiffBoundary * _getSDiffBoundary(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////

    static const tri_tets UNKNOWN_TRI_NEIGHBORS;
    static const tet_tets UNKNOWN_TET_NEIGHBORS;

private:
    /// Build pBars, pBarsN, pTri_bars from pTris.
    void buildBarData();

    ///////////////////////// DATA: VERTICES ///////////////////////////////
    ///
    /// The total number of vertices in the mesh
    index_t                                pVertsN{0};
    /// The vertices by x,y,z coordinates
    std::vector<point3d>                pVerts;

    /////////////////////////// DATA: BARS /////////////////////////////////
    ///
    /// The total number of 1D 'bars' in the mesh
    index_t                                pBarsN{0};
    /// The bars by the two vertices index
    std::vector<bar_verts>              pBars;

    /// The surface diffusion boundary a bar belongs to
    std::vector<SDiffBoundary *>  pBar_sdiffboundaries;

    /// The 2 triangle neighbours of each bar (by index) in the context of
    /// surface diffusion boundaries
    std::vector<bar_tris>				pBar_tri_neighbours;

    ///////////////////////// DATA: TRIANGLES //////////////////////////////
    ///
    /// The total number of triangles in the mesh
    index_t                         pTrisN{0};
    /// The triangles by vertices index
    std::vector<tri_verts>              pTris;
    // The bars of the triangle
    std::vector<tri_bars>               pTri_bars;
    /// The areas of the triangles
    std::vector<double>                 pTri_areas;
    /// The triangle barycenters
    std::vector<point3d>                pTri_barycs;
    /// The triangle normals
    std::vector<point3d>                pTri_norms;
    /// The patch a triangle belongs to
    std::vector<TmPatch *> pTri_patches;

    /// The diffusion boundary a triangle belongs to
    std::vector<DiffBoundary *> pTri_diffboundaries;

    /// The tetrahedron neighbours of each triangle (by index)
    std::vector<tri_tets>               pTri_tet_neighbours;

    ///////////////////////// DATA: TETRAHEDRA /////////////////////////////
    ///
    /// The total number of tetrahedron in the mesh
    index_t                         pTetsN{0};
    /// The tetrahedron by vertices index
    std::vector<tet_verts>              pTets;
    /// The volume of the tetrahedron
    std::vector<double>                 pTet_vols;
    /// The barycenters of the tetrahedra
    std::vector<point3d>                pTet_barycenters;
    /// The compartment a tetrahedron belongs to
    std::vector<TmComp  *> pTet_comps;
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
    MembPMap       pMembs;
    std::map<std::string, DiffBoundary *> pDiffBoundaries;
    std::map<std::string, SDiffBoundary *> pSDiffBoundaries;

    ////////////////////////// ROI Dataset /////////////////////////////////
    //std::map<std::string, ROISet>                       mROI;

public:
    RegionOfInterest rois;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace tetmesh
} // namespace steps

#endif

// STEPS_TETMESH_TETMESH_HPP
// END

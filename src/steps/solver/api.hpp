/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

#include <map>
#include <string>

#include "fwd.hpp"
#include "geom/fwd.hpp"
#include "model/fwd.hpp"
#include "rng/fwd.hpp"
#include "util/vocabulary.hpp"


namespace steps::model {
// Forward declarations
struct SubunitStateFilter;
}  // namespace steps::model

////////////////////////////////////////////////////////////////////////////////

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////
/// API class for a solver.
///
/// The API class provides all APIs for each solver.
///
/// \warning Not every API in this class is implemented in the solver.
///          If a API is not implemented in a solver,
///          STEPS will throw an error message to the user.
/// \warning Methods start with underscore are not exposed to Python.
///////////////////////////////////////////////////////////////////////////////
class API {
  public:
    // Constants for describing E-Field solver choices
    enum EF_solver {
        EF_NONE = 0,     // must be zero for API compatibility
        EF_DEFAULT = 1,  // must be one for API compatibility
        EF_DV_BDSYS,
        EF_DV_PETSC,
    };

    /// Constructor
    ///
    /// \param m Reference to the model.
    /// \param g Reference to the geometry container.
    /// \param r Reference to the random number generator.
    API(model::Model& m, wm::Geom& g, const rng::RNGptr& r);

    /// Destructor
    ///
    virtual ~API();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Return the solver's name.
    virtual std::string getSolverName() const = 0;
    /// Return the solver's description.
    virtual std::string getSolverDesc() const = 0;
    /// Return the solver's author.
    virtual std::string getSolverAuthors() const = 0;
    /// Return the solver's email.
    virtual std::string getSolverEmail() const = 0;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    /// checkpoint simulator state to a file
    virtual void checkpoint(std::string const& file_name) = 0;

    /// restore simulator state from a file
    virtual void restore(std::string const& file_name) = 0;

    /// Reset the solver.
    virtual void reset() = 0;

    /// Advance the simulation until endtime (given in seconds) is reached.
    /// The endtime must be larger or equal to the current simulation time.
    ///
    /// \param endtime Time to end the solver.
    virtual void run(double endtime) = 0;

    /// Advance the solver a given amount of time.
    ///
    /// \param adv Time to advance the solver (in seconds)
    virtual void advance(double adv);

    /// Advance the simulation for one 'step'. In stochastic solvers this is
    /// one 'realization' of the Gillespie SSA (one reaction 'event'). In
    /// numerical solvers (currently Wmrk4) this is one time-step, with the
    /// stepsize defined with the setDT method.
    virtual void step();

    /// Set DT of the numerical solver.
    ///
    /// \param dt Dt.
    virtual void setRk4DT(double dt);

    /// Set DT of the solver. Included for backwards compatibility but
    /// replaced by setRk4DT
    ///
    /// \param dt Dt.
    virtual void setDT(double dt);

    /// Set the stepsize for membrane potential solver (default 1e-6s).
    /// This is the time for each voltage calculation step. The SSA will
    /// run until passing this stepsize, so in fact each membrane potential
    /// time step will vary slightly around the dt so as to be aligned with the
    /// SSA.
    ///
    /// \param dt EField DT (in seconds)
    virtual void setEfieldDT(double efdt);

    virtual void setNSteps(uint nsteps);

    virtual void setTime(double time);

    /// Set the simulation temperature. Currently, this will only
    /// influence the GHK flux rate, so will only influence simulations
    /// including membrane potential calculation.
    ///
    /// \param temp EField temperature (in Kelvin)
    virtual void setTemp(double temp);

    // Set the default vesicle dt. The actual dt used in the simulation may be lower
    // depending on the presence of Rafts
    virtual void setVesicleDT(double dt);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////
    /// Returns the current simulation time in seconds.
    virtual double getTime() const = 0;

    /// Return the DT of the numerical solver
    virtual double getRk4DT() const;

    /// Return the DT
    /// Replaced by getRk4DT, but included for backwards compatability
    virtual double getDT() const;

    /// Return the stepsize for membrane potential solver (in seconds)
    virtual double getEfieldDT() const;

    /// Return the simulation temperature (in Kelvin)
    virtual double getTemp() const;

    /// Returns the total propensity of the current simulation state
    /// (the total propensity multiplied by an infinitesimally small
    /// time dt gives the probability that a reaction will occur in that dt).
    /// For Tetexact this includes the propensity from the extension of the SSA
    /// for diffusive flux between tetrahedral elements in the mesh.
    virtual double getA0() const;

    /// Return the number of 'realizations' of the SSA, the number of reaction
    /// (and diffusion) events in stochastic solvers.
    virtual uint getNSteps() const;

    // The actual dt used in the simulation, which may be lower than the default dt
    // depending on the presence of Rafts
    virtual double getVesicleDT() const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Returns the volume of compartment c (in m^3).
    ///
    /// \param c Name of the compartment.
    double getCompVol(std::string const& c) const;

    /// Sets the volume of compartment c.
    ///
    /// NOTE: this method may throw an exception if this does not make sense
    /// for a solver (e.g. a tetrahedral mesh-based solver).
    /// \param c Name of the compartment.
    /// \param vol Volume of the compartment.
    void setCompVol(std::string const& c, double vol);

    /// Returns the number of molecules of species s in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    double getCompSpecCount(std::string const& c, std::string const& s) const;

    /// Sets the number of molecules of species s in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all tetrahedrons in the compartment (i.e. a uniform distribution).
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    /// \param n Number of molecules of the species.
    void setCompSpecCount(std::string const& c, std::string const& s, double n);

    /// Returns the amount (in mols) of species s in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is computed as
    /// the sum of the amounts in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    double getCompSpecAmount(std::string const& c, std::string const& s) const;

    /// Set the amount (in mols) of species s in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all tetrahedrons in the compartment (i.e. a uniform distribution).
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    /// \param a Amount of the species.
    void setCompSpecAmount(std::string const& c, std::string const& s, double a);

    /// Returns the concentration (in molar units) of species s in compartment
    /// c.
    ///
    /// NOTE: in a mesh-based simulation, the overall concentration in a
    /// compartment is computed by taking the volume-weighted sum of the
    /// concentration in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    double getCompSpecConc(std::string const& c, std::string const& s) const;

    /// Sets the concentration (in molar units) of species s in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, this method changes the
    /// concentration to the same value in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    /// \param conc Concentration of the species.
    void setCompSpecConc(std::string const& c, std::string const& s, double conc);

    /// Returns whether the concentration of species s in compartment c
    /// remains constant over time (unless changed explicitly).
    ///
    /// NOTE: in a mesh-based simulation, this method will only return true
    /// only if the species has been clamped in all tetrahedrons of the
    /// compartment. \param c Name of the compartment. \param s Name of the
    /// species.
    bool getCompSpecClamped(std::string const& c, std::string const& s) const;

    /// Turns clamping of species s in compartment c on or off.
    ///
    /// NOTE: in a mesh based simulation, this method turns clamping on/off
    /// in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the species.
    /// \param b Flag to trun clamping of species on / off.
    void setCompSpecClamped(std::string const& c, std::string const& s, bool b);

    /// Returns the number of Complexes s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    double getCompComplexCount(std::string const& c,
                               std::string const& s,
                               const std::vector<std::vector<model::SubunitStateFilter>>& f) const;

    /// Sets the number of Complexes s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all voxels in the compartment (i.e. a uniform distribution).
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param i Init state for the complex.
    /// \param n Number of molecules of the species.
    void setCompComplexCount(std::string const& c,
                             std::string const& s,
                             const std::vector<std::vector<model::SubunitStateFilter>>& i,
                             double n);

    /// Returns the amount (in mols) of complex s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is computed as
    /// the sum of the amounts in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    double getCompComplexAmount(std::string const& c,
                                std::string const& s,
                                const std::vector<std::vector<model::SubunitStateFilter>>& f) const;

    /// Set the amount (in mols) of complex s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all voxels in the compartment (i.e. a uniform distribution).
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param i Init state for the complex.
    /// \param a Amount of the complex.
    void setCompComplexAmount(std::string const& c,
                              std::string const& s,
                              const std::vector<std::vector<model::SubunitStateFilter>>& i,
                              double a);

    /// Returns the concentration (in molar units) of complex s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the overall concentration in a
    /// compartment is computed by taking the volume-weighted sum of the
    /// concentration in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    double getCompComplexConc(std::string const& c,
                              std::string const& s,
                              const std::vector<std::vector<model::SubunitStateFilter>>& f) const;

    /// Sets the concentration (in molar units) of complex s matching filter f in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, this method changes the
    /// concentration to the same value in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param i Init state for the complex.
    /// \param conc Concentration of the complex.
    void setCompComplexConc(std::string const& c,
                            std::string const& s,
                            const std::vector<std::vector<model::SubunitStateFilter>>& i,
                            double conc);

    /// Returns the number of SubUnits in state m for Complexes s matching filter f in compartment
    /// c.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    //  \param m index of the subunitstate
    double getCompComplexSUSCount(std::string const& c,
                                  std::string const& s,
                                  const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                  uint m) const;

    /// Returns the concentration of SubUnits in state m for Complexes s matching filter f in
    /// compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the overall concentration in a
    /// compartment is computed by taking the volume-weighted sum of the
    /// concentration in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    //  \param m index of the subunitstate
    double getCompComplexSUSConc(std::string const& c,
                                 std::string const& s,
                                 const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                 uint m) const;

    /// Returns the amount (in mols) of SubUnits in state m for Complexes s matching filter f in
    /// compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is computed as
    /// the sum of the amounts in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    //  \param m index of the subunitstate
    double getCompComplexSUSAmount(std::string const& c,
                                   std::string const& s,
                                   const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                   uint m) const;

    // Returns the macroscopic reaction constant of reaction r in
    // compartment c.
    // Note: in a mesh-based simulation, the value is computed as the
    // volume-weighted sum of the reaction constants in all voxels of the
    // compartment.
    double getCompReacK(std::string const& c, std::string const& r) const;

    /// Sets the macroscopic reaction constant of reaction r in compartment c
    /// (units vary according to the order of the reaction).
    ///
    /// NOTE: in a mesh-based simulation, this method changes the reaction
    /// constant equally in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of te reaction.
    /// \param kf Reaction constant.
    void setCompReacK(std::string const& c, std::string const& r, double kf);

    /// Returns whether reaction r in compartment c is active or not
    ///
    /// NOTE: in a mesh-based simulation, this method returns false only when
    /// the reaction has been inactivated in all tetrahedrons.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    bool getCompReacActive(std::string const& c, std::string const& r) const;

    /// Activate or inactivate reaction r in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, activation/inactivation of a reaction
    /// turns it on/off in all tetrahedrons.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    /// \param a Flag to activate or deactivate the reaction.
    void setCompReacActive(std::string const& c, std::string const& r, bool a);

    /// Returns the diffusion constant of diffusion rule d in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param d Name of the diffusion.
    double getCompDiffD(std::string const& c, std::string const& d) const;

    /// Set the diffusion constant of diffusion rule d in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param d Name of the diffusion.
    /// \param dcst Rate constant of the diffusion.
    void setCompDiffD(std::string const& c, std::string const& d, double dcst);

    /// Returns whether diffusion rule d in compartment c is active or not.
    ///
    /// \param c Name of the compartment.
    /// \param d Name of the diffusion.
    bool getCompDiffActive(std::string const& c, std::string const& d) const;

    /// Activate or deactivate diffusion rule d in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param d Name of the diffusion.
    /// \param act Flag to activate or deactivate the diffusion.
    void setCompDiffActive(std::string const& c, std::string const& d, bool act);

    ////////////////////////////////////////////////////////////////////////

    /// Returns c_mu, the mesoscopic reaction constant of reaction r in
    /// compartment c.
    ///
    /// NOTE: in a mesh-based simulation, the mesoscopic reaction constant is
    /// computed as the sum of the mesoscopic constants in all tetrahedrons of
    /// the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    double getCompReacC(std::string const& c, std::string const& r) const;

    /// Returns h_mu, the distinct number of ways in which reaction r can
    /// occur in compartment c, by computing the product of its reactants.
    ///
    /// NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
    /// over all tetrahedrons of the compartment. This can become a very large
    /// value.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    double getCompReacH(std::string const& c, std::string const& r) const;

    /// Returns the propensity, a_mu, of reaction r in compartment c.
    /// The propensity value gives the probability per unit time that this
    /// reaction will occur in the current state.
    ///
    /// NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
    /// a_mu in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    double getCompReacA(std::string const& c, std::string const& r) const;

    /// Returns the extent of reaction r in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, returns the sum of the reaction
    /// extents in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    unsigned long long getCompReacExtent(std::string const& c, std::string const& r) const;

    /// Resets the extent of reaction r in compartment c to zero.
    ///
    /// NOTE: in a mesh-based simulation, resets the extents of the reaction
    /// in all tetrahedrons of the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    void resetCompReacExtent(std::string const& c, std::string const& r);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      COMPARTMENT, VESICLE-RELATED
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of vesicles v in compartment c
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    uint getCompVesicleCount(std::string const& c, std::string const& v) const;

    /// Sets the number of vesicles v in compartment c
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param n Number of vesicles
    void setCompVesicleCount(std::string const& c, std::string const& v, uint n);

    /// Adds one vesicle v to compartment c
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    vesicle_individual_id addCompVesicle(std::string const& c, std::string const& v);

    /// Deletes individual vesicle of type v with unique index ves_unique_index
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    void deleteSingleVesicle(std::string const& v, vesicle_individual_id ves_unique_index);

    /// Returns the count of 'link' species ls on vesicle of type v with unique index
    /// vesicle_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param ls Name of the link species
    uint getSingleVesicleSurfaceLinkSpecCount(std::string const& v,
                                              vesicle_individual_id ves_unique_index,
                                              std::string const& ls) const;

    /// Returns the indices of 'link' species ls on vesicle of type v with unique index
    /// vesicle_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param ls Name of the link species
    std::vector<linkspec_individual_id> getSingleVesicleSurfaceLinkSpecIndices(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& ls) const;

    /// Returns the indices of point species s on vesicle of type v with unique index
    /// vesicle_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the link species
    std::vector<pointspec_individual_id> getSingleVesicleSurfaceSpecIndices(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& s) const;


    /// Returns all the unique indices of all vesicles currently in the simulation
    ///
    std::vector<vesicle_individual_id> getAllVesicleIndices() const;

    /// Returns all the unique indices of all vesicles currently on a path
    ///
    std::vector<vesicle_individual_id> getAllVesicleIndicesOnPath() const;

    /// Returns all the unique indices of vesicles v currently in compartment c
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    std::vector<vesicle_individual_id> getCompVesicleIndices(std::string const& c,
                                                             std::string const& v) const;

    /// Get the compartment of vesicle of type v and unique index ves_unique_index
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    std::string getSingleVesicleCompartment(std::string const& v,
                                            vesicle_individual_id ves_unique_index) const;

    /// Returns the position of vesicle of type v with unique index
    /// ves_unique_index (in cartesian coordinates)
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    std::vector<double> getSingleVesiclePos(std::string const& v,
                                            vesicle_individual_id ves_unique_index) const;

    /// Set the position of vesicle of type v with unique index
    /// ves_unique_index to pos (cartesian coordinates)
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param pos Position of the vesicle in cartesian corrdinates
    void setSingleVesiclePos(std::string const& v,
                             vesicle_individual_id ves_unique_index,
                             const std::vector<double>& pos,
                             bool force = false);

    /// Returns the surface count of species s on vesicles v in compartment c.
    /// Return is a map, vesicle_unique_index : count of s
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param s Name of the species
    std::map<vesicle_individual_id, uint> getCompVesicleSurfaceSpecCountDict(
        std::string const& c,
        std::string const& v,
        std::string const& s) const;

    /// Returns the summed surface count of species s on vesicles v in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param s Name of the species
    uint getCompVesicleSurfaceSpecCount(std::string const& c,
                                        std::string const& v,
                                        std::string const& s) const;

    /// Returns the summed inner count of species s on vesicles v in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param s Name of the species
    uint getCompVesicleInnerSpecCount(std::string const& c,
                                      std::string const& v,
                                      std::string const& s) const;

    /// Get the surface count of species s on vesicle of type v and unique
    /// index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    uint getSingleVesicleSurfaceSpecCount(std::string const& v,
                                          vesicle_individual_id ves_unique_index,
                                          std::string const& s) const;

    /// Get the inner count of species s on vesicle of type v and unique
    /// index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    uint getSingleVesicleInnerSpecCount(std::string const& v,
                                        vesicle_individual_id ves_unique_index,
                                        std::string const& s) const;

    /// Set the surface count of species s on vesicle of type v and unique
    /// index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    /// \param count The number of the surface species
    void setSingleVesicleSurfaceSpecCount(std::string const& v,
                                          vesicle_individual_id ves_unique_index,
                                          std::string const& s,
                                          uint count);

    /// Get the cartesian coordinates of species s on vesicle of type v and
    /// unique index ves_unique_index. Position is absolute,
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    std::vector<std::vector<double>> getSingleVesicleSurfaceSpecPos(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& s);

    /// Set the count of the inner species, that is inside the vesicle volume,
    /// on vesicle of type v and unique index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    /// \param count The number of the inner species
    void setSingleVesicleInnerSpecCount(std::string const& v,
                                        vesicle_individual_id ves_unique_index,
                                        std::string const& s,
                                        uint count);

    /// Get the spherical coordinates of species s on vesicle of type v and
    /// unique index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    std::vector<std::vector<double>> getSingleVesicleSurfaceSpecPosSpherical(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& s) const;

    /// Set the spherical coordinates of species s on vesicle of type v and
    /// unique index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param s Name of the species
    /// \param pos_spherical Positions of the molecule in spherical coordinates
    /// (relative to the vesicle)
    void setSingleVesicleSurfaceSpecPosSpherical(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& s,
        const std::vector<std::vector<double>>& pos_spherical);

    /// Get the spherical coordinates of species s with unique id ps_unique_id.
    ///
    /// \param s Name of the species
    /// \param ps_unique_id Unique index of the individual point species
    std::vector<double> getSingleSpecPosSpherical(std::string const& s,
                                                  pointspec_individual_id ps_unique_id) const;

    /// Returns the count of 'link' species ls on vesicles v in compartment c.
    /// Return is a map, vesicle_unique_index : count of ls
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param ls Name of the link species
    std::map<vesicle_individual_id, uint> getCompVesicleSurfaceLinkSpecCountDict(
        std::string const& c,
        std::string const& v,
        std::string const& ls) const;

    /// Returns the summed count of 'link' species ls on vesicles v in compartment c.
    ///
    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param ls Name of the link species
    uint getCompVesicleSurfaceLinkSpecCount(std::string const& c,
                                            std::string const& v,
                                            std::string const& ls) const;

    /// Returns the positions of 'link' species ls on vesicle of type v and
    /// unique index ves_unique_index.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    /// \param ls Name of the link species
    std::vector<std::vector<double>> getSingleVesicleSurfaceLinkSpecPos(
        std::string const& v,
        vesicle_individual_id ves_unique_index,
        std::string const& ls);

    /// Returns the position of the link species with unique index ls_unique_id.
    ///
    /// \param ls_unique_id Unique index of the individual link species
    std::vector<double> getSingleLinkSpecPos(linkspec_individual_id ls_unique_id) const;

    /// Returns the unique index of the link species linked to the link species
    /// with unique index ls_unique_id.
    ///
    /// \param ls_unique_id Unique index of the individual link species
    linkspec_individual_id getSingleLinkSpecLinkedTo(linkspec_individual_id ls_unique_id) const;

    /// Returns the unique index of the vesicle that contains the link species
    /// with unique index ls_unique_id.
    ///
    /// \param ls_unique_id Unique index of the individual link species
    vesicle_individual_id getSingleLinkSpecVes(linkspec_individual_id ls_unique_id) const;

    /// Get the 'immobility' of vesicle of type v and unique index
    /// ves_unique_index. All non-zero numbers mean vesicle is
    /// immobile.
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    uint getSingleVesicleImmobility(std::string const& v,
                                    vesicle_individual_id ves_unique_index) const;

    /// Get the indexes of tetrahedrons that overlap vesicle of type v and
    /// unique index ves_unique_index
    ///
    /// \param v Name of the vesicle.
    /// \param ves_unique_index Unique index of the individual vesicle
    std::vector<tetrahedron_global_id> getSingleVesicleOverlapTets(
        std::string const& v,
        vesicle_individual_id ves_unique_index) const;

    /// Set the diffusion rate per tetrahedron of vesicles of type v. Vesicles
    /// will use this diffusion rate when vesicle centre is in this tet.
    ///
    /// \param tidx Tetrahrdron index
    /// \param v Name of the vesicle.
    /// \param dcst Diffusion coefficient
    void setTetVesicleDcst(tetrahedron_global_id tidx, std::string const& v, double dcst);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      VESICLE-RELATED, BUT NOT COMPARTMENT-SPECIFIC
    ////////////////////////////////////////////////////////////////////////

    ///  Set the diffusion rate of 'link' species ls on vesicles of type v.
    ///
    /// \param v Name of the vesicle.
    /// \param ls Name of the link species
    /// \param dcst Diffusion coefficient
    void setVesicleSurfaceLinkSpecSDiffD(std::string const& v, std::string const& ls, double dcst);

    /// \param c Name of the compartment.
    /// \param v Name of the vesicle.
    /// \param s Name of the species
    /// \param d Diffusion rate of the species in the vesicle surface
    // void setCompVesicleSpecDiffD(std::string const & c, std::string const &
    // v, std::string const & s, double d) const;

    /// Set the reaction rate of vesicle surface reaction. Not
    /// compartment-specific.
    ///
    /// \param vsr Name of the vesicle surface reaction
    /// \param kf Rate
    void setVesSReacK(std::string const& vsr, double kf);

    /// Get the reaction extent of vesicle surface reaction. Not
    /// compartment-specific.
    ///
    /// \param vsr Name of the vesicle surface reaction
    uint getVesSReacExtent(std::string const& vsr) const;

    /// Set rate of exocytosis events. Not compartment-specific.
    ///
    /// \param exo Name of the vesicle exocytosis
    /// \param kf Rate
    void setExocytosisK(std::string const& exo, double kf);

    /// Get the extent of vesicle exocytosis. Not compartment-specific.
    ///
    /// \param exo Name of the vesicle exocytosis
    uint getExocytosisExtent(std::string const& exo) const;

    /// Get the vesicle exocytosis events that happened since last call. Not compartment-specific.
    ///
    /// \param exo Name of the vesicle exocytosis
    std::vector<ExocytosisEvent> getExocytosisEvents(std::string const& exo);

    /// Get the extent of raft endocytosis. Not patch-specific.
    ///
    /// \param rendo Name of the raft endocytosis
    uint getRaftEndocytosisExtent(std::string const& rendo) const;

    /// Get the raft endocytosis events that happened since last call. Not patch-specific.
    ///
    /// \param rendo Name of the raft endocytosis
    std::vector<RaftEndocytosisEvent> getRaftEndocytosisEvents(std::string const& rendo);

    /// Set the rate of raft endocytosis. Not patch-specific.
    ///
    /// \param rendo Name of the raft endocytosis
    /// \param kcst Rate of raft endocytosis
    void setRaftEndocytosisK(std::string const& rendo, double kcst);

    /// Add a 'diffusion group' for vesicles of type v. Vesicles will diffuse
    /// freely amongst this group (if they border each other).
    ///
    /// \param v Name of the vesicle.
    /// \param comps List of compartment names.
    void addVesicleDiffusionGroup(std::string const& v, const std::vector<std::string>& comps);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      VESICLE PATHS
    ////////////////////////////////////////////////////////////////////////

    /// Create a path
    ///
    /// \param path Name of the path
    virtual void createPath(std::string const& path);

    /// Add a point to a path
    ///
    /// \param path Name of the path
    /// \param point_idx An index for the point, positive integer
    /// \param position Position of the point in cartesian coordinates
    virtual void addPathPoint(std::string const& path,
                              uint point_idx,
                              const std::vector<double>& position);

    /// Create branching in the path from a point
    ///
    /// \param path Name of the path
    /// \param sourcepoint_idx An index for the source point, positive integer
    /// \param dest_points_indxs Map of destination points with weight
    virtual void addPathBranch(std::string const& path,
                               uint sourcepoint_idx,
                               const std::map<uint, double>& destpoints_indxs);

    /// Return a structure describing all the paths added so far
    //
    //  Only rank 0 should call this method
    virtual std::map<std::string,
                     std::map<uint, std::pair<std::vector<double>, std::map<uint, double>>>>
    getAllPaths() const;

    /// Add a vesicle to this path. This means a vesicle of this type can
    /// interact with this path upon overlapping it
    ///
    /// \param path Name of the path
    /// \param ves Name of the vesicle
    /// \param speed Speed of the vesicle on this path in m/s
    /// \param spec_deps Optional species dependencies, in map of species names
    /// to number of the species required
    void addPathVesicle(std::string const& path,
                        std::string const& ves,
                        double speed,
                        const std::map<std::string, uint>& spec_deps,
                        const std::vector<double>& stoch_stepsize);

    ////////////////////////////////////////////////////////////////////////

    /// Returns the extent of complex reaction r in compartment c.
    ///
    /// NOTE: in a mesh-based simulation, returns the sum of the reaction
    /// extents in all voxels of the compartment.
    /// \param c Name of the compartment.
    /// \param r Name of the reaction.
    unsigned long long getCompComplexReacExtent(std::string const& c, std::string const& r) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Returns the volume of a tetrahedron (in m^3).
    ///
    /// \param tidx Index of the tetrahedron.
    double getTetVol(tetrahedron_global_id tidx) const;

    /// Set the volume of a tetrahedron (in m^3).
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param vol Volume of the tetrahedron.
    void setTetVol(tetrahedron_global_id tidx, double vol);

    /// Returns whether species s is defined in a tetrahedron
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    bool getTetSpecDefined(tetrahedron_global_id tidx, std::string const& s) const;

    /// Returns the number of molecules of species s in a tetrahedron
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    double getTetSpecCount(tetrahedron_global_id tidx, std::string const& s) const;

    /// Sets the number of molecules of species s in a tetrahedron
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    /// \param n Number of molecules of the species.
    void setTetSpecCount(tetrahedron_global_id tidx, std::string const& s, double n);

    /// Returns the amount (in mols) of species s in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    double getTetSpecAmount(tetrahedron_global_id tidx, std::string const& s) const;

    /// Sets the amount (in mols) of species s in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    /// \param m Amount of the species.
    void setTetSpecAmount(tetrahedron_global_id tidx, std::string const& s, double m);

    /// Returns the concentration (in molar units) of species s in a
    /// tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    double getTetSpecConc(tetrahedron_global_id tidx, std::string const& s) const;

    /// Sets the concentration (in molar units) of species s in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    /// \param c Concentration of the species.
    void setTetSpecConc(tetrahedron_global_id tidx, std::string const& s, double c);

    /// Returns whether the concentration of species s in a tetrahedron
    /// remains constant over time (unless changed explicitly).
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    bool getTetSpecClamped(tetrahedron_global_id tidx, std::string const& s) const;

    /// Sets clamping of species s in a tetrahedron on or off.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param s Name of the species.
    /// \param buf Flag to turn the clamping of species on or off.
    void setTetSpecClamped(tetrahedron_global_id tidx, std::string const& s, bool buf);

    /// Returns the macroscopic reaction constant of reaction r in a
    /// tetrahedron (units vary with order of reaction).
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param r Name of the reaction.
    double getTetReacK(tetrahedron_global_id tidx, std::string const& r) const;

    /// Sets the macroscopic reaction constant of reaction r in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param r Name of the reaction.
    /// \param kf Rate constant of the reaction.
    void setTetReacK(tetrahedron_global_id tidx, std::string const& r, double kf);

    /// Returns whether reaction r in a tetrahedron is active or not
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param r Name of the reaction.
    bool getTetReacActive(tetrahedron_global_id tidx, std::string const& r) const;

    /// Activates/deactivates reaction r in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param r Name of the reaction.
    /// \param act Flag to activate or deactivate the reaction.
    void setTetReacActive(tetrahedron_global_id tidx, std::string const& r, bool act);

    /// Returns the diffusion constant of diffusion rule d in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param d Name of the deffusion.
    /// \param direction_tet Tetrahedron index which specifies diffusion
    /// direction.
    double getTetDiffD(tetrahedron_global_id tidx,
                       std::string const& d,
                       tetrahedron_global_id direction_tet) const;

    /// Sets the diffusion constant of diffusion rule d in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param d Name of the diffusion.
    /// \param dk Rate constant of the diffusion.
    /// \param direction_tet Tetrahedron index which the diffusion towards.
    void setTetDiffD(tetrahedron_global_id tidx,
                     std::string const& d,
                     double dk,
                     tetrahedron_global_id direction_tet);

    /// Returns whether diffusion rule d in a tetrahedron is active or not.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param d Name of the diffusion.
    bool getTetDiffActive(tetrahedron_global_id tidx, std::string const& d) const;

    /// Activates/deactivates diffusion rule d in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param d Name of the diffusion.
    /// \param act Flag to activate / deactivate the diffusion.
    void setTetDiffActive(tetrahedron_global_id tidx, std::string const& d, bool act);

    ////////////////////////////////////////////////////////////////////////

    /// Returns c_mu, the mesoscopic reaction constant of reaction r in
    /// a tetrahedron
    ///
    /// \param tidx Index of the diffusion.
    /// \param r Name of the reaction.
    double getTetReacC(tetrahedron_global_id tidx, std::string const& r) const;

    /// Returns h_mu, the distinct number of ways in which reaction r can
    /// occur in a tetrahedron, by computing the product of its reactants.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param Name of the reaction.
    double getTetReacH(tetrahedron_global_id tidx, std::string const& r) const;

    /// Returns the propensity, a_mu, of reaction r in a tetrahedron.
    /// The propensity value gives the probability per unit time that this
    /// reaction will occur in the current state.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param r Name of the reaction.
    double getTetReacA(tetrahedron_global_id tidx, std::string const& r) const;

    /// Returns the propensity, a_mu of diffusion rule d in a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param d Name of the diffusion.
    double getTetDiffA(tetrahedron_global_id tidx, std::string const& d) const;

    ////////////////////////////////////////////////////////////////////////

    /// Returns the potential of tetrahedron in Volts.
    ///
    /// \param tidx Index of the tetrahedron.
    double getTetV(tetrahedron_global_id tidx) const;

    /// Set the potential of tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron
    /// \param V Potential in volts.
    void setTetV(tetrahedron_global_id tidx, double v);

    /// Returns whether the potential of tetrahedron is clamped over time
    /// (unless changed explicitly)
    ///
    /// \param tidx Index of the tetrahedron
    bool getTetVClamped(tetrahedron_global_id tidx) const;

    /// Sets voltage clamp in tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron
    /// \param cl Flag to turn the clamping on or off.
    void setTetVClamped(tetrahedron_global_id tidx, bool cl);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    /// Returns the area of patch p (in m^2)
    ///
    /// \param p Name of the patch.
    double getPatchArea(std::string const& p) const;

    /// Sets the area of patch p.
    ///
    /// NOTE: this method may throw an exception if this does not make sense
    /// for a solver (e.g. a tetrahedral mesh-based solver).
    /// \param p Name of the patch.
    /// \param area Area of the patch.
    void setPatchArea(std::string const& p, double area);

    /// Returns the number of molecules of species s in patch p.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all triangles of the patch.
    /// \param p Name of the path.
    /// \param s Name of the species.
    double getPatchSpecCount(std::string const& p, std::string const& s) const;

    /// Sets the number of molecules of species s in patch p.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all triangles in the patch.
    /// \param p Name of the patch.
    /// \param s Name of the species.
    /// \param n Number of molecules of species.
    void setPatchSpecCount(std::string const& p, std::string const& s, double n);

    /// Returns the amount (in mols) of species s in patch p.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is computed as
    /// the sum of the amounts in all triangles of the patch.
    /// \param p Name of the patch.
    /// \param s Name of the species.
    double getPatchSpecAmount(std::string const& p, std::string const& s) const;

    /// Sets the amount (in mols) of species s in patch p.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all triangles in the patch.
    /// \param p Name of the patch.
    /// \param s Name of the species.
    /// \param a Amount of the species.
    void setPatchSpecAmount(std::string const& p, std::string const& s, double a);

    /// Returns whether the count of species s in patch p remains constant.
    /// over time (unless changed explicitly).
    ///
    /// NOTE: this method will only return true if the species has been
    /// clamped in all triangles in the patch.
    /// \param p Name of the patch.
    /// \param s Name of the species.
    bool getPatchSpecClamped(std::string const& p, std::string const& s) const;

    /// Turns clamping of species in patch p on or off.
    ///
    /// NOTE: in a mesh-based simulation, this method turns clamping on/off
    /// in all triangles in the patch.
    /// \param p Name of the patch.
    /// \param s Name of the species.
    /// \param buf Flag to turn clamping of species on /off.
    void setPatchSpecClamped(std::string const& p, std::string const& s, bool buf);

    /// Returns the number of Complexes s matching filter f in patch c.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all voxels of the patch.
    /// \param c Name of the patch.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    double getPatchComplexCount(std::string const& c,
                                std::string const& s,
                                const std::vector<std::vector<model::SubunitStateFilter>>& f) const;

    /// Sets the number of Complexes s matching filter f in patch c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all voxels in the patch (i.e. a uniform distribution).
    /// \param c Name of the patch.
    /// \param s Name of the complex.
    /// \param i Init state for the complex.
    /// \param n Number of molecules of the species.
    void setPatchComplexCount(std::string const& c,
                              std::string const& s,
                              const std::vector<std::vector<model::SubunitStateFilter>>& i,
                              double n);

    /// Returns the amount (in mols) of complex s matching filter f in patch c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is computed as
    /// the sum of the amounts in all voxels of the patch.
    /// \param c Name of the patch.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    double getPatchComplexAmount(
        std::string const& c,
        std::string const& s,
        const std::vector<std::vector<model::SubunitStateFilter>>& f) const;

    /// Set the amount (in mols) of complex s matching filter f in patch c.
    ///
    /// NOTE: in a mesh-based simulation, the total amount is equally divided
    /// over all voxels in the patch (i.e. a uniform distribution).
    /// \param c Name of the patch.
    /// \param s Name of the complex.
    /// \param i Init state for the complex.
    /// \param a Amount of the complex.
    void setPatchComplexAmount(std::string const& c,
                               std::string const& s,
                               const std::vector<std::vector<model::SubunitStateFilter>>& i,
                               double a);

    /// Returns the number of SubUnits in state m for Complexes s matching filter f in patch
    /// c.
    ///
    /// NOTE: in a mesh-based simulation, the total count is computed as
    /// the sum of the counts in all voxels of the patch.
    /// \param c Name of the patch.
    /// \param s Name of the complex.
    /// \param f Filter for the complex.
    //  \param m index of the subunitstate
    double getPatchComplexSUSCount(std::string const& c,
                                   std::string const& s,
                                   const std::vector<std::vector<model::SubunitStateFilter>>& f,
                                   uint m) const;

    /// Returns the macroscopic reaction constant of surface reaction r
    /// in patch p.
    ///
    /// NOTE: in a mesh-based simulation, the value is computed as the
    /// area-weighted sum of the reaction constants in all triangles of
    /// the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    double getPatchSReacK(std::string const& p, std::string const& r) const;

    /// Sets the macroscopic reaction constant of surface reaction r
    /// in patch p.
    ///
    /// NOTE: in a mesh-based simulation this method changes the reaction
    /// constant equally in all triangles of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    /// \param kf Rate constant of the reaction.
    void setPatchSReacK(std::string const& p, std::string const& r, double kf);

    /// Returns whether surface reaction r in patch p is active or not.
    ///
    /// NOTE: in a mesh-based simulation, only returns false when the
    /// reaction has been inactivated in all triangles.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    bool getPatchSReacActive(std::string const& p, std::string const& r) const;

    /// Activate or inactivate surface reaction r in patch p.
    ///
    /// NOTE: in a mesh-based simulation, activation/inactivation of a
    /// surface reaction turns it on/off in all triangles.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    /// \param a Flag to activate / deactivate the reaction.
    void setPatchSReacActive(std::string const& p, std::string const& r, bool a);

    ////////////////////////////////////////////////////////////////////////

    /// Returns c_mu, the mesoscopic reaction constant of surface reaction r
    /// in patch p.
    ///
    /// NOTE: in a mesh_based simulation, the mesoscopic reaction constant
    /// is computed as the sum of the mesoscopic reaction constants from all
    /// triangles of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reacton.
    double getPatchSReacC(std::string const& p, std::string const& r) const;

    /// Returns h_mu, the distinct number of ways in which a surface reaction
    /// r can occur in patch p, by computing the product of its reactants.
    ///
    /// NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
    /// over all triangles triangles of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    double getPatchSReacH(std::string const& p, std::string const& r) const;

    /// Returns the propensity, a_mu of surface reaction r in patch p.
    ///
    /// This propensity value gives the probability per unit time that this
    /// surface reaction will occur in the current state.
    ///
    /// NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
    /// a_mu in all triangles in the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    double getPatchSReacA(std::string const& p, std::string const& r) const;

    /// Returns the extent of surface reaction r in patch p.
    ///
    /// NOTE: in a mesh-based simulation, returns the sum of the extents in
    /// all triangles of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    unsigned long long getPatchSReacExtent(std::string const& p, std::string const& r) const;

    /// Resets the extent of surface reaction r in patch p to zero.
    ///
    /// NOTE: in a mesh-based simulation, resets the extents of the
    /// surface reaction in all triangles of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    void resetPatchSReacExtent(std::string const& p, std::string const& r);

    /// Returns the extent of complex surface reaction r in patch p.
    ///
    /// NOTE: in a mesh-based simulation, returns the sum of the reaction
    /// extents in all voxels of the patch.
    /// \param p Name of the patch.
    /// \param r Name of the reaction.
    unsigned long long getPatchComplexSReacExtent(std::string const& p, std::string const& r) const;

    /// Returns whether voltage-dependent surface reaction vsr in patch p is
    /// active or not.
    ///
    /// NOTE: only returns false when the voltage-dependent surface
    /// reaction has been inactivated in all triangles.
    /// \param p Name of the patch.
    /// \param vsr Name of the voltage-dependent surface reaction.
    bool getPatchVDepSReacActive(std::string const& p, std::string const& vsr) const;

    /// Activate or inactivate voltage-dependent surface reaction vsr in patch
    /// p.
    ///
    /// NOTE: activation/inactivation of a voltage-dependent
    /// surface reaction turns it on/off in all triangles.
    /// \param p Name of the patch.
    /// \param vsr Name of the voltage-dependent surface reaction.
    /// \param a Flag to activate / deactivate the reaction.
    void setPatchVDepSReacActive(std::string const& p, std::string const& vsr, bool a);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      PATCH, RAFT-RELATED
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of rafts r in patch p
    ///
    /// \param p Name of the patch.
    /// \param r Name of the raft.
    uint getPatchRaftCount(std::string const& p, std::string const& r) const;

    /// Set the number of rafts r in patch p
    ///
    /// \param p Name of the patch.
    /// \param r Name of the raft.
    /// \param n Number of rafts
    void setPatchRaftCount(std::string const& p, std::string const& r, uint n);

    /// Returns the position of raft of type r with unique index
    /// raft_unique_index (in cartesian coordinates)
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    std::vector<double> getSingleRaftPos(std::string const& r,
                                         raft_individual_id raft_unique_index) const;

    /// Returns the count of species s in raft r in patch p. Return is a map,
    /// raft_unique_index : count of s
    ///
    /// \param p Name of the patch.
    /// \param r Name of the raft.
    /// \param s Name of the species
    std::map<solver::raft_individual_id, uint> getPatchRaftSpecCountDict(
        std::string const& p,
        std::string const& r,
        std::string const& s) const;
    /// Returns the summed count of species s in raft r in patch p.
    ///
    /// \param p Name of the patch.
    /// \param r Name of the raft.
    /// \param s Name of the species
    uint getPatchRaftSpecCount(std::string const& p,
                               std::string const& r,
                               std::string const& s) const;

    /// Get the count of species s on raft of type r and unique index
    /// raft_unique_index.
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param s Name of the species
    uint getSingleRaftSpecCount(std::string const& r,
                                raft_individual_id raft_unique_index,
                                std::string const& s) const;

    /// Set the count of species s on raft of type r and unique index
    /// raft_unique_index.
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param s Name of the species
    /// \param count The number of the species
    void setSingleRaftSpecCount(std::string const& r,
                                raft_individual_id raft_unique_index,
                                std::string const& s,
                                uint count);

    /// Get the 'immobility' of raft of type r and unique index
    /// raft_unique_index. All non-zero numbers mean raft is
    /// immobile.
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    uint getSingleRaftImmobility(std::string const& r, raft_individual_id raft_unique_index) const;

    /// Get the raftendocytosis rate of a raft
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param rendo Name of the raft endocytosis
    double getSingleRaftRaftEndocytosisK(std::string const& r,
                                         raft_individual_id raft_unique_index,
                                         std::string const& rendo) const;

    /// Set the raftendocytosis rate of a raft
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param rendo Name of the raft endocytosis
    /// \param k Rate
    void setSingleRaftRaftEndocytosisK(std::string const& r,
                                       raft_individual_id raft_unique_index,
                                       std::string const& rendo,
                                       double k);

    /// Activate or de-activate a raft surface reaction on a raft
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param rsreac Name of the raft surface reaction
    /// \param active Whether the raft surface reaction should be active
    void setSingleRaftSReacActive(std::string const& r,
                                  raft_individual_id raft_unique_index,
                                  std::string const& rsreac,
                                  bool active);

    /// returns whether raft surface reaction on a raft is active or not
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft.
    /// \param rsreac Name of the raft surface reaction
    bool getSingleRaftSReacActive(std::string const& r,
                                  raft_individual_id raft_unique_index,
                                  std::string const& rsreac) const;

    /// Returns all the unique indices of rafts r currently in patch p
    ///
    /// \param p Name of the patch.
    /// \param r Name of the raft.
    std::vector<raft_individual_id> getPatchRaftIndices(std::string const& p,
                                                        std::string const& r) const;

    /// Get the patch of raft of type r and unique index raft_unique_index
    ///
    /// \param r Name of the raft.
    /// \param raft_unique_index Unique index of the individual raft
    std::string getSingleRaftPatch(std::string const& r,
                                   raft_individual_id raft_unique_index) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      PATCH, ENDOCYTOSIS-RELATED
    ////////////////////////////////////////////////////////////////////////

    /// Activate or de-activate endocytosis reaction in an endocytic zone.
    ///
    /// \param p Name of the patch containing the zone
    /// \param zone Name of the zone
    /// \param endo Name of the endocytosis
    /// \param active Whether the endocytosis should be active
    void setPatchEndocyticZoneEndocytosisActive(std::string const& p,
                                                std::string const& zone,
                                                std::string const& endo,
                                                bool active);

    /// Set rate of endocytosis reaction in an endocytic zone.
    ///
    /// \param p Name of the patch containing the zone
    /// \param zone Name of the zone
    /// \param endo Name of the endocytosis
    /// \param k Rate of endocytosis in the zone
    void setPatchEndocyticZoneEndocytosisK(std::string const& p,
                                           std::string const& zone,
                                           std::string const& endo,
                                           double k);

    /// Returns the endocytosis extent in an endocytic zone
    ///
    /// \param p Name of the patch containing the zone
    /// \param zone Name of the zone
    /// \param endo Name of the endocytosis.
    uint getPatchEndocyticZoneEndocytosisExtent(std::string const& p,
                                                std::string const& zone,
                                                std::string const& endo) const;

    /// Returns the endocytosis events that happened since last call in an endocytic zone
    ///
    /// \param p Name of the patch containing the zone
    /// \param zone Name of the zone
    /// \param endo Name of the endocytosis.
    std::vector<EndocytosisEvent> getPatchEndocyticZoneEndocytosisEvents(
        std::string const& p,
        std::string const& zone,
        std::string const& endo) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    /// Activate or inactivate diffusion across a diffusion boundary for a
    /// species.
    ///
    /// \param db Name of the diffusion boundary.
    /// \param s Name of the species.
    /// \param act Bool to activate (true) or inactivate (false) diffusion.
    void setDiffBoundarySpecDiffusionActive(std::string const& db, std::string const& s, bool act);

    /// Returns whether diffusion is active across a diffusion boundary for a
    /// species.
    ///
    /// \param db Name of the diffusion boundary.
    /// \param s Name of the species.
    bool getDiffBoundarySpecDiffusionActive(std::string const& db, std::string const& s) const;

    /// Set the diffusion constant across a diffusion boundary.
    ///
    /// \param db Name of the diffusion boundary.
    /// \param s Name of the species.
    /// \param dcst diffusion constant.
    /// \param direction_comp ID of the compartment which the diffusion towards
    /// to.
    void setDiffBoundarySpecDcst(std::string const& db,
                                 std::string const& s,
                                 double dcst,
                                 std::string const& direction_comp = "");

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      SURFACE DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    /// Activate or inactivate diffusion across a surface diffusion boundary
    /// for a species.
    ///
    /// \param sdb Name of the surface diffusion boundary.
    /// \param s Name of the species.
    /// \param act Bool to activate (true) or inactivate (false) diffusion.
    void setSDiffBoundarySpecDiffusionActive(std::string const& sdb,
                                             std::string const& s,
                                             bool act);

    /// Returns whether diffusion is active across a surface diffusion boundary
    /// for a species.
    ///
    /// \param sdb Name of the surface diffusion boundary.
    /// \param s Name of the species.
    bool getSDiffBoundarySpecDiffusionActive(std::string const& sdb, std::string const& s) const;

    /// Set the diffusion constant across a surface diffusion boundary.
    ///
    /// \param sdb Name of the surface diffusion boundary.
    /// \param s Name of the species.
    /// \param dcst diffusion constant.
    /// \param direction_patch ID of the patch which the diffusion is towards to.
    void setSDiffBoundarySpecDcst(std::string const& sdb,
                                  std::string const& s,
                                  double dcst,
                                  std::string const& direction_patch = "");

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Returns the area of the triangle (in m^2).
    ///
    /// \param tidx Index of the triangle.
    double getTriArea(triangle_global_id tidx) const;

    /// Set the area (in m^2) of the triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param area Area of teh triangle.
    void setTriArea(triangle_global_id tidx, double area);

    /// Returns whether species s is defined in a triangle
    ///
    /// \param tidx Index of the triangle.
    /// \param s Name of the species.
    bool getTriSpecDefined(triangle_global_id tidx, std::string const& s) const;

    /// Returns the number of molecules of species s in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param s Name of the species.
    double getTriSpecCount(triangle_global_id tidx, std::string const& s) const;

    /// Sets the number of molecules of species s in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param s Name of the species.
    /// \param n Number of molecules of the species.
    void setTriSpecCount(triangle_global_id tidx, std::string const& s, double n);

    /// Returns the amount (in mols) of species s in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param s Name of the species.
    double getTriSpecAmount(triangle_global_id tidx, std::string const& s) const;

    /// Sets the amount (in mols) of species s in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param s Name of the species.
    /// \param m Amount of the species.
    void setTriSpecAmount(triangle_global_id tidx, std::string const& s, double m);

    /// Returns whether the number of molecules of species s in a triangle
    /// remains constant over time (unless changed explicitly)
    ///
    /// \param tidx Index of the triangle.
    /// \param s name of the species.
    bool getTriSpecClamped(triangle_global_id tidx, std::string const& s) const;

    /// Sets clamping of species s in a triangle on or off.
    ///
    /// \param tidx Index of the triangle.
    /// \param s name of the species.
    /// \param buf Flag to set clamping of species on /off.
    void setTriSpecClamped(triangle_global_id tidx, std::string const& s, bool buf);

    /// Returns the macroscopic reaction constant of surface reaction r
    /// in a triangle (units vary with order of reaction).
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    double getTriSReacK(triangle_global_id tidx, std::string const& r);

    /// Sets the macroscopic reaction constant of surface reaction r in
    /// a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    /// \param kf Rate constant of the reaction.
    void setTriSReacK(triangle_global_id tidx, std::string const& r, double kf);

    /// Returns whether surface reaction r in a triangle is active or not.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    bool getTriSReacActive(triangle_global_id tidx, std::string const& r);

    /// Activates/inactivates surface reaction r in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    /// \param act Flag to activate / deactivate the reaction.
    void setTriSReacActive(triangle_global_id tidx, std::string const& r, bool act);

    ////////////////////////////////////////////////////////////////////////

    /// Returns c_mu, the mesoscopic reaction constant of surface reaction r
    /// in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    double getTriSReacC(triangle_global_id tidx, std::string const& r);

    /// Returns h_mu, the distinct number of ways in which surface reaction r
    /// can occur in a triangle, by computing the product of it's reactants.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    double getTriSReacH(triangle_global_id tidx, std::string const& r);

    /// Returns the propensity, a_mu, of surface reaction r in a triangle.
    /// The propensity value gives the probability per unit time that this
    /// surface reaction will occur in the current state.
    ///
    /// \param tidx Index of the triangle.
    /// \param r name of the reaction.
    double getTriSReacA(triangle_global_id tidx, std::string const& r);

    /// outdate function
    double getTriDiffD(triangle_global_id tidx,
                       std::string const& d,
                       uint direction_tri = std::numeric_limits<uint>::max());

    /// Returns the diffusion constant of diffusion rule d in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param d Name of the diffusion.
    /// \param direction_tet Triangle index which specifies diffusion direction.
    double getTriSDiffD(triangle_global_id tidx,
                        std::string const& d,
                        triangle_global_id direction_tri);

    /// outdated function
    void setTriDiffD(triangle_global_id tidx,
                     std::string const& d,
                     double dk,
                     triangle_global_id direction_tri);

    /// Sets the diffusion constant of diffusion rule d on a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param d Name of the diffusion.
    /// \param dk Rate constant of the diffusion.
    /// \param direction_tri Triangle index which the diffusion towards
    void setTriSDiffD(triangle_global_id tidx,
                      std::string const& d,
                      double dk,
                      triangle_global_id direction_tri);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      TRIANGULAR SURFACE ELEMENTS: Raft/Vesicle related
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of rafts in a triangle
    ///
    /// \param tidx Index of the triangle.
    /// \param r Name of the raft.
    uint getTriRaftCount(triangle_global_id tidx, std::string const& r) const;

    /// Returns the number of rafts in a triangle
    ///
    /// \param tidx Index of the triangle.
    /// \param r Name of the raft.
    /// \param n Number of rafts
    void setTriRaftCount(triangle_global_id tidx, std::string const& r, uint n);

    /// Add one raft to a triangle
    ///
    /// \param tidx Index of the triangle.
    /// \param r Name of the raft.
    raft_individual_id addTriRaft(triangle_global_id tidx, std::string const& r);

    /// Returns whether exocytotic reaction in a triangle is active or not.
    ///
    /// \param tidx Index of the triangle.
    /// \param er name of the exocytotic reaction.
    bool getTriExocytosisActive(triangle_global_id tidx, std::string const& er) const;

    /// Activates/inactivates exocytotic reaction in a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param er name of the exocytotic reaction.
    /// \param act Flag to activate / deactivate the exocytotic reaction.
    void setTriExocytosisActive(triangle_global_id tidx, std::string const& er, bool act);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS:
    //      TRIANGULAR SURFACE ELEMENTS: Efield-related
    ////////////////////////////////////////////////////////////////////////

    /// Returns the potential of triangle in Volts.
    ///
    /// \param tidx Index of the triangle.
    double getTriV(triangle_global_id tidx) const;

    /// Set the potential of triangle.
    ///
    /// \param tidx Index of the triangle
    /// \param V Potential in volts.
    void setTriV(triangle_global_id tidx, double v);

    /// Returns whether the potential of triangle is clamped over time
    /// (unless changed explicitly).
    ///
    /// \param tidx Index of the triangle
    bool getTriVClamped(triangle_global_id tidx) const;

    /// Sets voltage clamp in triangle.
    ///
    /// \param tidx Index of the triangle
    /// \param cl Flag to turn the clamping on or off.
    void setTriVClamped(triangle_global_id tidx, bool cl);

    /// Gets the reversal potential of ohmic current of triangle in volts.
    ///
    /// \param tidx Index of the triangle.
    /// \param oc name of the ohmic current
    double getTriOhmicErev(triangle_global_id tidx, std::string const& oc) const;

    /// Sets the reversal potential of ohmic current of triangle in volts.
    ///
    /// \param tidx Index of the triangle.
    /// \param oc name of the ohmic current
    /// \param erev reversal potential
    void setTriOhmicErev(triangle_global_id tidx, std::string const& oc, double erev);

    /// Returns the ohmic current of triangle in amperes.
    ///
    /// \param tidx Index of the triangle.
    double getTriOhmicI(triangle_global_id tidx) const;

    /// Returns the ohmic current of triangle in amperes.
    ///
    /// \param tidx Index of the triangle.
    /// \param oc name of the ohmic current
    double getTriOhmicI(triangle_global_id tidx, std::string const& oc) const;

    /// Returns the GHK current of triangle in amperes.
    ///
    /// \param tidx Index of the triangle.
    double getTriGHKI(triangle_global_id tidx) const;

    /// Returns the GHK current of triangle in amperes.
    ///
    /// \param tidx Index of the triangle.
    /// \param ghk name of the ghk current
    double getTriGHKI(triangle_global_id tidx, std::string const& ghk) const;

    /// Returns the current of a triangle in amperes from the last EField
    /// calculation step.
    ///
    /// \param tidx Index of the triangle
    double getTriI(triangle_global_id tidx) const;

    /// Gets current injection to triangle.
    /// \param tidx Index of the triangle
    double getTriIClamp(triangle_global_id tidx) const;

    /// Sets current injection to triangle.
    /// Will be assumed to be constant for one EField DT
    /// \param tidx Index of the triangle
    /// \param I Current in amperes.
    void setTriIClamp(triangle_global_id tidx, double i);

    /// Returns the macroscopic reaction constant of voltage-dependent surface reaction vsr
    /// in a triangle (units vary with order of reaction).
    ///
    /// \param tidx Index of the triangle.
    /// \param vsr name of the voltage-dependent surface reaction.
    double getTriVDepSReacK(triangle_global_id tidx, std::string const& vsr) const;

    /// Returns whether voltage-dependent surface reaction vsr in a triangle is
    /// active or not.
    ///
    /// \param tidx Index of the triangle.
    /// \param vsr name of the voltage-dependent surface reaction.
    bool getTriVDepSReacActive(triangle_global_id tidx, std::string const& vsr) const;

    /// Activates/inactivates voltage-dependent surface reaction vsr in a
    /// triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param vsr name of the voltage-dependent surface reaction.
    /// \param act Flag to activate / deactivate the reaction.
    void setTriVDepSReacActive(triangle_global_id tidx, std::string const& vsr, bool act);

    /// Set the specific capacitance of a triangle surface element.
    ///
    /// \param tidx Index of the triangle surface element
    /// \param cm Specific membrane capacitance (farad / m^2)
    void setTriCapac(triangle_global_id tidx, double cm);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VERTICES ELEMENTS (all EField-related)
    ////////////////////////////////////////////////////////////////////////

    /// Returns the potential of vertex in Volts.
    ///
    /// \param vidx Index of the vertex.
    double getVertV(vertex_id_t vidx) const;

    /// Set the potential of vertex.
    ///
    /// \param vidx Index of the vertex
    /// \param V Potential in volts.
    void setVertV(vertex_id_t vidx, double v);

    /// Returns whether the potential of vertex is clamped over time
    /// (unless changed explicitly).
    ///
    /// \param vidx Index of the vertex
    bool getVertVClamped(vertex_id_t vidx) const;

    /// Sets voltage clamp in vertex.
    ///
    /// \param vidx Index of the vertex
    /// \param cl Flag to turn the clamping on or off.
    void setVertVClamped(vertex_id_t vidx, bool cl);

    /// Gets current injection to vertex.
    /// \param vidx Index of the vertex
    double getVertIClamp(vertex_id_t vidx) const;

    /// Sets current injection to vertex.
    /// Will be assumed to be constant for one EField DT
    /// \param vidx Index of the vertex
    /// \param I Current in amperes.
    void setVertIClamp(vertex_id_t vidx, double i);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      MEMBRANES (all EField-related)
    ////////////////////////////////////////////////////////////////////////

    /// Set the electric potential of the membrane, including all nodes
    /// in the conduction volume.
    /// \param m Name of the membrane
    /// \param v Potential (volts)
    void setMembPotential(std::string const& m, double v);

    /// Set the specific membrane capacitance of the membrane
    /// \param m Name of the membrane
    /// \param cm Specific membrane capacitance (farad / m^2)
    void setMembCapac(std::string const& m, double cm);

    /// Set the bulk electrical resistivity of the section of the mesh
    /// representing the volume conductor
    ///
    /// \param m Name of the membrane
    /// \param ro Electrical resistivity (ohm.m)
    void setMembVolRes(std::string const& m, double ro);

    /// Set the resistivity of the membrane
    /// \param m Name of the membrane
    /// \param ro membrane resistivity (ohm.m^2)
    /// \param vrev Reversal potential (Volts)
    void setMembRes(std::string const& m, double ro, double vrev);
    std::pair<double, double> getMembRes(std::string const& m);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // DATA RECORDING:
    ////////////////////////////////////////////////////////////////////////
    /// Count the number of compartments
    uint getNComps() const;

    /// Count the number of patches
    uint getNPatches() const;

    /// Get a compartment's name by its index in the solver
    std::string getCompName(comp_global_id c_idx) const;

    /// Get a patch's name by its index in the solver
    std::string getPatchName(patch_global_id p_idx) const;

    /// Get number of species in a compartment
    uint getNCompSpecs(comp_global_id c_idx) const;

    /// Get number of species in a patch
    uint getNPatchSpecs(patch_global_id p_idx) const;

    /// Get species name of a compartment
    std::string getCompSpecName(comp_global_id c_idx, spec_local_id s_idx) const;

    /// Get species name of a patch
    std::string getPatchSpecName(patch_global_id p_idx, spec_local_id s_idx) const;

    ////////////////////////////////////////////////////////////////////////
    // Batch Data Access
    ////////////////////////////////////////////////////////////////////////

    /// Get species counts of a list of tetrahedrons
    virtual std::vector<double> getBatchTetSpecCounts(const std::vector<index_t>& tets,
                                                      std::string const& s) const;

    /// Get species concentrations of a list of tetrahedrons
    virtual std::vector<double> getBatchTetSpecConcs(const std::vector<index_t>& tets,
                                                     std::string const& s) const;

    /// Set species concentrations of a list of tetrahedrons
    virtual void setBatchTetSpecConcs(const std::vector<index_t>& tets,
                                      std::string const& s,
                                      const std::vector<double>& concs);

    /// Get species counts of a list of triangles
    virtual std::vector<double> getBatchTriSpecCounts(const std::vector<index_t>& tris,
                                                      std::string const& s) const;

    /// Get species counts of a list of tetrahedrons
    virtual void getBatchTetSpecCountsNP(const index_t* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size) const;

    /// Get species concentrations of a list of tetrahedrons
    virtual void getBatchTetSpecConcsNP(const index_t* indices,
                                        size_t input_size,
                                        std::string const& s,
                                        double* counts,
                                        size_t output_size) const;

    /// Set species concentrations of a list of tetrahedrons
    virtual void setBatchTetSpecConcsNP(const index_t* indices,
                                        size_t ntets,
                                        std::string const& s,
                                        const double* concs,
                                        size_t output_size);

    /// Get species counts of a list of triangles
    virtual void getBatchTriSpecCountsNP(const index_t* indices,
                                         size_t input_size,
                                         std::string const& s,
                                         double* counts,
                                         size_t output_size) const;

    ////////////////////////////////////////////////////////////////////////
    // ROI Data Access
    ////////////////////////////////////////////////////////////////////////

    /// Get species counts of a list of tetrahedrons
    virtual std::vector<double> getROITetSpecCounts(const std::string& ROI_id,
                                                    std::string const& s) const;

    /// Get species counts of a list of triangles
    virtual std::vector<double> getROITriSpecCounts(const std::string& ROI_id,
                                                    std::string const& s) const;

    /// Get species counts of a list of tetrahedrons
    virtual void getROITetSpecCountsNP(const std::string& ROI_id,
                                       std::string const& s,
                                       double* counts,
                                       size_t output_size) const;

    /// Get species counts of a list of triangles
    virtual void getROITriSpecCountsNP(const std::string& ROI_id,
                                       std::string const& s,
                                       double* counts,
                                       size_t output_size) const;

    /// Get the volume of a ROI.
    virtual double getROIVol(const std::string& ROI_id) const;

    /// Get the area of a ROI.
    virtual double getROIArea(const std::string& ROI_id) const;

    /// Get the count of a species in a ROI.
    virtual double getROISpecCount(const std::string& ROI_id, std::string const& s) const;

    /// Set the count of a species in a ROI.
    virtual void setROISpecCount(const std::string& ROI_id, std::string const& s, double count);

    /// Get the amount of a species in a ROI.
    virtual double getROISpecAmount(const std::string& ROI_id, std::string const& s) const;

    /// Set the amount of a species in a ROI.
    virtual void setROISpecAmount(const std::string& ROI_id, std::string const& s, double amount);

    /// Get the concentration of a species in a ROI.
    virtual double getROISpecConc(const std::string& ROI_id, std::string const& s) const;

    /// Set the concentration of a species in a ROI.
    virtual void setROISpecConc(const std::string& ROI_id, std::string const& s, double conc);

    /// Set a species in a ROI to be clamped or not. The count of species s in
    /// the ROI is clamped if b is True, not clamped if b is False.
    virtual void setROISpecClamped(const std::string& ROI_id, std::string const& s, bool b);

    /// Sets the macroscopic reaction constant of reaction with identifier
    /// string r in a ROI with identifier string ROI_id to kf. The unit of the
    /// reaction constant depends on the order of the reaction.
    ///
    /// Note: The default value still comes from the steps.model description,
    /// so calling reset() will return the reaction constant to that value.
    virtual void setROIReacK(const std::string& ROI_id, std::string const& r, double kf);

    /// Sets the macroscopic reaction constant of surface reaction with
    /// identifier string sr in a ROI with identifier string ROI_id to kf. The
    /// unit of the reaction constant depends on the order of the reaction.
    ///
    /// Note: The default value still comes from the steps.model description,
    /// so calling reset() will return the reaction constant to that value.
    virtual void setROISReacK(const std::string& ROI_id, std::string const& sr, double kf);

    /// Sets the macroscopic diffusion constant of diffusion with identifier
    /// string d in a ROI with identifier string ROI_id to dk.
    ///
    /// Note: The default value still comes from the steps.model description,
    /// so calling reset() will return the diffusion constant to that value.
    virtual void setROIDiffD(const std::string& ROI_id, std::string const& d, double dk);

    /// Set reaction r in a ROI to be active or not.
    virtual void setROIReacActive(const std::string& ROI_id, std::string const& r, bool a);

    /// Set surface reaction sr in a ROI to be active or not.
    virtual void setROISReacActive(const std::string& ROI_id, std::string const& sr, bool a);

    /// Set diffusion d in a ROI to be active or not.
    virtual void setROIDiffActive(const std::string& ROI_id, std::string const& d, bool act);

    /// Set voltage dependent surface reaction vsr in a ROI to be active or
    /// not.
    virtual void setROIVDepSReacActive(const std::string& ROI_id, std::string const& vsr, bool a);

    /// Return the extent of reaction with identifier string r in ROI with
    /// identifier string ROI_id, that is the number of times the reaction has
    /// occurred up to the current simulation time.
    virtual unsigned long long getROIReacExtent(const std::string& ROI_id,
                                                std::string const& r) const;

    /// Reset the extent of reaction with identifier string r in ROI with
    /// identifier string ROI_id, that is the number of times the reaction has
    /// occurred up to the current simulation time, to 0.
    virtual void resetROIReacExtent(const std::string& ROI_id, std::string const& r);

    /// Return the extent of surface reaction with identifier string sr in ROI
    /// with identifier string ROI_id, that is the number of times the reaction
    /// has occurred up to the current simulation time.
    virtual unsigned long long getROISReacExtent(const std::string& ROI_id,
                                                 std::string const& sr) const;

    /// Reset the extent of surface reaction with identifier string r in ROI
    /// with identifier string ROI_id, that is the number of times the reaction
    /// has occurred up to the current simulation time, to 0.
    virtual void resetROISReacExtent(const std::string& ROI_id, std::string const& sr);

    /// Return the extent of diffusion with identifier string d in ROI with
    /// identifier string ROI_id, that is the number of times the diffusion has
    /// occurred up to the current simulation time.
    virtual unsigned long long getROIDiffExtent(const std::string& ROI_id,
                                                std::string const& d) const;

    /// Reset the extent of diffusion with identifier string d in ROI with
    /// identifier string ROI_id, that is the number of times the diffusion has
    /// occurred up to the current simulation time, to 0.
    virtual void resetROIDiffExtent(const std::string& ROI_id, std::string const& s);

    ////////////////////////////////////////////////////////////////////////
    // DEPRECATED SPEC METHODS
    ////////////////////////////////////////////////////////////////////////

    virtual std::vector<double> getBatchTetCounts(const std::vector<index_t>& tets,
                                                  std::string const& s) const;

    virtual std::vector<double> getBatchTetConcs(const std::vector<index_t>& tets,
                                                 std::string const& s) const;

    virtual void setBatchTetConcs(const std::vector<index_t>& tets,
                                  std::string const& s,
                                  const std::vector<double>& concs);

    virtual std::vector<double> getBatchTriCounts(const std::vector<index_t>& tris,
                                                  std::string const& s) const;

    virtual void getBatchTetCountsNP(const index_t* indices,
                                     size_t input_size,
                                     std::string const& s,
                                     double* counts,
                                     size_t output_size) const;

    virtual void getBatchTetConcsNP(const index_t* indices,
                                    size_t input_size,
                                    std::string const& s,
                                    double* counts,
                                    size_t output_size) const;

    virtual void setBatchTetConcsNP(const index_t* indices,
                                    size_t ntets,
                                    std::string const& s,
                                    const double* concs,
                                    size_t output_size);

    virtual void getBatchTriCountsNP(const index_t* indices,
                                     size_t input_size,
                                     std::string const& s,
                                     double* counts,
                                     size_t output_size) const;

    virtual std::vector<double> getROITetCounts(const std::string& ROI_id,
                                                std::string const& s) const;

    virtual std::vector<double> getROITriCounts(const std::string& ROI_id,
                                                std::string const& s) const;

    virtual void getROITetCountsNP(const std::string& ROI_id,
                                   std::string const& s,
                                   double* counts,
                                   size_t output_size) const;

    virtual void getROITriCountsNP(const std::string& ROI_id,
                                   std::string const& s,
                                   double* counts,
                                   size_t output_size) const;

    virtual double getROICount(const std::string& ROI_id, std::string const& s) const;

    virtual void setROICount(const std::string& ROI_id, std::string const& s, double count);

    virtual double getROIAmount(const std::string& ROI_id, std::string const& s) const;

    virtual void setROIAmount(const std::string& ROI_id, std::string const& s, double amount);

    virtual double getROIConc(const std::string& ROI_id, std::string const& s) const;

    virtual void setROIConc(const std::string& ROI_id, std::string const& s, double conc);

    double getCompCount(std::string const& c, std::string const& s) const;

    void setCompCount(std::string const& c, std::string const& s, double n);

    double getCompAmount(std::string const& c, std::string const& s) const;

    void setCompAmount(std::string const& c, std::string const& s, double a);

    double getCompConc(std::string const& c, std::string const& s) const;

    void setCompConc(std::string const& c, std::string const& s, double conc);

    double getTetCount(tetrahedron_global_id tidx, std::string const& s) const;

    void setTetCount(tetrahedron_global_id tidx, std::string const& s, double n);

    double getTetAmount(tetrahedron_global_id tidx, std::string const& s) const;

    void setTetAmount(tetrahedron_global_id tidx, std::string const& s, double m);

    double getTetConc(tetrahedron_global_id tidx, std::string const& s) const;

    void setTetConc(tetrahedron_global_id tidx, std::string const& s, double c);

    double getPatchCount(std::string const& p, std::string const& s) const;

    void setPatchCount(std::string const& p, std::string const& s, double n);

    double getPatchAmount(std::string const& p, std::string const& s) const;

    void setPatchAmount(std::string const& p, std::string const& s, double a);

    double getTriCount(triangle_global_id tidx, std::string const& s) const;

    void setTriCount(triangle_global_id tidx, std::string const& s, double n);

    double getTriAmount(triangle_global_id tidx, std::string const& s) const;

    void setTriAmount(triangle_global_id tidx, std::string const& s, double m);

    virtual void setROIClamped(const std::string& ROI_id, std::string const& s, bool b);

    bool getCompClamped(std::string const& c, std::string const& s) const;

    void setCompClamped(std::string const& c, std::string const& s, bool b);

    bool getTetClamped(tetrahedron_global_id tidx, std::string const& s) const;

    void setTetClamped(tetrahedron_global_id tidx, std::string const& s, bool buf);

    bool getPatchClamped(std::string const& p, std::string const& s) const;

    void setPatchClamped(std::string const& p, std::string const& s, bool buf);

    bool getTriClamped(triangle_global_id tidx, std::string const& s) const;

    void setTriClamped(triangle_global_id tidx, std::string const& s, bool buf);

    void setDiffBoundaryDiffusionActive(std::string const& db, std::string const& s, bool act);

    bool getDiffBoundaryDiffusionActive(std::string const& db, std::string const& s) const;

    void setDiffBoundaryDcst(std::string const& db,
                             std::string const& s,
                             double dcst,
                             std::string const& direction_comp = "");

    void setSDiffBoundaryDiffusionActive(std::string const& sdb, std::string const& s, bool act);

    bool getSDiffBoundaryDiffusionActive(std::string const& sdb, std::string const& s) const;

    void setSDiffBoundaryDcst(std::string const& sdb,
                              std::string const& s,
                              double dcst,
                              std::string const& direction_patch = "");


  protected:
    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    void checkpoint(std::ostream& cp_file) const;
    void restore(std::istream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    virtual double _getCompVol(comp_global_id cidx) const = 0;
    virtual void _setCompVol(comp_global_id cidx, double vol);

    virtual double _getCompSpecCount(comp_global_id cidx, spec_global_id sidx) const = 0;
    virtual void _setCompSpecCount(comp_global_id cidx, spec_global_id sidx, double n) = 0;

    virtual double _getCompSpecAmount(comp_global_id cidx, spec_global_id sidx) const = 0;
    virtual void _setCompSpecAmount(comp_global_id cidx, spec_global_id sidx, double a) = 0;

    virtual double _getCompSpecConc(comp_global_id cidx, spec_global_id sidx) const = 0;
    virtual void _setCompSpecConc(comp_global_id cidx, spec_global_id sidx, double c) = 0;

    virtual bool _getCompSpecClamped(comp_global_id cidx, spec_global_id sidx) const = 0;
    virtual void _setCompSpecClamped(comp_global_id cidx, spec_global_id sidx, bool b) = 0;

    virtual double _getCompReacK(comp_global_id cidx, reac_global_id ridx) const = 0;
    virtual void _setCompReacK(comp_global_id cidx, reac_global_id ridx, double kf) = 0;

    virtual double _getCompComplexCount(
        solver::comp_global_id cidx,
        solver::complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f)
        const;
    virtual void _setCompComplexCount(comp_global_id cidx,
                                      complex_global_id sidx,
                                      const util::strongid_vector<complex_substate_id, uint>& i,
                                      double n);

    virtual double _getCompComplexAmount(
        comp_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f)
        const;
    virtual void _setCompComplexAmount(comp_global_id cidx,
                                       complex_global_id sidx,
                                       const util::strongid_vector<complex_substate_id, uint>& i,
                                       double a);

    virtual double _getCompComplexConc(
        comp_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f)
        const;
    virtual void _setCompComplexConc(comp_global_id cidx,
                                     complex_global_id sidx,
                                     const util::strongid_vector<complex_substate_id, uint>& i,
                                     double c);

    virtual bool _getCompReacActive(comp_global_id cidx, reac_global_id ridx) const = 0;
    virtual void _setCompReacActive(comp_global_id cidx, reac_global_id ridx, bool a) = 0;

    virtual double _getCompComplexSUSCount(
        comp_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f,
        complex_substate_id m) const;

    virtual double _getCompComplexSUSConc(
        comp_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f,
        complex_substate_id m) const;

    virtual double _getCompComplexSUSAmount(
        comp_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f,
        complex_substate_id m) const;

    virtual double _getCompDiffD(comp_global_id cidx, diff_global_id didx) const;
    virtual void _setCompDiffD(comp_global_id cidx, diff_global_id didx, double dcst);

    virtual bool _getCompDiffActive(comp_global_id cidx, diff_global_id didx) const;
    virtual void _setCompDiffActive(comp_global_id cidx, diff_global_id didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getCompReacH(comp_global_id cidx, reac_global_id ridx) const;
    virtual double _getCompReacC(comp_global_id cidx, reac_global_id ridx) const;
    virtual long double _getCompReacA(comp_global_id cidx, reac_global_id ridx) const;

    virtual unsigned long long _getCompReacExtent(comp_global_id cidx, reac_global_id ridx) const;
    virtual void _resetCompReacExtent(comp_global_id cidx, reac_global_id ridx);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VESICLE-RELATED
    ////////////////////////////////////////////////////////////////////////

    virtual uint _getCompVesicleCount(comp_global_id cidx, vesicle_global_id vidx) const;
    virtual void _setCompVesicleCount(comp_global_id cidx, vesicle_global_id vidx, uint n);

    virtual vesicle_individual_id _addCompVesicle(comp_global_id cidx, vesicle_global_id vidx);
    virtual void _deleteSingleVesicle(vesicle_global_id vidx,
                                      vesicle_individual_id ves_unique_index);

    virtual uint _getSingleVesicleSurfaceLinkSpecCount(vesicle_global_id vidx,
                                                       vesicle_individual_id ves_unique_index,
                                                       linkspec_global_id lsidx) const;

    virtual std::vector<linkspec_individual_id> _getSingleVesicleSurfaceLinkSpecIndices(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        linkspec_global_id lsidx) const;

    virtual std::vector<pointspec_individual_id> _getSingleVesicleSurfaceSpecIndices(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        spec_global_id sidx) const;

    virtual std::vector<vesicle_individual_id> _getAllVesicleIndices() const;

    virtual std::vector<vesicle_individual_id> _getAllVesicleIndicesOnPath() const;

    virtual std::vector<vesicle_individual_id> _getCompVesicleIndices(comp_global_id cidx,
                                                                      vesicle_global_id vidx) const;

    virtual comp_global_id _getSingleVesicleCompartment(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index) const;

    virtual std::vector<double> _getSingleVesiclePos(vesicle_global_id vidx,
                                                     vesicle_individual_id ves_unique_index) const;

    virtual void _setSingleVesiclePos(solver::vesicle_global_id vidx,
                                      solver::vesicle_individual_id ves_unique_index,
                                      const std::vector<double>& pos,
                                      bool force = false);

    virtual std::map<vesicle_individual_id, uint> _getCompVesicleSurfaceSpecCountMap(
        comp_global_id cidx,
        vesicle_global_id vidx,
        spec_global_id sidx) const;

    virtual uint _getCompVesicleSurfaceSpecCount(comp_global_id cidx,
                                                 vesicle_global_id vidx,
                                                 spec_global_id sidx) const;

    virtual uint _getCompVesicleInnerSpecCount(comp_global_id cidx,
                                               vesicle_global_id vidx,
                                               spec_global_id sidx) const;

    virtual uint _getSingleVesicleSurfaceSpecCount(vesicle_global_id vidx,
                                                   vesicle_individual_id ves_unique_index,
                                                   spec_global_id sidx) const;

    virtual uint _getSingleVesicleInnerSpecCount(vesicle_global_id vidx,
                                                 vesicle_individual_id ves_unique_index,
                                                 spec_global_id sidx) const;

    virtual void _setSingleVesicleSurfaceSpecCount(vesicle_global_id vidx,
                                                   vesicle_individual_id ves_unique_index,
                                                   spec_global_id sidx,
                                                   uint c);

    virtual std::vector<std::vector<double>> _getSingleVesicleSurfaceSpecPos(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        spec_global_id sidx);

    virtual std::vector<std::vector<double>> _getSingleVesicleSurfaceSpecPosSpherical(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        spec_global_id sidx) const;

    virtual void _setSingleVesicleSurfaceSpecPosSpherical(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        spec_global_id sidx,
        const std::vector<std::vector<double>>& pos_spherical);

    virtual std::vector<double> _getSingleSpecPosSpherical(
        spec_global_id sidx,
        pointspec_individual_id ps_unique_id) const;

    virtual void _setSingleVesicleInnerSpecCount(vesicle_global_id vidx,
                                                 vesicle_individual_id ves_unique_index,
                                                 spec_global_id sidx,
                                                 uint c);

    virtual std::map<vesicle_individual_id, uint> _getCompVesicleSurfaceLinkSpecCountMap(
        comp_global_id cidx,
        vesicle_global_id vidx,
        linkspec_global_id lsidx) const;

    virtual uint _getCompVesicleSurfaceLinkSpecCount(comp_global_id cidx,
                                                     vesicle_global_id vidx,
                                                     linkspec_global_id lsidx) const;

    virtual std::vector<std::vector<double>> _getSingleVesicleSurfaceLinkSpecPos(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index,
        linkspec_global_id lsidx);

    virtual std::vector<double> _getSingleLinkSpecPos(linkspec_individual_id ls_unique_id) const;

    virtual linkspec_individual_id _getSingleLinkSpecLinkedTo(
        linkspec_individual_id ls_unique_id) const;

    virtual vesicle_individual_id _getSingleLinkSpecVes(linkspec_individual_id ls_unique_id) const;

    virtual uint _getSingleVesicleImmobility(vesicle_global_id vidx,
                                             vesicle_individual_id ves_unique_index) const;

    virtual std::vector<tetrahedron_global_id> _getSingleVesicleOverlapTets(
        vesicle_global_id vidx,
        vesicle_individual_id ves_unique_index) const;

    virtual void _setTetVesicleDcst(tetrahedron_global_id tidx,
                                    vesicle_global_id vidx,
                                    double dcst);

    //////////////// NOT COMPARTMENT-SPECIFIC ///////////////////////////////

    virtual void _setVesicleSurfaceLinkSpecSDiffD(vesicle_global_id vidx,
                                                  linkspec_global_id lsidx,
                                                  double dcst);

    virtual void _setVesSReacK(vessreac_global_id vsridx, double kf);
    virtual uint _getVesSReacExtent(vessreac_global_id vsridx) const;

    virtual void _setExocytosisK(exocytosis_global_id exoidx, double kf);

    virtual uint _getExocytosisExtent(exocytosis_global_id exoidx) const;

    virtual std::vector<ExocytosisEvent> _getExocytosisEvents(exocytosis_global_id exoidx);

    virtual uint _getRaftEndocytosisExtent(raftendocytosis_global_id rendoidx) const;

    virtual std::vector<RaftEndocytosisEvent> _getRaftEndocytosisEvents(
        raftendocytosis_global_id rendoidx);

    virtual void _setRaftEndocytosisK(solver::raftendocytosis_global_id rendoidx, double kcst);

    // virtual void _setCompVesicleSpecDiffD(comp_global_id cidx,
    // vesicle_global_id vidx, spec_global_id sidx, double d) const;

    virtual void _addVesicleDiffusionGroup(vesicle_global_id vidx,
                                           const std::vector<comp_global_id>& comp_indices);

    virtual void _addPathVesicle(std::string const& path_name,
                                 vesicle_global_id vidx,
                                 double speed,
                                 const std::map<spec_global_id, uint>& spec_deps,
                                 const std::vector<double>& stoch_stepsize);

    ////////////////////////////////////////////////////////////////////////

    virtual unsigned long long _getCompComplexReacExtent(comp_global_id cidx,
                                                         complexreac_global_id ridx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    virtual double _getTetVol(tetrahedron_global_id tidx) const;
    virtual void _setTetVol(tetrahedron_global_id tidx, double vol);

    virtual bool _getTetSpecDefined(tetrahedron_global_id tidx, spec_global_id sidx) const;

    virtual double _getTetSpecCount(tetrahedron_global_id tidx, spec_global_id sidx) const;
    virtual void _setTetSpecCount(tetrahedron_global_id tidx, spec_global_id sidx, double n);

    virtual double _getTetSpecAmount(tetrahedron_global_id tidx, spec_global_id sidx) const;
    virtual void _setTetSpecAmount(tetrahedron_global_id tidx, spec_global_id sidx, double m);

    virtual double _getTetSpecConc(tetrahedron_global_id tidx, spec_global_id sidx) const;
    virtual void _setTetSpecConc(tetrahedron_global_id tidx, spec_global_id sidx, double c);

    virtual bool _getTetSpecClamped(tetrahedron_global_id tidx, spec_global_id sidx) const;
    virtual void _setTetSpecClamped(tetrahedron_global_id tidx, spec_global_id sidx, bool buf);

    virtual double _getTetReacK(tetrahedron_global_id tidx, reac_global_id ridx) const;
    virtual void _setTetReacK(tetrahedron_global_id tidx, reac_global_id ridx, double kf);

    virtual bool _getTetReacActive(tetrahedron_global_id tidx, reac_global_id ridx) const;
    virtual void _setTetReacActive(tetrahedron_global_id tidx, reac_global_id ridx, bool act);

    virtual double _getTetDiffD(tetrahedron_global_id tidx,
                                diff_global_id didx,
                                tetrahedron_global_id direction_tet = {}) const;
    virtual void _setTetDiffD(tetrahedron_global_id tidx,
                              diff_global_id didx,
                              double dk,
                              tetrahedron_global_id direction_tet = {});

    virtual bool _getTetDiffActive(tetrahedron_global_id tidx, diff_global_id didx) const;
    virtual void _setTetDiffActive(tetrahedron_global_id tidx, diff_global_id didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getTetReacH(tetrahedron_global_id tidx, reac_global_id ridx) const;
    virtual double _getTetReacC(tetrahedron_global_id tidx, reac_global_id ridx) const;
    virtual double _getTetReacA(tetrahedron_global_id tidx, reac_global_id ridx) const;

    virtual double _getTetDiffA(tetrahedron_global_id tidx, diff_global_id didx) const;

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    virtual double _getTetV(tetrahedron_global_id tidx) const;
    virtual void _setTetV(tetrahedron_global_id tidx, double v);
    virtual bool _getTetVClamped(tetrahedron_global_id tidx) const;
    virtual void _setTetVClamped(tetrahedron_global_id tidx, bool cl);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    virtual double _getPatchArea(patch_global_id pidx) const = 0;
    virtual void _setPatchArea(patch_global_id pidx, double area);

    virtual double _getPatchSpecCount(patch_global_id pidx, spec_global_id sidx) const = 0;
    virtual void _setPatchSpecCount(patch_global_id pidx, spec_global_id sidx, double n) = 0;

    virtual double _getPatchSpecAmount(patch_global_id pidx, spec_global_id sidx) const = 0;
    virtual void _setPatchSpecAmount(patch_global_id pidx, spec_global_id sidx, double a) = 0;

    virtual bool _getPatchSpecClamped(patch_global_id pidx, spec_global_id sidx) const = 0;
    virtual void _setPatchSpecClamped(patch_global_id pidx, spec_global_id sidx, bool buf) = 0;

    virtual double _getPatchComplexCount(
        patch_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f)
        const;
    virtual void _setPatchComplexCount(patch_global_id cidx,
                                       complex_global_id sidx,
                                       const util::strongid_vector<complex_substate_id, uint>& f,
                                       double n);

    virtual double _getPatchComplexAmount(
        patch_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f)
        const;
    virtual void _setPatchComplexAmount(patch_global_id cidx,
                                        complex_global_id sidx,
                                        const util::strongid_vector<complex_substate_id, uint>& f,
                                        double a);

    virtual double _getPatchComplexSUSCount(
        patch_global_id cidx,
        complex_global_id sidx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& f,
        complex_substate_id m) const;

    virtual double _getPatchSReacK(patch_global_id pidx, sreac_global_id sridx) const = 0;
    virtual void _setPatchSReacK(patch_global_id pidx, sreac_global_id sridx, double kf) = 0;

    virtual bool _getPatchSReacActive(patch_global_id pidx, sreac_global_id ridx) const = 0;
    virtual void _setPatchSReacActive(patch_global_id pidx, sreac_global_id ridx, bool a) = 0;

    virtual bool _getPatchVDepSReacActive(patch_global_id pidx, vdepsreac_global_id vsridx) const;
    virtual void _setPatchVDepSReacActive(patch_global_id pidx, vdepsreac_global_id vsridx, bool a);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getPatchSReacH(patch_global_id pidx, sreac_global_id ridx) const;
    virtual double _getPatchSReacC(patch_global_id pidx, sreac_global_id ridx) const;
    virtual double _getPatchSReacA(patch_global_id pidx, sreac_global_id ridx) const;

    virtual unsigned long long _getPatchSReacExtent(patch_global_id pidx,
                                                    sreac_global_id ridx) const;
    virtual void _resetPatchSReacExtent(patch_global_id pidx, sreac_global_id ridx);

    virtual unsigned long long _getPatchComplexSReacExtent(patch_global_id pidx,
                                                           complexsreac_global_id ridx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      RAFT/ENDOCYTOSIS -RELATED
    ////////////////////////////////////////////////////////////////////////

    virtual uint _getPatchRaftCount(patch_global_id pidx, raft_global_id ridx) const;
    virtual void _setPatchRaftCount(patch_global_id pidx, raft_global_id ridx, uint n);

    virtual std::vector<double> _getSingleRaftPos(raft_global_id ridx,
                                                  raft_individual_id raft_unique_index) const;

    virtual std::map<solver::raft_individual_id, uint>
    _getPatchRaftSpecCountMap(patch_global_id pidx, raft_global_id ridx, spec_global_id sidx) const;

    virtual uint _getPatchRaftSpecCount(patch_global_id pidx,
                                        raft_global_id ridx,
                                        spec_global_id sidx) const;

    virtual uint _getSingleRaftSpecCount(raft_global_id ridx,
                                         raft_individual_id raft_unique_index,
                                         spec_global_id sidx) const;

    virtual void _setSingleRaftSpecCount(raft_global_id ridx,
                                         raft_individual_id raft_unique_index,
                                         spec_global_id sidx,
                                         uint c);

    virtual uint _getSingleRaftImmobility(raft_global_id ridx,
                                          raft_individual_id raft_unique_index) const;

    virtual double _getSingleRaftRaftEndocytosisK(raft_global_id ridx,
                                                  raft_individual_id raft_unique_index,
                                                  raftendocytosis_global_id rendoidx) const;

    virtual void _setSingleRaftRaftEndocytosisK(raft_global_id ridx,
                                                raft_individual_id raft_unique_index,
                                                raftendocytosis_global_id rendoidx,
                                                double k);

    virtual void _setSingleRaftSReacActive(raft_global_id ridx,
                                           raft_individual_id raft_unique_index,
                                           raftsreac_global_id rsreacidx,
                                           bool active);

    virtual bool _getSingleRaftSReacActive(raft_global_id ridx,
                                           raft_individual_id raft_unique_index,
                                           raftsreac_global_id rsreacidx) const;

    virtual std::vector<raft_individual_id> _getPatchRaftIndices(patch_global_id pidx,
                                                                 raft_global_id ridx) const;

    virtual patch_global_id _getSingleRaftPatch(raft_global_id ridx,
                                                raft_individual_id raft_unique_index) const;

    virtual void _setPatchEndocyticZoneEndocytosisActive(patch_global_id pidx,
                                                         std::string const& zone,
                                                         endocytosis_global_id endogidx,
                                                         bool active);

    virtual void _setPatchEndocyticZoneEndocytosisK(patch_global_id pidx,
                                                    std::string const& zone,
                                                    endocytosis_global_id endogidx,
                                                    double k);

    virtual uint _getPatchEndocyticZoneEndocytosisExtent(patch_global_id pidx,
                                                         std::string const& zone,
                                                         endocytosis_global_id endogidx) const;

    virtual std::vector<EndocytosisEvent> _getPatchEndocyticZoneEndocytosisEvents(
        patch_global_id pidx,
        std::string const& zone,
        endocytosis_global_id endogidx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    virtual void _setDiffBoundarySpecDiffusionActive(diffboundary_global_id dbidx,
                                                     spec_global_id sidx,
                                                     bool act);
    virtual bool _getDiffBoundarySpecDiffusionActive(diffboundary_global_id dbidx,
                                                     spec_global_id sidx) const;
    virtual void _setDiffBoundarySpecDcst(diffboundary_global_id dbidx,
                                          spec_global_id sidx,
                                          double dcst,
                                          comp_global_id direction_comp = {});

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      SURFACE DIFFUSION BOUNDARIES
    ////////////////////////////////////////////////////////////////////////

    virtual void _setSDiffBoundarySpecDiffusionActive(sdiffboundary_global_id sdbidx,
                                                      spec_global_id sidx,
                                                      bool act);
    virtual bool _getSDiffBoundarySpecDiffusionActive(sdiffboundary_global_id sdbidx,
                                                      spec_global_id sidx) const;
    virtual void _setSDiffBoundarySpecDcst(sdiffboundary_global_id sdbidx,
                                           spec_global_id sidx,
                                           double dcst,
                                           patch_global_id direction_patch = {});

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    virtual double _getTriArea(triangle_global_id tidx) const;
    virtual void _setTriArea(triangle_global_id tidx, double area);

    virtual bool _getTriSpecDefined(triangle_global_id tidx, spec_global_id sidx) const;

    virtual double _getTriSpecCount(triangle_global_id tidx, spec_global_id sidx) const;
    virtual void _setTriSpecCount(triangle_global_id tidx, spec_global_id sidx, double n);

    virtual double _getTriSpecAmount(triangle_global_id tidx, spec_global_id sidx) const;
    virtual void _setTriSpecAmount(triangle_global_id tidx, spec_global_id sidx, double m);

    virtual bool _getTriSpecClamped(triangle_global_id tidx, spec_global_id sidx) const;
    virtual void _setTriSpecClamped(triangle_global_id tidx, spec_global_id sidx, bool buf);

    virtual double _getTriSReacK(triangle_global_id tidx, sreac_global_id ridx) const;
    virtual void _setTriSReacK(triangle_global_id tidx, sreac_global_id ridx, double kf);

    virtual double _getTriSDiffD(triangle_global_id tidx,
                                 surfdiff_global_id didx,
                                 triangle_global_id direction_tri) const;

    virtual void _setTriSDiffD(triangle_global_id, surfdiff_global_id, double, triangle_global_id);

    virtual bool _getTriSReacActive(triangle_global_id tidx, sreac_global_id ridx) const;
    virtual void _setTriSReacActive(triangle_global_id tidx, sreac_global_id ridx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getTriSReacH(triangle_global_id tidx, sreac_global_id ridx) const;
    virtual double _getTriSReacC(triangle_global_id tidx, sreac_global_id ridx) const;
    virtual double _getTriSReacA(triangle_global_id tidx, sreac_global_id ridx) const;

    ////////////////////////// RAFT/VESICLE -RELATED ////////////////////////////

    virtual uint _getTriRaftCount(triangle_global_id tidx, raft_global_id ridx) const;
    virtual void _setTriRaftCount(triangle_global_id tidx, raft_global_id ridx, uint n);
    virtual solver::raft_individual_id _addTriRaft(triangle_global_id tidx, raft_global_id ridx);

    virtual bool _getTriExocytosisActive(triangle_global_id tidx, exocytosis_global_id eidx) const;
    virtual void _setTriExocytosisActive(triangle_global_id tidx,
                                         exocytosis_global_id eidx,
                                         bool act);

    ////////////////////////// ADDED FOR EFIELD ////////////////////////////

    virtual double _getTriV(triangle_global_id tidx) const;
    virtual void _setTriV(triangle_global_id tidx, double v);
    virtual bool _getTriVClamped(triangle_global_id tidx) const;
    virtual void _setTriVClamped(triangle_global_id tidx, bool cl);

    virtual double _getTriOhmicErev(triangle_global_id tidx, ohmiccurr_global_id ocidx) const;
    virtual void _setTriOhmicErev(triangle_global_id tidx, ohmiccurr_global_id ocidx, double erev);
    virtual double _getTriOhmicI(triangle_global_id tidx) const;
    virtual double _getTriOhmicI(triangle_global_id tidx, ohmiccurr_global_id ocidx) const;
    virtual double _getTriI(triangle_global_id tidx) const;

    virtual double _getTriIClamp(triangle_global_id tidx) const;
    virtual void _setTriIClamp(triangle_global_id tidx, double i);

    virtual double _getTriGHKI(triangle_global_id tidx) const;
    virtual double _getTriGHKI(triangle_global_id tidx, ghkcurr_global_id ocidx) const;

    virtual double _getTriVDepSReacK(triangle_global_id tidx, vdepsreac_global_id vsridx) const;
    virtual bool _getTriVDepSReacActive(triangle_global_id tidx, vdepsreac_global_id vsridx) const;
    virtual void _setTriVDepSReacActive(triangle_global_id tidx,
                                        vdepsreac_global_id vsridx,
                                        bool act);

    virtual void _setTriCapac(triangle_global_id tidx, double cm);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      VERTICES ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    virtual double _getVertV(vertex_id_t vidx) const;
    virtual void _setVertV(vertex_id_t vidx, double v);

    virtual bool _getVertVClamped(vertex_id_t vidx) const;
    virtual void _setVertVClamped(vertex_id_t vidx, bool cl);

    virtual double _getVertIClamp(vertex_id_t vidx) const;
    virtual void _setVertIClamp(vertex_id_t vidx, double i);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROL:
    //      MEMBRANES
    ////////////////////////////////////////////////////////////////////////

    virtual void _setMembPotential(membrane_global_id midx, double v);
    virtual void _setMembCapac(membrane_global_id midx, double cm);
    virtual void _setMembVolRes(membrane_global_id midx, double ro);
    virtual void _setMembRes(membrane_global_id midx, double ro, double vrev);
    virtual std::pair<double, double> _getMembRes(membrane_global_id midx) const;

    ////////////////////////////////////////////////////////////////////////

  public:
    /// Return a reference of the Model object.
    inline model::Model& model() const noexcept {
        return pModel;
    }

    /// Return a reference of the Geom object.
    inline wm::Geom& geom() const noexcept {
        return pGeom;
    }

    /// Return a reference of the RNG object
    inline const rng::RNGptr& rng() const noexcept {
        return pRNG;
    }

    /// Return a reference of the Statedef object.
    inline const solver::Statedef& statedef() const noexcept {
        return *pStatedef;
    }

    inline solver::Statedef& statedef() noexcept {
        return *pStatedef;
    }

  private:
    model::Model& pModel;

    wm::Geom& pGeom;

    const rng::RNGptr pRNG;

    std::unique_ptr<Statedef> pStatedef;
};

std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>
_convertComplexFilters(const std::vector<std::vector<model::SubunitStateFilter>>& f);

util::strongid_vector<complex_substate_id, uint> _convertComplexState(
    const std::vector<std::vector<model::SubunitStateFilter>>& f);

}  // namespace steps::solver

////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SOLVER_API_HPP
#define STEPS_SOLVER_API_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/geom/geom.hpp>
#include <steps/model/model.hpp>
#include <steps/rng/rng.hpp>


////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forward declarations
class Statedef;

////////////////////////////////////////////////////////////////////////////////
/// API class for a solver.
///
/// The API class provides all APIs for each solver.
///
/// \warning Not every API in this class is implemented in the solver.
///          If a API is not implemented in a solver, 
///          STEPS will throw an error message to the user.
///////////////////////////////////////////////////////////////////////////////
class API
{

public:

    /// Constructor
    ///
    /// \param m Pointer to the model.
    /// \param g Pointer to the geometr container.
    /// \param r Pointer to the random number generator.
    API(steps::model::Model * m, steps::wm::Geom * g, steps::rng::RNG * r);

    /// Destructor
    ///
    virtual ~API(void);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Return the solver's name.
    virtual std::string getSolverName(void) const = 0;
    /// Return the solver's description.
    virtual std::string getSolverDesc(void) const = 0;
    /// Return the solver's author.
    virtual std::string getSolverAuthors(void) const = 0;
    /// Return the solver's email.
    virtual std::string getSolverEmail(void) const = 0;

    ////////////////////////////////////////////////////////////////////////
    // SIMULATION CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    /// Save the solver's state in a file.
    ///
    /// \param filename The file name to save the solver's state.
    virtual void saveState(std::string const & filename);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER CONTROLS
    ////////////////////////////////////////////////////////////////////////

    /// Reset the solver.
    virtual void reset(void) = 0;

    /// Run the solver until a given end time.
    ///
    /// \param endtime Time to end the solver.
    virtual void run(double endtime) = 0;

    /// Run the solver for a step.
    virtual void step(void);

    /// Set DT of the solver.
    ///
    /// \param dt Dt.
    /// \todo Ask Iain about this.
    virtual void setDT(double dt);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      GENERAL
    ////////////////////////////////////////////////////////////////////////
    /// Return the time.
    virtual double getTime(void) const = 0;
    /// Return the DT
    /// \todo ask iain
    virtual double getDT(void) const;
    /// Return the A0
    /// \todo ask iain
    virtual double getA0(void) const;
    /// Return the number of steps.
    virtual uint getNSteps(void) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Returns the volume of compartment c (in m^3).
    ///
    /// \param c Name of the compartment.
    double getCompVol(std::string const & c) const;

    // Sets the volume of compartment c.
    // NOTE: this method may throw an exception if this does not make sense
    // for a solver (e.g. a tetrahedral mesh-based solver).
    void setCompVol(std::string const & c, double vol);

    // Returns the number of molecules of species s in compartment c.
    // NOTE: in a mesh-based simulation, the total count is computed as
    // the sum of the counts in all voxels of the compartment.
    double getCompCount(std::string const & c, std::string const & s) const;

    // Sets the number of molecules of species s in compartment c.
    // NOTE: in a mesh-based simulation, the total amount is equally divided
    // over all voxels in the compartment (i.e. a uniform distribution).
    void setCompCount(std::string const & c, std::string const & s, double n);

    // Returns the amount (in mols) of species s in compartment c.
    // NOTE: in a mesh-based simulation, the total amount is computed as
    // the sum of the amounts in all voxels of the compartment.
    double getCompAmount(std::string const & c, std::string const & s) const;

    // Set the amount (in mols) of species s in compartment c.
    // NOTE: in a mesh-based simulation, the total amount is equally divided
    // over all voxels in the compartment (i.e. a uniform distribution).
    void setCompAmount(std::string const & c, std::string const & s, double a);

    // Returns the concentration (in molar units) of species s in compartment c.
    // NOTE: in a mesh-based simulation, the overall concentration in a
    // compartment is computed by taking the volume-weighted sum of the
    // concentration in all voxels of the compartment.
    double getCompConc(std::string const & c, std::string const & s) const;

    // Sets the concentration (in molar units) of species s in compartment c.
    // NOTE: in a mesh-based simulation, this method changes the
    // concentration to the same value in all voxels of the compartment.
    void setCompConc(std::string const & c, std::string const & s, double conc);

    // Returns whether the concentration of species s in compartment c
    // remains constant over time (unless changed explicitly).
    // NOTE: in a mesh-based simulation, this method will only return true
    // only if the species has been clamped in all voxels of the compartment.
    bool getCompClamped(std::string const & c, std::string const & s) const;

    // Turns clamping of species s in compartment c on or off.
    // NOTE: in a mesh based simulation, this method turns clamping on/off
    // in all voxels of the compartment.
    void setCompClamped(std::string const & c, std::string const & s, bool b);

    // Returns the macroscopic reaction constant of reaction r in
    // compartment c.
    // Note: in a mesh-based simulation, the value is computed as the
    // volume-weighted sum of the reaction constants in all voxels of the
    // compartment.
    double getCompReacK(std::string const & c, std::string const & r) const;

    // Sets the macroscopic reaction constant of reaction r in compartment c
    // (units vary according to the order of the reaction).
    // NOTE: in a mesh-based simulation, this method changes the reaction
    // constant equally in all voxels of the compartment.
    void setCompReacK(std::string const & c, std::string const & r, double kf);

    // Returns whether reaction r in compartment c is active or not
    // NOTE: in a mesh-based simulation, this method returns false only when
    // the reaction has been inactivated in all voxels.
    bool getCompReacActive(std::string const & c, std::string const & r) const;

    // Activate or inactivate reaction r in compartment c.
    // NOTE: in a mesh-based simulation, activation/inactivation of a reaction
    // turns it on/off in all voxels.
    void setCompReacActive(std::string const & c, std::string const & r, bool a);

    // Returns the diffusion constant of diffusion rule d in compartment c.
    double getCompDiffD(std::string const & c, std::string const & d) const;

    // Set the diffusion constant of diffusion rule d in compartment c.
    void setCompDiffD(std::string const & c, std::string const & d, double dcst);

    // Returns whether diffusion rule d in compartment c is active or not.
    bool getCompDiffActive(std::string const & c, std::string const & d) const;

    // Activate or deactivate diffusion rule d in compartment c.
    void setCompDiffActive(std::string const & c, std::string const & d, bool act);

    ////////////////////////////////////////////////////////////////////////

    // Returns c_mu, the mesoscopic reaction constant of reaction r in
    // compartment c.
    // NOTE: in a mesh-based simulation, the mesoscopic reaction constant is
    // computed as the sum of the mesoscopic constants in all voxels of the
    // compartment.
    double getCompReacC(std::string const & c, std::string const & r) const;

    // Returns h_mu, the distinct number of ways in which reaction r can
    // occur in compartment c, by computing the product of its reactants.
    // NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
    // over all voxels of the compartment. This can become a very large value.
    double getCompReacH(std::string const & c, std::string const & r) const;

    // Returns the propensity, a_mu, of reaction r in compartment c.
    // The propensity value gives the probability per unit time that this
    // reaction will occur in the current state.
    // NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
    // a_mu in all voxels of the compartment.
    double getCompReacA(std::string const & c, std::string const & r) const;

    // Returns the extent of reaction r in compartment c.
    // NOTE: in a mesh-based simulation, returns the sum of the reaction
    // extents in all voxels of the compartment.
    uint getCompReacExtent(std::string const & c, std::string const & r) const;

    // Resets the extent of reaction r in compartment c to zero.
    // NOTE: in a mesh-based simulation, resets the extents of the reaction
    // in all voxels of the compartment.
    void resetCompReacExtent(std::string const & c, std::string const & r);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    // Returns the volume of a tetrahedron (in m^3).
    double getTetVol(uint tidx) const;
																							////////////
    void setTetVol(uint tidx, double vol);

    // Returns the number of molecules of species s in a tetrahedral volume
    // element (voxel).
    double getTetCount(uint tidx, std::string const & s) const;

    // Sets the number of molecules of species s in a voxel.
    void setTetCount(uint tidx, std::string const & s, double n);

    // Returns the amount (in mols) of species s in a voxel
    double getTetAmount(uint tidx, std::string const & s) const;

    // Sets the amount (in mols) of species s in a voxel.
    void setTetAmount(uint tidx, std::string const & s, double m);

    // Returns the concentration (in molar units) of species s in a voxel..
    double getTetConc(uint tidx, std::string const & s) const;

    // Sets the concentration (in molar units) of species s in a voxel.
    void setTetConc(uint tidx, std::string const & s, double c);

    // Returns whether the concentration of species s in a voxel
    // remains constant over time (unless changed explicitly).
    bool getTetClamped(uint tidx, std::string const & s) const;

    // Sets clamping of species s in a voxel on or off.
    void setTetClamped(uint tidx, std::string const & s, bool buf);

    // Returns the macroscopic reaction constant of reaction r in a voxel
    // (units vary with order of reaction).
    double getTetReacK(uint tidx, std::string const & r) const;

    // Sets the macroscopic reaction constant of reaction r in a voxel.
    void setTetReacK(uint tidx, std::string const & r, double kf);

    // Returns whether reaction r in a voxel is active or not
    bool getTetReacActive(uint tidx, std::string const & r) const;

    // Activates/inactivates reaction r in a voxel.
    void setTetReacActive(uint tidx, std::string const & r, bool act);

    // Returns the diffusion constant of diffusion rule d in a voxel.
    double getTetDiffD(uint tidx, std::string const & d) const;

    // Sets the diffusion constant of diffusion rule d in a voxel.
    void setTetDiffD(uint tidx, std::string const & d, double dk);

    // Returns whether diffusion rule d in a voxel is active or not.
    bool getTetDiffActive(uint tidx, std::string const & d) const;

    // Activates/deactivates diffusion rule d in a voxel.
    void setTetDiffActive(uint tidx, std::string const & d, bool act);

    ////////////////////////////////////////////////////////////////////////

    // Returns c_mu, the mesoscopic reaction constant of reaction r in
    // a voxel
    double getTetReacC(uint tidx, std::string const & r) const;

    // Returns h_mu, the distinct number of ways in which reaction r can
    // occur in a voxel, by computing the product of its reactants.
    double getTetReacH(uint tidx, std::string const & r) const;

    // Returns the propensity, a_mu, of reaction r in a voxel.
    // The propensity value gives the probability per unit time that this
    // reaction will occur in the current state.
    double getTetReacA(uint tidx, std::string const & r) const;

    // Returns the propensity, a_mu of diffusion rule d in a voxel.
    double getTetDiffA(uint tidx, std::string const & d) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    // Returns the area of patch p (in m^2)
    double getPatchArea(std::string const & p) const;

    // Sets the area of patch p.
    // NOTE: this method may throw an exception if this does not make sense
    // for a solver (e.g. a tetrahedral mesh-based solver).
    void setPatchArea(std::string const & p, double area);

    // Returns the number of molecules of species s in patch p.
    // NOTE: in a mesh-based simulation, the total count is computed as
    // the sum of the counts in all triangles of the patch.
    double getPatchCount(std::string const & p, std::string const & s) const;

    // Sets the number of molecules of species s in patch p.
    // NOTE: in a mesh-based simulation, the total amount is equally divided
    // over all triangles in the patch.
    void setPatchCount(std::string const & p, std::string const & s, double n);

    // Returns the amount (in mols) of species s in patch p.
    // NOTE: in a mesh-based simulation, the total amount is computed as
    // the sum of the amounts in all triangles of the patch.
    double getPatchAmount(std::string const & p, std::string const & s) const;

    // Sets the amount (in mols) of species s in patch p.
    // NOTE: in a mesh-based simulation, the total amount is equally divided
    // over all triangles in the patch.
    void setPatchAmount(std::string const & p, std::string const & s, double a);

    // Returns whether the count of species s in patch p remains constant
    // over time (unless changed explicitly).
    // NOTE: this method will only return true if the species has been
    // clamped in all triangles in the patch.
    bool getPatchClamped(std::string const & p, std::string const & s) const;

    // Turns clamping of species in patch p on or off.
    // NOTE: in a mesh-based simulation, this method turns clamping on/off
    // in all triangles in the patch.
    void setPatchClamped(std::string const & p, std::string const & s, bool buf);

    // Returns the macroscopic reaction constant of surface reaction r
    // in patch p.
    // NOTE: in a mesh-based simulation, the value is computed as the
    // area-weighted sum of the reaction constants in all triangles of
    // the patch.
    double getPatchSReacK(std::string const & p, std::string const & r) const;

    // Sets the macroscopic reaction constant of surface reaction r
    // in patch p.
    // NOTE: in a mesh-based simulation this method changes the reaction
    // constant equally in all triangles of the patch.
    void setPatchSReacK(std::string const & p, std::string const & r, double kf);

    // Returns whether surface reaction r in patch p is active or not.
    // NOTE: in a mesh-based simulation, only returns false when the
    // reaction has been inactivated in all triangles.
    bool getPatchSReacActive(std::string const & p, std::string const & r) const;

    // Activate or inactivate surface reaction r in patch p.
    // NOTE: in a mesh-based simulation, activation/inactivation of a
    // surface reaction turns it on/off in all triangles.
    void setPatchSReacActive(std::string const & p, std::string const & r, bool a);

    ////////////////////////////////////////////////////////////////////////

    // Returns c_mu, the mesoscopic reaction constant of surface reaction r
    // in patch p.
    // NOTE: in a mesh_based simulation, the mesoscopic reaction constant
    // is computed as the sum of the mesoscopic reaction constants from all
    // triangles of the patch.
    double getPatchSReacC(std::string const & p, std::string const & r) const;

    // Returns h_mu, the distinct number of ways in which a surface reaction
    // r can occur in patch p, by computing the product of its reactants.
    // NOTE: in a mesh-based simulation, it returns the sum of the h_mu's
    // over all triangles triangles of the patch.
    double getPatchSReacH(std::string const & p, std::string const & r) const;

    // Returns the propensity, a_mu of surface reaction r in patch p.
    // This propensity value gives the probability per unit time that this
    // surface reaction will occur in the current state.
    // NOTE: in a mesh-based simulation, a_mu is computed as the sum of the
    // a_mu in all triangles in the patch.
    double getPatchSReacA(std::string const & p, std::string const & r) const;

    // Returns the extent of surface reaction r in patch p.
    // NOTE: in a mesh-based simulation, returns the sum of the extents in
    // all triangles of the patch.
    uint getPatchSReacExtent(std::string const & p, std::string const & r) const;

    // Resets the extent of surface reaction r in patch p to zero.
    // NOTE: in a mesh-based simulation, resets the extents of the
    // surface reaction in all triangles of the patch.
    void resetPatchSReacExtent(std::string const & p, std::string const & r);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    // Returns the area of the triangle (in m^2).
    double getTriArea(uint tidx) const;

    //
    void setTriArea(uint tidx, double area);

    // Returns the number of molecules of species s in a triangle.
    double getTriCount(uint tidx, std::string const & s) const;

    // Sets the number of molecules of species s in a triangle.
    void setTriCount(uint tidx, std::string const & s, double n);

    // Returns whether the number of molecules of species s in a triangle
    // remains constant over time (unless changed explicitly)
    bool getTriClamped(uint tidx, std::string const & s) const;

    // Sets clamping of species s in a triangle on or off.
    void setTriClamped(uint tidx, std::string const & s, bool buf);

    // Returns the macroscopic reaction constant of surface reaction r
    // in a triangle (units vary with order of reaction).
    double getTriSReacK(uint tidx, std::string const & r) const;

    // Sets the macroscopic reaction constant of surface reaction r in
    // a triangle.
    void setTriSReacK(uint tidx, std::string const & r, double kf);

    // Returns whether surface reaction r in a triangle is active or not.
    bool getTriSReacActive(uint tidx, std::string const & r) const;

    // Activates/inactivates surface reaction r in a triangle.
    void setTriSReacActive(uint tidx, std::string const & r, bool act);

    ////////////////////////////////////////////////////////////////////////

    // Returns c_mu, the mesoscopic reaction constant of surface reaction r
    // in a triangle.
    double getTriSReacC(uint tidx, std::string const & r) const;

    // Returns h_mu, the distinct number of ways in which surface reaction r
    // can occur in a triangle, by computing the product of it's reactants.
    double getTriSReacH(uint tidx, std::string const & r) const;

    // Returns the propensity, a_mu, of surface reaction r in a triangle.
    // The propensity value gives the probability per unit time that this
    // surface reaction will occur in the current state.
    double getTriSReacA(uint tidx, std::string const & r) const;

    ////////////////////////////////////////////////////////////////////////

protected:

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    virtual double _getCompVol(uint cidx) const = 0;
    virtual void _setCompVol(uint cidx, double vol);

    virtual double _getCompCount(uint cidx, uint sidx) const = 0;
    virtual void _setCompCount(uint cidx, uint sidx, double n) = 0;

    virtual double _getCompAmount(uint cidx, uint sidx) const = 0;
    virtual void _setCompAmount(uint cidx, uint sidx, double a) = 0;

    virtual double _getCompConc(uint cidx, uint sidx) const = 0;
    virtual void _setCompConc(uint cidx, uint sidx, double c) = 0;

    virtual bool _getCompClamped(uint cidx, uint sidx) const = 0;
    virtual void _setCompClamped(uint cidx, uint sidx, bool b) = 0;

    virtual double _getCompReacK(uint cidx, uint ridx) const = 0;
    virtual void _setCompReacK(uint cidx, uint ridx, double kf) = 0;

    virtual bool _getCompReacActive(uint cidx, uint ridx) const = 0;
    virtual void _setCompReacActive(uint cidx, uint ridx, bool a) = 0;

    virtual double _getCompDiffD(uint cidx, uint didx) const;
    virtual void _setCompDiffD(uint cidx, uint didx, double dcst);

    virtual bool _getCompDiffActive(uint cidx, uint didx) const;
    virtual void _setCompDiffActive(uint cidx, uint didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getCompReacH(uint cidx, uint ridx) const;
    virtual double _getCompReacC(uint cidx, uint ridx) const;
    virtual double _getCompReacA(uint cidx, uint ridx) const;

    virtual uint _getCompReacExtent(uint cidx, uint ridx) const;
    virtual void _resetCompReacExtent(uint cidx, uint ridx);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TETRAHEDRAL VOLUME ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    virtual double _getTetVol(uint tidx) const;
    virtual void _setTetVol(uint tidx, double vol);

    virtual double _getTetCount(uint tidx, uint sidx) const;
    virtual void _setTetCount(uint tidx, uint sidx, double n);

    virtual double _getTetAmount(uint tidx, uint sidx) const;
    virtual void _setTetAmount(uint tidx, uint sidx, double m);

    virtual double _getTetConc(uint tidx, uint sidx) const;
    virtual void _setTetConc(uint tidx, uint sidx, double c);

    virtual bool _getTetClamped(uint tidx, uint sidx) const;
    virtual void _setTetClamped(uint tidx, uint sidx, bool buf);

    virtual double _getTetReacK(uint tidx, uint ridx) const;
    virtual void _setTetReacK(uint tidx, uint ridx, double kf);

    virtual bool _getTetReacActive(uint tidx, uint ridx) const;
    virtual void _setTetReacActive(uint tidx, uint ridx, bool act);

    virtual double _getTetDiffD(uint tidx, uint didx) const;
    virtual void _setTetDiffD(uint tidx, uint didx, double dk);

    virtual bool _getTetDiffActive(uint tidx, uint didx) const;
    virtual void _setTetDiffActive(uint tidx, uint didx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getTetReacH(uint tidx, uint ridx) const;
    virtual double _getTetReacC(uint tidx, uint ridx) const;
    virtual double _getTetReacA(uint tidx, uint ridx) const;

    virtual double _getTetDiffA(uint tidx, uint didx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      PATCH
    ////////////////////////////////////////////////////////////////////////

    virtual double _getPatchArea(uint pidx) const = 0;
    virtual void _setPatchArea(uint pidx, double area);

    virtual double _getPatchCount(uint pidx, uint sidx) const = 0;
    virtual void _setPatchCount(uint pidx, uint sidx, double n) = 0;

    virtual double _getPatchAmount(uint pidx, uint sidx) const = 0;
    virtual void _setPatchAmount(uint pidx, uint sidx, double a) = 0;

    virtual bool _getPatchClamped(uint pidx, uint sidx) const = 0;
    virtual void _setPatchClamped(uint pidx, uint sidx, bool buf) = 0;

    virtual double _getPatchSReacK(uint pidx, uint ridx) const = 0;
    virtual void _setPatchSReacK(uint pidx, uint ridx, double kf) = 0;

    virtual bool _getPatchSReacActive(uint pidx, uint ridx) const = 0;
    virtual void _setPatchSReacActive(uint pidx, uint ridx, bool a) = 0;

    ////////////////////////////////////////////////////////////////////////

    virtual double _getPatchSReacH(uint pidx, uint ridx) const;
    virtual double _getPatchSReacC(uint pidx, uint ridx) const;
    virtual double _getPatchSReacA(uint pidx, uint ridx) const;

    virtual uint _getPatchSReacExtent(uint pidx, uint ridx) const;
    virtual void _resetPatchSReacExtent(uint pidx, uint ridx);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER STATE ACCESS:
    //      TRIANGULAR SURFACE ELEMENTS
    ////////////////////////////////////////////////////////////////////////

    virtual double _getTriArea(uint tidx) const;
    virtual void _setTriArea(uint tidx, double area);

    virtual double _getTriCount(uint tidx, uint sidx) const;
    virtual void _setTriCount(uint tidx, uint sidx, double n);

    virtual bool _getTriClamped(uint tidx, uint sidx) const;
    virtual void _setTriClamped(uint tidx, uint sidx, bool buf);

    virtual double _getTriSReacK(uint tidx, uint ridx) const;
    virtual void _setTriSReacK(uint tidx, uint ridx, double kf);

    virtual bool _getTriSReacActive(uint tidx, uint ridx) const;
    virtual void _setTriSReacActive(uint tidx, uint ridx, bool act);

    ////////////////////////////////////////////////////////////////////////

    virtual double _getTriSReacH(uint tidx, uint ridx) const;
    virtual double _getTriSReacC(uint tidx, uint ridx) const;
    virtual double _getTriSReacA(uint tidx, uint ridx) const;

    ////////////////////////////////////////////////////////////////////////

    steps::model::Model * model(void) const
    { return pModel; }

    steps::wm::Geom * geom(void) const
    { return pGeom; }

    steps::rng::RNG * rng(void) const
    { return pRNG; }

    steps::solver::Statedef * statedef(void) const
    { return pStatedef; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    steps::model::Model *               pModel;

    steps::wm::Geom *                   pGeom;

    steps::rng::RNG *                   pRNG;

    Statedef *                          pStatedef;						/////////

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_API_HPP

// END

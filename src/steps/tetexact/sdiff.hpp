/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_TETEXACT_SDIFF_HPP
#define STEPS_TETEXACT_SDIFF_HPP 1


// Standard library & STL headers.
#include <array>
#include <map>
#include <string>
#include <vector>
#include <fstream>

// STEPS headers.
#include "util/common.h"
#include "kproc.hpp"
#include "tetexact.hpp"
#include "math/constants.hpp"
#include "solver/diffdef.hpp"
#include "solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace tetexact{

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tri;

////////////////////////////////////////////////////////////////////////////////

class SDiff
: public steps::tetexact::KProc
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    SDiff(steps::solver::Diffdef * sdef, steps::tetexact::Tri * tri);
    ~SDiff() override;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file) override;

    /// restore data
    void restore(std::fstream & cp_file) override;

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Diffdef * def() const noexcept
    { return pSDiffdef; }

    double dcst(int direction = -1);
    void setDcst(double d);
    void setDirectionDcst(int direction, double dcst);

    void setupDeps() override;

    bool depSpecTet(uint gidx, steps::tetexact::WmVol * tet) override;
    bool depSpecTri(uint gidx, steps::tetexact::Tri * tri) override;

    void reset() override;
    double rate(steps::tetexact::Tetexact * solver = nullptr) override;

    std::vector<KProc*> const & apply(const rng::RNGptr &rng, double dt, double simtime) override;

    uint updVecSize() const override;

    ////////////////////////////////////////////////////////////////////////

    void setSDiffBndActive(uint i, bool active);

    bool getSDiffBndActive(uint i) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                lidxTri;
    steps::solver::Diffdef              * pSDiffdef;
    steps::tetexact::Tri                * pTri;
    std::vector<KProc*>                 pUpdVec[3];

    // Storing the species local index for each neighbouring tri: Needed
    // because neighbours may belong to different patches for
    // diffusion boundary for surfaces
    std::array<solver::lidxT, 3>  pNeighbPatchLidx{{solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED, solver::LIDX_UNDEFINED}};

    /// Properly scaled diffusivity constant.
    double                              pScaledDcst{};
    // Compartmental dcst. Stored for convenience
    double                              pDcst{};
    std::map<uint, double>              directionalDcsts;

    std::array<double, 2>  				pCDFSelector{0.0, 0.0};

    // A flag to see if the species can move between compartments
    std::array<bool, 3>     pSDiffBndActive{false, false, false};

    // Flags to store if a direction is a diffusion boundary direction
    std::array<bool, 3>     pSDiffBndDirection{false, false, false};

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif // STEPS_TETEXACT_SDIFF_HPP

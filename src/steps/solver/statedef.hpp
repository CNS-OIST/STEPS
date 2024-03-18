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

#include <iosfwd>
#include <string>

#include "geom/fwd.hpp"
#include "model/fwd.hpp"
#include "rng/fwd.hpp"
#include "solver/fwd.hpp"

namespace steps::solver {

/// Defined State
class Statedef {
  public:
    enum PoolFlags { CLAMPED_POOLFLAG = 1 };
    static const uint PoolFlagDefault = 0;

    enum ReacFlags { INACTIVE_REACFLAG = 1, KCONST_REACFLAG = 2 };
    static const uint ReacFlagDefault = INACTIVE_REACFLAG | KCONST_REACFLAG;

    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param m Reference to the model object.
    /// \param g Reference to the geometry container.
    /// \param r Reference to the random number generator.
    Statedef(model::Model& m, wm::Geom& g, const rng::RNGptr& r);

    Statedef(const Statedef&) = delete;
    Statedef& operator=(const Statedef&) = delete;

    membrane_global_id getMembIdx(std::string const& m) const;

    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return reference to Compdef object specified by global index argument.
    ///
    /// \param gidx Global index of the Compdef object.
    Compdef& compdef(comp_global_id gidx) const;

    /// Return the total number of compartments in the simulation state.
    inline uint countComps() const noexcept {
        return pCompdefs.size();
    }

    /// Return the global index of compartment identified by string argument.
    ///
    /// \param c Name of the compartment.
    /// \exception Throw exception if geometry does not contain comp with this
    /// identifier.
    comp_global_id getCompIdx(std::string const& c) const;

    /// Return the global index of compartment identified by  object argument.
    ///
    /// \param comp Reference to the Comp object..
    /// \exception Throw exception if geometry does not contain comp with this
    /// identifier.
    comp_global_id getCompIdx(const wm::Comp& comp) const;

    inline const auto& comps() const noexcept {
        return pCompdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return reference to Patchdef object specified by global index argument.
    ///
    /// \param gidx Global index of the patch.
    Patchdef& patchdef(patch_global_id gidx) const;

    /// Return the total number of patches in the simulation state.
    inline uint countPatches() const noexcept {
        return pPatchdefs.size();
    }

    /// Return the global index of patch identified by string argument.
    ///
    /// \param p Name of the patch.
    /// \exception Throw exception if geometry does not contain patch with this
    /// identifier.
    patch_global_id getPatchIdx(std::string const& p) const;

    /// Return the global index of patch identified by string argument.
    ///
    /// \param patch Reference to the patch.
    /// \exception Throw exception if geometry does not contain patch with this
    /// identifier.
    patch_global_id getPatchIdx(const wm::Patch& patch) const;

    inline const auto& patches() const noexcept {
        return pPatchdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species in the simulation state.
    inline uint countSpecs() const noexcept {
        return pSpecdefs.size();
    }

    /// Return reference to Specdef object specified by global index argument.
    ///
    /// \param gidx Global index of the species.
    Specdef& specdef(spec_global_id gidx) const;

    /// Return the global index of species identified by string argument.
    ///
    /// \param s Name of the species.
    /// \exception Throw exception if model does not contain species with this
    /// identifier.
    spec_global_id getSpecIdx(std::string const& s) const;
    spec_global_id getSpecIdx(const model::Spec& spec) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEXES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of complexes in the simulation state.
    inline uint countComplexes() const noexcept {
        return pComplexdefs.size();
    }

    /// Return pointer to Complexdef object specified by global index argument.
    ///
    /// \param gidx Global index of the complex.
    Complexdef& complexdef(complex_global_id gidx) const;

    /// Return the global index of complexes identified by string argument.
    ///
    /// \param s Name of the complex.
    /// \exception Throw exception if model does not contain complexes with this identifier.
    complex_global_id getComplexIdx(std::string const& s) const;
    complex_global_id getComplexIdx(steps::model::Complex* spec) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of reactions in the simulation state.
    inline uint countReacs() const noexcept {
        return pReacdefs.size();
    }

    /// Return reference to Reacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the reaction.
    Reacdef& reacdef(reac_global_id gidx) const;

    /// Return the global index of reac identified by string argument.
    ///
    /// \param r Name of the reaction.
    /// \exception Throw exception if model does not contain reac with this
    /// identifier.
    reac_global_id getReacIdx(std::string const& r) const;

    /// Return the global index of reac identified by object argument.
    ///
    /// \param reac Reference to the reaction object.
    /// \exception Throw exception if model does not contain reac with this
    /// identifier.
    reac_global_id getReacIdx(const model::Reac& reac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of complex reactions in the simulation state.
    inline uint countComplexReacs() const noexcept {
        return pComplexReacdefs.size();
    }

    /// Return pointer to ComplexReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the complex reaction.
    ComplexReacdef& complexreacdef(complexreac_global_id gidx) const;

    /// Return the global index of complex reac identified by string argument.
    ///
    /// \param r Name of the complex reaction.
    /// \exception Throw exception if model does not contain complex reac with this identifier.
    complexreac_global_id getComplexReacIdx(std::string const& r) const;

    /// Return the global index of complex reac identified by object argument.
    ///
    /// \param reac Pointer to the complex reaction object.
    /// \exception Throw exception if model does not contain complex reac with this identifier.
    complexreac_global_id getComplexReacIdx(model::ComplexReac* reac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface reactions in the simulation state.
    inline uint countSReacs() const noexcept {
        return pSReacdefs.size();
    }

    /// Return reference to SReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    SReacdef& sreacdef(sreac_global_id gidx) const;

    /// Return the global index of surface reaction identified by string argument.
    ///
    /// \param sr Name of the surface reaction.
    /// \exception Throw exception if model does not contain sreac with this
    /// identifier.
    sreac_global_id getSReacIdx(std::string const& sr) const;
    /// Return the global index of surface reaction identified by object argument.
    ///
    /// \param sreac Reference to the surface reaction object..
    /// \exception Throw exception if model does not contain sreac with this
    /// identifier.
    sreac_global_id getSReacIdx(const model::SReac& sreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of complex reactions in the simulation state.
    inline uint countComplexSReacs() const noexcept {
        return pComplexSReacdefs.size();
    }

    /// Return pointer to ComplexSReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the complex reaction.
    ComplexSReacdef& complexsreacdef(complexsreac_global_id gidx) const;

    /// Return the global index of complex reac identified by string argument.
    ///
    /// \param r Name of the complex reaction.
    /// \exception Throw exception if model does not contain complex reac with this identifier.
    complexsreac_global_id getComplexSReacIdx(std::string const& r) const;

    /// Return the global index of complex reac identified by object argument.
    ///
    /// \param reac Pointer to the complex reaction object.
    /// \exception Throw exception if model does not contain complex reac with this identifier.
    complexsreac_global_id getComplexSReacIdx(model::ComplexSReac* reac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion rules in the simulation state.
    inline uint countDiffs() const noexcept {
        return pDiffdefs.size();
    }

    /// Return reference to Diffdef object specified by global index argument.
    ///
    /// \param gidx Global index of the diffusion.
    Diffdef& diffdef(diff_global_id gidx) const;

    /// Return the global index of diffusion identified by string argument.
    ///
    /// \param d Name of the diffusion.
    /// \exception Throw exception if model does not contain diff with this
    /// identifier.
    diff_global_id getDiffIdx(std::string const& d) const;

    /// Return the global index of diffusion identified by object argument.
    ///
    /// \param diff Reference to the diffusion object.
    /// \exception Throw exception if model does not contain this diff.
    diff_global_id getDiffIdx(const model::Diff& diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules in the simulation
    /// state.
    inline uint countSurfDiffs() const noexcept {
        return pSurfDiffdefs.size();
    }

    /// Return reference to SurfDiffdef object specified by global index argument.
    ///
    /// \param gidx Global index of the surface diffusion.
    SurfDiffdef& surfdiffdef(surfdiff_global_id gidx) const;

    /// Return the global index of surface diffusion identified by string
    /// argument.
    ///
    /// \param d Name of the surface diffusion.
    /// \exception Throw exception if model does not contain diff with this
    /// identifier.
    surfdiff_global_id getSurfDiffIdx(std::string const& d) const;

    /// Return the global index of surface diffusion identified by object
    /// argument.
    ///
    /// \param diff Reference to the surface diffusion object.
    /// \exception Throw exception if model does not contain this diff.
    surfdiff_global_id getSurfDiffIdx(const model::Diff& diff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of voltage-dependent transitions in the simulation
    /// state.
    inline uint countVDepSReacs() const noexcept {
        return pVDepSReacdefs.size();
    }

    /// Return reference to VDepSReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the voltage-dependent reaction.
    VDepSReacdef& vdepsreacdef(vdepsreac_global_id gidx) const;

    /// Return the global index of voltage-dependent reaction identified by string
    /// argument.
    ///
    /// \param vdt Name of the voltage-dependent reaction.
    /// \exception Throw exception if model does not contain vdepsreac with this
    /// identifier.
    vdepsreac_global_id getVDepSReacIdx(std::string const& vdt) const;

    /// Return the global index of voltage-dependent reaction identified by object
    /// argument.
    ///
    /// \param sreac Reference to the voltage-dependent reaction object..
    /// \exception Throw exception if model does not contain this vdepsreac.
    vdepsreac_global_id getVDepSReacIdx(const model::VDepSReac& vdepsreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of ohmic currents in the simulation state.
    inline uint countOhmicCurrs() const noexcept {
        return pOhmicCurrdefs.size();
    }

    /// Return reference to OhmicCurr object specified by global index argument.
    ///
    /// \param gidx Global index of the ohmic current.
    OhmicCurrdef& ohmiccurrdef(ohmiccurr_global_id gidx) const;

    /// Return the global index of ohmic current identified by string argument.
    ///
    /// \param oc Name of the ohmic current.
    /// \exception Throw exception if model does not contain ohmic current with
    /// this identifier.
    ohmiccurr_global_id getOhmicCurrIdx(std::string const& oc) const;

    /// Return the global index of ohmic current identified by object argument.
    ///
    /// \param ocurr Reference to the ohmic current object..
    /// \exception Throw exception if model does not contain this ohmiccurr.
    ohmiccurr_global_id getOhmicCurrIdx(const model::OhmicCurr& ohmiccurr) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of ghk currents in the simulation state.
    inline uint countGHKcurrs() const noexcept {
        return pGHKcurrdefs.size();
    }

    /// Return reference to GHKcurr object specified by global index argument.
    ///
    /// \param gidx Global index of the ghk current.
    GHKcurrdef& ghkcurrdef(ghkcurr_global_id gidx) const;

    /// Return the global index of ghk current identified by string argument.
    ///
    /// \param ghk Name of the ghk current.
    /// \exception Throw exception if model does not contain ghk current with this
    /// identifier.
    ghkcurr_global_id getGHKcurrIdx(std::string const& ghk) const;

    /// Return the global index of ghk current identified by object argument.
    ///
    /// \param ghkcurr Reference to the ghk current object..
    /// \exception Throw exception if model does not contain this ghkcurr.
    ghkcurr_global_id getGHKcurrIdx(const model::GHKcurr& ghkcurr) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion boundaries in the simulation state.
    inline uint countDiffBoundaries() const noexcept {
        return pDiffBoundarydefs.size();
    }

    /// Return reference to DiffBoundarydef object specified by global index
    /// argument.
    ///
    /// \param gidx Global index of the diffusion boundary.
    DiffBoundarydef& diffboundarydef(diffboundary_global_id gidx) const;

    /// Return the global index of diffusion boundary identified by string
    /// argument.
    ///
    /// \param d Name of the diffusion boundary.
    /// \exception Throw exception if geometry does not contain diff boundary with
    /// this identifier.
    diffboundary_global_id getDiffBoundaryIdx(std::string const& d) const;

    /// Return the global index of diffusion boundary identified by object
    /// argument.
    ///
    /// \param diff Reference to the diffusion boundary object.
    /// \exception Throw exception if geoemtry does not contain diff boundary with
    /// this identifier.
    diffboundary_global_id getDiffBoundaryIdx(const tetmesh::DiffBoundary& diffb) const;

    inline const auto& diffBoundaries() const noexcept {
        return pDiffBoundarydefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion boundaries in the simulation
    /// state.
    inline uint countSDiffBoundaries() const noexcept {
        return pSDiffBoundarydefs.size();
    }

    /// Return reference to SDiffBoundarydef object specified by global index
    /// argument.
    ///
    /// \param gidx Global index of the surface diffusion boundary.
    SDiffBoundarydef& sdiffboundarydef(sdiffboundary_global_id gidx) const;

    /// Return the global index of surface diffusion boundary identified by string
    /// argument.
    ///
    /// \param d Name of the surface diffusion boundary.
    /// \exception Throw exception if geometry does not contain sdiff boundary
    /// with this identifier.
    sdiffboundary_global_id getSDiffBoundaryIdx(std::string const& d) const;

    /// Return the global index of surface diffusion boundary identified by object
    /// argument.
    ///
    /// \param sdiffb Reference to the surface diffusion boundary object.
    /// \exception Throw exception if geometry does not contain sdiff boundary
    /// with this identifier.
    sdiffboundary_global_id getSDiffBoundaryIdx(const tetmesh::SDiffBoundary& sdiffb) const;

    inline const auto& sdiffBoundaries() const noexcept {
        return pSDiffBoundarydefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: BRIDGE SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species in the simulation state.
    inline uint countLinkSpecs() const noexcept {
        return pLinkSpecdefs.size();
    }

    /// Return reference to LinkSpecdef object specified by global index argument.
    ///
    /// \param gidx Global index of the Linkspecies.
    LinkSpecdef& linkspecdef(linkspec_global_id gidx) const;

    /// Return the global index of Linkspecies identified by string argument.
    ///
    /// \param s Name of the Linkspecies.
    /// \exception Throw exception if model does not contain Linkspecies with this
    /// identifier.
    linkspec_global_id getLinkSpecIdx(std::string const& s) const;
    linkspec_global_id getLinkSpecIdx(const model::LinkSpec& spec) const;

    inline const auto& linkspecs() const noexcept {
        return pLinkSpecdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicles in the simulation state.
    inline uint countVesicles() const noexcept {
        return pVesicledefs.size();
    }

    /// Return reference to Vesicledef object specified by global index argument.
    ///
    /// \param gidx Global index of the vesicle.
    Vesicledef& vesicledef(vesicle_global_id gidx) const;

    vesicle_global_id getVesicleIdx(std::string const& v) const;
    vesicle_global_id getVesicleIdx(const model::Vesicle& vesicle) const;

    inline const auto& vesicles() const noexcept {
        return pVesicledefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of rafts in the simulation state.
    inline uint countRafts() const noexcept {
        return pRaftdefs.size();
    }

    /// Return reference to Raftdef object specified by global index argument.
    ///
    /// \param gidx Global index of the raft.
    Raftdef& raftdef(raft_global_id gidx) const;

    raft_global_id getRaftIdx(std::string const& r) const;
    raft_global_id getRaftIdx(const model::Raft& raft) const;

    inline const auto& rafts() const noexcept {
        return pRaftdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE BINDING REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicle binding reactions in the simulation
    /// state.
    inline uint countVesBinds() const noexcept {
        return pVesBinddefs.size();
    }

    /// Return reference to VesBinddef object specified by global index argument.
    ///
    /// \param gidx Global index of the vesicle binding reaction.
    VesBinddef& vesbinddef(vesbind_global_id gidx) const;

    /// Return the global index of vesicle binding reaction identified by string
    /// argument.
    ///
    /// \param vb Name of the vesicle binding reaction.
    /// \exception Throw exception if model does not contain vesbind with this
    /// identifier.
    vesbind_global_id getVesBindIdx(std::string const& vb) const;

    /// Return the global index of vesicle binding reaction identified by object
    /// argument.
    ///
    /// \param vesbind Reference to the vesicle binding reaction object.
    /// \exception Throw exception if model does not contain vesbind with this
    /// identifier.
    vesbind_global_id getVesBindIdx(const model::VesBind& vesbind) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE UNBINDING REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicle unbinding reactions in the simulation
    /// state.
    inline uint countVesUnbinds() const noexcept {
        return pVesUnbinddefs.size();
    }

    /// Return reference to VesUnbinddef object specified by global index argument.
    ///
    /// \param gidx Global index of the vesicle unbinding reaction.
    VesUnbinddef& vesunbinddef(vesunbind_global_id gidx) const;

    /// Return the global index of vesicle unbinding reaction identified by string
    /// argument.
    ///
    /// \param vb Name of the vesicle unbinding reaction.
    /// \exception Throw exception if model does not contain vesunbind with this
    /// identifier.
    vesunbind_global_id getVesUnbindIdx(std::string const& vb) const;

    /// Return the global index of vesicle unbinding reaction identified by object
    /// argument.
    ///
    /// \param vesbind Reference to the vesicle unbinding reaction object.
    /// \exception Throw exception if model does not contain vesunbind with this
    /// identifier.
    vesunbind_global_id getVesUnbindIdx(const model::VesUnbind& vesunbind) const;

    inline const auto& vesunbinds() const noexcept {
        return pVesUnbinddefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicle surface reactions in the simulation
    /// state.
    inline uint countVesSReacs() const noexcept {
        return pVesSReacdefs.size();
    }

    /// Return reference to VesSReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the vesicle surface reaction.
    VesSReacdef& vessreacdef(vessreac_global_id gidx) const;

    /// Return the global index of vesicle surface reaction identified by string
    /// argument.
    ///
    /// \param vsr Name of the vesicle surface reaction.
    /// \exception Throw exception if model does not contain vessreac with this
    /// identifier.
    vessreac_global_id getVesSReacIdx(std::string const& vsr) const;

    /// Return the global index of vesicle surface reaction identified by object
    /// argument.
    ///
    /// \param vessreac Reference to the surface reaction object..
    /// \exception Throw exception if model does not contain vessreac with this
    /// identifier.
    vessreac_global_id getVesSReacIdx(const model::VesSReac& vessreac) const;

    inline const auto& vessreacs() const noexcept {
        return pVesSReacdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft surface reactions in the simulation state.
    inline uint countRaftSReacs() const noexcept {
        return pRaftSReacdefs.size();
    }

    /// Return reference to RaftSReacdef object specified by global index argument.
    ///
    /// \param gidx Global index of the raft surface reaction.
    RaftSReacdef& raftsreacdef(raftsreac_global_id gidx) const;

    /// Return the global index of raft surface reaction identified by string
    /// argument.
    ///
    /// \param rsr Name of the raft surface reaction.
    /// \exception Throw exception if model does not contain vessreac with this
    /// identifier.
    raftsreac_global_id getRaftSReacIdx(std::string const& rsr) const;

    /// Return the global index of raft surface reaction identified by object
    /// argument.
    ///
    /// \param raftsreac Reference to the surface reaction object..
    /// \exception Throw exception if model does not contain raftsreac with this
    /// identifier.
    raftsreac_global_id getRaftSReacIdx(const model::RaftSReac& raftsreac) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT ENDOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of endocytotic reactions in the simulation state.
    inline uint countRaftEndocytosis() const noexcept {
        return pRaftEndocytosisdefs.size();
    }

    /// Return reference to RaftEndocytosisdef object specified by global index
    /// argument.
    ///
    /// \param gidx Global index of the endocytotic reaction.
    RaftEndocytosisdef& raftendocytosisdef(raftendocytosis_global_id gidx) const;

    /// Return the global index of endocytotic reaction identified by string
    /// argument.
    ///
    /// \param endo Name of the endocytotic reaction.
    /// \exception Throw exception if model does not contain endocytosis with this
    /// identifier.
    raftendocytosis_global_id getRaftEndocytosisIdx(std::string const& endo) const;

    /// Return the global index of endocytotic reaction identified by object
    /// argument.
    ///
    /// \param endocyt Reference to the endocytotic reaction object..
    /// \exception Throw exception if model does not contain endocytosis with this
    /// identifier.
    raftendocytosis_global_id getRaftEndocytosisIdx(const model::RaftEndocytosis& endocyt) const;

    inline const auto& raftendocytosiss() const noexcept {
        return pRaftEndocytosisdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT GENESIS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft geneses in the simulation state.
    inline uint countRaftGens() const noexcept {
        return pRaftGendefs.size();
    }

    /// Return reference to pRaftGendefs object specified by global index argument.
    ///
    /// \param gidx Global index of the raft genesis.
    RaftGendef& raftgendef(raftgen_global_id gidx) const;

    /// Return the global index of raft genesis identified by string argument.
    ///
    /// \param endo Name of the raft genesis.
    /// \exception Throw exception if model does not contain raft genesis with
    /// this identifier.
    raftgen_global_id getRaftGenIdx(std::string const& raftgen) const;

    /// Return the global index of raft genesis identified by object argument.
    ///
    /// \param endocyt Reference to the raft genesis object..
    /// \exception Throw exception if model does not contain raft genesis with
    /// this identifier.
    raftgen_global_id getRaftGenIdx(const model::RaftGen& raftgen) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT DISSOLUTION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft dissolutions in the simulation state.
    inline uint countRaftDiss() const noexcept {
        return pRaftDisdefs.size();
    }

    /// Return reference to pRaftDisdefs object specified by global index argument.
    ///
    /// \param gidx Global index of the raft dissolution.
    RaftDisdef& raftdisdef(raftdis_global_id gidx) const;

    /// Return the global index of raft dissolution identified by string argument.
    ///
    /// \param endo Name of the raft dissolution.
    /// \exception Throw exception if model does not contain raft dissolution with
    /// this identifier.
    raftdis_global_id getRaftDisIdx(std::string const& raftdis) const;

    /// Return the global index of raft dissolution identified by object argument.
    ///
    /// \param endocyt Reference to the raft dissolution object..
    /// \exception Throw exception if model does not contain raft dissolution with
    /// this identifier.
    raftdis_global_id getRaftDisIdx(const model::RaftDis& raftdis) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: ENDOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of endocytotic reactions in the simulation state.
    inline uint countEndocytosis() const noexcept {
        return pEndocytosisdefs.size();
    }

    /// Return reference to Endocytosisdef object specified by global index
    /// argument.
    ///
    /// \param gidx Global index of the endocytotic reaction.
    Endocytosisdef& endocytosisdef(endocytosis_global_id gidx) const;

    /// Return the global index of endocytotic reaction identified by string
    /// argument.
    ///
    /// \param endo Name of the endocytotic reaction.
    /// \exception Throw exception if model does not contain endocytosis with this
    /// identifier.
    endocytosis_global_id getEndocytosisIdx(std::string const& endo) const;

    /// Return the global index of endocytotic reaction identified by object
    /// argument.
    ///
    /// \param endocyt Reference to the endocytotic reaction object..
    /// \exception Throw exception if model does not contain endocytosis with this
    /// identifier.
    endocytosis_global_id getEndocytosisIdx(const model::Endocytosis& endocyt) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EXOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of exocytotic reactions in the simulation state.
    inline uint countExocytosis() const noexcept {
        return pExocytosisdefs.size();
    }

    /// Return reference to Exocytosisdef object specified by global index argument.
    ///
    /// \param gidx Global index of the exocytotic reaction.
    Exocytosisdef& exocytosisdef(exocytosis_global_id gidx) const;

    /// Return the global index of exocytotic reaction identified by string
    /// argument.
    ///
    /// \param exo Name of the exocytotic reaction.
    /// \exception Throw exception if model does not contain exocytosis with this
    /// identifier.
    exocytosis_global_id getExocytosisIdx(std::string const& exo) const;

    /// Return the global index of exocytotic reaction identified by object
    /// argument.
    ///
    /// \param exocyt Reference to the exocytotic reaction object..
    /// \exception Throw exception if model does not contain exocytosis with this
    /// identifier.
    exocytosis_global_id getExocytosisIdx(const model::Exocytosis& exocyt) const;

    inline const auto& exocytosiss() const noexcept {
        return pExocytosisdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules in the simulation
    /// state.
    inline uint countVesSDiffs() const noexcept {
        return pVesSDiffdefs.size();
    }

    /// Return reference to VesSDiffdefs object specified by global index argument.
    ///
    /// \param gidx Global index of the vesicle surface diffusion.
    VesSDiffdef& vessdiffdef(vessdiff_global_id gidx) const;

    /// Return the global index of vesicle surface diffusion identified by string
    /// argument.
    ///
    /// \param vsd Name of the vesicle surface diffusion.
    /// \exception Throw exception if model does not contain diff with this
    /// identifier.
    vessdiff_global_id getVesSDiffIdx(std::string const& vsd) const;

    /// Return the global index of vesicle surface diffusion identified by object
    /// argument.
    ///
    /// \param vsdiff Reference to the vesicle surface diffusion object.
    /// \exception Throw exception if model does not contain this diff.
    vessdiff_global_id getVesSDiffIdx(const model::VesSDiff& vsdiff) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STATE
    ////////////////////////////////////////////////////////////////////////

    /// Return the current simulation time.
    inline double time() const noexcept {
        return pTime;
    }

    /// Return the model object.
    inline model::Model& model() const noexcept {
        return pModel;
    }

    /// Return the geometry object.
    inline wm::Geom& geom() const noexcept {
        return pGeom;
    }

    /// Return the random number generator object.
    inline const rng::RNGptr& rng() const noexcept {
        return pRNG;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: STATE
    ////////////////////////////////////////////////////////////////////////

    /// Set the current simulation time.
    ///
    /// \param t Simulation time.
    void setTime(double t);

    /// Increase the current simulation time.
    ///
    /// \param dt Discrete time.
    void incTime(double dt);

    /// Reset the simulation time to 0s.
    inline void resetTime() noexcept {
        pTime = 0.0;
    }

    /// Increase the time step.
    ///
    /// \param Time step to be increased.
    void incNSteps(uint i = 1);

    /// Reset the time step to 0.
    inline void resetNSteps() noexcept {
        pNSteps = 0;
    }

    /// Return current simulation time step.
    inline uint nsteps() const noexcept {
        return pNSteps;
    }

    inline void setNSteps(uint nsteps) noexcept {
        pNSteps = nsteps;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    model::Model& pModel;
    wm::Geom& pGeom;
    const rng::RNGptr pRNG;

    double pTime;

    uint pNSteps;

    util::strongid_vector<spec_global_id, std::unique_ptr<Specdef>> pSpecdefs;
    util::strongid_vector<complex_global_id, std::unique_ptr<Complexdef>> pComplexdefs;
    util::strongid_vector<chan_global_id, std::unique_ptr<Chandef>> pChandefs;
    util::strongid_vector<comp_global_id, std::unique_ptr<Compdef>> pCompdefs;
    util::strongid_vector<patch_global_id, std::unique_ptr<Patchdef>> pPatchdefs;
    util::strongid_vector<reac_global_id, std::unique_ptr<Reacdef>> pReacdefs;
    util::strongid_vector<complexreac_global_id, std::unique_ptr<ComplexReacdef>> pComplexReacdefs;
    util::strongid_vector<sreac_global_id, std::unique_ptr<SReacdef>> pSReacdefs;
    util::strongid_vector<complexsreac_global_id, std::unique_ptr<ComplexSReacdef>>
        pComplexSReacdefs;
    util::strongid_vector<diff_global_id, std::unique_ptr<Diffdef>> pDiffdefs;
    util::strongid_vector<surfdiff_global_id, std::unique_ptr<SurfDiffdef>> pSurfDiffdefs;
    util::strongid_vector<vdepsreac_global_id, std::unique_ptr<VDepSReacdef>> pVDepSReacdefs;
    util::strongid_vector<ohmiccurr_global_id, std::unique_ptr<OhmicCurrdef>> pOhmicCurrdefs;
    util::strongid_vector<ghkcurr_global_id, std::unique_ptr<GHKcurrdef>> pGHKcurrdefs;
    util::strongid_vector<diffboundary_global_id, std::unique_ptr<DiffBoundarydef>>
        pDiffBoundarydefs;
    util::strongid_vector<sdiffboundary_global_id, std::unique_ptr<SDiffBoundarydef>>
        pSDiffBoundarydefs;

    util::strongid_vector<linkspec_global_id, std::unique_ptr<LinkSpecdef>> pLinkSpecdefs;
    util::strongid_vector<vesicle_global_id, std::unique_ptr<Vesicledef>> pVesicledefs;
    util::strongid_vector<raft_global_id, std::unique_ptr<Raftdef>> pRaftdefs;
    util::strongid_vector<endocytosis_global_id, std::unique_ptr<Endocytosisdef>> pEndocytosisdefs;
    util::strongid_vector<exocytosis_global_id, std::unique_ptr<Exocytosisdef>> pExocytosisdefs;
    util::strongid_vector<vesbind_global_id, std::unique_ptr<VesBinddef>> pVesBinddefs;
    util::strongid_vector<vesunbind_global_id, std::unique_ptr<VesUnbinddef>> pVesUnbinddefs;
    util::strongid_vector<vessdiff_global_id, std::unique_ptr<VesSDiffdef>> pVesSDiffdefs;
    util::strongid_vector<vessreac_global_id, std::unique_ptr<VesSReacdef>> pVesSReacdefs;
    util::strongid_vector<raftsreac_global_id, std::unique_ptr<RaftSReacdef>> pRaftSReacdefs;
    util::strongid_vector<raftendocytosis_global_id, std::unique_ptr<RaftEndocytosisdef>>
        pRaftEndocytosisdefs;
    util::strongid_vector<raftgen_global_id, std::unique_ptr<RaftGendef>> pRaftGendefs;
    util::strongid_vector<raftdis_global_id, std::unique_ptr<RaftDisdef>> pRaftDisdefs;
};

}  // namespace steps::solver

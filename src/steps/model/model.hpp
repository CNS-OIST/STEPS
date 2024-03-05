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
#include "solver/fwd.hpp"

namespace steps::model {

/// Top-level container for the objects in a kinetic model.
///
/// A model::Model object is parent to the following objects:
/// <UL>
/// <LI>model::Spec
/// <LI>model::Volsys
/// <LI>model::Surfsys
/// <LI>model::Vesicle
/// <LI>model::VesSurfsys
/// <LI>model::Raft
/// <LI>model::Raftsys
/// </UL>
/// \sa Spec, Volsys, Surfsys, Vesicle, VesSurfsys, Raft, Raftsys
/// \warning Methods start with an underscore are not exposed to Python.
///

class Model {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////
    Model() = default;
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;


    /// Destructor
    ~Model();

    // Model * deepcopy();

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: SPECIES (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////
    /// Return a species with name id.
    ///
    /// \param id ID of the species.
    /// \return Reference to the species.
    Spec& getSpec(std::string const& id) const;

    /// Delete a species with name id.
    ///
    /// \param id ID of the species.
    void delSpec(std::string const& id) const;

    /// Return a list of all species in the Model object.
    ///
    /// \return List of pointers to the species in the Model object.
    std::vector<Spec*> getAllSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: LINK SPECIES (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////
    /// Return a link species with name id.
    ///
    /// \param id ID of the link species.
    /// \return Reference to the link species.
    LinkSpec& getLinkSpec(std::string const& id) const;

    /// Delete a link species with name id.
    ///
    /// \param id ID of the link species.
    void delLinkSpec(std::string const& id) const;

    /// Return a list of all link species in the Model object.
    ///
    /// \return List of pointers to the link species in the Model object.
    std::vector<LinkSpec*> getAllLinkSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: CHANNELS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////
    /// Return a channel with name id.
    ///
    /// \param id ID of the channel.
    /// \return Reference to the channel.
    Chan& getChan(std::string const& id) const;

    /// Return a list of all channels in the Model object.
    ///
    /// \return List of pointers to the channels in the Model object.
    std::vector<Chan*> getAllChans() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: VESICLES (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a vesicle with name id.
    ///
    /// \param id ID of the vesicle.
    /// \return Reference to the vesicle.
    Vesicle& getVesicle(std::string const& id) const;

    /// Return a list of all vesicles in the Model object.
    ///
    /// \return List of pointers to the vesicles in the Model object.
    std::vector<Vesicle*> getAllVesicles() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: RAFTS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a raft with name id.
    ///
    /// \param id ID of the raft.
    /// \return Reference to the raft.
    Raft& getRaft(std::string const& id) const;

    /// Return a list of all rafts in the Model object.
    ///
    /// \return List of pointers to the rafts in the Model object.
    std::vector<Raft*> getAllRafts() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: VOLSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a volume system with name id.
    ///
    /// \param id ID of the volume system.
    /// \return Reference to the volume system.
    Volsys& getVolsys(std::string const& id) const;

    /// Delete a volume system with name id.
    ///
    /// \param id ID of the volume system.
    void delVolsys(std::string const& id) const;

    /// Return a list of all volume systems in the Model object.
    ///
    /// \return List of pointers to the volume systems in the Model object.
    std::vector<Volsys*> getAllVolsyss() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: SURFSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a surface system with name id.
    ///
    /// \param id ID of the surface system.
    /// \return Reference to the surface system.
    Surfsys& getSurfsys(std::string const& id) const;

    /// Delete a surface system with name id.
    ///
    /// \param id ID of the surface system.
    void delSurfsys(std::string const& id) const;

    /// Return a list of all surface systems in the Model object.
    ///
    /// \return List of pointers to the surface systems in the Model object.
    std::vector<Surfsys*> getAllSurfsyss() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: VESICLE SURFSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a vesicle surface system with name id.
    ///
    /// \param id ID of the vesicle surface system.
    /// \return Reference to the vesicle surface system.
    VesSurfsys& getVesSurfsys(std::string const& id) const;

    /// Delete a vesicle surface system with name id.
    ///
    /// \param id ID of the vesicle surface system.
    void delVesSurfsys(std::string const& id) const;

    /// Return a list of all vesiclev surface systems in the Model object.
    ///
    /// \return List of pointers to the vesicle surface systems in the Model
    /// object.
    std::vector<VesSurfsys*> getAllVesSurfsyss() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS: RAFTSYS (EXPOSED TO PYTHON)
    ////////////////////////////////////////////////////////////////////////

    /// Return a raft system with name id.
    ///
    /// \param id ID of the raft system.
    /// \return Reference to the raft system.
    Raftsys& getRaftsys(std::string const& id) const;

    /// Delete a raft system with name id.
    ///
    /// \param id ID of the raft system.
    void delRaftsys(std::string const& id) const;

    /// Return a list of all raft systems in the Model object.
    ///
    /// \return List of pointers to the raft systems in the Model object.
    std::vector<Raftsys*> getAllRaftsyss() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the species in the Model object.
    ///
    /// \return Number of species.
    inline uint _countSpecs() const noexcept {
        return pSpecs.size();
    }

    /// Return a species with index gidx.
    ///
    /// \param gidx Index of the species.
    /// \return Reference to the species.
    Spec& _getSpec(solver::spec_global_id gidx) const;

    /// Count the link species in the Model object.
    ///
    /// \return Number of link species.
    inline uint _countLinkSpecs() const noexcept {
        return pLinkSpecs.size();
    }

    /// Return a link species with index gidx.
    ///
    /// \param gidx Index of the link species.
    /// \return Reference to the link species.
    LinkSpec& _getLinkSpec(solver::linkspec_global_id gidx) const;

    /// Count the vesicles in the Model object.
    ///
    /// \return Number of vesicles.
    inline uint _countVesicles() const noexcept {
        return pVesicles.size();
    }

    /// Return a vesicle with index gidx.
    ///
    /// \param gidx Index of the vesicle.
    /// \return Reference to the vesicle.
    Vesicle& _getVesicle(solver::vesicle_global_id gidx) const;

    /// Count the rafts in the Model object.
    ///
    /// \return Number of rafts.
    inline uint _countRafts() const noexcept {
        return pRafts.size();
    }

    /// Return a raft with index gidx.
    ///
    /// \param gidx Index of the raft.
    /// \return Reference to the raft.
    Raft& _getRaft(solver::raft_global_id gidx) const;

    /// Count the complexes in the Model object.
    ///
    /// \return Number of complexes.
    inline uint _countComplexes() const noexcept {
        return pComplexes.size();
    }

    /// Return a complex with index gidx.
    ///
    /// \param gidx Index of the species.
    /// \return Pointer to the complex.
    Complex& _getComplex(solver::complex_global_id gidx) const;

    /// Count the channels in the Model object.
    ///
    /// \return Number of channels.
    inline uint _countChans() const noexcept {
        return pChans.size();
    }

    /// Return a channel with index gidx.
    ///
    /// \param gidx Index of the channel.
    /// \return Reference to the channel.
    Chan& _getChan(solver::chan_global_id gidx) const;

    /// Count the reactions in the Model object.
    ///
    /// \return Number of reactions.
    uint _countReacs() const;

    /// Count the complex reactions in the Model object.
    ///
    /// \return Number of complex reactions.
    uint _countComplexReacs() const;

    /// Return a reaction with index gidx
    ///
    /// \param gidx Index of the reaction.
    /// \return Reference to the reaction.
    Reac& _getReac(solver::reac_global_id gidx) const;

    /// Return a complex reaction with index gidx
    ///
    /// \param gidx Index of the complex reaction.
    /// \return Pointer to the conplex reaction.
    ComplexReac& _getComplexReac(solver::complexreac_global_id gidx) const;

    /// Count the surface reactions in the Model object.
    ///
    /// \return Number of surface reactions.
    uint _countSReacs() const;

    /// Count the complex surface reactions in the Model object.
    ///
    /// \return Number of complex reactions.
    uint _countComplexSReacs() const;

    /// Return a surface with index gidx.
    ///
    /// \param gidx Index of the surface reaction.
    /// \return Reference to the surface reaction.
    SReac& _getSReac(solver::sreac_global_id gidx) const;

    /// Return a complex surface reaction with index gidx
    ///
    /// \param gidx Index of the complex reaction.
    /// \return Pointer to the conplex reaction.
    ComplexSReac& _getComplexSReac(solver::complexsreac_global_id gidx) const;

    /// Count the raft geneses in the Model object.
    ///
    /// \return Number of raft geneses
    uint _countRaftGeneses() const;

    /// Return a raft genesis with index gidx.
    ///
    /// \param gidx Index of the raft genesis.
    /// \return Reference to the raft genesis.
    RaftGen& _getRaftGen(solver::raftgen_global_id gidx) const;

    /// Count the raft dissolutions in the Model object.
    ///
    /// \return Number of raft dissolutions
    uint _countRaftDiss() const;

    /// Return a raft dissolution with index gidx.
    ///
    /// \param gidx Index of the raft dissolution.
    /// \return Reference to the raft dissolution.
    RaftDis& _getRaftDis(solver::raftdis_global_id gidx) const;

    /// Count the endocytotic reactions in the Model object.
    ///
    /// \return Number of endocytotic reactions.
    uint _countEndocytosis() const;

    /// Return a endocytotic reaction with index gidx.
    ///
    /// \param gidx Index of the endocytotic reaction.
    /// \return Reference to the endocytotic reaction.
    Endocytosis& _getEndocytosis(solver::endocytosis_global_id gidx) const;

    /// Count the raft endocytotic reactions in the Model object.
    ///
    /// \return Number of raft endocytotic reactions.
    uint _countRaftEndocytosis() const;

    /// Return a raft endocytotic reaction with index gidx.
    ///
    /// \param gidx Index of the raft endocytotic reaction.
    /// \return Reference to the raft endocytotic reaction.
    RaftEndocytosis& _getRaftEndocytosis(solver::raftendocytosis_global_id gidx) const;

    /// Count the exocytotic reactions in the Model object.
    ///
    /// \return Number of exocytotic reactions.
    uint _countExocytosis() const;

    /// Return a exocytotic reaction with index gidx.
    ///
    /// \param gidx Index of the exocytotic reaction.
    /// \return Reference to the exocytotic reaction.
    Exocytosis& _getExocytosis(solver::exocytosis_global_id gidx) const;

    /// Return a voltage-dependent reaction with index gidx.
    ///
    /// \param gidx Index of the voltage-dependent reaction.
    /// \return Reference to the voltage-dependent reaction.
    VDepSReac& _getVDepSReac(solver::vdepsreac_global_id gidx) const;

    /// Count the voltage-dependent reactions in the Model object.
    ///
    /// \return Number of voltage-dependent reactions.
    uint _countVDepSReacs() const;

    /// Count the ohmic currents in the Model object.
    ///
    /// \return Number of ohmic currents.
    uint _countOhmicCurrs() const;

    /// Return an ohmic current with index gidx.
    ///
    /// \param gidx Index of the ohmic current.
    /// \return Reference to the ohmic current.
    OhmicCurr& _getOhmicCurr(solver::ohmiccurr_global_id gidx) const;

    /// Count the ghk currents in the Model object.
    ///
    /// \return Number of ghk currents.
    uint _countGHKcurrs() const;

    /// Return an ghk current with index gidx.
    ///
    /// \param gidx Index of the ghk current.
    /// \return Reference to the ghk current.
    GHKcurr& _getGHKcurr(solver::ghkcurr_global_id gidx) const;

    /// Count the volume diffusions in the Model object.
    ///
    /// \return Number of volume diffusions.
    uint _countVDiffs() const;

    /// Return a volume diffusion with index gidx.
    ///
    /// \param gidx Index of the volume diffusion.
    /// \return Reference to the volume diffusion.
    Diff& _getVDiff(solver::diff_global_id gidx) const;

    /// Count the surface diffusions in the Model object.
    ///
    /// \return Number of surface diffusions.
    uint _countSDiffs() const;

    /// Return a surface diffusion with index gidx.
    ///
    /// \param gidx Index of the surface diffusion.
    /// \return Reference to the surface diffusion.
    Diff& _getSDiff(solver::surfdiff_global_id gidx) const;

    /// Count the vesicle binding reactions in the Model object.
    ///
    /// \return Number of vesicle binding reactions.
    uint _countVesBinds() const;

    /// Return a vesicle binding reaction with index gidx.
    ///
    /// \param gidx Index of the vesicle binding reaction.
    /// \return Reference to the vesicle binding reaction.
    VesBind& _getVesBind(solver::vesbind_global_id gidx) const;

    /// Count the vesicle unbinding reactions in the Model object.
    ///
    /// \return Number of vesicle unbinding reactions.
    uint _countVesUnbinds() const;

    /// Return a vesicle unbinding reaction with index gidx.
    ///
    /// \param gidx Index of the vesicle unbinding reaction.
    /// \return Reference to the vesicle unbinding reaction.
    VesUnbind& _getVesUnbind(solver::vesunbind_global_id gidx) const;

    /// Count the vesicle surface diffusions in the Model object.
    ///
    /// \return Number of vesicle surface diffusions.
    uint _countVesSDiffs() const;

    /// Return a vesicle surface diffusion with index gidx.
    ///
    /// \param gidx Index of the vesicle surface diffusion.
    /// \return Reference to the surface diffusion.
    VesSDiff& _getVesSDiff(solver::vessdiff_global_id gidx) const;

    /// Count the vesicle surface reactions in the Model object.
    ///
    /// \return Number of vesicle surface reactions.
    uint _countVesSReacs() const;

    /// Return a vesicle surface reaction with index gidx.
    ///
    /// \param gidx Index of the vesicle surface reaction.
    /// \return Reference to the vesicle surface reaction.
    VesSReac& _getVesSReac(solver::vessreac_global_id gidx) const;

    /// Count the raft surface reactions in the Model object.
    ///
    /// \return Number of raft surface reactions.
    uint _countRaftSReacs() const;

    /// Return a raft surface reaction with index gidx.
    ///
    /// \param gidx Index of the raft surface reaction.
    /// \return Reference to the raft surface reaction.
    RaftSReac& _getRaftSReac(solver::raftsreac_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a species id is occupied.
    ///
    /// \param id ID of the species.
    void _checkSpecID(std::string const& id) const;

    /// Change the id of a species from o to n.
    ///
    /// \param o Old id of the species.
    /// \param n New id of the species.
    void _handleSpecIDChange(std::string const& o, std::string const& n);

    /// Add a species to the Model.
    ///
    /// \param spec Reference to the species being added.
    void _handleSpecAdd(Spec& spec);

    /// Delete a species in the Model.
    ///
    /// \param spec Reference to the species being deleted.
    void _handleSpecDel(Spec& spec);

    /// Check if a link species id is occupied.
    ///
    /// \param id ID of the link species.
    void _checkLinkSpecID(std::string const& id) const;

    /// Change the id of a link species from o to n.
    ///
    /// \param o Old id of the link species.
    /// \param n New id of the link species.
    void _handleLinkSpecIDChange(std::string const& o, std::string const& n);

    /// Add a link species to the Model.
    ///
    /// \param spec Reference to the link species being added.
    void _handleLinkSpecAdd(LinkSpec& spec);

    /// Delete a link species in the Model.
    ///
    /// \param spec Reference to the link species being deleted.
    void _handleLinkSpecDel(LinkSpec& spec);

    /// Check if a vesicle id is occupied.
    ///
    /// \param id ID of the vesicle.
    void _checkVesicleID(std::string const& id) const;

    /// Change the id of a vesicle from o to n.
    ///
    /// \param o Old id of the vesicle.
    /// \param n New id of the vesicle.
    void _handleVesicleIDChange(std::string const& o, std::string const& n);

    /// Add a vesicle to the Model.
    ///
    /// \param spec Reference to the vesicle being added.
    void _handleVesicleAdd(Vesicle& vesicle);

    // Not fully supported yet
    /// Delete a vesicle in the Model.
    ///
    /// \param vesicle Reference to the vesicle being deleted.
    void _handleVesicleDel(Vesicle& vesicle);

    /// Check if a raft id is occupied.
    ///
    /// \param id ID of the raft.
    void _checkRaftID(std::string const& id) const;

    /// Change the id of a raft from o to n.
    ///
    /// \param o Old id of the raft.
    /// \param n New id of the raft.
    void _handleRaftIDChange(std::string const& o, std::string const& n);

    /// Add a raft to the Model.
    ///
    /// \param spec Reference to the raft being added.
    void _handleRaftAdd(Raft& raft);

    // Not fully supported yet
    /// Delete a raft in the Model.
    ///
    /// \param vesicle Reference to the raft being deleted.
    void _handleRaftDel(Raft& raft);

    /// Add a complex to the Model.
    ///
    /// \param cmplx Pointer to the complex being added.
    void _handleComplexAdd(Complex& cmplx);

    /// Delete a species in the Model.
    ///
    /// \param spec Pointer to the species being deleted.
    void _handleComplexDel(Complex& cmplx);


    /// Check if a channel id is occupied.
    ///
    /// \param id ID of the channel.
    void _checkChanID(std::string const& id) const;

    /// Change the id of a channel from o to n.
    ///
    /// \param o Old id of the channel.
    /// \param n New id of the channel.
    void _handleChanIDChange(std::string const& o, std::string const& n);

    /// Add a channel to the Model.
    ///
    /// \param spec Reference to the channel being added.
    void _handleChanAdd(Chan& chan);

    /// Delete a channel in the Model.
    ///
    /// \param chan Reference to the channel being deleted.
    void _handleChanDel(Chan& chan);

    /// Check if a volume system id is occupied.
    ///
    /// \param id ID of the volume system.
    void _checkVolsysID(std::string const& id) const;

    /// Change the id of a volume system from o to n.
    ///
    /// \param o Old id of the volume system.
    /// \param n New id of the volume system.
    void _handleVolsysIDChange(std::string const& o, std::string const& n);

    /// Add a volume system to the Model.
    ///
    /// \param volsys Reference to the volume system being added.
    void _handleVolsysAdd(Volsys& volsys);

    /// Delete a volume system in the Model.
    ///
    /// \param volsys Reference to the volume system being deleted.
    void _handleVolsysDel(Volsys& volsys);

    /// Check if a surface system id is occupied.
    ///
    /// \param id ID of the surface system.
    void _checkSurfsysID(std::string const& id) const;

    /// Change the id of a surface system from o to n.
    ///
    /// \param o Old id of the surface system.
    /// \param n New id of the surface system.
    void _handleSurfsysIDChange(std::string const& o, std::string const& n);

    /// Add a surface system to the Model.
    ///
    /// \param surfsys Reference to the surface system being added.
    void _handleSurfsysAdd(Surfsys& surfsys);

    /// Delete a surface system in the Model.
    ///
    /// \param surfsys Reference to the surface system being deleted.
    void _handleSurfsysDel(Surfsys& surfsys);

    /// Check if a vesicle surface system id is occupied.
    ///
    /// \param id ID of the vesicle surface system.
    void _checkVesSurfsysID(std::string const& id) const;

    /// Add a vesicle surface system to the Model.
    ///
    /// \param vessurfsys Reference to the vesicle surface system being added.
    void _handleVesSurfsysAdd(VesSurfsys& vessurfsys);

    /// Change the id of a vesicle surface system from o to n.
    ///
    /// \param o Old id of the vesicle surface system.
    /// \param n New id of the vesicle surface system.
    void _handleVesSurfsysIDChange(std::string const& o, std::string const& n);

    /// Delete a vesicle surface system in the Model.
    ///
    /// \param vessurfsys Reference to the vesicle surface system being deleted.
    void _handleVesSurfsysDel(VesSurfsys& vessurfsys);

    /// Check if a raft system id is occupied.
    ///
    /// \param id ID of the raft system.
    void _checkRaftsysID(std::string const& id) const;

    /// Add a raft system to the Model.
    ///
    /// \param raftsys Reference to the raft system being added.
    void _handleRaftsysAdd(Raftsys& raftsys);

    /// Change the id of a raft system from o to n.
    ///
    /// \param o Old id of the raft system.
    /// \param n New id of the raft system.
    void _handleRaftsysIDChange(std::string const& o, std::string const& n);

    /// Delete a raft system in the Model.
    ///
    /// \param raftsys Reference to the raft system being deleted.
    void _handleRaftsysDel(Raftsys& raftsys);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::map<std::string, Spec*> pSpecs;
    std::map<std::string, LinkSpec*> pLinkSpecs;
    std::map<std::string, Vesicle*> pVesicles;
    std::map<std::string, Raft*> pRafts;
    std::map<std::string, Chan*> pChans;
    std::map<std::string, Volsys*> pVolsys;
    std::map<std::string, Surfsys*> pSurfsys;
    std::map<std::string, VesSurfsys*> pVesSurfsys;
    std::map<std::string, Raftsys*> pRaftsys;
    std::map<std::string, Complex*> pComplexes;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model

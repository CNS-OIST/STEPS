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


#include "model/complexevents.hpp"

#include "util/error.hpp"

namespace steps::model {

const std::array<ComplexLocation, 3> AllPatchLocations = {ComplexLocation::PATCH_IN,
                                                          ComplexLocation::PATCH_SURF,
                                                          ComplexLocation::PATCH_OUT};

ComplexEvent::ComplexEvent(std::string const& cmplxId)
    : pcomplexId(cmplxId) {}

ComplexCreateEvent::ComplexCreateEvent(std::string const& cmplxId, const std::vector<uint>& in)
    : ComplexEvent::ComplexEvent(cmplxId)
    , pinit(in) {}

ComplexLHSEvent::ComplexLHSEvent(std::string const& cmplxId,
                                 const std::vector<std::vector<SubunitStateFilter>>& filts)
    : ComplexEvent::ComplexEvent(cmplxId) {
    pfilters.reserve(filts.size());
    for (const auto& filt: filts) {
        pfilters.emplace_back(filt);
    }
}

ComplexUpdateEvent::ComplexUpdateEvent(std::string const& cmplxId,
                                       const std::vector<std::vector<SubunitStateFilter>>& filts,
                                       const std::vector<uint>& reac,
                                       const std::vector<ComplexUpdate>& upd,
                                       ComplexLocation destLoc)
    : ComplexLHSEvent::ComplexLHSEvent(cmplxId, filts)
    , preactants(reac)
    , pupdates(upd)
    , pdestLoc(destLoc) {}

ComplexDeleteEvent::ComplexDeleteEvent(std::string const& cmplxId,
                                       const std::vector<std::vector<SubunitStateFilter>>& filts)
    : ComplexLHSEvent::ComplexLHSEvent(cmplxId, filts) {}


}  // namespace steps::model

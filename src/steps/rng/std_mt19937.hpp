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

#include <random>

#include "rng/rng.hpp"

namespace steps::rng {

class STDMT19937: public RNG {
  public:
    STDMT19937(unsigned int bufsize);
    ~STDMT19937() override;

    void checkpoint(std::ostream& cp_file) const override;

    void restore(std::istream& cp_file) override;

  protected:
    void concreteInitialize(unsigned long seed) override;

    void concreteFillBuffer() override;

    std::mt19937 rng_;
};

}  // namespace steps::rng

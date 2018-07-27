/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

     This file is part of MSL GUI.

     MSL GUI is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     MSL GUI is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with MSL GUI.  If not, see <https://www.gnu.org/licenses/>.
*/


#include "dm_aggframe.h"


namespace ldomain
{
namespace poisson_bvp
{


using FrameTL::AggregatedFrame;
using FrameTL::SimpleEllipticEquation;


template <class IBASIS>
AggregatedFrame<IBASIS,2>* DM_AggFrame::Problem<IBASIS>::createWaveletSystem(int jmax)
{
    return createAggFrame_on_LDomain<IBASIS>(jmax, this->lastInput_.overlap, 1);
}

template <class IBASIS>
SimpleEllipticEquation<IBASIS,2>*
DM_AggFrame::Problem<IBASIS>::createDiscretizedEquation(const MathTL::PoissonBVP<2>& rawProblem,
                                                        const AggregatedFrame<IBASIS,2>& frame)
{
    return new SimpleEllipticEquation<IBASIS,2>(&rawProblem, &frame, this->lastInput_.jmax);
}

}
}


// explicit template instantiation:
template class Instantiate<ldomain::poisson_bvp::DM_AggFrame>;

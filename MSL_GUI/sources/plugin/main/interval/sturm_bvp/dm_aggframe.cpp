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

namespace interval
{
namespace sturm_bvp
{


using FrameTL::AggregatedFrame;
using FrameTL::EllipticEquation;


template <class IBASIS>
AggregatedFrame<IBASIS,1>* DM_AggFrame::Problem<IBASIS>::createWaveletSystem(int jmax)
{
    return createAggFrame_on_Interval<IBASIS>(jmax, this->lastInput_.overlap, 1);
}

template <class IBASIS>
EllipticEquation<IBASIS,1>*
DM_AggFrame::Problem<IBASIS>::createDiscretizedEquation(const MathTL::EllipticBVP<1>& rawProblem,
                                                        const AggregatedFrame<IBASIS,1>& frame)
{
    return new EllipticEquation<IBASIS,1>(&rawProblem, &frame, this->lastInput_.jmax);
}

}
}


// explicit template instantiation:
template class Instantiate<interval::sturm_bvp::DM_AggFrame>;

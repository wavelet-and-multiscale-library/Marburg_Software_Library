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


#ifndef INTERVAL_STURM_BVP_DM_AGGFRAME_H
#define INTERVAL_STURM_BVP_DM_AGGFRAME_H


#include "discr/generic_aggframe_discretization_module.h"

#include "FrameTL/elliptic_equation.h"



namespace interval
{
namespace sturm_bvp
{

class DM_AggFrame : public GenericAggFrameDiscretizationModule<10>
{
public:

    template <class IBASIS>
    class Problem
        : public AggFrameDiscretizedProblem<FrameTL::EllipticEquation<IBASIS,1>, MathTL::EllipticBVP<1> >
    {
        FrameTL::AggregatedFrame<IBASIS,1>* createWaveletSystem(int jmax) override;

        FrameTL::EllipticEquation<IBASIS,1>*
        createDiscretizedEquation(const MathTL::EllipticBVP<1>& rawProblem,
                                  const FrameTL::AggregatedFrame<IBASIS,1>& frame) override;
    };
};


}
}

#endif // INTERVAL_STURM_BVP_DM_AGGFRAME_H

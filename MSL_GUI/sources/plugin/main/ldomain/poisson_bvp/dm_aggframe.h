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


#ifndef LDOMAIN_POISSON_BVP_DM_AGGFRAME_H
#define LDOMAIN_POISSON_BVP_DM_AGGFRAME_H


#include "discr/generic_aggframe_discretization_module.h"

#include "FrameTL/simple_elliptic_equation.h"



namespace ldomain
{
namespace poisson_bvp
{

class DM_AggFrame : public GenericAggFrameDiscretizationModule<6>
{
public:

    template <class IBASIS>
    class Problem
        : public AggFrameDiscretizedProblem<FrameTL::SimpleEllipticEquation<IBASIS,2>, MathTL::PoissonBVP<2> >
    {
        FrameTL::AggregatedFrame<IBASIS,2>* createWaveletSystem(int jmax) override;

        FrameTL::SimpleEllipticEquation<IBASIS,2>*
        createDiscretizedEquation(const MathTL::PoissonBVP<2>& rawProblem,
                                  const FrameTL::AggregatedFrame<IBASIS,2>& frame) override;
    };
};


}
}

#endif // LDOMAIN_POISSON_BVP_DM_AGGFRAME_H

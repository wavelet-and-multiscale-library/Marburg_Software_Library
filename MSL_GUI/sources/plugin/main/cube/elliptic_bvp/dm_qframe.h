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


#ifndef CUBE_ELLIPTIC_BVP_DM_QFRAME_H
#define CUBE_ELLIPTIC_BVP_DM_QFRAME_H


#define DYADIC
#define FRAME


#include "discr/generic_qframe_discretization_module.h"

#include "WaveletTL/galerkin/tframe_equation.h"


#undef DYADIC
#undef FRAME


namespace cube
{
namespace elliptic_bvp
{

class DM_QFrame : public GenericQFrameDiscretizationModule<6>
{
public:

    template <class IFRAME>
    class Problem
        : public QFrameDiscretizedCubeProblem<WaveletTL::TensorFrameEquation<IFRAME,2>, MathTL::EllipticBVP<2> >
    {
        WaveletTL::TensorFrame<IFRAME,2>* createWaveletSystem(int jmax) override;

        WaveletTL::TensorFrameEquation<IFRAME, 2>*
        createDiscretizedEquation(const MathTL::EllipticBVP<2>& rawProblem,
                                  const WaveletTL::TensorFrame<IFRAME,2>& frame) override;
    };
};


}
}

#endif // CUBE_ELLIPTIC_BVP_DM_QFRAME_H

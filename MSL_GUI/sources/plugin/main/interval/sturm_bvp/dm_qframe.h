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


#ifndef INTERVAL_STURM_BVP_DM_QFRAME_H
#define INTERVAL_STURM_BVP_DM_QFRAME_H


#define FRAME
#define DYADIC


#include "discr/generic_qframe_discretization_module.h"

#include "sturm_qframe_equation.h"


#undef DYADIC
#undef FRAME


namespace interval
{
namespace sturm_bvp
{

class DM_QFrame : public GenericQFrameDiscretizationModule<10, IntervalQFramesList_1D>
{
public:

    template<class IBASIS>
    class Problem
       : public QFrameDiscretizedIntervalProblem<WaveletTL::SturmQFrameEquation<IBASIS>, MathTL::SimpleSturmBVP>
    {
        IBASIS* createWaveletSystem(int jmax) override;

        WaveletTL::SturmQFrameEquation<IBASIS>*
        createDiscretizedEquation(const MathTL::SimpleSturmBVP& rawProblem,
                                  const IBASIS& basis) override;
    };

};

}
}

#endif // INTERVAL_STURM_BVP_DM_QFRAME_H

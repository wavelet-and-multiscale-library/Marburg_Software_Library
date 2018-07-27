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


#include "dm_basis.h"

namespace interval
{
namespace sturm_bvp
{


template<class IBASIS>
IBASIS* DM_Basis::Problem<IBASIS>::createWaveletSystem(int jmax)
{
    bool bc_left = true;
    bool bc_right = true;

    IBASIS* basis = new IBASIS(bc_left, bc_right);
    basis->set_jmax(jmax);

    return basis;
}



template<class IBASIS>
WaveletTL::SturmEquation<IBASIS>*
DM_Basis::Problem<IBASIS>::createDiscretizedEquation(const MathTL::SimpleSturmBVP& rawProblem,
                                                     const IBASIS& basis)
{
    return new WaveletTL::SturmEquation<IBASIS>(rawProblem, basis, true);
}


}
}


// explicit template instantiation:
template class Instantiate<interval::sturm_bvp::DM_Basis>;

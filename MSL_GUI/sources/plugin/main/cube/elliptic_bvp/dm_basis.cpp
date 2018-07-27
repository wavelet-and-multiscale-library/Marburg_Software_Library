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

namespace cube
{
namespace elliptic_bvp
{

template <class IBASIS>
WaveletTL::CubeBasis<IBASIS,2>* DM_Basis::Problem<IBASIS>::createWaveletSystem(int jmax)
{
    MathTL::FixedArray1D<bool,4> bc;
    bc[0] = bc[1] = bc[2] = bc[3] = true;   // homogeneous Dirichlet boundary conditions

    WaveletTL::CubeBasis<IBASIS,2>* basis = new WaveletTL::CubeBasis<IBASIS,2>(bc);
    basis->set_jmax(jmax);

    return basis;
}

template <class IBASIS>
WaveletTL::CubeEquation<IBASIS,2>*
DM_Basis::Problem<IBASIS>::createDiscretizedEquation(const MathTL::EllipticBVP<2>& rawProblem,
                                                     const WaveletTL::CubeBasis<IBASIS,2>& basis)
{
    return new WaveletTL::CubeEquation<IBASIS,2>(&rawProblem, basis);
}

}
}


// explicit template instantiation:
template class Instantiate<cube::elliptic_bvp::DM_Basis>;

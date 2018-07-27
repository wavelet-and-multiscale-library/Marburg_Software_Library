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


#ifndef INTERVAL_BIHARMONIC_BVP_DM_BASIS_H
#define INTERVAL_BIHARMONIC_BVP_DM_BASIS_H


#include "discr/generic_basis_discretization_module.h"

#include "WaveletTL/galerkin/biharmonic_equation.h"



namespace interval
{
namespace biharmonic_bvp
{


class DM_Basis : public GenericBasisDiscretizationModule<10, TypeList< PBasis<3,3>, PBasis<4,4> >,
                                                             TypeList<methods::basis::CDD2> >
{
public:

    template<class IBASIS>
    class Problem
       : public BasisDiscretizedProblem< WaveletTL::BiharmonicEquation1D<IBASIS>, MathTL::Function<1> >
    {
        IBASIS* createWaveletSystem(int jmax) override;

        WaveletTL::BiharmonicEquation1D<IBASIS>*
        createDiscretizedEquation(const MathTL::Function<1>& rhs,
                                  const IBASIS& basis) override;
    };

};

}
}

#endif // INTERVAL_BIHARMONIC_BVP_DM_BASIS_H

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


#ifndef METHODS_BASIS_CDD2_H
#define METHODS_BASIS_CDD2_H


#include "implementation_cdd2.h"
#include "methods/method_base.h"


namespace methods
{
namespace basis
{

class CDD2 : public MethodBase<CDD2, 5, 4>
{
public:
    static QString name()
    {
        return QStringLiteral("CDD2");
    }


    static void addParametersTo(ParameterWidget& widget)
    {

    }


    template<class CPROBLEM>
    static void solve(CPROBLEM& cproblem, const GuiInputData& input,
                      MathTL::InfiniteVector<double, typename CPROBLEM::Index>& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        const double nu = cproblem.norm_Ainv() * cproblem.F_norm();
        CDD2_GuiSolve(cproblem, nu, input.epsilon, coeffs, logger, input.jmax, WaveletTL::CDD1);
    }
};


}
}



#endif // METHODS_BASIS_CDD2_H

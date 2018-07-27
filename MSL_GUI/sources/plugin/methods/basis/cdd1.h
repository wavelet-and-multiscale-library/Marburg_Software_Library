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


#ifndef METHODS_BASIS_CDD1_H
#define METHODS_BASIS_CDD1_H

#define _WAVELETTL_CDD1_VERBOSITY 1

#include "WaveletTL/adaptive/cdd1.h"

#include "methods/method_base.h"


namespace methods
{
namespace basis
{

class CDD1 : public MethodBase<CDD1, 5, 3>
{
public:
    static QString name()
    {
        return QStringLiteral("CDD1");
    }


    static void addParametersTo(ParameterWidget& widget)
    {
        widget.addBoolParameter("Coarsening", false);
    }


    template<class CPROBLEM>
    static void solve(CPROBLEM& cproblem, const GuiInputData& input,
                      MathTL::InfiniteVector<double, typename CPROBLEM::Index>& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        bool coarsening = boolParameter("Coarsening");
        WaveletTL::CDD1_SOLVE(cproblem, input.epsilon, coeffs, logger, input.jmax, coarsening);
    }
};


}

}

#endif // METHODS_BASIS_CDD1_H

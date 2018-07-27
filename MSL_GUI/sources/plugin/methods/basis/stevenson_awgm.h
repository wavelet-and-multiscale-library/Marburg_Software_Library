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


#ifndef METHODS_BASIS_STEVENSON_AWGM_H
#define METHODS_BASIS_STEVENSON_AWGM_H


#include "WaveletTL/adaptive/stevenson_AWGM.h"

#include "methods/method_base.h"


namespace methods
{
namespace basis
{

class StevensonAWGM : public MethodBase<StevensonAWGM, 5, 3>
{
public:
    static QString name()
    {
        return QStringLiteral("Stevenson AWGM");
    }


    static void addParametersTo(ParameterWidget& widget)
    {
        widget.addDoubleParameter("alpha", 0.9, 0.05, 0.95, QChar(0x03B1), 0.05);
        widget.addDoubleParameter("omega", 0.01, 0.01, 0.99, QChar(0x03C9), 0.01);
        widget.addDoubleParameter("gamma", 0.01, 0.01, 0.99, QChar(0x03B3), 0.01);
        widget.addDoubleParameter("theta", 0.4, 0.1, 99.9, QChar(0x03B8), 0.1);
    }


    template<class CPROBLEM>
    static void solve(CPROBLEM& cproblem, const GuiInputData& input,
                      MathTL::InfiniteVector<double, typename CPROBLEM::Index>& coeffs,
                      MathTL::AbstractConvergenceLogger& logger)
    {
        WaveletTL::AWGM_SOLVE(cproblem, input.epsilon, coeffs, input.jmax, logger,
                              doubleParameter("alpha"),
                              doubleParameter("omega"),
                              doubleParameter("gamma"),
                              doubleParameter("theta"));
    }
};


}
}

#endif // METHODS_BASIS_STEVENSON_AWGM_H


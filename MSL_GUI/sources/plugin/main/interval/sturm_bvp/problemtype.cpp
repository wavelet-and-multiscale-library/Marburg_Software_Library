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


#include "problemtype.h"

namespace interval
{
namespace sturm_bvp
{

void ProblemType::setupNamesAndDescription()
{
    problemDescription_ = QString("Sturm boundary value problem on the interval (0,1) with homogeneous Dirichlet boundary conditions:\n"
                                  "\n"
                                  "Find y(t) with\n"
                                  "\n"
                                  "  -(p%1y')'(t) + q(t)%1y(t) = g(t),   0 %2 t %2 1\n"
                                  "\n"
                                  "and y(0) = y(1) = 0.")
                          .arg(QChar(0x22C5), QChar(0x2264));
    functionIdentifiers_ = {"p", "q", "g"};
    functionVariables_ = {"t"};
}



void ProblemType::setupExamples()
{

}


}
}

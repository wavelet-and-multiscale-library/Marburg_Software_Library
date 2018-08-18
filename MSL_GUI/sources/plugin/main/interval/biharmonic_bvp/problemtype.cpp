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
namespace biharmonic_bvp
{

void ProblemType::setupNamesAndDescription()
{
    problemDescription_ = QString("Biharmonic boundary value problem on the interval (0,1) with homogeneous Dirichlet boundary conditions: "
                                  "Find u(t) with\n"
                                  "\n"
                                  "    -u%1%2%3(t) = g(t),   0 %4 t %4 1    and u(0) = u(1) = u'(0) = u'(1) = 0.")
                          .arg(QChar(0x207D), QChar(0x2074), QChar(0x207E), QChar(0x2264));
    functionIdentifiers_ = {"g"};
    functionVariables_ = {"t"};
}



class RHS_1
        : public Function<1, double>
{
public:
    double value(const Point<1>& p, const unsigned int component = 0) const
    {
        return  384;
    }

    void vector_value(const Point<1> &p, Vector<double>& values) const
    {
        values[0] = value(p);
    }

};



void ProblemType::setupExamples()
{
    std::array<QString, 1> definition1 = { "384 with solution u(t) = 16*(t⁴-2t³+t²)" };
    addExampleProblem("Simple polynomial example", definition1, new RHS_1);
}


}
}

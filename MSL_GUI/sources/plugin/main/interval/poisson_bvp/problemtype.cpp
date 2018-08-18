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

#include "WaveletTL/galerkin/TestProblem.h"

namespace interval
{
namespace poisson_bvp
{

void ProblemType::setupNamesAndDescription()
{
    problemDescription_ = QString("Poisson boundary value problem on the interval (0,1) with homogeneous Dirichlet boundary conditions: "
                                  "Find y(t) with\n"
                                  "\n"
                                  "    -y''(t) = g(t),   0 %1 t %1 1    and y(0) = y(1) = 0.")
                          .arg(QChar(0x2264));
    functionIdentifiers_ = {"g"};
    functionVariables_ = {"t"};
}



class RHS_BBCCDDU1
        : public Function<1, double>
{
public:
    double value(const Point<1>& p, const unsigned int component = 0) const
    {
        return (200.0-(200.0*p[0]-100.0)*(200.0*p[0]-100.0)) * exp(-100.0*(p[0]-0.5)*(p[0]-0.5));
    }

    void vector_value(const Point<1> &p, Vector<double>& values) const
    {
        values[0] = value(p);
    }

};



class RHS_BBCCDDU2
        : public Function<1, double>
{
public:
    double value(const Point<1>& p, const unsigned int component = 0) const
    {
        return  -100*exp(5*p[0])*(1-(exp(5*p[0])-1)/(exp(5.)-1))/(exp(5.)-1)+200*exp(10*p[0]) /
                ((exp(5.)-1)*(exp(5.)-1))+100*(exp(5*p[0])-1)*exp(5*p[0])/((exp(5.)-1)*(exp(5.)-1));
    }

    void vector_value(const Point<1> &p, Vector<double>& values) const
    {
        values[0] = value(p);
    }

};



void ProblemType::setupExamples()
{
//    std::array<QString, 1> definitionProblem1 = {"2"};
//    addExampleProblem("Test Problem 1", definitionProblem1, new TestProblem<1>);

    std::array<QString, 1> definitionBBCCDDU1 = {"-d²/dt² [exp(-100(t-0.5)²)]"};
    addExampleProblem("1D example 1 from [BBCCDDU]", definitionBBCCDDU1, new TestProblem<2>, new PoissonBVP<1>(new RHS_BBCCDDU1));

    std::array<QString, 1> definitionBBCCDDU2 = {"-d²/dt² [4(exp(5t)-1)/(exp(5)-1) * (1-(exp(5t)-1)/(exp(5)-1))]"};
    addExampleProblem("1D example 2 from [BBCCDDU]", definitionBBCCDDU2, new TestProblem<3>, new PoissonBVP<1>(new RHS_BBCCDDU2));


}


}
}

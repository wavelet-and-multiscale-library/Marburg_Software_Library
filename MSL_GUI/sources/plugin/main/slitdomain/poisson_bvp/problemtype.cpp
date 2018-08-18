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

#include "MathTL/numerics/corner_singularity.h"

namespace slitdomain
{
namespace poisson_bvp
{

void ProblemType::setupNamesAndDescription()
{
    problemDescription_ = QString("Poisson boundary value problem on the slit domain %1 = (-1,1)² \\ {0}x(0,1) with homogeneous Dirichlet boundary conditions: "
                                  "Find u(x,y) with\n"
                                  "\n"
                                  "    -%2u(x,y) = f(x,y),   (x,y) %3 %1    and u(x,y) = 0 on %4%1.")
                          .arg(QChar(0x03A9), QChar(0x2206), QChar(0x2208), QChar(0x2202));
    functionIdentifiers_ = {"f"};
    functionVariables_ = {"x", "y"};
}



class RHS_TestProblem1
        : public MathTL::Function<2,double>
{
public:
    RHS_TestProblem1() : csrhs_(Point<2>(0,0), 0.5, 2.0) { }

    double value(const MathTL::Point<2>& p, const unsigned int component = 0) const
    {
        return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])+5*csrhs_.value(p);
    }

    void vector_value(const MathTL::Point<2>& p, MathTL::Vector<double>& values) const
    {
        values[0] = value(p);
    }

private:
    MathTL::CornerSingularityRHS csrhs_;
};




void ProblemType::setupExamples()
{
    std::array<QString, 1> rhs_problem1 = {QString("-%1(sin(%2x)sin(%2y) + 5s(x,y)) with corner singularity s").arg(QChar(0x2206), QChar(0x03C0))};
    addExampleProblem("Sin-problem with Corner Singularity", rhs_problem1, new MathTL::PoissonBVP<2>(new RHS_TestProblem1));
}

}
}

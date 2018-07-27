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

#include "WaveletTL/galerkin/test_rhs_cube.h"

namespace cube
{
namespace poisson_bvp
{

void ProblemType::setupNamesAndDescription()
{
    problemDescription_ = QString("Poisson boundary value problem on the cube %1 = (0,1)² with homogeneous Dirichlet boundary conditions:\n"
                                  "\n"
                                  "Find u(x,y) with\n"
                                  "\n"
                                  "  -%2u(x,y) = f(x,y),   (x,y) %3 %1\n"
                                  "\n"
                                  "and u(x,y) = 0 on %4%1.")
                          .arg(QChar(0x03A9), QChar(0x2206), QChar(0x2208), QChar(0x2202));
    functionIdentifiers_ = {"f"};
    functionVariables_ = {"x", "y"};
}



void ProblemType::setupExamples()
{
    std::array<QString, 1> rhs_problem1 = {QString("-%1(x(1-x)y(1-y))").arg(QChar(0x2206))};
    MathTL::Function<2>* rhs1 = new TestRHS<1>();
    addExampleProblem("Test Problem Cube 1", rhs_problem1, new MathTL::PoissonBVP<2>(rhs1));

    std::array<QString, 1> rhs_problem2 = {QString("-%1 exp(-50*((x-0.5)²+(y-0.5)²))").arg(QChar(0x2206))};
    MathTL::Function<2>* rhs2 = new TestRHS<2>();
    addExampleProblem("Test Problem Cube 2", rhs_problem2, new MathTL::PoissonBVP<2>(rhs2));

    std::array<QString, 1> rhs_problem3 = {QString("-%1(x(1-x)²y²(1-y))").arg(QChar(0x2206))};
    MathTL::Function<2>* rhs3 = new TestRHS<3>();
    addExampleProblem("Test Problem Cube 3", rhs_problem3, new MathTL::PoissonBVP<2>(rhs3));

}

}
}


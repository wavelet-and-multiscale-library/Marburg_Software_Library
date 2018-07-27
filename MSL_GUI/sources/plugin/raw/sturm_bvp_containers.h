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


#ifndef STURM_BVP_CONTAINERS_H
#define STURM_BVP_CONTAINERS_H

#include "raw_problem_container.h"
#include "parsed/parsed_sturm_bvp.h"



class Poisson1DContainer : public RawProblemContainer<1, MathTL::SimpleSturmBVP, Parsed1DFunction>
{
    MathTL::SimpleSturmBVP* createParsedProblem(Parsed1DFunction* g_parsed) override;
};




class SturmBvpContainer : public RawProblemContainer<1, MathTL::SimpleSturmBVP,
                                                        Parsed1DFunction,
                                                        Parsed1DFunction,
                                                        Parsed1DFunction >
{
    MathTL::SimpleSturmBVP* createParsedProblem(Parsed1DFunction* p_parsed,
                                                Parsed1DFunction* q_parsed,
                                                Parsed1DFunction* g_parsed) override;
};

#endif // STURM_BVP_CONTAINERS_H

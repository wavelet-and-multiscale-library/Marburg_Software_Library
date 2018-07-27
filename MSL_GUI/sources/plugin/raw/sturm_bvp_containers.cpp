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


#include "sturm_bvp_containers.h"


MathTL::SimpleSturmBVP* Poisson1DContainer::createParsedProblem(Parsed1DFunction* g_parsed)
{
    return new Parsed1DPoissonBVP(g_parsed);
}



MathTL::SimpleSturmBVP* SturmBvpContainer::createParsedProblem(Parsed1DFunction* p_parsed,
                                                               Parsed1DFunction* q_parsed,
                                                               Parsed1DFunction* g_parsed)
{
    return new ParsedSturmBVP(p_parsed, q_parsed, g_parsed);
}

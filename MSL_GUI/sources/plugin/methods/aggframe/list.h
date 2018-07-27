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


#ifndef METHODS_AGGFRAME_LIST_H
#define METHODS_AGGFRAME_LIST_H

#include "misc/typelist.h"

#include "steepest_descent_solver.h"
#include "richardson_solver.h"
#include "additive_schwarz_solver.h"
#include "multiplicative_schwarz_solver.h"

namespace methods
{
namespace aggframe
{

typedef TypeList< RichardsonSolver,
                  SteepestDescentSolver,
                  AdditiveSchwarzSolver,
                  MultiplicativeSchwarzSolver> List;

}
}


#endif // METHODS_AGGFRAME_LIST_H

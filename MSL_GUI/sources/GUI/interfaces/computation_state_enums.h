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


#ifndef COMPUTATION_STATE_ENUMS_H
#define COMPUTATION_STATE_ENUMS_H

namespace CoeffComputationState
{
enum Enum
{
    RUNNING,
    ERROR_PROBLEM_CREATION,
    ERROR_NO_SOLUTION_CREATED,
    ERROR_COEFF_COMPUTATION,
    ABORTED,
    COMPLETE,
};
}


namespace PlotComputationState
{
enum Enum
{
    NOT_STARTED_YET,
    RUNNING,
    ERROR,
    ABORTED,
    COMPLETE
};
}


#endif // COMPUTATION_STATE_ENUMS_H

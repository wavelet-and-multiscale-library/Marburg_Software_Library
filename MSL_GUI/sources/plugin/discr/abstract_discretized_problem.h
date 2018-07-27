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


#ifndef ABSTRACT_DISCRETIZED_PROBLEM_H
#define ABSTRACT_DISCRETIZED_PROBLEM_H

class QString;
struct GuiInputData;
class AbstractSolution;
class AbstractGuiCommunicator;



class AbstractDiscretizedProblem
{
public:
    virtual ~AbstractDiscretizedProblem() {}

    virtual AbstractSolution* computeSolution(const QString& method, double epsilon, int jmax,
                                              const AbstractGuiCommunicator* communicator) = 0;

    /* The member function "canBeReusedFor" determines whether the given problem
     * can be reused for a new computation (defined by "newInput") without or
     * only minor modifications. In that case, a call to "updateTo(newInput)"
     * performs these modifications. */
    virtual bool canBeReusedFor(const GuiInputData& newInput) const = 0;
    virtual void updateTo(const GuiInputData& newInput) = 0;

    virtual double norm_A() = 0;
    virtual double norm_Ainv() = 0;
};

#endif // ABSTRACT_DISCRETIZED_PROBLEM_H

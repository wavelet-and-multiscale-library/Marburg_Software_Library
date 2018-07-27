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


#ifndef DISCRETIZATION_MODULE_GUIDATA_H
#define DISCRETIZATION_MODULE_GUIDATA_H

#include <QStringList>

struct DiscretizationModuleGuiData
{
    QStringList basis1DList;
    QStringList methodList;

    unsigned int jMaxStandard;

    std::vector<int> epsilonStandardExponents;
    // epsilonStandardExponents.at(i) is supposed to be the standard exponent (to base 10)
    // for epsilon belonging to the i-th method in methodList
};

#endif // DISCRETIZATION_MODULE_GUIDATA_H

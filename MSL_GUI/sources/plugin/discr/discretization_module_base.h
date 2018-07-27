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


#ifndef DISCRETIZATION_MODULE_BASE_H
#define DISCRETIZATION_MODULE_BASE_H


// TODO: find a better place for the following line:
#define _WAVELETTL_USE_TFRAME 1


#include "misc/typelist.h"
#include "misc/instantiate.h"
#include "discretization_module_guidata.h"
#include "misc/convert_basis_to_qstring.h"
#include "generic_discretized_problem.h"
#include "solution/generic_solution.h"


template<unsigned int jmax_standard, class LIST_1D_BASES, class METHOD_LIST>
class DiscretizationModuleBase;


template<unsigned int jmax_standard, class ... BASES1D, class ... METHODS>
class DiscretizationModuleBase< jmax_standard, TypeList<BASES1D...>, TypeList<METHODS...> >
{
public:

    typedef TypeList<BASES1D...> List_1D_Bases;

    static DiscretizationModuleGuiData getGuiData(unsigned int problemDimension);
};



/*##################################################################################################
    Implementation
##################################################################################################*/


template<unsigned int jmax_standard, class ... BASES1D, class ... METHODS>
DiscretizationModuleGuiData
DiscretizationModuleBase< jmax_standard, TypeList<BASES1D...>, TypeList<METHODS...> >::getGuiData(unsigned int problemDimension)
{
    DiscretizationModuleGuiData guiData;
    guiData.jMaxStandard = jmax_standard;

    (void) std::initializer_list<int>{
                          (guiData.basis1DList.append(Convert<BASES1D>::toQString()), 0)...};
    (void) std::initializer_list<int>{(guiData.methodList.append(METHODS::name()), 0)...};

    if (problemDimension == 1) {
        guiData.epsilonStandardExponents = {METHODS::standardEpsilonExponent1D...};
    }
    else {
        guiData.epsilonStandardExponents = {METHODS::standardEpsilonExponent2D...};
    }

    return guiData;
}


#endif // DISCRETIZATION_MODULE_BASE_H

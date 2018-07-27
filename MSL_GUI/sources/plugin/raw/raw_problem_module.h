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


#ifndef RAW_PROBLEM_MODULE_H
#define RAW_PROBLEM_MODULE_H

#include <array>
#include <QString>



//template<unsigned int DIM, unsigned int N_FUNCTIONS, class RP_CONTAINER1, class ... RP_CONTAINERS>
template<class RP_CONTAINER1, class ... RP_CONTAINERS>
class RawProblemModule : private RP_CONTAINER1, private RP_CONTAINERS...
{
public:

    static const unsigned int dimension = RP_CONTAINER1::dimension;
    static const std::size_t n_functions = RP_CONTAINER1::n_functions;

    template<class REQUESTED_RP>
    const REQUESTED_RP& getExampleProblem(unsigned int problemIndex) const
    {
        return ContainerFor<REQUESTED_RP>::getExampleProblem(problemIndex);
    }


    template<class REQUESTED_RP>
    const REQUESTED_RP& getParsedProblem(const std::array<QString,n_functions>& funcDefinitions)
    {
        return ContainerFor<REQUESTED_RP>::getParsedProblem(funcDefinitions);
    }


    void setupParsedProblemAndFunctions(const std::array<QString,dimension>& funcVariables)
    {
        RP_CONTAINER1::setupParsedProblemAndFunctions(funcVariables);
        (void) std::initializer_list<int>{
                              (RP_CONTAINERS::setupParsedProblemAndFunctions(funcVariables), 0)...};
    }


    void addExampleProblem(const typename RP_CONTAINER1::RawProblemType* exampleFirstType,
                           const typename RP_CONTAINERS::RawProblemType* ... exampleOtherTypes)
    {
        RP_CONTAINER1::addExampleProblem(exampleFirstType);
        (void) std::initializer_list<int>{
                                       (RP_CONTAINERS::addExampleProblem(exampleOtherTypes), 0)...};
    }

private:

    template<class REQUESTED_RP, bool is_base, class ... CONTAINERS>
    struct ContainerTypeFor;

    //specializations:

    template<class REQUESTED_RP, class CONTAINER1, class ... CONTAINERS>
    struct ContainerTypeFor<REQUESTED_RP, true, CONTAINER1, CONTAINERS...>
    {
        typedef CONTAINER1 type;
    };

    template<class REQUESTED_RP, class CONTAINER1, class CONTAINER2, class ... CONTAINERS>
    struct ContainerTypeFor<REQUESTED_RP, false, CONTAINER1, CONTAINER2, CONTAINERS...>
    {
        typedef typename ContainerTypeFor<REQUESTED_RP,
        std::is_base_of<REQUESTED_RP, typename CONTAINER2::RawProblemType>::value,
        CONTAINER2, CONTAINERS...>::type
        type;
    };


    template<class REQUESTED_RP>
    using ContainerFor = struct ContainerTypeFor<REQUESTED_RP,
                       std::is_base_of<REQUESTED_RP, typename RP_CONTAINER1::RawProblemType>::value,
                       RP_CONTAINER1, RP_CONTAINERS...>::type;

};

#endif // RAW_PROBLEM_MODULE_H

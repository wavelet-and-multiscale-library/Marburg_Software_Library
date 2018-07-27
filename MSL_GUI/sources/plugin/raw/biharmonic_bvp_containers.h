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


#ifndef BIHARMONIC_BVP_CONTAINERS_H
#define BIHARMONIC_BVP_CONTAINERS_H

#include "raw_problem_container.h"
#include "parsed/parsed_function.h"

#include "MathTL/numerics/bvp.h"


template<unsigned int DIM>
class BiharmonicBvpContainer : public RawProblemContainer< DIM, MathTL::BiharmonicBVP<DIM>,
                                                           ParsedFunction<DIM> >
{
    MathTL::BiharmonicBVP<DIM>* createParsedProblem(ParsedFunction<DIM>* f_parsed) override
    {
        return new MathTL::BiharmonicBVP<DIM>(f_parsed);
    }
};


// To be used for FrameTL::SimpleBiharmonicEquation(const Functional<IBASIS,DIM>* rhs, ...)
// We postpone the construction of the Functional via the constructor
// Functional(const Function<DIM>* g, const AggregatedFrame<IBASIS, DIM>* frame) to the method
// createDiscretizedEquation(const MathTL::Function<DIM>& rawProblem, const FrameTL::AggregatedFrame<IBASIS,DIM>& frame)
// of the corresponding discretization module.
template<unsigned int DIM>
class SimpleBiharmonicBvpContainer : public RawProblemContainer< DIM, MathTL::Function<DIM>,
                                                           ParsedFunction<DIM> >
{
    MathTL::Function<DIM>* createParsedProblem(ParsedFunction<DIM>* f_parsed) override
    {
        return f_parsed;
    }
};


#endif // BIHARMONIC_BVP_CONTAINERS_H

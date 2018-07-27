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


#ifndef PARSED_FUNCTION_H
#define PARSED_FUNCTION_H

#include "parsedfunction_base.h"
#include "MathTL/utils/function.h"


template<unsigned int DIM>
class ParsedFunction : public ParsedFunctionBase, public MathTL::Function<DIM>
{
public:
    ParsedFunction(const std::array<QString, DIM>& funcVariables,
                   const unsigned int n_components = 1,
                   const double initial_time = 0.0);

    double value(const MathTL::Point<DIM>& p, const unsigned int component = 0) const override;

    void vector_value(const MathTL::Point<DIM>& p, MathTL::Vector<double>& values) const override;

private:
    mutable MathTL::Point<DIM> evalPoint_;
};



/*##################################################################################################
    Implementation
##################################################################################################*/


template<unsigned int DIM>
ParsedFunction<DIM>::ParsedFunction(const std::array<QString, DIM>& funcVariables,
                                    const unsigned int n_components,
                                    const double initial_time)
    : ParsedFunctionBase(), MathTL::Function<DIM>(n_components, initial_time)
{
    for (int i = 0; i < DIM; ++i)
        defineVariable(funcVariables[i], &evalPoint_[i]);
}



template<unsigned int DIM>
double ParsedFunction<DIM>::value(const MathTL::Point<DIM>& p, const unsigned int component) const
{
//    assert(component == MathTL::Function<DIM>::n_components-1);
    evalPoint_ = p;
    return parser_.Eval();
}



template<unsigned int DIM>
void ParsedFunction<DIM>::vector_value(const MathTL::Point<DIM>& p, MathTL::Vector<double>& values)
const
{
    evalPoint_ = p;
    int number_of_expressions;
    double* v = parser_.Eval(number_of_expressions);
//    assert(values.size() >= number_of_expressions);

    for (int i = 0; i < number_of_expressions; ++i)
        values[i] = v[i];
}



#endif // PARSED_FUNCTION_H

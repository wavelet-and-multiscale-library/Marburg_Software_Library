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


#ifndef PARSED_EXPRESSION_CHECKER_H
#define PARSED_EXPRESSION_CHECKER_H

#include "parsed_function.h"
#include "misc/string_conversion.h"



template<unsigned int DIM, unsigned int N_COMPONENTS = 1>
class ParsedExpressionChecker : private ParsedFunction<DIM>
{
public:
    ParsedExpressionChecker(const std::array<QString,DIM>& funcVariables)
        : ParsedFunction<DIM>(funcVariables, N_COMPONENTS)
    {

    }

    bool detectsSyntaxErrors(const QString& expression, QString& out_errors);
};



/*##################################################################################################
    Implementation
##################################################################################################*/


template<unsigned int DIM, unsigned int N_COMPONENTS>
bool ParsedExpressionChecker<DIM, N_COMPONENTS>::detectsSyntaxErrors(const QString& expression,
                                                                     QString& out_errors)
{
    try
    {
        this->setNewDefinition(expression);
        this->parser_.Eval();
    }
    catch(const mu::ParserError& e)
    {
        out_errors = qStringFromMuString(e.GetMsg());
        return true;
    }

    if (this->parser_.GetNumResults() != int (N_COMPONENTS))
    {
        QString error = QString("Wrong number of components!\n"
                                "Number of components in definition is: %1\n"
                                "Number of components should be: %2\n"
                                "(components of vector-valued functions are separated by \'%3\')")
                        .arg(this->parser_.GetNumResults()).arg(N_COMPONENTS)
                        .arg(qStringFromMuChar(this->parser_.GetArgSep()));
        out_errors = error;
        return true;
    }
    else
    {
        return false;
    }
}


#endif // PARSED_EXPRESSION_CHECKER_H

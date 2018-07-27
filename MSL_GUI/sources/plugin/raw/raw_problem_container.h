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


#ifndef RAW_PROBLEM_CONTAINER_H
#define RAW_PROBLEM_CONTAINER_H

#include <QStringList>

#include "parsed/parsedfunction_base.h"


template<unsigned int DIM, class RAW_PROBLEM, class ... PARSED_FUNCTIONS>
class RawProblemContainer
{
public:

    typedef RAW_PROBLEM RawProblemType;

    static const unsigned int dimension = DIM;
    static const std::size_t n_functions = sizeof...(PARSED_FUNCTIONS);


    virtual ~RawProblemContainer() {}


    const RAW_PROBLEM& getExampleProblem(unsigned int problemIndex) const
    {
        return *(examples_.at(problemIndex));
    }


    const RAW_PROBLEM& getParsedProblem(const std::array<QString, n_functions>& funcDefinitions)
    {
        for (int i = 0; i < n_functions; ++i)
            (parsedProblemFunctions_[i])->setNewDefinition(funcDefinitions[i]);

        return *parsedProblem_;
    }


    void addExampleProblem(const RAW_PROBLEM* example)
    {
        examples_.push_back(example);
    }


    void setupParsedProblemAndFunctions(const std::array<QString, DIM>& funcVariables)
    {
        setupHelper(new PARSED_FUNCTIONS(funcVariables)...);
    }


private:

    void setupHelper(PARSED_FUNCTIONS* ... parsedFunctions)
    {
        parsedProblem_ = createParsedProblem(parsedFunctions...);
        parsedProblemFunctions_ = {(ParsedFunctionBase*) parsedFunctions...};
        //parsedProblemFunctions_ = std::initializer_list<ParsedFunctionBase*>{parsedFunctions...};
    }

    virtual RAW_PROBLEM* createParsedProblem(PARSED_FUNCTIONS* ... parsedFunctions) = 0;

    std::vector<const RAW_PROBLEM*> examples_;

    const RAW_PROBLEM* parsedProblem_;

    std::array<ParsedFunctionBase*, n_functions> parsedProblemFunctions_;
//    std::vector<ParsedFunctionBase*> parsedProblemFunctions_;
};


#endif // RAW_PROBLEM_CONTAINER_H

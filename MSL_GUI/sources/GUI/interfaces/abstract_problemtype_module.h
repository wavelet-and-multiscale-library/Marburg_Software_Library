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


#ifndef ABSTRACT_PROBLEMTYPE_MODULE_H
#define ABSTRACT_PROBLEMTYPE_MODULE_H


class QString;
class QStringList;
class AbstractDiscretizedProblem;
struct GuiInputData;

#include <vector>


class AbstractProblemTypeModule
{
public:
    virtual ~AbstractProblemTypeModule() {}

    virtual void initialize(const QString& domainTag, const QString& problemTypeTag) = 0;

    virtual const QString& getDomainTag() const = 0;
    virtual const QString& getProblemTypeTag() const = 0;

    virtual const QString& getProblemTypeDescription() const = 0;
    virtual const QStringList& getFunctionNames() const = 0;

    virtual const QStringList& getExampleNames() const = 0;
    virtual const std::vector<QStringList>& getExampleDefinitions() const = 0;
    virtual const QStringList& getDiscretizationTypeList() const = 0;

    virtual const QStringList& get1DBasisList(int discretizationIndex) const = 0;
    virtual const QStringList& getMethodList(int discretizationIndex) const = 0;
    virtual unsigned int getJmaxStandard(int discretizationIndex) const = 0;
    virtual int getEpsilonStandardExponent(int discretizationIndex, int methodIndex) const = 0;

    virtual bool detectsErrorsInProblemDefinition(const QString& problemName,
                                                  const QStringList& functionDefs,
                                                  QString& out_errors) const = 0;

    virtual AbstractDiscretizedProblem* createDiscretizedProblem(const GuiInputData& input) = 0;
};

#endif // ABSTRACT_PROBLEMTYPE_MODULE_H

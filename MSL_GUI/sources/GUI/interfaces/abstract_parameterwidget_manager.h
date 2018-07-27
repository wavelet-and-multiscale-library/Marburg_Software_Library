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


#ifndef ABSTRACT_PARAMETERWIDGET_MANAGER_H
#define ABSTRACT_PARAMETERWIDGET_MANAGER_H


class QString;
class QStackedWidget;

class AbstractParameterWidgetManager
{
public:
    virtual ~AbstractParameterWidgetManager() {}

    virtual void addAllWidgetsTo(QStackedWidget* stackedWidget) = 0;
    virtual void displayWidgetForMethod(const QString& methodName) const = 0;

    virtual QString getParamStringForMethod(const QString& methodName) const = 0;
};

#endif // ABSTRACT_PARAMETERWIDGET_MANAGER_H

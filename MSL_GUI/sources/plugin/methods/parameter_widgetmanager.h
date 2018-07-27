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


#ifndef METHODS_PARAMETER_WIDGETMANAGER_H
#define METHODS_PARAMETER_WIDGETMANAGER_H


#include <map>
#include "GUI/interfaces/abstract_parameterwidget_manager.h"
#include "misc/typelist.h"


namespace methods
{

class ParameterWidget;


class ParameterWidgetManager : public AbstractParameterWidgetManager
{
public:
    ParameterWidgetManager();

    void addAllWidgetsTo(QStackedWidget* stackedWidget) override;
    void displayWidgetForMethod(const QString& methodName) const override;

    QString getParamStringForMethod(const QString& methodName) const override;

private:
    template<class METHOD>
    void addWidgetForMethod();

    template<class ... METHODS>
    void addWidgetsOnList(const TypeList<METHODS...>& list);

    std::map<QString, ParameterWidget*> widgetMap_; // maps method names to parameter widgets
    QStackedWidget* stackedWidget_;
};

}

#endif // METHODS_PARAMETER_WIDGETMANAGER_H

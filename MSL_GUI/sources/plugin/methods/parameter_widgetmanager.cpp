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


#include "parameter_widgetmanager.h"
#include <QStackedWidget>


#include "basis/list.h"
#include "qframe/list.h"
#include "aggframe/list.h"


namespace methods
{

ParameterWidgetManager::ParameterWidgetManager()
{
    addWidgetsOnList(basis::List());
    addWidgetsOnList(qframe::List());
    addWidgetsOnList(aggframe::List());
}



void ParameterWidgetManager::addAllWidgetsTo(QStackedWidget* stackedWidget)
{
    stackedWidget_ = stackedWidget;

    for (const auto& pair : widgetMap_) {
       stackedWidget->addWidget(pair.second);
    }
}



void ParameterWidgetManager::displayWidgetForMethod(const QString& methodName) const
{
    ParameterWidget* widget = widgetMap_.at(methodName);
    stackedWidget_->setCurrentWidget(widget);
}



QString ParameterWidgetManager::getParamStringForMethod(const QString &methodName) const
{
    return widgetMap_.at(methodName)->getParamString();
}



template<class METHOD>
void ParameterWidgetManager::addWidgetForMethod()
{
    ParameterWidget* widget = METHOD::widget();
    widgetMap_[METHOD::name()] = widget;
}



template<class ... METHODS>
void ParameterWidgetManager::addWidgetsOnList(const TypeList<METHODS...>& list)
{
    (void) list;
     // A trick for running through all METHODS without using recursion:
    (void) std::initializer_list<int>{(addWidgetForMethod<METHODS>(), 0)...};
}

}

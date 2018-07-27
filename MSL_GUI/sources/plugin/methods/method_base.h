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


#ifndef METHODS_METHOD_BASE_H
#define METHODS_METHOD_BASE_H

#include "parameter_widget.h"
#include "GUI/interfaces/gui_inputdata.h"

#include "MathTL/utils/convergence_logger.h"

namespace methods
{

template <class METHOD, int EPS_EXPONENT_1D, int EPS_EXPONENT_2D>
class MethodBase
{
public:
    static const int standardEpsilonExponent1D = EPS_EXPONENT_1D;
    static const int standardEpsilonExponent2D = EPS_EXPONENT_2D;


    static ParameterWidget* widget()
    {
        static ParameterWidget* widget(new ParameterWidget(METHOD::name()));
        // ownership of widget is later transferred to a QStackedWidget by ParameterWidgetManager.
        static bool widgetHasBeenSetUp = false;

        if (widgetHasBeenSetUp)
        {
            return widget;
        }
        else
        {
            METHOD::addParametersTo(*widget);
            widget->setTheLayout();
            widgetHasBeenSetUp = true;
            return widget;
        }
    }

    static double doubleParameter(const QString& identifier)
    {
        return widget()->getDoubleParameter(identifier);
    }

    static int intParameter(const QString& identifier)
    {
        return widget()->getIntParameter(identifier);
    }

    static bool boolParameter(const QString& identifier)
    {
        return widget()->getBoolParameter(identifier);
    }

    static QString stringParameter(const QString& identifier)
    {
        return widget()->getStringParameter(identifier);
    }

// Disallow creating an instance of this class:
private:
    MethodBase() {}
};


}

#endif // METHODS_METHOD_BASE_H

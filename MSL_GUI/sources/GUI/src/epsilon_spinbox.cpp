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


#include "epsilon_spinbox.h"

EpsilonSpinBox::EpsilonSpinBox(QWidget* parent)
    : QSpinBox(parent)
{

}


QString EpsilonSpinBox::textFromValue(int val) const
{
    QString exponent;

    if (val >= 10)
        exponent = QChar(0x00B9);

    switch (val%10)
    {
    case 1:
        exponent += QChar(0x00B9);
        break;
    case 2:
        exponent += QChar(0x00B2);
        break;
    case 3:
        exponent += QChar(0x00B3);
        break;
    case 4:
        exponent += QChar(0x2074);
        break;
    case 5:
        exponent += QChar(0x2075);
        break;
    case 6:
        exponent += QChar(0x2076);
        break;
    case 7:
        exponent += QChar(0x2077);
        break;
    case 8:
        exponent += QChar(0x2078);
        break;
    case 9:
        exponent += QChar(0x2079);
        break;
    case 0:
        exponent += QChar(0x2070);
        break;
    }

    return QString("10%1%2").arg(QChar(0x207B), exponent);
}

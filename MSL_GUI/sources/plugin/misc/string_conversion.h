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


#ifndef STRING_CONVERSION_H
#define STRING_CONVERSION_H

#include <QString>

#if defined(_UNICODE) // In this case mu::string_type is defined as std::wstring

inline std::wstring muStringFromQString(const QString& qstr)
{
    return qstr.toStdWString();
}


inline QString qStringFromMuString(const std::wstring& str)
{
    return QString::fromStdWString(str);
}


inline QString qStringFromMuChar(wchar_t c)
{
    return QString::fromWCharArray(&c,1);
}

#else   // Otherwise mu::string_type is defined as std::string

inline std::string muStringFromQString(const QString& qstr)
{
    return qstr.toStdString();
}


inline QString qStringFromMuString(const std::string& str)
{
    return QString::fromStdString(str);
}


inline QString qStringFromMuChar(char c)
{
    return QString(QChar(c));
}

#endif // _UNICODE


inline const char* cStringFromQString(const QString& qstr)
{
    QByteArray ba = qstr.toLocal8Bit();
    return ba.data();
}


#endif // STRING_CONVERSION_H

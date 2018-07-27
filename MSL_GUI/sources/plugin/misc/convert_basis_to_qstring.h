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


#ifndef CONVERT_BASIS_TO_QSTRING_H
#define CONVERT_BASIS_TO_QSTRING_H


#include <QString>

#include "WaveletTL/interval/ds_bio.h"


namespace WaveletTL {

template <int d, int dT>
class PBasis;

template <int d, int dT, DSBiorthogonalizationMethod BIO>
class DSBasis;

template <int d, int dT>
class PQFrame;

}



template<class C>
struct Convert;


template<int d, int dT>
struct Convert< WaveletTL::PBasis<d, dT> >
{
    static const QString& toQString()
    {
        static const QString str = QString("Primbs Basis (d=%1, d%3=%2)").arg(d).arg(dT).arg(QChar(0x0303));
        return str;
    }
};


template<int d, int dT, WaveletTL::DSBiorthogonalizationMethod BIO>
struct Convert< WaveletTL::DSBasis<d, dT, BIO> >
{
    static const QString& toQString()
    {
        static const QString str = QString("DS Basis (d=%1, d%3=%2)").arg(d).arg(dT).arg(QChar(0x0303));
        return str;
    }
};


template<int d, int dT>
struct Convert< WaveletTL::PQFrame<d, dT> >
{
    static const QString& toQString()
    {
        static const QString str = QString("Primbs Quarklet Frame (d=%1, d%3=%2)").arg(d).arg(dT).arg(QChar(0x0303));
        return str;
    }
};


#endif // CONVERT_BASIS_TO_QSTRING_H

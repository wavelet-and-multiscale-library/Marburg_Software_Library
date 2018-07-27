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


#ifndef INTERVAL_BASES_LIST_H
#define INTERVAL_BASES_LIST_H

#include "typelist.h"

#include "WaveletTL/interval/p_basis.h"
#include "WaveletTL/interval/ds_basis.h"

using WaveletTL::PBasis;
using WaveletTL::DSBasis;


typedef TypeList< PBasis<3,3>, PBasis<4,4>, PBasis<2,2>,
                  DSBasis<3,3>, DSBasis<4,4>, DSBasis<2,2> > IntervalBasesList_all;


typedef TypeList< PBasis<3,3>, PBasis<4,4>, PBasis<2,2> > IntervalBasesList_Primbs;


#endif // INTERVAL_BASES_LIST_H

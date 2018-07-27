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


#include "parsed_sturm_bvp.h"

ParsedSturmBVP::ParsedSturmBVP(const Parsed1DFunction* p,
                               const Parsed1DFunction* q,
                               const Parsed1DFunction* g)
    : p_(p), q_(q), g_(g)
{

}



double ParsedSturmBVP::p(const double t) const
{
    return p_->value(t);
}



double ParsedSturmBVP::p_prime(const double t) const
{
    return 0.0;     // TODO: implement this properly
}



double ParsedSturmBVP::q(const double t) const
{
    return q_->value(t);
}



double ParsedSturmBVP::g(const double t) const
{
    return g_->value(t);
}



bool ParsedSturmBVP::bc_left() const
{
    return true;    /* boundary conditions left: true */
}



bool ParsedSturmBVP::bc_right() const
{
    return true;    /* boundary conditions right: true */
}




// implementation for class Parsed1DPoissonBVP:

Parsed1DPoissonBVP::Parsed1DPoissonBVP(const Parsed1DFunction* g)
    : ParsedSturmBVP(nullptr, nullptr, g)
{

}



double Parsed1DPoissonBVP::p(const double t) const
{
    return 1.0;
}



double Parsed1DPoissonBVP::p_prime(const double t) const
{
    return 0.0;
}



double Parsed1DPoissonBVP::q(const double t) const
{
    return 0.0;
}


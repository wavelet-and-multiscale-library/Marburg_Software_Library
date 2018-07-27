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


#ifndef PARSEDSTURMBVP_H
#define PARSEDSTURMBVP_H

#include "MathTL/numerics/sturm_bvp.h"

#include "parsed_1d_function.h"


class ParsedSturmBVP : public MathTL::SimpleSturmBVP
{
public:
    ParsedSturmBVP(const Parsed1DFunction* p,
                   const Parsed1DFunction* q,
                   const Parsed1DFunction* g);


    double p(const double t) const override;

    double p_prime(const double t) const override;

    double q(const double t) const override;

    double g(const double t) const override;

    bool bc_left() const override;

    bool bc_right() const override;


private:
    //! diffusion coefficient
    const Parsed1DFunction* p_;

    //! reaction coefficient
    const Parsed1DFunction* q_;

    //! right-hand side
    const Parsed1DFunction* g_;
};




class Parsed1DPoissonBVP : public ParsedSturmBVP
{
public:
    Parsed1DPoissonBVP(const Parsed1DFunction* g);

    double p(const double t) const override;

    double p_prime(const double t) const override;

    double q(const double t) const override;
};





#endif // PARSEDSTURMBVP_H

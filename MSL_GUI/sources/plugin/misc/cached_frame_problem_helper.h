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


#ifndef CACHED_FRAME_PROBLEM_HELPER_H
#define CACHED_FRAME_PROBLEM_HELPER_H

#include "WaveletTL/galerkin/cached_problem.h"

template <class PROBLEM>
class CachedFrameProblemHelper
{
public:

    CachedFrameProblemHelper(const PROBLEM* P,
                             const double normA = 0.0,
                             const double normAinv = 0.0)
        : cachedProblem_(P, normA, normAinv),
          cachedProblemLocal_(P, normA, normAinv)
    {

    }


    const WaveletTL::CachedProblem<PROBLEM>& cachedProblem()
    {
        lastUsedCachedProblemLocal_ = false;
        return cachedProblem_;
    }


    const WaveletTL::CachedProblemLocal<PROBLEM>& cachedProblemLocal()
    {
        lastUsedCachedProblemLocal_ = true;
        return cachedProblemLocal_;
    }


    void set_normA(double new_normA)
    {
        cachedProblem_.set_normA(new_normA);
        cachedProblemLocal_.set_normA(new_normA);
    }


    void set_normAinv(double new_normAinv)
    {
        cachedProblem_.set_normAinv(new_normAinv);
        cachedProblemLocal_.set_normAinv(new_normAinv);
    }


    bool lastUsedCachedProblemLocal() const
    {
        return lastUsedCachedProblemLocal_;
    }


    double norm_A()
    {
        if (lastUsedCachedProblemLocal_)
            return cachedProblemLocal_.norm_A();
        else
            return cachedProblem_.norm_A();
    }


    double norm_Ainv()
    {
        if (lastUsedCachedProblemLocal_)
            return cachedProblemLocal_.norm_Ainv();
        else
            return cachedProblem_.norm_Ainv();
    }


private:

    WaveletTL::CachedProblem<PROBLEM> cachedProblem_;
    WaveletTL::CachedProblemLocal<PROBLEM> cachedProblemLocal_;

    bool lastUsedCachedProblemLocal_;

};

#endif // CACHED_FRAME_PROBLEM_HELPER_H

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


#ifndef GENERIC_BASIS_DISCRETIZATION_MODULE_H
#define GENERIC_BASIS_DISCRETIZATION_MODULE_H


#include "discretization_module_base.h"
#include "misc/interval_bases_list.h"
#include "methods/basis/list.h"

#include "WaveletTL/galerkin/cached_problem.h"


template<unsigned int jmax_standard, class LIST_1D_BASES = IntervalBasesList_all, class METHOD_LIST = methods::basis::List>
class GenericBasisDiscretizationModule
        : public DiscretizationModuleBase<jmax_standard, LIST_1D_BASES, METHOD_LIST>
{

public:

    static QString getDiscretizationTypeName()
    {
        return QStringLiteral("Wavelet Riesz Basis");
    }


protected:

    template<class EQUATION>
    using SimpleSolution = GenericSolution< InfiniteVector<double, typename EQUATION::Index>,
                                            SampledMapping<EQUATION::space_dimension>,
                                            typename EQUATION::WaveletBasis >;

    template<class EQUATION, class RAW_PROBLEM>
    using BasisDiscretizedProblem
    = GenericDiscretizedProblem<WaveletTL::CachedProblem, EQUATION, RAW_PROBLEM,
                                SimpleSolution<EQUATION>, METHOD_LIST>;


    // For more complex domains like the L-domain:

    template<class EQUATION>
    using MultiPatchSolution = GenericSolution< InfiniteVector<double, typename EQUATION::Index>,
                                              Array1D<SampledMapping<EQUATION::space_dimension> >,
                                              typename EQUATION::WaveletBasis >;

    template<class EQUATION, class RAW_PROBLEM>
    using MultiPatchBasisDiscretizedProblem
    = GenericDiscretizedProblem<WaveletTL::CachedProblem, EQUATION, RAW_PROBLEM,
                                MultiPatchSolution<EQUATION>, METHOD_LIST>;
};


#endif // GENERIC_BASIS_DISCRETIZATION_MODULE_H

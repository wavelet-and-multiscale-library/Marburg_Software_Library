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


#ifndef GENERIC_QFRAME_DISCRETIZATION_MODULE_H
#define GENERIC_QFRAME_DISCRETIZATION_MODULE_H


#include "discretization_module_base.h"
#include "WaveletTL/galerkin/cached_quarklet_ldomain_problem.h"
#include "WaveletTL/galerkin/cached_quarklet_recring_problem.h"
#include "WaveletTL/galerkin/cached_quarklet_slitdomain_problem.h"
#include "WaveletTL/galerkin/cached_quarklet_problem.h"
#include "WaveletTL/galerkin/cached_quarklet_tproblem.h"

#include "misc/interval_qframes_list.h"
#include "methods/qframe/list.h"


template<unsigned int jmax_standard, class LIST_1D_BASES = IntervalQFramesList, class METHOD_LIST = methods::qframe::List>
class GenericQFrameDiscretizationModule
        : public DiscretizationModuleBase<jmax_standard, LIST_1D_BASES, METHOD_LIST>
{

public:

    static QString getDiscretizationTypeName()
    {
        return QStringLiteral("Quarklet Frame");
    }


protected:

    template<class EQUATION>
    using QFrameSolution = GenericSolution<InfiniteVector<double, typename EQUATION::Index>,
                                          Array1D<SampledMapping<EQUATION::space_dimension> >,
                                          typename EQUATION::Frame>;

    template<class EQUATION, class RAW_PROBLEM>
    using QFrameDiscretizedLDomainProblem
    = GenericDiscretizedProblem<WaveletTL::CachedQuarkletLDomainProblem, EQUATION, RAW_PROBLEM,
                                QFrameSolution<EQUATION>, METHOD_LIST>;

    template<class EQUATION, class RAW_PROBLEM>
    using QFrameDiscretizedRecRingProblem
    = GenericDiscretizedProblem<WaveletTL::CachedQuarkletRecRingProblem, EQUATION, RAW_PROBLEM,
                                QFrameSolution<EQUATION>, METHOD_LIST>;

    template<class EQUATION, class RAW_PROBLEM>
    using QFrameDiscretizedSlitDomainProblem
    = GenericDiscretizedProblem<WaveletTL::CachedQuarkletSlitDomainProblem, EQUATION, RAW_PROBLEM,
                                QFrameSolution<EQUATION>, METHOD_LIST>;

    template<class EQUATION>
    using SimpleQFrameSolution = GenericSolution<InfiniteVector<double, typename EQUATION::Index>,
                                                 SampledMapping<EQUATION::space_dimension>,
                                                 typename EQUATION::WaveletBasis>;

    template<class EQUATION, class RAW_PROBLEM>
    using QFrameDiscretizedIntervalProblem
    = GenericDiscretizedProblem<WaveletTL::CachedQuarkletProblem, EQUATION, RAW_PROBLEM,
                                SimpleQFrameSolution<EQUATION>, METHOD_LIST>;

    template<class EQUATION, class RAW_PROBLEM>
    using QFrameDiscretizedCubeProblem
    = GenericDiscretizedProblem<WaveletTL::CachedQuarkletTProblem, EQUATION, RAW_PROBLEM,
                                SimpleQFrameSolution<EQUATION>, METHOD_LIST>;
};



// To avoid memory leaks:
template <class FRAMECLASS>
class QFrameWithOwnership : public FRAMECLASS
{
    typedef typename FRAMECLASS::IntervalFrame IFrame;
public:
    QFrameWithOwnership(IFrame* frame1d, IFrame* frame1d_11,
                        IFrame* frame1d_01, IFrame* frame1d_10)
        : FRAMECLASS(frame1d, frame1d_11, frame1d_01, frame1d_10),
          frame1d_(frame1d),
          frame1d_11_(frame1d_11),
          frame1d_01_(frame1d_01),
          frame1d_10_(frame1d_10) { }

private:
    std::unique_ptr<IFrame> frame1d_, frame1d_11_, frame1d_01_, frame1d_10_;
};



template <class FRAMECLASS>
FRAMECLASS* createQFrame(int jmax, int pmax)
{
    typedef typename FRAMECLASS::IntervalFrame IFrame;

    IFrame* frame1d = new IFrame(false,false);
    frame1d->set_jpmax(jmax-frame1d->j0(),pmax);
    IFrame* frame1d_11 = new IFrame(true,true);
    frame1d_11->set_jpmax(jmax-frame1d->j0(),pmax);
    IFrame* frame1d_01 = new IFrame(false,true);
    frame1d_01->set_jpmax(jmax-frame1d->j0(),pmax);
    IFrame* frame1d_10 = new IFrame(true,false);
    frame1d_10->set_jpmax(jmax-frame1d->j0(),pmax);

    FRAMECLASS* frame = new QFrameWithOwnership<FRAMECLASS>(frame1d, frame1d_11, frame1d_01, frame1d_10);
    frame->set_jpmax(jmax,pmax);

    return frame;
}


#endif // GENERIC_QFRAME_DISCRETIZATION_MODULE_H

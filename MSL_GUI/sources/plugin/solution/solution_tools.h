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


#ifndef SOLUTION_TOOLS_H
#define SOLUTION_TOOLS_H

#include "WaveletTL/interval/p_evaluate.h"
#include "WaveletTL/interval/ds_evaluate.h"
#include "WaveletTL/interval/pq_evaluate.h"
#include "WaveletTL/cube/cube_evaluate.h"

#include "WaveletTL/Ldomain/ldomain_frame_evaluate.h"
#include "WaveletTL/recring/recring_frame_evaluate.h"
#include "WaveletTL/slitdomain/slitdomain_frame_evaluate.h"
#include "WaveletTL/cube/tframe_evaluate.h"

#include "FrameTL/frame_evaluate.h"

#include "FrameTL/aggregated_frame.h"
#include "WaveletTL/interval/i_indexplot.h"
#include "WaveletTL/interval/i_q_indexplot.h"
#include "WaveletTL/cube/cube_indexplot.h"
#include "WaveletTL/Ldomain/ldomain_frame_indexplot.h"
#include "WaveletTL/recring/recring_frame_indexplot.h"
#include "WaveletTL/slitdomain/slitdomain_frame_indexplot.h"
#include "WaveletTL/cube/tframe_indexplot.h"

#include "misc/cached_frame_problem_helper.h"


namespace solution_tools
{

using MathTL::Array1D;
using MathTL::InfiniteVector;
using MathTL::InfiniteDiagonalMatrix;
using FrameTL::AggregatedFrame;
using WaveletTL::CubeBasis;
using WaveletTL::LDomainFrame;
using WaveletTL::RecRingFrame;


template <class WBASIS, class COEFF_VECTOR, class PLOT_DATA>
void evaluate(const WBASIS& basis, const COEFF_VECTOR& scaledCoeffs, int resolution, PLOT_DATA& plotData)
{
    plotData = WaveletTL::evaluate(basis, scaledCoeffs, true, resolution);
}



template <class IBASIS, unsigned int DIM, class INDEX>
void evaluate(const AggregatedFrame<IBASIS,DIM>& frame, const Array1D< InfiniteVector<double, INDEX> >& scaledCoeffs, int resolution, Array1D< SampledMapping<DIM> >& plotData)
{
    FrameTL::EvaluateFrame<IBASIS, DIM, DIM> evalObj;
    plotData = evalObj.evaluate(frame, scaledCoeffs[frame.n_p()], true, resolution);
}



template <class INDEX, class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
void adjustCoeffVector(Array1D<InfiniteVector<double, INDEX> >& coeffs,
                       const AggregatedFrame<IBASIS, DIM_d, DIM_m>& frame)
{
    coeffs.resize(frame.n_p() + 1);
}



template <class INDEX, class WBASIS>
void adjustCoeffVector(const InfiniteVector<double, INDEX>& coeffs, const WBASIS& basis)
{
    // just a dummy (since no adjustment is required)
    (void) coeffs;
    (void) basis;
}



template <class INDEX>
void scaleCoefficients(InfiniteVector<double, INDEX>& coeffs,
                       const InfiniteDiagonalMatrix<double, INDEX>& D, const int k)
{
    coeffs.scale(&D, k);
}



template <class INDEX>
void scaleCoefficients(Array1D<InfiniteVector<double, INDEX> >& coeffs,
                       const InfiniteDiagonalMatrix<double, INDEX>& D, const int k)
{
    for (int i = 0; i < coeffs.size(); i++)
        coeffs[i].scale(&D, k);
}



template <class INDEX, class PROBLEM>
void scaleCoefficients(Array1D<InfiniteVector<double, INDEX> >& coeffs,
                       CachedFrameProblemHelper<PROBLEM>& problem, const int k)
{
    if (problem.lastUsedCachedProblemLocal())
    {
        scaleCoefficients(coeffs, problem.cachedProblemLocal(), k);
    }
    else
    {
        scaleCoefficients(coeffs, problem.cachedProblem(), k);
    }
}



// For IntervalIndex-bases:
template <class IBASIS>
void plotIndices(const IBASIS* basis,
                 const InfiniteVector<double, WaveletTL::IntervalIndex<IBASIS> >& coeffs,
                 const int jmax, std::ostream& os)
{
    WaveletTL::plot_indices(basis, coeffs, jmax, os, "jet", true, true);
}


// For IntervalIndex2-bases:
template <class IBASIS>
void plotIndices(const IBASIS* basis,
                 const InfiniteVector<double, WaveletTL::IntervalIndex2<IBASIS> >& coeffs,
                 const int jmax, std::ostream& os)
{
    WaveletTL::plot_indices2(basis, coeffs, jmax, os, "jet", true, true);
}


// For IntervalQIndex-bases:
template <class IBASIS>
void plotIndices(const IBASIS* basis,
                 const InfiniteVector<double, WaveletTL::IntervalQIndex<IBASIS> >& coeffs,
                 const int jmax, std::ostream& os)
{
    int pmax = basis->get_pmax_();

    int rows = int (ceil(sqrt(pmax+1)));
    int columns = rows - 1;
    if (rows*columns < pmax + 1) columns++;

    for (int p = 0; p <= pmax; p++)
    {
        os << "subplot(" << rows << "," << columns << "," << p+1 << ")" << std::endl;
        WaveletTL::plot_indices_iq(basis, coeffs, jmax, os, p, "jet", true, false);     // from i_q_indexplot.h
        os << "title('coeffs on level p=" << p <<  "')\n" << std::endl;
    }
}

// For IntervalQIndex2-bases:
template <class IBASIS>
void plotIndices(const IBASIS* basis,
                 const InfiniteVector<double, WaveletTL::IntervalQIndex2<IBASIS> >& coeffs,
                 const int jmax, std::ostream& os)
{
    int pmax = basis->get_pmax_();

    int rows = int (ceil(sqrt(pmax+1)));
    int columns = rows - 1;
    if (rows*columns < pmax + 1) columns++;

    for (int p = 0; p <= pmax; p++)
    {
        os << "subplot(" << rows << "," << columns << "," << p+1 << ")" << std::endl;
        WaveletTL::plot_indices_iq2(basis, coeffs, jmax, os, p, "jet", true, false);     // from i_q_indexplot.h
        os << "title('coeffs on level p=" << p <<  "')\n" << std::endl;
    }
}


// For cube bases:
template<class IBASIS, unsigned int DIM>
void plotIndices(const CubeBasis<IBASIS, DIM>* basis,
                 const InfiniteVector<double, typename CubeBasis<IBASIS, DIM>::Index>& coeffs,
                 const int jmax, std::ostream& os)
{
    WaveletTL::plot_indices_cube(basis, coeffs, jmax, os, "jet", true, true);
}




// For LDomain-frame:
template<class IFRAME>
void plotIndices(const LDomainFrame<IFRAME>* frame,
                 const InfiniteVector<double, typename LDomainFrame<IFRAME>::Index>& coeffs,
                 const int jmax, std::ostream& os)
{
    MultiIndex<int, 2> jstart;// start=basis1.j0();
    MultiIndex<int, 2> estart;
    for (typename InfiniteVector<double, typename LDomainFrame<IFRAME>::Index>::const_iterator it(coeffs.begin()); it != coeffs.end(); ++it)
    {
        typename LDomainFrame<IFRAME>::Index lambda = it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart))
        {
            jstart=lambda.j();
            estart=lambda.e();
            WaveletTL::plot_indices_ldomain(frame, coeffs, os, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            os <<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
}



// For RecRing-frame:
template<class IFRAME>
void plotIndices(const RecRingFrame<IFRAME>* frame,
                 const InfiniteVector<double, typename RecRingFrame<IFRAME>::Index>& coeffs,
                 const int jmax, std::ostream& os)
{
    MultiIndex<int, 2> jstart;// start=basis1.j0();
    MultiIndex<int, 2> estart;
    for (typename InfiniteVector<double, typename RecRingFrame<IFRAME>::Index>::const_iterator it(coeffs.begin()); it != coeffs.end(); ++it)
    {
        typename RecRingFrame<IFRAME>::Index lambda = it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart))
        {
            jstart=lambda.j();
            estart=lambda.e();
            WaveletTL::plot_indices_recring(frame, coeffs, os, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            os <<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
}



// For SlitDomain-frame:
template<class IFRAME>
void plotIndices(const SlitDomainFrame<IFRAME>* frame,
                 const InfiniteVector<double, typename SlitDomainFrame<IFRAME>::Index>& coeffs,
                 const int jmax, std::ostream& os)
{
    MultiIndex<int, 2> jstart;// start=basis1.j0();
    MultiIndex<int, 2> estart;
    for (typename InfiniteVector<double, typename SlitDomainFrame<IFRAME>::Index>::const_iterator it(coeffs.begin()); it != coeffs.end(); ++it)
    {
        typename SlitDomainFrame<IFRAME>::Index lambda = it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart))
        {
            jstart=lambda.j();
            estart=lambda.e();
            WaveletTL::plot_indices_slitdomain(frame, coeffs, os, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            os <<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
}



// For tensor frames:
template<class IFRAME, unsigned int DIM>
void plotIndices(const TensorFrame<IFRAME, DIM>* frame,
                 const InfiniteVector<double, typename TensorFrame<IFRAME, DIM>::Index>& coeffs,
                 const int jmax, std::ostream& os)
{
    WaveletTL::plot_indices_tframe(frame, coeffs, jmax, os);
}



// For AggregatedFrame<IBASIS,1>:
template <class IBASIS, class INDEX>
void plotIndices(const AggregatedFrame<IBASIS,1>* frame,
                 const Array1D< InfiniteVector<double, INDEX> >& coeffs,
                 const int jmax, std::ostream& os)
{
    typedef typename IBASIS::Index BasisIndex;
    int n_p = frame->n_p();

    std::vector< InfiniteVector<double, BasisIndex> > indices(n_p);

    typename InfiniteVector<double, INDEX>::const_iterator it = coeffs[n_p].begin();
    for (; it!= coeffs[n_p].end(); ++it)
    {
        INDEX ind(it.index());
        //cout << "level = " << ind.j() << endl;
        indices[ind.p()].set_coefficient(BasisIndex(ind.j(), ind.e()[0], ind.k()[0], frame->bases()[0]->bases()[0]), *it);
    }

    int rows = int (ceil(sqrt(n_p)));
    int columns = rows - 1;
    if (rows*columns < n_p) columns++;

    for (int p = 0; p < n_p; p++)
    {
        os << "subplot(" << rows << "," << columns << "," << p+1 << ")" << std::endl;
        plotIndices(frame->bases()[p]->bases()[0], indices[p], jmax, os);
        os << "title('Coefficients on patch " << p+1 <<  "')\n" << std::endl;
    }
}




// For AggregatedFrame<IBASIS,2>:
template <class IBASIS, class INDEX>
void plotIndices(const AggregatedFrame<IBASIS,2>* frame,
                 const Array1D< InfiniteVector<double, INDEX> >& coeffs,
                 const int jmax, std::ostream& os)
{
    typedef CubeIndex<IBASIS,2,MappedCubeBasis<IBASIS,2,2> > CIndex;

    Array1D< InfiniteVector<double, CIndex> > approximations_cube(frame->n_p());

    // convert indices to CubeIndices
    for (int i = 0; i < frame->n_p(); i++)
    {
        MappedCubeBasis<IBASIS,2,2>* mapped_basis = frame->bases()[i];

        for (typename InfiniteVector<double, INDEX>::const_iterator it = coeffs[frame->n_p()].begin(),
             itend = coeffs[frame->n_p()].end(); it != itend; ++it)
        {
            if (it.index().p() == i)
            {
                approximations_cube[i].set_coefficient(CIndex(it.index().j(),it.index().e(),it.index().k(), mapped_basis),*it);
            }
        }

        if (i > 0)
        {
            os << "\nfigure\n" << std::endl;
        }
        plot_indices_cube(frame->bases()[i], approximations_cube[i], jmax, os, "jet", false, true);
        os << "str = 'Coefficients on patch " << i+1 << "';" << std::endl;
        os << "annotation('textbox', [0 0.9 1 0.1], 'String', str, 'FontWeight', 'bold', 'FontSize', 12, 'EdgeColor', 'none', 'HorizontalAlignment', 'center')" << std::endl;
    }
}

}

#endif // SOLUTION_TOOLS_H

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


#ifndef GENERIC_AGGFRAME_DISCRETIZATION_MODULE_H
#define GENERIC_AGGFRAME_DISCRETIZATION_MODULE_H


#include "discretization_module_base.h"
#include "misc/cached_frame_problem_helper.h"
#include "FrameTL/aggregated_frame.h"

#include "misc/interval_bases_list.h"
#include "methods/aggframe/list.h"


template<unsigned int jmax_standard, class LIST_1D_BASES = IntervalBasesList_Primbs, class METHOD_LIST = methods::aggframe::List>
class GenericAggFrameDiscretizationModule
        : public DiscretizationModuleBase<jmax_standard, LIST_1D_BASES, METHOD_LIST>
{

public:

    static QString getDiscretizationTypeName()
    {
        return QStringLiteral("Aggregated Wavelet Frame");
    }


protected:

    template<class EQUATION>
    using AggFrameSolution = GenericSolution<Array1D<InfiniteVector<double, typename EQUATION::Index> >,
                                             Array1D<SampledMapping<EQUATION::space_dimension> >,
                                             typename EQUATION::Frame>;

    template<class EQUATION, class RAW_PROBLEM>
    using AggFrameDiscretizedProblem = GenericDiscretizedProblem<CachedFrameProblemHelper, EQUATION,
                                                                 RAW_PROBLEM, AggFrameSolution<EQUATION>,
                                                                 METHOD_LIST>;
};



// To avoid memory leaks:
template <class IBASIS, unsigned int DIM>
class AggFrameWithOwnership : public FrameTL::AggregatedFrame<IBASIS, DIM>
{
public:
    AggFrameWithOwnership(Atlas<DIM>* atlas,
                          const Array1D<FixedArray1D<int,2*DIM> >& bc,
                          const int jmax)
        : FrameTL::AggregatedFrame<IBASIS,DIM>(atlas, bc, jmax),
          atlas_(atlas),
          chart1_((atlas->charts())[0]),
          chart2_((atlas->charts())[1]) { }

private:
    std::unique_ptr< MathTL::Atlas<DIM> > atlas_;
    std::unique_ptr< MathTL::Chart<DIM> > chart1_, chart2_;
};



template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 1>* createAggFrame_on_Interval(int jmax, double overlap, int bcOrder);

template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 2>* createAggFrame_on_LDomain(int jmax, double overlap, int bcOrder);

template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 2>* createAggFrame_on_RecRing(int jmax, int bcOrder);



/*##################################################################################################
    Implementation
##################################################################################################*/


template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 1>* createAggFrame_on_Interval(int jmax, double overlap, int bcOrder)
{
    const int DIM = 1;
    const int d  = IBASIS::primal_polynomial_degree();

    const double patch_width = 0.5 + overlap/2.0;

    Matrix<double> A(DIM,DIM);
    A(0,0) = patch_width;
    Point<1> b;
    b[0] = 0.0;
    AffineLinearMapping<1>* chart1 = new AffineLinearMapping<1>(A,b);

    Matrix<double> A2(DIM,DIM);
    A2(0,0) = patch_width;
    Point<1> b2;
    b2[0] = 1.0-patch_width;
    AffineLinearMapping<1>* chart2 = new AffineLinearMapping<1>(A2,b2);

    Array1D<Chart<DIM,DIM>* > charts(2);
    charts[0] = chart1;
    charts[1] = chart2;

    SymmetricMatrix<bool> adj(2);
    adj(0,0) = 1;
    adj(1,1) = 1;
    adj(1,0) = 1;
    adj(0,1) = 1;

    //to specify primal boundary the conditions
    Array1D<FixedArray1D<int,2*DIM> > bc(2);

    if (bcOrder == 1)
    {
        //primal boundary conditions for first patch: all Dirichlet
        FixedArray1D<int,2*DIM> bound_1;
        bound_1[0] = 1;
        bound_1[1] = d-1;

        bc[0] = bound_1;

        //primal boundary conditions for second patch: all Dirichlet
        FixedArray1D<int,2*DIM> bound_2;
        bound_2[0] = d-1;
        bound_2[1] = 1;

        bc[1] = bound_2;
    }
    else if(bcOrder == 2)
    {
        //primal boundary conditions for first patch: all Dirichlet
        FixedArray1D<int,2*DIM> bound_1;
        bound_1[0] = 2;
        bound_1[1] = 2;

        bc[0] = bound_1;

        //primal boundary conditions for second patch: all Dirichlet
        FixedArray1D<int,2*DIM> bound_2;
        bound_2[0] = 2;
        bound_2[1] = 2;

        bc[1] = bound_2;
    }

    Atlas<1>* atlas = new Atlas<1>(charts, adj);
    return new AggFrameWithOwnership<IBASIS,1>(atlas, bc, jmax);
}



template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 2>* createAggFrame_on_LDomain(int jmax, double overlap, int bcOrder)
{
    const int DIM = 2;

    //##############################
    Matrix<double> A(DIM,DIM);
    // A(0,0) = 2.0;
    A(0,0) = 1.0 + overlap;
    A(1,1) = 1.0;
    Point<2> b;
    //b[0] = -1.;
    b[0] = -overlap;
    b[1] = -1.0;
    AffineLinearMapping<2>* chart1 = new AffineLinearMapping<2>(A,b);

    Matrix<double> A2(DIM,DIM);
    A2(0,0) = 1.0;
    A2(1,1) = 2.0;
    Point<2> b2;
    b2[0] = -1.0;
    b2[1] = -1.0;
    AffineLinearMapping<2>* chart2 = new AffineLinearMapping<2>(A2,b2);


    //##############################
    Array1D<Chart<2>* > charts(2);
    charts[0] = chart1;
    charts[1] = chart2;

    SymmetricMatrix<bool> adj(2);
    adj(0,0) = 1;
    adj(1,1) = 1;
    adj(1,0) = 1;
    adj(0,1) = 1;

    //to specify the primal boundary conditions
    Array1D<FixedArray1D<int,2*DIM> > bc(2);

    //primal boundary conditions for first patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_1;
    bound_1[0] = bcOrder;
    bound_1[1] = bcOrder;
    bound_1[2] = bcOrder;
    bound_1[3] = bcOrder; // d-1

    bc[0] = bound_1;

    //primal boundary conditions for second patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_2;
    bound_2[0] = bcOrder;
    bound_2[1] = bcOrder; // d-1
    bound_2[2] = bcOrder;
    bound_2[3] = bcOrder;

    bc[1] = bound_2;

    Atlas<2>* atlas = new Atlas<2>(charts, adj);
    return new AggFrameWithOwnership<IBASIS,2>(atlas, bc, jmax);
}



template <class IBASIS>
FrameTL::AggregatedFrame<IBASIS, 2>* createAggFrame_on_RecRing(int jmax, int bcOrder)
{
    const int DIM = 2;

    //##############################
    Matrix<double> A(DIM,DIM);
    A(0,0) = 3.0;
    A(1,1) = 1.0;
    Point<2> b;
    b[0] = -1.0;
    b[1] = -1.0;
    AffineLinearMapping<2>* chart1 = new AffineLinearMapping<2>(A,b);

    Matrix<double> A2(DIM,DIM);
    A2(0,0) = 1.0;
    A2(1,1) = 3.0;
    Point<2> b2;
    b2[0] =  1.0;
    b2[1] = -1.0;
    AffineLinearMapping<2>* chart2 = new AffineLinearMapping<2>(A2,b2);

    Matrix<double> A3(DIM,DIM);
    A3(0,0) = 3.0;
    A3(1,1) = 1.0;
    Point<2> b3;
    b3[0] = -1.0;
    b3[1] =  1.0;
    AffineLinearMapping<2>* chart3 = new AffineLinearMapping<2>(A3,b3);

    Matrix<double> A4(DIM,DIM);
    A4(0,0) = 1.0;
    A4(1,1) = 3.0;
    Point<2> b4;
    b4[0] = -1.0;
    b4[1] = -1.0;
    AffineLinearMapping<2>* chart4 = new AffineLinearMapping<2>(A4,b4);

    //##############################
    Array1D<Chart<2>* > charts(4);
    charts[0] = chart1;
    charts[1] = chart2;
    charts[2] = chart3;
    charts[3] = chart4;

    SymmetricMatrix<bool> adj(4);

    // patch 0
    adj(0,0) = 1;
    adj(0,1) = 1;
    adj(0,3) = 1;

    // patch 1
    adj(1,0) = 1;
    adj(1,1) = 1;
    adj(1,2) = 1;

    // patch 2
    adj(2,1) = 1;
    adj(2,2) = 1;
    adj(2,3) = 1;

    // patch 3
    adj(3,0) = 1;
    adj(3,2) = 1;
    adj(3,3) = 1;


    //to specify primal boundary the conditions
    Array1D<FixedArray1D<int,2*DIM> > bc(4);

    //primal boundary conditions for first patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_1;
    bound_1[0] = bcOrder;//2
    bound_1[1] = bcOrder;
    bound_1[2] = bcOrder;
    bound_1[3] = bcOrder;//2;

    bc[0] = bound_1;

    //primal boundary conditions for second patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_2;
    bound_2[0] = bcOrder;
    bound_2[1] = bcOrder;//2;
    bound_2[2] = bcOrder;
    bound_2[3] = bcOrder;

    bc[1] = bound_2;

    //primal boundary conditions for second patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_3;
    bound_3[0] = bcOrder;
    bound_3[1] = bcOrder;//2;
    bound_3[2] = bcOrder;
    bound_3[3] = bcOrder;

    bc[2] = bound_3;

    //primal boundary conditions for second patch: all Dirichlet
    FixedArray1D<int,2*DIM> bound_4;
    bound_4[0] = bcOrder;
    bound_4[1] = bcOrder;//2;
    bound_4[2] = bcOrder;
    bound_4[3] = bcOrder;

    bc[3] = bound_4;

    Atlas<2>* atlas = new Atlas<2>(charts, adj);
    return new AggFrameWithOwnership<IBASIS,2>(atlas, bc, jmax);
}


#endif // GENERIC_AGGFRAME_DISCRETIZATION_MODULE_H

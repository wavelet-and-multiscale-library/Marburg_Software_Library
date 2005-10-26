#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <numerics/corner_singularity.h>
#include "elliptic_equation.h"
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <frame_evaluate.h>


using std::cout;
using std::endl;

using FrameTL::EllipticEquation;
using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::CornerSingularityRHS;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;

using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;
using namespace FLAT;

int main()
{
  
  cout << "Testing class EllipticEquation..." << endl;
  
  const int DIM = 2;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.0;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = .0;
  b[1] = .0;
  AffineLinearMapping<2> affineP(A,b);
  //##############################

  //##############################
  LinearBezierMapping bezierP(Point<2>(-1.,-1),Point<2>(-1.,1),
			      Point<2>(0.,-1.), Point<2>(0,1.));
  //##############################
  Array1D<Chart<DIM,DIM>* > charts(1);
  charts[0] = &affineP;
  //charts[1] = &affineP;
  //charts[1] = &bezierP;
  

  SymmetricMatrix<bool> adj(1);
  adj(0,0) = 1;
  //adj(1,1) = 1;
  //adj(1,0) = 1;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(1);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;
  bound_2[2] = 1;
  bound_2[3] = 1;

  //bc[1] = bound_2;

//to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;
 
  //bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  typedef DSBasis<2,2> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;

  //finally a frame can be constructed
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT);

  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;

  CornerSingularityRHS singRhs(origin, 0.5 * M_PI, 1.5);
  
  //PoissonBVP<DIM> poisson(&singRhs);
  PoissonBVP<DIM> poisson(&const_fun);

  EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame);

  FrameIndex<Basis1D,2,2> ind = FrameTL::first_wavelet<Basis1D,2,2,Frame2D>(&frame,3);

   double max_val = 0.0;
   double tmp = 0.0;
   int c = 0;
   int d = 0;

   cout.precision(12);
#if 0
   SparseMatrix<double> stiff(225,225);
   int i = 0, j;
   for (FrameIndex<Basis1D,2,2> ind1 = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, 3);
	ind1 <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, 3); ++ind1)
     {
       j=0;
       cout << "i = " << i << endl;
       for (FrameIndex<Basis1D,2,2> ind2 = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, 3);
	    ind2 <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, 3); ++ind2)
	 {
	   //cout << i << " " << j << endl;
	   double scale_fac = 1.0 / (discrete_poisson.D(ind1)*discrete_poisson.D(ind2));
	   stiff.set_entry(i,j,discrete_poisson.a(ind1,ind2,2));
	   //int level = 
	   double dd = stiff.get_entry(i,j);
 	   if (dd != 0.0) {
 	     cout << ind1 << " " << ind2 << endl;
 	     cout << "a results in: " << dd << endl;
 	   }
	   j++;
	 }
       i++;
     }
#endif

   InfiniteVector<double,Frame2D::Index> rhs;
   Vector<double> v(225);

#if 0
   int k = 0;
   for (FrameIndex<Basis1D,2,2> ind = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, 3);
	ind <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, 3); ++ind) 
     {
       
       if (ind.p() != 0)
	 continue;
       double g = discrete_poisson.f(ind);
       rhs.set_coefficient(ind, g);
       v[k] = g;
       cout << "k=" << k << endl;
       k++;
       if (c == 0)
	 max_val = discrete_poisson.f(ind);
       else {
	 if (tmp > max_val)
	   max_val = tmp;
       }
       c++;
       //       if (fabs(1.0/(1<<ind.j())*tmp) > 1.0e-14) {
       // 	cout << "evaluation of discrete right hand side at  " 
       // 	     << ind << " results in: "
       // 	     << ldexp(1.0,-ind.j()) * tmp << endl;
       // 	d++;
       //       }
       
     }
   cout << "maximum = " << max_val << endl;
#endif
   cout << rhs << endl;

   discrete_poisson.rescale(rhs,-1);
   cout << "++++++++++++++++++++++" << endl;
   cout << rhs << endl;

   InfiniteVector<double,Frame2D::Index> xk;
   Vector<double> xk_(225);
   unsigned int iter = 0;


//    MathTL::Richardson<Vector<double>, SparseMatrix<double> >
//      (stiff, v, xk_, 0.5, 0.001, 300, iter);

//    cout << xk_ << endl;

//   WaveletTL::evaluate<Basis1D,2>(*(frame.bases()[0]),
// 			       const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
// 			       const bool primal,
// 			       const int resolution);
 



   FrameIndex<Basis1D,2,2> index = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, 3); 
   SampledMapping<2> expansion_out = FrameTL::evaluate(frame,index,1,5);
   

   std::ofstream ofs("wav_out.m");
   expansion_out.matlab_output(ofs);
   ofs.close();

   //Array1D<SampledMapping<2> > expansion_out =  FrameTL::evaluate(frame,rhs,1,5);


  return 0;
}

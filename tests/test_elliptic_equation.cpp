
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <numerics/corner_singularity.h>
#include "elliptic_equation.h"


using std::cout;
using std::endl;

using FrameTL::EllipticEquation;
using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::CornerSingularityRHS;

using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  
  cout << "Testing class EllipticEquation..." << endl;
  
  const int DIM = 2;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 1.0;
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
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP;
  //charts[1] = &bezierP;
  

  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

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

  bc[1] = bound_2;

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
 
  bcT[1] = bound_4;

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

  EllipticEquation<Basis1D,DIM> discrete_poisson(poisson, &frame);

  FrameIndex<Basis1D,2,2> ind = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame,4);

  double max_val = 0.0;
  double tmp = 0.0;
  int c = 0;
  int d = 0;

  for (FrameIndex<Basis1D,2,2> ind = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, 3);
       ind <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, 4); ++ind) 
    {
      
      if (ind.p() != 0)
	continue;

      tmp = discrete_poisson.f(ind);

      if (c == 0)
	max_val = discrete_poisson.f(ind);
      else {
	if (tmp > max_val)
	  max_val = tmp;
      }
      c++;
      if (fabs(1.0/(1<<ind.j())*tmp) > 1.0e-14) {
	cout << "evaluation of discrete right hand side at  " 
	     << ind << " results in: "
	     << ldexp(1.0,-ind.j()) * tmp << endl;
	d++;
      }
      
    }
  cout << "d= " << d << endl;
  cout << "maximum = " << max_val << endl;

  return 0;
}

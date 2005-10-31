#include <map>
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
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin_frame_utils.h>
#include <frame_support.h>


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
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  
  cout << "Testing class EllipticEquation..." << endl;
  
  const int DIM = 2;

  typedef DSBasis<2,2> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;

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
  LinearBezierMapping bezierP(Point<2>(-1.4,-1),Point<2>(-1.,1),
			      Point<2>(0.,-1.), Point<2>(0,1.));
  //##############################
  Array1D<Chart<DIM,DIM>* > charts(1);
  charts[0] = &affineP;
  //charts[1] = &affineP;
  charts[0] = &bezierP;
  

  SymmetricMatrix<bool> adj(1);
  adj(0,0) = 1;
  // adj(1,1) = 1;
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
  Array1D<FixedArray1D<int,2*DIM> > bcT(1);

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

   double tmp = 0.0;
   int c = 0;
   int d = 0;

   cout.precision(12);

   InfiniteVector<double,Frame2D::Index> rhs;
   Vector<double> v(225);

   FixedArray1D<bool,4> bcn;
   bcn[0] = bcn[1] = true;
   bcn[2] = bcn[3] = true;
   Basis basis(bcn);
 
   FrameIndex<Basis1D,2,2> index = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0()); 

   SampledMapping<2> wav_out = FrameTL::evaluate(frame,index,1,5);

   std::ofstream ofs("wav_out.m");
   wav_out.matlab_output(ofs);
   ofs.close();

   Array1D<SampledMapping<2> > expansion_out = FrameTL::evaluate(frame,rhs,1,5);

   std::ofstream ofs2("expan_out.m");
   expansion_out[0].matlab_output(ofs2);
   ofs2.close();

   std::ofstream ofs3("full_expan_out.m");
   matlab_output(ofs3,expansion_out);
   ofs3.close();
   
   //###############################################
   // testing correctness of routine for
   // computation of right hand side
   
   
   cout << "++++++++++++++++++++++" << endl;
   cout << rhs << endl;

   Array1D<SampledMapping<2> > S = FrameTL::evaluate<Basis1D,2>(frame, rhs, false, 5);// expand in dual basis
   std::ofstream ofs4("rhs_out.m");
   matlab_output(ofs4,S);
   ofs4.close();
   //###############################################   

   typedef Frame2D::Index Index;
   set<Index> Lambda;
   for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
      lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, frame.j0()); ++lambda)
     Lambda.insert(lambda);

   cout << "setting up full stiffness matrix..." << endl;
   SparseMatrix<double> stiff;
   FrameTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff);
   cout << "setting up full right hand side..." << endl;
   Vector<double> rh;
   FrameTL::setup_righthand_side(discrete_poisson, Lambda, rh);

   cout << "perfrming CG algorithm to solve projected problem..." << endl;
   Vector<double> xk(Lambda.size()), err(Lambda.size()); xk = 0;
   unsigned int iter= 0;
   CG(stiff, rh, xk, 0.0001, 200, iter);
   
   cout << "performing output..." << endl;

   InfiniteVector<double,Frame2D::Index> u;
   unsigned int i = 0;
   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, xk[i]);

   discrete_poisson.rescale(u,-1);

   Array1D<SampledMapping<2> > U = FrameTL::evaluate<Basis1D,2>(frame, u, true, 6);// expand in dual basis

   std::ofstream ofs5("approx_solution_out.m");
   matlab_output(ofs5,U);
   ofs5.close();

   //########## testing runtime type information functionality #############
   //FrameTL::intersect_supports<Basis1D,2,2>(frame, index, index);
   FrameIndex<Basis1D,2,2> index2 = FrameTL::last_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0()); 
   //discrete_poisson.a(index,index2,3);
   //#######################################################################

  return 0;
}

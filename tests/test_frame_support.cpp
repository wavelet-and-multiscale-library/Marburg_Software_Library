
#include <map>
#include <fstream>
#include <iostream>
#include <time.h> 
#include <list>
#include <set>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <numerics/corner_singularity.h>
#include <elliptic_equation.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/cached_problem.h>
#include <frame_support.h>
#include <frame_index.h>
//#include <utils/multiindex.h>

using std::cout;
using std::endl;

using FrameTL::EvaluateFrame;
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

  cout << "Testing Frame support..." << endl;
  
  const int DIM = 2;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<3,3> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

  //! The index type.
  typedef MultiIndex<int,DIM> type_type;
    
  //! The translation index type.
  typedef MultiIndex<int,DIM> translation_type;

  typedef Frame2D::Support SuppType;
  typedef Basis::Support SuppTypeB;




  // #####################################################################################
  // We set up a distorted L-shaped domain similar to the one in Figure 5.4 right in
  // Manuel's PhD thesis.
  // #####################################################################################
  const double t = 2./3.;
  const double theta0 = 0.5-t*0.25;
  const double omega = 1.5+2*t*0.25;
  
  const double alpha = tan(t*M_PI*0.25);
  cout << "alpha = " << alpha << endl;

//  LinearBezierMapping bezierP(Point<2>(-1.,-1.), Point<2>(-1.,-alpha),
//			      Point<2>(1.,-1.), Point<2>(1.,alpha));

  
//  LinearBezierMapping bezierP2(Point<2>(-1.,-1.), Point<2>(-1.,1.),
//			       Point<2>(-alpha,-1.), Point<2>(alpha,1.));

  // L- shaped:
   LinearBezierMapping bezierP(Point<2>(-1.,-1.),Point<2>(-1.,1.),
  			      Point<2>(0.,-1.), Point<2>(0.,1.));
  
  LinearBezierMapping bezierP2(Point<2>(-1.,-1.),Point<2>(-1.,0.),
 			       Point<2>(1.,-1.), Point<2>(1.,0.));

  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &bezierP;
  charts[1] = &bezierP2;
  // #####################################################################################

  // setup the adjacency relation of the patches
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
  // to specify the primal boundary conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

  // primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;//2

  bc[0] = bound_1;

  // primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;//2
  bound_2[2] = 1;
  bound_2[3] = 1;

  bc[1] = bound_2;

  // to specify the dual boundary conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  // dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  // dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;

  
  bcT[1] = bound_4;

  // create the atlas
  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  // finally, a frame can be constructed
  // AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, jmax);
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, 9);
 

  cout << "Anzahl Patches  =  " << frame.n_p() << endl;
  
  type_type e(1,1);
  translation_type k(3,3);

  Index lambda(&frame, 3, e,0, k);

  cout << "Test  =  " <<  lambda << endl;

  list<Index> intersecting;

  clock_t tstart, tend;
  double time;
  tstart = clock();

  cout << "Calculating support" << endl;


  SuppType Supp;
  SuppTypeB SuppB;


  support<Basis1D, 2, 2>(frame, lambda, Supp);
  support_on_cube<Basis1D, 2, 2>(frame, lambda, SuppB);
  
  cout << "The support of lambda is  =  (" <<  Supp.a[0] << "," << Supp.b[0] << ") X (" <<  Supp.a[1] << "," <<  Supp.b[1] << ")"<< endl;
  cout << "The support of lambda on the cube is  =  (" <<  SuppB.a[0] << "," << SuppB.b[0] << ") X (" <<  SuppB.a[1] << "," <<  SuppB.b[1] << ")"<< endl;
  
//  for(int i = 0; i<100; i++)
//  {
    intersecting_wavelets<Basis1D, 2, 2>(frame, lambda, 3, false, intersecting);
//  }

  tend = clock();

 for(list<Index>::iterator it=intersecting.begin();
	 it != intersecting.end(); ++it) {
  support<Basis1D, 2, 2>(frame, *it, Supp);
  support_on_cube<Basis1D, 2, 2>(frame, *it, SuppB);
  cout << "Anzahl Patches  =  " <<  *it << "  The support on cube is  =  (" <<  SuppB.a[0] << "," << SuppB.b[0] << ") X (" <<  SuppB.a[1] << "," <<  SuppB.b[1] << ")" << endl;
}


  time = (double)(tend-tstart)/CLOCKS_PER_SEC;

  cout << "Time needed for calculating support  =  " <<  time << endl;

  //EvaluateFrame<Basis1D,2,2> evalObj;

  return 0;
}

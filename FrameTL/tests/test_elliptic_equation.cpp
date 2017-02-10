#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

#define TWO_D

#include <map>
#include <fstream>
#include <iostream>
#include <time.h> 
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
  
  cout << "Testing class EllipticEquation..." << endl;
  
  const int DIM = 2;
  const int jmax = 3;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<3,3> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

  EvaluateFrame<Basis1D,2,2> evalObj;

  // #####################################################################################
  // We set up a distorted L-shaped domain similar to the one in Figure 5.4 right in
  // Manuel's PhD thesis.
  // #####################################################################################
  const double t = 2./3.;
  const double theta0 = 0.5-t*0.25;
  const double omega = 1.5+2*t*0.25;
  
  const double alpha = tan(t*M_PI*0.25);
  cout << "alpha = " << alpha << endl;

  LinearBezierMapping bezierP(Point<2>(-1.,-1.), Point<2>(-1.,-alpha),
			      Point<2>(1.,-1.), Point<2>(1.,alpha));

  
  LinearBezierMapping bezierP2(Point<2>(-1.,-1.), Point<2>(-1.,1.),
			       Point<2>(-alpha,-1.), Point<2>(alpha,1.));

  // L- shaped:
//   LinearBezierMapping bezierP(Point<2>(-1.,-1.),Point<2>(-1.,1.),
//  			      Point<2>(0.,-1.), Point<2>(0.,1.));
  
//   LinearBezierMapping bezierP2(Point<2>(-1.,-1.),Point<2>(-1.,0.),
// 			       Point<2>(1.,-1.), Point<2>(1.,0.));

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
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);

  // #####################################################################################
  // Setting up a constant function which can serve as the right-hand side.
  // #####################################################################################
  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;
  // #####################################################################################

  CornerSingularityRHS singRhs(origin, theta0, omega);
  CornerSingularity sing2D(origin, theta0, omega);

  // // L-shaped:
  //  CornerSingularityRHS singRhs(origin, 0.5, 1.5);
  // CornerSingularity sing2D(origin, 0.5, 1.5);
  
  PoissonBVP<DIM> poisson(&singRhs);
  //PoissonBVP<DIM> poisson(&const_fun);

  EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0048, 1.0/0.01);

  // #####################################################################################
  // Plotting exact solution and right-hand side.
  // #####################################################################################
  //double tmp = 0.0;
  //int c = 0;
  //int d = 0;
  
  cout.precision(12);
  
  InfiniteVector<double,Frame2D::Index> rhs;
  Vector<double> v(225);
  
  FixedArray1D<bool,4> bcn;
  bcn[0] = bcn[1] = true;
  bcn[2] = bcn[3] = true;
  Basis basis(bcn);
  

  Point<2> x1;
  x1[0] = -1;
  x1[1] = -1;
  Point<2> x2;
  x2[0] = 1;
  x2[1] = 1;

  Grid<2> grid(x1,x2,128);
  
  SampledMapping<2> rhsOut (grid, singRhs);
  std::ofstream ofs_rhs("sing_rhs_out.m");
  rhsOut.matlab_output(ofs_rhs);
  ofs_rhs.close();
  
  SampledMapping<2> singout (grid, sing2D);
  std::ofstream ofs_sing("sing_out.m");
  singout.matlab_output(ofs_sing);
  ofs_sing.close();
  // #####################################################################################
  
#if 1

  // #####################################################################################
  // We set up a full stiffness matrix and right-hand side up to a maximal level
  // and try to solve the projected problem.
  // #####################################################################################

  // setup Galerkin index set
  set<Index> Lambda;
  for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
    //cout << lambda << endl;
  }
  
  // set up full right-hand side
  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(discrete_poisson, Lambda, rh);
  //cout << rh << endl;

  // set up full stiffness matrix
  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;

  clock_t tstart, tend;
  double time;
  tstart = clock();

  WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  
  // write the stiffness matrix into matab readable output
  stiff.matlab_output("stiff_2D_out", "stiff",1);  
  
  
  //unsigned int iter= 0;
  Vector<double> x(Lambda.size()); x = 1;
  //double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  double lmax = 1;
  cout << "lmax = " << lmax << endl;

  x = 1;
  //double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  
  cout << "performing iterative scheme to solve projected problem..." << endl;
  Vector<double> xk(Lambda.size()), err(Lambda.size()); xk = 0;


  //CG(stiff, rh, xk, 1.0e-6, 100, iter);
  //cout << "CG iterations needed: "  << iter << endl;
  //Richardson(stiff, rh, xk, 2. / lmax - 0.01, 0.0001, 2000, iter);
  double alpha_n = 2. / lmax - 0.001;
  
  // simple non-adaptive steepest descent iteration
  Vector<double> resid(xk.size());
  Vector<double> help(xk.size());
  for (int i = 0; i < 500;i++) {
    stiff.apply(xk,help);
    resid = rh - help;
    cout << "i = " << i << " " << sqrt(resid*resid) << endl;
    stiff.apply(resid,help);
    alpha_n = (resid * resid) * (1.0 / (resid * help));
    resid *= alpha_n;
    xk = xk + resid;
  }

//   // check parts of the matrix for symmetry
//   for (int i = 0; i < 450; i++) 
//     for (int j = 0; j < 450; j++) {
//       if (! (fabs(stiff.get_entry(i,j) -  stiff.get_entry(j,i)) < 1.0e-13)) {
// 	cout << stiff.get_entry(i,j) << endl;
// 	cout << stiff.get_entry(j,i) << endl;
// 	cout << "i = " << i << " j = " << j << endl;
// 	cout << "#######################" << endl;
// 	abort();
//       }
//     }


  // #####################################################################################

  // #####################################################################################
  // We sample the generated approximate solution.
  // #####################################################################################
  cout << "performing output..." << endl;
  
  InfiniteVector<double,Frame2D::Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, xk[i]);
  
  // we have to apply the inverse diagonal scaling matrix D^{-1}
  u.scale(&discrete_poisson,-1);
  
  // sample the frame expansion given by the coefficients u
  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u, true, 6);//expand in primal basis
  
  // write matlab m-file
  std::ofstream ofs5("approx_solution_out.m");
  matlab_output(ofs5,U);
  ofs5.close();
  // #####################################################################################
 

  // #####################################################################################
  // The rest of the code in this file can be neglected.
  // #####################################################################################

#if 0
  char filename1[50];
  u.clear();
  Lambda.clear();
  int l = 0;
  for (Index lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, frame.j0()); ++lambda) {

    if (l != 193) {
      l++;
      continue;
    }

    Index mu = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
   
    
    cout << "result = " << discrete_poisson.a(lambda,mu) << endl;

    typedef WaveletTL::CubeBasis<Basis1D,DIM> CUBEBASIS;
    
    typedef CUBEBASIS::Index CubeIndex;
    
    CUBEBASIS::Support supp_lambda;
    
    WaveletTL::support<Basis1D,DIM>(*frame.bases()[lambda.p()], 
				    CubeIndex(lambda.j(),
					      lambda.e(),
					      lambda.k(),
					      frame.bases()[lambda.p()]),
				    supp_lambda);
    



    double dx = 1./(1 << supp_lambda.j);
//     cout << "l= " << l << " lambda:" << endl;
//     for (unsigned int i = 0; i < 2; i++) {
//       cout << "a_lambda = " << supp_lambda.a[i]*dx <<  " b= " << supp_lambda.b[i]*dx << " p = " << lambda.p() << endl;
//       cout << "a_mu = " << supp_mu.a[i]*dx <<  " b= " << supp_mu.b[i]*dx << " p = " << mu.p() << endl;
//     }    
    
//     if (l != 14) {
//       l++;
//       continue;
//     }

    cout << lambda << endl;
    cout << "##################" << endl;
    u.set_coefficient(lambda, 1);
    
    Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u, true, 7);//expand in primal basis
    sprintf(filename1, "%s%d%s", "wav_", l, ".m");
    
    
    std::ofstream ofs5(filename1);
    matlab_output(ofs5,U);
    ofs5.close();

    u.clear();
    ++l;

  }
#endif
  //cout << "lmin = " << lmin << endl;
  //########## testing runtime type information functionality #############
  //FrameTL::intersect_supports<Basis1D,2,2>(frame, index, index);
  //FrameIndex<Basis1D,2,2> index2 = FrameTL::last_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0()); 
  //discrete_poisson.a(index,index2,3);
  //#######################################################################
#endif
 
  cout << "  ... done, time needed: " << time << " seconds" << endl;
  
#if 0
  MultiIndex<unsigned int, 2> e1;
  e1[0] = 0;
  e1[1] = 0;
  MultiIndex<int, 2> k1;
  k1[0] = 1;
  k1[1] = 1;
  
  MultiIndex<unsigned int, 2> e2;
  e2[0] = 0;
  e2[1] = 0;
  MultiIndex<int, 2> k2;
  k2[0] = 1;
  k2[1] = 1;
   
  unsigned int p1 = 1, p2 = 0;
  int j2 = 3;
   

  FrameIndex<Basis1D,2,2> la(&frame,j2,e1,p1,k1);
  FrameIndex<Basis1D,2,2> mu(&frame,j2,e2,p2,k2);
  cout << "val  " << discrete_poisson.a(la,mu,2) << endl;
#endif

  //    std::list<Index> intersecting;
  //    FrameTL::intersecting_wavelets<Basis1D,2,2>(frame, la, 4, false, intersecting);

  //    cout << intersecting.size() << endl;

  //    for (std::list<Index>::const_iterator  it = intersecting.begin();
  // 	it != intersecting.end(); ++it) {
  //      cout << *it << endl;
  //    }


   return 0;
}

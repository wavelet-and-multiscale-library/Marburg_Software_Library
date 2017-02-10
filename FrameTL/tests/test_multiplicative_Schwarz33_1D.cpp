#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define _WAVELETTL_CDD1_VERBOSITY 0

#define OVERLAP 0.7
#define JMAX 20

#define PRIMALORDER 3
#define DUALORDER   3

//#define PRECOMP
//#define COMPUTECONSTANTS

#define SPARSE
//#define FULL
#define ONE_D

#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/spline_basis.h>
//#include <elliptic_equation.h>
#include <simple_elliptic_equation.h>
#include <biharmonic_equation.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#include <frame_support.h>
#include <frame_index.h>
//#include <steepest_descent.h>
#include <adaptive_multiplicative_Schwarz.h>
//#include <multiplicative_Schwarz.h>
//#include <additive_Schwarz.h>
//#include <steepest_descent_basis.h>
//#include <richardson_CDD2.h>
#include <galerkin/cached_problem.h>
#include <utils/plot_tools.h>
#include <interval/i_indexplot.h>

using std::cout;
using std::endl;

using FrameTL::FrameIndex;
//using FrameTL::EllipticEquation;
using FrameTL::SimpleEllipticEquation;
using FrameTL::BiharmonicEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;
using WaveletTL::CachedProblem;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D_RHS
  : public Function<1, VALUE>
{
public:
  Singularity1D_RHS() {};
  virtual ~Singularity1D_RHS() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    return  -100*exp(5*p[0])*(1-(exp(5*p[0])-1)/(exp(5.)-1))/(exp(5.)-1)+200*exp(10*p[0]) / 
      ((exp(5.)-1)*(exp(5.)-1))+100*(exp(5*p[0])-1)*exp(5*p[0])/((exp(5.)-1)*(exp(5.)-1));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};

/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D
  : public Function<1, VALUE>
{
public:
  Singularity1D() {};
  virtual ~Singularity1D() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    double res = 1.0 / (exp(5.) - 1.0);
    res = (exp(5.*p[0]) - 1.0) * res;
    return  (4.0 * res * (1.0 - res));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};


int main()
{


  cout << "testing steepest descent in 1D..." << endl;
  
  const int DIM = 1;

  const int jmax = JMAX;
  
  const int d = PRIMALORDER, dT = DUALORDER;

  //typedef SplineBasis<d,dT,P_construction> Basis1D;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;


  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = OVERLAP;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = OVERLAP;
  Point<1> b2;
  b2[0] = 1-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);

  //######### Three Patches ######  
  //   Matrix<double> A(DIM,DIM);
  //   A(0,0) = 0.8;
  //   Point<1> b;
  //   b[0] = 0.;
  //   AffineLinearMapping<1> affineP(A,b);
  
  //   Matrix<double> A2(DIM,DIM);
  //   A2(0,0) = 0.002;
  //   Point<1> b2;
  //   b2[0] = 0.5-0.001;
  //   AffineLinearMapping<1> affineP2(A2,b2);
  
  //   Matrix<double> A3(DIM,DIM);
  //   A3(0,0) = 0.8;
  //   Point<1> b3;
  //   b3[0] = 1-A3.get_entry(0,0);
  //   AffineLinearMapping<1> affineP3(A3,b3);

  //##############################  

  //   FixedArray1D<double,1> A3;
  //   A3[0] = 0.75;
  //   SimpleAffineLinearMapping<1> simlpeaffine1(A3,b);
  
  //   FixedArray1D<double,1> A4;
  //   A4[0] = 0.75;
  //   SimpleAffineLinearMapping<1> simlpeaffine2(A4,b2);

  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;
  //charts[2] = &affineP3;
  
  //charts[0] = &simlpeaffine1;
  //charts[1] = &simlpeaffine2;


  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;

  //   adj(2,2) = 1;
  //   adj(2,0) = 1;
  //   adj(2,1) = 1;
  //   adj(0,2) = 1;
  //   adj(1,2) = 1;

  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);
  
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

  // //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);
   
  //   //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
   
  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
 
  bcT[1] = bound_4;

  Atlas<DIM,DIM> interval(charts,adj);  
  cout << interval << endl;
   
  //finally a frame can be constructed
  //Frame1D frame(&Lshaped, bc, bcT, jmax);
  Frame1D frame(&interval, bc, jmax);
   
  EvaluateFrame<Basis1D,1,1> evalObj;


//   // ########## plotting local primal/dual frame elements #########
  
//   for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
//        lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, jmax); ++lambda) {
//     if (lambda.p() != 1)
//       continue;

//     cout << lambda << endl;
//     //SampledMapping<1> dual = evalObj.evaluate(frame, lambda, true, 10);
//     SampledMapping<1> dual = evalObj.evaluate(frame, lambda, true, 10);
//     char fname[50];
//     //sprintf(fname, "%s%d%s", "dual_wavelet_", lambda.number(), ".m");
//     sprintf(fname, "%s%d%s", "primal_wavelet_", lambda.number(), ".m");
//     std::ofstream of(fname);
//     dual.matlab_output(of);
//     of.close();
//   }
 
//   // #######################################################




  Vector<double> value(1);
  value[0] = 384;
  ConstantFunction<DIM> const_fun(value);

  //  Singularity1D_RHS<double> sing1D;
  //  Singularity1D<double> exactSolution1D;

  Singularity1D_RHS_2<double> sing1D;
  Singularity1D_2<double> exactSolution1D;
  //PolySolBiharmonic<double> exactSolution1D;
  
  //PoissonBVP<DIM> poisson(&const_fun);
  PoissonBVP<DIM> poisson(&sing1D);
  //BiharmonicBVP<DIM> biharmonic(&const_fun);

  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  //BiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&biharmonic, &frame, jmax, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.19, 1.0/0.146);

  CachedProblemLocal<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
  discrete_poisson.set_norm_A(1.0);
  discrete_poisson.set_Ainv(1.0);

  const double epsilon = 1.0e-5;

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  cout << "precompute the energy norm of the exact solution ..." << endl;
  Singularity1D_2_prime<double> exact1D_der;
  InfiniteVector<double, Index> full_vector_zero;
  double H1err = error_H_scale_interval<Basis1D>(1, frame, full_vector_zero, exact1D_der);
  cout << "H1err = " << H1err << endl;

  // ##########################################################################################
  // estimate extremal eigenvalues of local stiffness matrices and largest eigenvalue
  // of whole stiffness matrix
  
#ifdef COMPUTECONSTANTS
  set<Index> Lambda_0;
  set<Index> Lambda_1;
  set<Index> Lambda;
  for (FrameIndex<Basis1D,DIM,DIM> lambda = FrameTL::first_generator<Basis1D,DIM,DIM,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,DIM,DIM,Frame1D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
    if (lambda.p() == 0)
      Lambda_0.insert(lambda);
    else {
      Lambda_1.insert(lambda);
    }
  }
  
  SparseMatrix<double> stiff;
  
  // starting vector for Power and Inverse Power Iteration
  Vector<double> x(Lambda_0.size()); x = 1;
  // number of iterations in Power and Inverse Power Iteration
  unsigned int iter= 0;  

  WaveletTL::setup_stiffness_matrix(problem, Lambda_0, stiff);

  cout << "computing smallest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 0 is " << lmin << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 0 is " << lmax << endl;

  WaveletTL::setup_stiffness_matrix(problem, Lambda_1, stiff);

  x.resize(Lambda_1.size()); x = 1;
  cout << "computing smallest eigenvalue of stiffness matrix on patch 1" << endl;
  lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 1 is " << lmin << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 1" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 1 is " << lmax << endl;

  x.resize(Lambda.size()); x = 1;
  WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  cout << "computing largest eigenvalue of whole stiffness matrix" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of whole stiffness matrix is " << lmax << endl;
  abort();

  // (d,dt) = (2,2), jmax = 11, jmin = 3:
  // patch 0: \lambda_{\min}^0 = 0.0749363, \lambda_{\max}^0 = 2.50189
  // patch 1: \lambda_{\min}^1 = 0.0749363, \lambda_{\max}^1 = 2.50189
  // whole domain: \lambda_{\max} = 3.37459

  // (d,dt) = (3,3), jmax = 11, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 0.01238440, \lambda_{\max}^0 = 2.75083
  // patch 1: \lambda_{\min}^1 = 0.01238440, \lambda_{\max}^1 = 2.75083
  // whole domain: \lambda_{\max} = 4.17833

  // (d,dt) = (3,3), jmax = 11, jmin = 3:
  // patch 0: \lambda_{\min}^0 = 0.0113602, \lambda_{\max}^0 = 2.75056
  // patch 1: \lambda_{\min}^1 = 0.0113651, \lambda_{\max}^1 = 2.75094
  // whole domain: \lambda_{\max} = 4.87718

  // (d,dt) = (4,4), jmax = 11, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 0.00053711, \lambda_{\max}^0 = 3.83605
  // patch 1: \lambda_{\min}^1 = 0.00405706, \lambda_{\max}^1 = 3.33542
  // whole domain: \lambda_{\max} = 5.62803

  // (d,dt) = (5,5), jmax = 11, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 1.63179e-05 , \lambda_{\max}^0 = 4.46068
  // patch 1: \lambda_{\min}^1 = 0.00201216, \lambda_{\max}^1 = 4.74052
  // whole domain: \lambda_{\max} = 7.37132


#endif
  // ##########################################################################################


  clock_t tstart, tend;
  double time;
  tstart = clock();

  //multiplicative_Schwarz_SOLVE(problem, epsilon, u_epsilon0, u_epsilon1, u_epsilon);
  //adaptive_multiplicative_Schwarz_SOLVE(problem, epsilon, approximations);
  MultSchw(problem, epsilon, approximations);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "adaptive multiplicative Schwarz done, generating output..." << endl;

 
  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 12);//expand in primal basis
  cout << "...finished plotting global approximate solution" << endl;
  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], exactSolution1D, 12);
  cout << "...finished plotting global error" << endl;

  char filename1[128];
  sprintf(filename1, "%s%d%s%d%s", "./ms_results33/approx1D_global_d", d, "_dT", dT, ".m");
  std::ofstream ofs(filename1);
  matlab_output(ofs,U);
  //uplot_output(ofs,U);
  ofs.close();

  char filename2[128];
  sprintf(filename2, "%s%d%s%d%s", "./ms_results33/error1D_global_d", d, "_dT", dT, ".m");
  std::ofstream ofs1(filename2);
  matlab_output(ofs1,Error);
  //gnuplot_output(ofs1,Error);
  ofs1.close();


  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[128];
    sprintf(filename3, "%s%d%s%d%s%d%s", "./ms_results33/approx1D_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

    U = evalObj.evaluate(frame, approximations[i], true, 12);//expand in primal basis
    std::ofstream ofsloc(filename3);
    matlab_output(ofsloc,U);
    //gnuplot_output(ofsloc,U);
    ofsloc.close();
  }

  typedef Basis1D::Index Index1D;

  FixedArray1D<InfiniteVector<double, Index1D>, 2> indices;
 
  InfiniteVector<double, Index>::const_iterator it = approximations[frame.n_p()].begin();
  for (; it!= approximations[frame.n_p()].end(); ++it) {
    //cout << *it << endl;
    Index ind(it.index());
    //cout << "level = " << ind.j() << endl;
    indices[ind.p()].set_coefficient(Index1D(ind.j(),ind.e()[0],ind.k()[0],
					     frame.bases()[0]->bases()[0]), *it);
  }

  std::ofstream ofs7("./ms_results33/indices_patch_0.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[0]->bases()[0], indices[0], JMAX, ofs7, "jet", true, -16);

  std::ofstream ofs8("./ms_results33/indices_patch_1.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[1]->bases()[0], indices[1], JMAX, ofs8, "jet", true, -16);

  return 0;
}

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1

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
#include <multiplicative_Schwarz.h>
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

  const int jmax = 9;
  
  const int d = 3, dT = 3;

  //typedef SplineBasis<d,dT,P_construction> Basis1D;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;


  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 0.9;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.9;
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
  bound_1[1] = 2;
  
  bc[0] = bound_1;
  
  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 2;
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
  //  CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.19, 1.0/0.146);


//   // (d,dT) = (3,5)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.6548, 1.0/0.146);
//   discrete_poisson.set_norm_A(3.6548);
//   optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);

  // (d,dT) = (3,5)
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.6548, 1.0/0.146);
  discrete_poisson.set_norm_A(3.6548);
  //optimistic guess:
  discrete_poisson.set_Ainv(1.0/0.146);

  //CachedProblem<BiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 4.19, 1.0/0.146);

  const double epsilon = 0.00005;

  InfiniteVector<double, Index> u_epsilon0;
  InfiniteVector<double, Index> u_epsilon1;
  InfiniteVector<double, Index> u_epsilon;
  

#if 0
  
  //la = (5,(1),0,(23))
  //nu = (4,(0),1,(4))

  //d=dt=3
  MultiIndex<unsigned int,1> e1;
  e1[0] = 1;
  MultiIndex<int,1> k1;
  k1[0] = 23;

  MultiIndex<unsigned int,1> e2;
  e2[0] = 0;
  MultiIndex<int,1> k2;
  k2[0] = 4;

  Index ind_1(&frame, 5, e1, 0, k1);
  Index ind_2(&frame, 4, e2, 1, k2);
 
  //d=dt=2
//   MultiIndex<unsigned int,1> e1;
//   e1[0] = 1;
//   MultiIndex<int,1> k1;
//   k1[0] = 11;

//   MultiIndex<unsigned int,1> e2;
//   e2[0] = 0;
//   MultiIndex<int,1> k2;
//   k2[0] = 4;

//   Index ind_1(&frame, 6, e1, 1, k1);
//   Index ind_2(&frame, 3, e2, 0, k2);

  map<double,double> integral_values;

  cout.precision(12);

  for (unsigned int i = 0; i < 300; i++) {
    double d = discrete_poisson.a_quad(ind_1, ind_2, 2, i+1);
    cout << "Result = " << d << endl;
    integral_values[i+1] = fabs(d);
  }

  std::ofstream os3("integral_values.m");
  matlab_output(integral_values,os3);
  os3.close();
#endif

  clock_t tstart, tend;
  double time;
  tstart = clock();
//   for (Index ind = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
//        ind <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, frame.j0()+2); ++ind)
//     {
//       cout << ind << endl;
//     }

  multiplicative_Schwarz_SOLVE(problem, epsilon, u_epsilon0, u_epsilon1, u_epsilon);
  //additive_Schwarz_SOLVE(problem, epsilon, u_epsilon);
  //steepest_descent_SOLVE(problem, epsilon, u_epsilon);
  //  steepest_descent_SOLVE_basis(problem, epsilon, u_epsilon);
  //richardson_SOLVE_CDD2(problem, epsilon, u_epsilon);
  //CDD1_SOLVE(problem, epsilon, u_epsilon, 7, CDD1);

  //steepest_descent_SOLVE(discrete_poisson, epsilon, u_epsilon);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "steepest descent done" << endl;

  u_epsilon0.scale(&discrete_poisson,-1);
  u_epsilon1.scale(&discrete_poisson,-1);
  u_epsilon.scale(&discrete_poisson,-1);

  EvaluateFrame<Basis1D,1,1> evalObj;

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u_epsilon0, true, 10);//expand in primal basis
  Array1D<SampledMapping<1> > U1 = evalObj.evaluate(frame, u_epsilon1, true, 10);//expand in primal basis
  Array1D<SampledMapping<1> > U2 = evalObj.evaluate(frame, u_epsilon, true, 10);//expand in primal basis
  cout << "...finished plotting approximate solution" << endl;
  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(frame, u_epsilon, exactSolution1D, 10);
  cout << "...finished plotting error" << endl;
  
  std::ofstream ofs5("approx_sol_mult_schw_1D_out0.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("approx_sol_mult_schw_1D_out1.m");
  matlab_output(ofs6,U1);
  ofs6.close();

  std::ofstream ofs7("approx_sol_mult_schw_1D_out_global.m");
  matlab_output(ofs7,U2);
  ofs6.close();


  std::ofstream ofs8("error_mult_schw_1D_out.m");
  matlab_output(ofs8,Error);
  ofs8.close();
  return 0;


}

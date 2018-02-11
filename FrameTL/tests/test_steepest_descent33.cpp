#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define PARALLEL 0

#define JMAX 18
#define ONE_D
#define _DIM 1

#define PRIMALORDER 3
#define DUALORDER   3

#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <elliptic_equation.h>
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
#include <steepest_descent.h>
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

  //typedef DSBasis<3,3> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  //typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;


  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 0.7;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.7;
  Point<1> b2;
  b2[0] = 1-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);


  FixedArray1D<double,1> A3;
  A3[0] = 0.75;
  SimpleAffineLinearMapping<1> simlpeaffine1(A3,b);
  
  FixedArray1D<double,1> A4;
  A4[0] = 0.75;
  SimpleAffineLinearMapping<1> simlpeaffine2(A4,b2);

  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;
  
  //charts[0] = &simlpeaffine1;
  //charts[1] = &simlpeaffine2;


  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
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
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
 
  bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //Frame1D frame(&Lshaped, bc, bcT, jmax);
  Frame1D frame(&Lshaped, bc, jmax);

//   Vector<double> value(1);
//   value[0] = 384;
//   ConstantFunction<DIM> const_fun(value);

  //  Singularity1D_RHS<double> sing1D;
  //  Singularity1D<double> exactSolution1D;

  Singularity1D_RHS_2<double> sing1D;
  Singularity1D_2<double> exactSolution1D;

  
  //PoissonBVP<DIM> poisson(&const_fun);
  PoissonBVP<DIM> poisson(&sing1D);
  //BiharmonicBVP<DIM> biharm(&const_fun);  
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
//   EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  //BiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&biharm, &frame, jmax);

  // (0,0.7) \cup (0.3,1) DSBasis
  
  //   // (d,dT) = (2,2)
  //   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.3743, 1.0/0.146);
  //   discrete_poisson.set_norm_A(3.3743);
  //   // optimistic guess:
  //   discrete_poisson.set_Ainv(1.0/0.146);
  
//     // (d,dT) = (3,5)
//     CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.4545, 1.0/0.146);
//     discrete_poisson.set_norm_A(4.4545);
//     // optimistic guess:
//     discrete_poisson.set_Ainv(1.0/0.146);

  // (d,dT) = (3,3)
    CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
    discrete_poisson.set_norm_A(1.0);
    // optimistic guess:
    discrete_poisson.set_Ainv(1.0);

//    CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0581, 1.0/0.146);
//    discrete_poisson.set_norm_A(5.0581);
//    // optimistic guess:
//    discrete_poisson.set_Ainv(1.0/0.146);

//   CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
//   // most optimistic guess:
//   discrete_poisson.set_norm_A(1.0);
//   discrete_poisson.set_Ainv(1.0);

  //CachedProblem<BiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 4.19, 1.0/0.146);
//   discrete_biharmonic.set_norm_A(5.0581);
//   // optimistic guess:
//   discrete_biharmonic.set_Ainv(1.0/0.146);



  
//     // (d,dT) = (4,6)
//     CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 7.1803, 1.0/0.146);
//     discrete_poisson.set_norm_A(7.1803);
//     // optimistic guess:
//     discrete_poisson.set_Ainv(1.0/0.146);
  
 
  // (0,0.7) \cup (0.3,1) PBasis
  //   // (d,dT) = (2,2)
  //   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.3697, 1.0/0.146);
  //   discrete_poisson.set_norm_A(3.3697);
  //   // optimistic guess:
  //   discrete_poisson.set_Ainv(1.0/0.146);
  
  //   // (d,dT) = (3,5)
  //   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.6548, 1.0/0.146);
  //   discrete_poisson.set_norm_A(3.6548);
  //   // optimistic guess:
  //   discrete_poisson.set_Ainv(1.0/0.146);

//     // (d,dT) = (3,3)
//     CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.8782, 1.0/0.146);
//     discrete_poisson.set_norm_A(4.8782);
//     // optimistic guess:
//     discrete_poisson.set_Ainv(1.0/0.146);

  
  //   // (d,dT) = (4,6)
  //   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.8748, 1.0/0.146);
  //   discrete_poisson.set_norm_A(5.8748);
  //   // optimistic guess:
  //   discrete_poisson.set_Ainv(1.0/0.146);

  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 2.13, 1.0/0.0038);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 2.47, 1.0/0.0751);
  //  CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.3076, 1.0/0.1);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.19, 1.0/0.146);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.032, 1.0/(1.0e-3*0.672));
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 3.032, 1.0/0.01);


  //interval case:

  // p basis d = 2 dt = 2
  //   discrete_poisson.set_norm_A(3.3076);
  //   discrete_poisson.set_Ainv(1.0/0.1);
  
  
  // p basis d = 3 dt = 3
  //   discrete_poisson.set_norm_A(4.19);
  //   discrete_poisson.set_Ainv(1.0/0.146);
  
  
  // d = 2 dt = 2
  //   discrete_poisson.set_norm_A(2.47);
  //   discrete_poisson.set_Ainv(1.0/0.0751);
  
  // d = 3 dt = 3
  //discrete_poisson.set_norm_A(2.13);
  //discrete_poisson.set_Ainv(1.0/0.0038);
  
  // d = 4 dt = 4
  //  discrete_pisson.set_norm_A(3.032);
  //  discrete_poisson.set_Ainv(1.0/(1.0e-3*0.672));
  
  
  const double epsilon = 1.0e-6;
  //InfiniteVector<double, Index> u_epsilon;


  clock_t tstart, tend;
  double time;
  tstart = clock();

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  steepest_descent_SOLVE(problem, epsilon, approximations);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "steepest descent done" << endl;

  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);
  //u_epsilon.scale(&discrete_poisson,-1);
  //u_epsilon.scale(&discrete_biharmonic,-1);
  

  EvaluateFrame<Basis1D,1,1> evalObj;

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 12);//expand in primal basis
  cout << "...finished plotting approximate solution" << endl;
  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], exactSolution1D, 12);
  cout << "...finished plotting error" << endl;
  
  std::ofstream ofs5("./sd_results33/approx_sol_steep_1D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("sd_results33/error_steep_1D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[128];
    sprintf(filename3, "%s%d%s%d%s%d%s", "sd_results33/approx1Dsteep_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

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
  
    //cout << log10(fabs(*it)) << endl;
  }

  std::ofstream ofs7("./sd_results33/indices_patch_0.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[0]->bases()[0], indices[0], JMAX, ofs7, "jet", true, -16);

  std::ofstream ofs8("./sd_results33/indices_patch_1.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[1]->bases()[0], indices[1], JMAX, ofs8, "jet", true, -16);
  // compute infinite vectors of 1D indices, one for each patch
  // and plot them

  return 0;

}

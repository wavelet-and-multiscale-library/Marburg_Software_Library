#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define PARALLEL 0

#define JMAX 12
#define ONE_D
#define _DIM 1

#define AGGREGATED
#undef BASIS
#undef FRAME



#ifdef BASIS
#define DYADIC
#endif
#ifdef FRAME
#define PMAX 0
#define DYADIC
#endif

#define PRIMALORDER 3
#define DUALORDER   3

#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#ifdef AGGREGATED
#include <elliptic_equation.h>
#include <simple_elliptic_equation.h>
#include <biharmonic_equation.h>
#endif
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <numerics/sturm_bvp.h>
#include <galerkin/sturm_equation.h>
#ifdef AGGREGATED
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#endif
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#ifdef AGGREGATED
#include <frame_support.h>
#include <frame_index.h>
#include <steepest_descent.h>
#include <additive_Schwarz.h>
#include <steepest_descent.h>
#include <richardson_CDD2.h>

#endif
#include <galerkin/cached_problem.h>
#include <galerkin/TestProblem.h>
#include <utils/plot_tools.h>
#include <interval/i_indexplot.h>
#include <interval/pq_frame.h>
#ifdef FRAME
#include <galerkin/cached_quarklet_problem.h>
#endif

using std::cout;
using std::endl;
#ifdef AGGREGATED
using FrameTL::FrameIndex;
//using FrameTL::EllipticEquation;
using FrameTL::SimpleEllipticEquation;
using FrameTL::BiharmonicEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
#endif
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
#ifdef AGGREGATED
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;
using WaveletTL::CachedProblem;
#endif

using namespace std;
#ifdef AGGREGATED
using namespace FrameTL;
#endif
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
  Basis1D basis(1,1);
  basis.set_jmax(jmax);
#ifdef AGGREGATED
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef Frame1D::Index Index;
#endif
#ifdef FRAME
  const int pmax = PMAX;
  typedef PQFrame<d,dT> Frame1D;
  Frame1D frame(true, true, false);
  
  frame.set_jpmax(jmax,pmax);
  typedef Frame1D::Index Index;
#endif
  //typedef CubeBasis<Basis1D,1> IntervalBasis;
  

#ifdef AGGREGATED
  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 1;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(1);
  charts[0] = &affineP;

  SymmetricMatrix<bool> adj(1);
  adj(0,0) = 1;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(1);
  
  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  
  bc[0] = bound_1;
  
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
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
#endif
#ifdef BASIS
  TestProblem<9> T;
  SturmEquation<Basis1D> discrete_poisson(T, basis);
  CachedProblem<SturmEquation<Basis1D> > problem(&discrete_poisson);
#endif
#ifdef FRAME
  TestProblem<9> T;
  SturmEquation<Frame1D> discrete_poisson(T, frame);
  CachedQuarkletProblem<SturmEquation<Frame1D> > problem(&discrete_poisson);
#endif
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
    
    
    
    cout.precision(12);
  
  
#ifdef AGGREGATED 
  set<Index> Lambda;
  for (int i=0; i<frame.degrees_of_freedom();i++) {
    Lambda.insert(*frame.get_wavelet(i));
        cout << *frame.get_wavelet(i) << endl;
  }
#endif
 
#ifdef FRAME 
set<Index> Lambda;
for (int i=0; i<frame.degrees_of_freedom();i++) {
  Lambda.insert(*frame.get_quarklet(i));
      cout << *frame.get_quarklet(i) << endl;
}
#endif
  
#ifdef BASIS
  set<Basis1D::Index> Lambda;
  for (int i=0; i<basis.degrees_of_freedom();i++) {
    Lambda.insert(*basis.get_wavelet(i));
        cout << *basis.get_wavelet(i) << endl;
  }
#endif  
//  
//  
//  
////  for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
////       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, jmax); ++lambda) {
////    Lambda.insert(lambda);
////    cout << lambda << endl;
////  }
//  
//  
  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(problem, Lambda, rh);
//  cout << rh << endl;
  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;
  
  clock_t tstart, tend;
  double time;
  tstart = clock();

//    WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff, false);
//    WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff, false);
//  WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff);
  

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  stiff.matlab_output("stiff_1D_out", "stiff",1); 



  

  return 0;

}

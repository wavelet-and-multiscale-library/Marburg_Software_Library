#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <simple_biharmonic_equation.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <adaptive/apply.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/cached_problem.h>
#include <frame_support.h>
#include <frame_index.h>


using std::cout;
using std::endl;


using FrameTL::FrameIndex;
using FrameTL::SimpleBiharmonicEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using MathTL::BiharmonicBVP;
using MathTL::ConstantFunction;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;



/*!
*/
template<class VALUE = double>
class Singularity1D_RHS_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_RHS_2() {};
  virtual ~Singularity1D_RHS_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    double x=p[0];
    return 16*M_PI*M_PI*M_PI*M_PI*cos(2*M_PI*x)-72;
    //return 0;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};

/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_2() {};
  virtual ~Singularity1D_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    //function for constant rhs
    //return 16*(p[0]*p[0]*p[0]*p[0]-2*p[0]*p[0]*p[0]+p[0]*p[0]);
    
    //Non-constant function with singularity in their first order
    double x=p[0];
    if(x<0.5)
      return cos(2*M_PI*x)+2*x*x*x-3*x*x*x*x-1;
    else return (cos(2*M_PI*x)-3*(1-x)*(1-x)*(1-x)*(1-x)+2*(1-x)*(1-x)*(1-x)-1);

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};


int main()
{
  cout << "Testing class BiharmonicEquation..." << endl;
  
  const int DIM = 1;
  int jmax=9;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<3,3> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;
  typedef Basis1D::Index IIndex;

  EvaluateFrame<Basis1D,1,1> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 0.7;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.7;
  Point<1> b2;
  b2[0] = 1.0-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);


  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;
  
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

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

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //Frame1D frame(&Lshaped, bc, bcT);
  Frame1D frame(&Lshaped, bc, jmax);

  Vector<double> value(1);
  value[0] = 384;
  ConstantFunction<DIM> const_fun(value);

  Singularity1D_2<double> exactSolution;
  Singularity1D_RHS_2<double> sing1D;
  
  BiharmonicBVP<DIM> biharmonic(&sing1D);

  SimpleBiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&biharmonic, &frame, jmax, TrivialAffine);
 
  CachedProblem<SimpleBiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 5, 1.0/0.146);

  cout.precision(12);
  
  //############### 1D galerkin scheme test ##################
#if 1

  set<Index> Lambda;
  for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
  }

  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(discrete_biharmonic, Lambda, rh);
  
  cout << rh << endl;
  
 //  Basis1D basis(3,3);

//   InfiniteVector<double, IIndex> coeff;
//   IIndex index(first_generator(&basis, basis.j0()));
//   for (int i = 0;; ++index, i++) {
//     coeff.set_coefficient(index, rh[i]);
//     if (index == last_generator(&basis, basis.j0())) break;
//   }

//   SampledMapping<1> res = evaluate(basis, coeff, false, 8);
  
//   std::ofstream ofs4("reproduced_function.m");
//   res.matlab_output(ofs4);
//   ofs4.close();

//   abort();

  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;
  
  clock_t tstart, tend;
  double time;
  tstart = clock();
  
  WaveletTL::setup_stiffness_matrix(discrete_biharmonic, Lambda, stiff);
  
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  stiff.matlab_output("stiff_out", "stiff",1);

//   unsigned int iter= 0;
//   Vector<double> x(Lambda.size()); x = 1;
//   double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
//   cout << "lmax = " << lmax << endl;

  cout << "performing iterative scheme to solve projected problem..." << endl;
  Vector<double> xk(Lambda.size()); xk = 0;


  double alpha_n = 0.07;
  Vector<double> resid(xk.size());
  Vector<double> help(xk.size());
  for (int i = 0; i <30000; i++) {
    stiff.apply(xk,help);
    resid = rh - help;
    cout << ".loop = " << i << " " << " res = " << sqrt((resid*resid)) << endl;
    stiff.apply(resid,help);
    alpha_n = (resid * resid) * (1.0 / (resid * help));
    //cout  << alpha_n << endl;
    resid *= alpha_n;
    //resid *= 0.3;
    xk = xk + resid;
  }



  //CG(stiff, rh, xk, 0.0001, 1000, iter);
  //Richardson(stiff, rh, xk, 2. / lmax, 0.0001, 1000, iter);
  //Richardson(stiff, rh, xk, 0.07, 0.0001, 2000, iter);  
  cout << "performing output..." << endl;
  
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, xk[i]);
  
  u.scale(&discrete_biharmonic,-1);

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u, true, 11);//expand in primal basis
  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(frame, u, exactSolution, 11);


 
  std::ofstream ofs5("approx_solution_1D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("error_1D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  // cout << "  ... done, time needed: " << time << " seconds" << endl;
   
#endif



   return 0;

}


#define ONE_D

#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <simple_elliptic_equation.h>
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
#include <algebra/vector_norms.h>
#include <utils/plot_tools.h>
//#include <steepest_descent.h>

using std::cout;
using std::endl;


using FrameTL::FrameIndex;
using FrameTL::SimpleEllipticEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
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
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
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
    if (0. <= p[0] && p[0] < 0.5)
      return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];

    if (0.5 <= p[0] && p[0] <= 1.0)
      return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);

    return 0.;

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};




int main()
{
  const int DIM = 1;
  const int jmax = 14;

  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<3,3> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;
  typedef Basis1D::Index IIndex;

  EvaluateFrame<Basis1D,1,1> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 1;//0.7
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.7;
  Point<1> b2;
  b2[0] = 1.0-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);

  FixedArray1D<double,1> A3;
  A3[0] = 0.75;
  SimpleAffineLinearMapping<1> simpleaffine1(A3,b);
  
  FixedArray1D<double,1> A4;
  A4[0] = 0.75;
  SimpleAffineLinearMapping<1> simpleaffine2(A4,b2);


  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(1);
  charts[0] = &affineP;
  //  charts[1] = &affineP2;

  //charts[0] = &simpleaffine1;
  //charts[1] = &simpleaffine2;
  
  SymmetricMatrix<bool> adj(1);
  adj(0,0) = 1;
//   adj(1,1) = 1;
//   adj(1,0) = 1;
//   adj(0,1) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(1);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 2;
  bound_2[1] = 1;

  //bc[1] = bound_2;

//to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(1);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
 
  //bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //Frame1D frame(&Lshaped, bc, bcT);
  Frame1D frame(&Lshaped, bc, jmax);

  Vector<double> value(1);
  value[0] = 1;
  ConstantFunction<DIM> const_fun(value);

  Singularity1D_2<double> exactSolution;
  Singularity1D_RHS_2<double> sing1D;
  
  //PoissonBVP<DIM> poisson(&const_fun);
  PoissonBVP<DIM> poisson(&sing1D);
  //IdentityBVP<DIM> trivial_bvp(&const_fun);
  //IdentityBVP<DIM> trivial_bvp(&exactSolution);

  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&trivial_bvp, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.19, 1.0/0.146);

  cout.precision(12);
  
  InfiniteVector<double,Frame1D::Index> rhs;
  Vector<double> v(225);
  
  set<Index> Lambda;
  for (FrameIndex<Basis1D,1,1> lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
  }
  
  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(discrete_poisson, Lambda, rh);
//   //cout << rh << endl;

  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;

 
  //WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  //cout << "NUMBER OF NONZERO ENTRIES = " << stiff.size() << endl;

  InfiniteVector<double, Index> w, f, w_res1, w_res2, w_eps, w_exact, error, sparse_f;
  
//   // copy rh to w
//   for (unsigned int i = 0; i < rh.size(); i++) {
//     if (rh[i] != 0.) {
//       Index ind = *problem.basis().get_wavelet(i);
//       f.set_coefficient(ind, rh[i]);
//     }
//   }
  
  cout << "number of nonzeros in right hand side " << w.size() << endl;


  
//   //APPLY_TEST(problem, f, 0.0001, w_res1, jmax, St04a);
//   APPLY(problem, f, 0.00001, w_res1, jmax, St04a);

//   APPLY_OPTIMIZED(problem, f, 0.00001, w_res2, jmax, St04a);

//   cout << "difference = " << l2_norm(w_res1-w_res2) << endl;
  
//  abort();
//   APPLY(problem, w, 0.0001, w_res2, jmax, CDD1);
//   cout << "difference = " << l2_norm(w_res1-w_res2) << endl;


//   cout << "size of result 1 = " << w_res1.size() << endl;
//   cout << "size of result 2 = " << w_res2.size() << endl;
  
  const double alpha = 0.2;
  w = alpha * f;
  double eta = 1.;
  int step = 1;
  map<double,double> asymptotic;
  map<double,double> alltime;
  double time = 0;
  for (int i = 0; i < 1000; i++) {
    
    //APPLY(problem, w, eta, w_res1, jmax, CDD1);
    cout << "eta = " << eta << endl;
    //APPLY_OPTIMIZED(problem, w, eta, w_res1, time, jmax, St04a);
    APPLY_TEST(problem, w, eta, w_res1, jmax, St04a);
    problem.RHS(eta,sparse_f);
    w_res2 = sparse_f - w_res1;

//     if (w.size() != 0)
//       asymptotic[log10( (double)w.size() )] = log10(l2_norm(w_res2));
//     std::ofstream os3("asymptotic_0205_07.m");
//     matlab_output(asymptotic,os3);
//     os3.close();

//     alltime[step] = time;
//     //std::ofstream os4("old_approach_time_0205_07.m");
//     std::ofstream os4("time_0205_07.m");
//     matlab_output(alltime,os4);
//     os4.close();
 
    cout << "residual norm in step " << step << "= " << l2_norm(w_res2) << endl;
    w += alpha * w_res2;
    w.COARSE(eta,w_res2);
    w = w_res2;
    cout << "size of result 1 = " << w_res1.size() << endl;
    eta *= 0.95;
    step++;
   
    
    
  }


  // old code
#if 0

  Vector<double> exact;
  //stiff.apply<Vector<double> >(rh, exact);
  APPLY(problem, w, 0., w_exact, jmax, St04a);

  map<double,double> list_of_errors_supp;
  map<double,double> list_of_errors_time;

  double s = 0.9999;

  clock_t tstart, tend;
  double time = 0;
  

  cout << "size of vector w = " << w.size() << endl;
  double eps = 2.;
  for (int i = 0; i < 100; i++) {
    cout << "i = " << i << " eps = " << eps << endl;
    tstart = clock();
    APPLY_COARSE(problem, w, s*eps, w_eps, (1-s)*eps, jmax, St04a);
    //APPLY(problem, w, eps, w_eps, 5, St04a);
    tend = clock();
    time += (double)(tend-tstart)/CLOCKS_PER_SEC;
        
//     Vector<double> w_eps_vec;
//     w_eps_vec.resize(rh.size());
    
//     // copy w_eps to w_eps
//     InfiniteVector<double, Index>::const_iterator it = w_eps.begin();
//     int k = 0;
//     for (; it != w_eps.end(); it++) {
//       w_eps_vec[k] = (*it);
//       k++;
//     }


//     Vector<double> error;
    
//     error = exact - w_eps_vec;
//     double d = l2_norm<Vector<double> >(error);
//     cout << "d = " << d << endl;
    
    list_of_errors_time[log10(time)] = log10(l2_norm(w_exact - w_eps));
    list_of_errors_supp[log10((double)w_eps.size())] = log10(l2_norm(w_exact - w_eps));
    
    eps *= 0.9;
    
  }
  
  std::ofstream os3("apply_errors_time.m");
  matlab_output(list_of_errors_time,os3);
  os3.close();
  
  std::ofstream os4("apply_errors_supp.m");
  matlab_output(list_of_errors_supp,os4);
  os4.close();
#endif


   return 0;

}

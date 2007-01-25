#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <biharmonic_equation.h>
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
using FrameTL::BiharmonicEquation;
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
    return 0;
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
    
    return 16*(p[0]*p[0]*p[0]*p[0]-2*p[0]*p[0]*p[0]+p[0]*p[0]);
    

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};


int main()
{
  cout << "Testing class BiharmonicEquation..." << endl;
  
  const int DIM = 1;
  int jmax=7;

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

  FixedArray1D<double,1> A3;
  A3[0] = 0.75;
  SimpleAffineLinearMapping<1> simpleaffine1(A3,b);
  
  FixedArray1D<double,1> A4;
  A4[0] = 0.75;
  SimpleAffineLinearMapping<1> simpleaffine2(A4,b2);


  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;

  //charts[0] = &simpleaffine1;
  //charts[1] = &simpleaffine2;
  
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
  
  BiharmonicBVP<DIM> biharmonic(&const_fun);
  //PoissonBVP<DIM> poisson(&sing1D);
  //IdentityBVP<DIM> trivial_bvp(&const_fun);
  //IdentityBVP<DIM> trivial_bvp(&exactSolution);

    BiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&biharmonic, &frame, jmax, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&trivial_bvp, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  
  CachedProblem<BiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 4.19, 1.0/0.146);

  // double tmp = 0.0;
//   int c = 0;
//   int d = 0
    ;
  
  cout.precision(12);
  
  //############### 1D galerkin scheme test ##################
#if 1

  //int z = 0;
  set<Index> Lambda;
  for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, jmax); ++lambda) {
    //cout << lambda << endl;
    //cout << z++ << endl;
    Lambda.insert(lambda);
  }

  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(discrete_biharmonic, Lambda, rh);
  
//   cout << rh << endl;
  
//   Basis1D basis(3,3);

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


//   for (int i = 0; i < 30 ; i++) 
//     for (int j = 0; j < 30 ; j++) {
//       if (! (fabs(stiff.get_entry(i,j) -  stiff.get_entry(j,i)) < 1.0e-13)) {
// 	cout << stiff.get_entry(i,j) << endl;
// 	cout << stiff.get_entry(j,i) << endl;
// 	cout << "i = " << i << " j = " << j << endl;
// 	//abort();
// 	cout << "#######################" << endl;
//       }
//     } 

  double alpha_n = 0.07;
  Vector<double> resid(xk.size());
  Vector<double> help(xk.size());
  for (int i = 0; i <3000; i++) {
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
	  
  
#if 0
  char filename1[50];
  u.clear();
  Lambda.clear();
  int l = 0;
  for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, frame.j0()+7); ++lambda) {


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
    cout << "l= " << l << " lambda:" << endl;
    for (unsigned int i = 0; i < 1; i++) {
      cout << "j = " << supp_lambda.j << endl;
      cout << "a = " << supp_lambda.a[i]*dx <<  " b= " << supp_lambda.b[i]*dx << " p = " << lambda.p() << endl;
    }    
    
//     if (l != 14) {
//       l++;
//       continue;
 //     }

    cout << lambda << endl;
    cout << "##################" << endl;
    u.set_coefficient(lambda, 1);
    
    Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u, true, 11);//expand in primal basis
    sprintf(filename1, "%s%d%s", "wav_1D_", l, ".m");
    
    
    std::ofstream ofs5(filename1);
    matlab_output(ofs5,U);
    ofs5.close();

    u.clear();
    ++l;

  }
#endif


   return 0;

}


#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include "elliptic_equation.h"
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
using FrameTL::EllipticEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;
using MathTL::IdentityBVP;
using MathTL::ConstantFunction;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


//   /*!
//     special function with steep gradients
//     near the right end of the interval
//   */
// template<class VALUE = double>
// class Singularity1D
//   : public Function<1, VALUE>
// {
// public:
//   Singularity1D() {};
//   virtual ~Singularity1D() {};
//   VALUE value(const Point<1>& p,
// 	      const unsigned int component = 0) const
//   {
//     return  -100*exp(5*p[0])*(1-(exp(5*p[0])-1)/(exp(5.)-1))/(exp(5.)-1)+200*exp(10*p[0]) / 
//       ((exp(5.)-1)*(exp(5.)-1))+100*(exp(5*p[0])-1)*exp(5*p[0])/((exp(5.)-1)*(exp(5.)-1));
//   }
  
//   void vector_value(const Point<1> &p,
// 		    Vector<VALUE>& values) const { ; }
  
// };


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
  
  cout << "Testing class EllipticEquation..." << endl;
  
  const int DIM = 1;
  const int d  = 3;
  const int dt = 5;

  typedef DSBasis<d,dt> Basis1D;
  //typedef PBasis<2,4> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
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
  bound_1[0] = 1;
  bound_1[1] = 1;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
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

  Atlas<DIM,DIM> interval(charts,adj);  
  cout << interval << endl;

  //finally a frame can be constructed
  //Frame1D frame(&interval, bc, bcT);
  Frame1D frame(&interval, bc, 9);

  Vector<double> value(1);
  value[0] = 1;
  ConstantFunction<DIM> const_fun(value);

  Singularity1D_2<double> exactSolution;
  Singularity1D_RHS_2<double> sing1D;
  
  // the constant function in the argument of the boundary
  // value problem is just a dummy.
  // this object is anyway only created to have access to the gramian 
  // matrix of the underlying frame.
  IdentityBVP<DIM> trivial_bvp(&const_fun);
  //PoissonBVP<DIM> trivial_bvp(&sing1D);

  EllipticEquation<Basis1D,DIM> discrete_poisson(&trivial_bvp, &frame, TrivialAffine);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 4.19, 1.0/0.146);

  cout.precision(12);
  
  //############### 1D galerkin scheme test ##################
#if 1

  int z = 0;
  set<Index> Lambda;
  for (Index lambda = FrameTL::first_generator<Basis1D,1,1,Frame1D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,1,1,Frame1D>(&frame, frame.j0()+4); ++lambda) {
    //cout << z++ << endl;
    Lambda.insert(lambda);
  }

  cout << "setting up full gramian matrix..." << endl;
  SparseMatrix<double> gramian;
  
  clock_t tstart, tend;
  double time;
  tstart = clock();
  
  WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, gramian);
  //WaveletTL::setup_stiffness_matrix(problem, Lambda, gramian);  
  // symmetry check
  for (int i = 0; i < Lambda.size() ; i++) 
    for (int j = 0; j < Lambda.size()  ; j++) {
      if (! (fabs(gramian.get_entry(i,j) -  gramian.get_entry(j,i)) < 1.0e-13)) {
	cout << gramian.get_entry(i,j) << endl;
	cout << gramian.get_entry(j,i) << endl;
	cout << "i = " << i << " j = " << j << endl;
	//abort();
	cout << "#######################" << endl;
      }
    } 

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  typedef SparseMatrix<double>::size_type size_type;


  //  gramian.matlab_output("gramian_out", "gramian",1);

  SparseMatrix<double> G;
  G.resize(Lambda.size(), Lambda.size());
  // ##########################################################################
  cout << "computing square of gramian explicitely..." << endl;
 
  //G = gramian;
  size_type row = 0;
  for (std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
       it1 != itend; ++it1, ++row) {
    std::list<size_type> indices;
    std::list<double> entries;
    
    size_type column = 0;
    for (std::set<Index>::const_iterator it2(Lambda.begin());
	 it2 != itend; ++it2, ++column) {
      InfiniteVector<double, size_type> row_of_gramian;
      InfiniteVector<double, size_type> col_of_gramian;
      gramian.get_row(row, row_of_gramian);
      gramian.get_row(column, col_of_gramian); // gramian is symmetric!
      double entry = row_of_gramian * col_of_gramian;
      if (entry != 0.) {
	indices.push_back(column);
	entries.push_back(entry);
      }
    }
    cout << "row " << row << " completed" << endl;
    G.set_row(row, indices, entries);
  }

  //G = gramian * gramian;
  cout << "done computing square of gramian..." << endl;
  // ##########################################################################


//   cout << "setting up full right hand side..." << endl;
//   WaveletTL::setup_righthand_side(discrete_poisson, Lambda, f);
//   //cout << f << endl;

  for (int i = 100; i < 200; i++) {
    
    InfiniteVector<double, size_type> ff;
    int num = i;
    gramian.get_row(num, ff);
    
    Vector<double> f;
    f.resize(Lambda.size());
    
    InfiniteVector<double, size_type>::const_iterator it = ff.begin();
    for (; it != ff.end(); ++it) {
      f[it.index()] = *it;
    }
  
    cout << "performing iterative scheme to solve projected problem..." << endl;
    Vector<double> xk(Lambda.size()); xk = 0;
    
    double alpha_n = 0.07;
    Vector<double> resid(xk.size());
    Vector<double> help(xk.size());
    for (int i = 0; i < 1000; i++) {
      G.apply(xk,help);
      resid = f - help;
      cout << ".loop = " << i << " " << " res = " << sqrt((resid*resid)) << endl;
      G.apply(resid,help);
      alpha_n = (resid * resid) * (1.0 / (resid * help));
      //cout  << alpha_n << endl;
      resid *= alpha_n;
      xk = xk + resid;
    }

    //   unsigned int iter = 0;
    //   CG(G, f, xk, 0.00001, 1022, iter);
    //   cout << iter << " iterations needed in CG" << endl; 
    //Richardson(G, f, xk, 2. / lmax, 0.0001, 1000, iter);
    //Richardson(G, f, xk, 0.07, 0.0001, 2000, iter);  
    cout << "performing output..." << endl;
    
    InfiniteVector<double,Index> u;
    unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      u.set_coefficient(*it, xk[i]);
    
    //u.scale(&discrete_poisson,-1);
    
    EvaluateFrame<Basis1D,1,1> evalObj;
    
    Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u, true, 11);//expand in primal basis
    
    char filename1[50];
    
    sprintf(filename1, "%s%d%s%d%d%s", "dual_frame_elem_1D_ds_", num, "_", d, dt, ".m");
    //sprintf(filename1, "%s%d%s", "dual_basis_elem_1D_ds_", num, ".m");
    
    std::ofstream ofs5(filename1);
    matlab_output(ofs5,U);
    ofs5.close();
  }
    cout << "  ... done, time needed: " << time << " seconds" << endl;
#endif
	  
  
   return 0;

}

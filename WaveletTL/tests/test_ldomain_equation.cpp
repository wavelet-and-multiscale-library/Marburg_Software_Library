#undef BASIS
#define FRAME

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 1
#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_index.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_frame_index.h>
#include <Ldomain/ldomain_frame_support.h>
#include <Ldomain/ldomain_frame_evaluate.h>


#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <galerkin/ldomain_equation.h>
#include <galerkin/ldomain_helmholtz_equation.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the L-shaped domain ..." << endl;

  const int d  = 2;
  const int dT = 2;
//  typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
  
#ifdef FRAME
  typedef PQFrame<d,dT> Frame1D;
  
  Frame1D frame1D;
  
  typedef LDomainFrame<Frame1D> LFrame;
  typedef LFrame::Index Index;
  typedef Index::polynomial_type polynomial_type;
  polynomial_type p;
  typedef Index::level_type level_type;
  LFrame frame(frame1D);
  frame.set_jpmax(5,3);
  
  InfiniteVector<double, Index> coeffs;
  Index lambda = frame.first_generator(frame.j0());
  typedef LFrame::Support Support;
  
  
  
  Support suppconst;
  suppconst.j[0]=0, suppconst.j[1]=0;
  for(int i=0;i<3;i++){
      suppconst.xmin[i]=0, suppconst.xmax[i]=0,suppconst.ymin[i]=0,suppconst.ymax[i]=0;
  }
  
  
  set<Index> Lambda;
  int counter = 0;
  for (Index lambda = frame.first_generator(frame.j0());counter<=3000;) {
      Support supp = suppconst;
      frame.support(lambda,supp);
      cout << "("<< supp.j[0] << "," << supp.j[1] << ", [" << supp.xmin[0]  << "," << supp.xmax[0]<< "]x["   << supp.ymin[0]<< "," << supp.ymax[0]<< "]"<< ", ["<< supp.xmin[1]<< "," << supp.xmax[1]<< "]x["<< supp.ymin[1]
           << ","<< supp.ymax[1]<< "]"<< ", ["<< supp.xmin[2]<< ","<< supp.xmax[2]<< "]x["<< supp.ymin[2]<< ","<< supp.ymax[2]<< "])"<< endl;
    Lambda.insert(lambda);
//    cout << lambda << ", " << lambda.number() << endl;
    if(lambda==frame.last_quarklet(5,p)){
            ++p;
            lambda=frame.first_generator(frame.j0(), p, counter+1);
        }
            else
                ++lambda;
    
    counter++;
//    if (lambda == frame.last_quarklet(frame.j0())) break;
  }
  cout << "Test" << endl;
  polynomial_type p2(2,3);
  level_type j2(2,3);
  Support supp;
  Index lambda2 = frame.first_generator(frame.j0(), p2);
  Index lambda3 = frame.first_quarklet(j2);
  intersect_supports(frame,lambda2,lambda3,supp);
  cout << "("<< supp.j[0] << "," << supp.j[1] << ", [" << supp.xmin[0]  << "," << supp.xmax[0]<< "]x["   << supp.ymin[0]<< "," << supp.ymax[0]<< "]"<< ", ["<< supp.xmin[1]<< "," << supp.xmax[1]<< "]x["<< supp.ymin[1]
       << ","<< supp.ymax[1]<< "]"<< ", ["<< supp.xmin[2]<< ","<< supp.xmax[2]<< "]x["<< supp.ymin[2]<< ","<< supp.ymax[2]<< "])"<< endl;
//  cout << lambda2 << ", " << lambda2.number() << endl;
  
#endif
#ifdef BASIS
   typedef PBasis<d,dT> Basis1D;
  
  Basis1D basis1D;
//  basis1D.orthogonalize_boundary_wavelets(); // cf. [Barsch]

  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;
  typedef LBasis::Support Support;
  Support supp;

  CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
  CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);

//#if 0
//  PoissonBVP<2> poisson(&f_sing);
//#else
  EigenRHS rhs;
  PoissonBVP<2> poisson(&rhs);
//#endif
//
  LBasis basis(basis1D);
//
//#if 1
  LDomainEquation<Basis1D> eq(&poisson, basis);
//#else
//  LDomainHelmholtzEquation<Basis1D> eq(&rhs, basis, 0.0);
//#endif
//  
  InfiniteVector<double, Index> coeffs;
  eq.RHS(1e-8, coeffs);

  set<Index> Lambda;
  for (Index lambda = eq.basis().first_generator(eq.basis().j0());; ++lambda) {
    basis.support(lambda,supp);
      cout << "("
         << supp.j
         << ", ["
         << supp.xmin[0]
         << ","
         << supp.xmax[0]
         << "]x["
         << supp.ymin[0]
         << ","
         << supp.ymax[0]
         << "]"
         << ", ["
         << supp.xmin[1]
         << ","
         << supp.xmax[1]
         << "]x["
         << supp.ymin[1]
         << ","
         << supp.ymax[1]
         << "]"
         << ", ["
         << supp.xmin[2]
         << ","
         << supp.xmax[2]
         << "]x["
         << supp.ymin[2]
         << ","
         << supp.ymax[2]
         << "])"
         << endl;
      
      Lambda.insert(lambda);
    cout << lambda << endl;
    if (lambda == eq.basis().last_wavelet(eq.basis().j0()+1)) break;
  }

   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
     cout << *it << endl;

  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();

  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
//#if 0
//  std::ofstream ofs("stiff_out.m");
//  ofs << "M=";
//  print_matrix(A,ofs);
//  ofs.close();
//#endif
//
//  tend = clock();
//  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//  cout << "  ... done, time needed: " << time << " seconds" << endl;
////   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;
//
  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- right hand side: " << b << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-8, 100, iterations);
//  
//// //   cout << "- solution coefficients: " << x;
//  cout << " with residual (infinity) norm ";
//  A.apply(x, err);
//  err -= b;
//  cout << linfty_norm(err) << endl;
//
//#if 1
  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);

  u.scale(&eq, -1);
  Array1D<SampledMapping<2> > s((eq.basis()).evaluate(u, 6));
//  std::ofstream u_Lambda_stream("u_lambda.m");
//  matlab_output(u_Lambda_stream, s);
//  u_Lambda_stream.close();
//  cout << "  ... done, see file 'u_lambda.m'" << endl;
//   EigenSolution u_Lambda;
//   s.add(-1.0, SampledMapping<2>(s, u_Lambda));
//   cout << "  ... done, pointwise error: " << row_sum_norm(s.values()) << endl;
//#endif
#endif
  return 0;
}

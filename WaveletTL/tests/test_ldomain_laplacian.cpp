#include <iostream>
#include <set>

#include <utils/function.h>
#include <algebra/sparse_matrix.h>
#include <interval/spline_basis.h>
#include <interval/ds_basis.h>
#include <numerics/iteratsolv.h>
#include <numerics/bezier.h>
#include <numerics/corner_singularity.h>
#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 0
#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_evaluate.h>
#include <Ldomain/ldomain_expansion.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <galerkin/ldomain_laplacian.h>
#include <galerkin/galerkin_utils.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainLaplacian ..." << endl;

//  const int d  = 2;
//  const int dT = 2;

// #if 0
//   typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
//   Basis1D basis1d;
// #else
//   typedef SplineBasis<d,dT,DS_construction> Basis1D;
//   Basis1D basis1d("bio5-energy",0,0,0,0);
// #endif
//   
//   typedef LDomainBasis<Basis1D> Basis;
//   Basis basis(basis1d);
// 
//   typedef Basis::Index Index;
// 
//   const int solution = 1;
//   Function<2> *uexact = 0, *f = 0;
//   switch(solution) {
//   case 1:
//     uexact = new PolySolution();
//     f      = new PolyRHS();
//     break;
//   case 2:
//     uexact = new EigenSolution();
//     f      = new EigenRHS();
//     break;
//   case 3:
//     uexact = new CornerSingularity   (Point<2>(0,0), 0.5, 1.5);
//     f      = new CornerSingularityRHS(Point<2>(0,0), 0.5, 1.5);
//     break;
//   default:
//     break;
//   }
// 
// #if 0
//   {
//     cout << "- plot exact solution into a file..." << endl;
//     const bool matlab = false; // false -> Octave output
//     Grid<2> grid(Point<2>(-1.0, 0.0), Point<2>(0.0, 1.0), 1<<5, 1<<5); // only patch 0
//     SampledMapping<2> u_sm(grid, *uexact);
//     std::ofstream u_fs("uexact.m");
//     
//     if (matlab)
//       u_sm.matlab_output(u_fs);
//     else
//       u_sm.octave_output(u_fs);
//     
//     if (matlab)
//       u_fs << "surf(x,y,z)" << endl;
//     else
//       u_fs << "mesh(x,y,z)" << endl;
//     
//     u_fs.close();
//     cout << "  ... done, see uexact.m" << endl;
//   }
// #endif
// 
//   const int jmax = basis.j0();
//   
//   typedef LDomainLaplacian<Basis1D> Problem;
//   Problem problem(basis, InfiniteVector<double, Index>());
// 
//   cout << "- expand right-hand side..." << endl;
//   InfiniteVector<double,Index> fcoeffs;
//   expand(f, basis, true, jmax, fcoeffs);
//   problem.set_rhs(fcoeffs);
//   cout << "  ... done!" << endl;
// 
//   //   cout << "- integrals of f against the primal wavelets:" << endl
//   //        << fcoeffs << endl;
// 
//   cout << "- set up index set of active wavelets..." << endl;
//   set<Index> Lambda;
//   for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
//     Lambda.insert(lambda);
//     if (lambda == basis.last_wavelet(jmax)) break;
//   }
//   cout << "  ... done!" << endl;
// 
//   cout << "- set up stiffness matrix..." << endl;
//   clock_t tstart, tend;
//   double time;
//   tstart = clock();
// 
//   SparseMatrix<double> A;
//   setup_stiffness_matrix(problem, Lambda, A);
//   tend = clock();
//   time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//   cout << "  ... done, time needed: " << time << " seconds" << endl;
//   //   cout << "- stiffness matrix A=" << endl << A << endl;
// #if 0
//   A.matlab_output("Ldomain_laplacian", "A", 1);
// #endif
// 
//   cout << "- set up right-hand side..." << endl;
//   tstart = clock();
//   Vector<double> b;
//   setup_righthand_side(problem, Lambda, b);
//   tend = clock();
//   time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//   cout << "  ... done, time needed: " << time << " seconds" << endl;
// //   cout << "- right hand side: " << b << endl;
// 
//   Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
//   unsigned int iterations;
//   CG(A, b, x, 1e-15, 200, iterations);
//   
//   //   cout << "- solution coefficients: " << x;
//   cout << "- residual (infinity) norm: ";
//   A.apply(x, err);
//   err -= b;
//   cout << linfty_norm(err) << endl;
//   
// #if 1
//   {
//     cout << "- plot point values of the solution:" << endl;
//     InfiniteVector<double,Index> u;
//     unsigned int i = 0;
//     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
//       u.set_coefficient(*it, x[i]);
//     u.scale(&problem, -1);
//     Array1D<SampledMapping<2> > s(basis.evaluate(u, 5));
//     std::ofstream u_Lambda_stream("u_lambda.m");
// //     octave_output(u_Lambda_stream, s);
//     matlab_output(u_Lambda_stream, s);
//     u_Lambda_stream.close();
//     cout << "  ... done, see file 'u_lambda.m'" << endl;
//  
//     for (int i = 0; i <= 2; i++) {
//       s[i].add(-1.0, SampledMapping<2>(s[i], *uexact));
//       cout << "  pointwise error on patch " << i << ": " << row_sum_norm(s[i].values()) << endl;
//     }
//   }
// #endif
// 
// // #if 0
// //   {
// //     cout << "- compute error:" << endl;
// //     InfiniteVector<double,Index> u;
// //     unsigned int i = 0;
// //     for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
// //       u.set_coefficient(*it, x[i]);
// //     u.scale(&problem, -1);
// //     const int N = 64;
// //     Array1D<SampledMapping<2> > s(evaluate(problem.basis(), u, N));
// //     std::ofstream u_Lambda_stream("u_lambda.m");
// //     octave_output(u_Lambda_stream, s);
// //     u_Lambda_stream.close();
// //     cout << "  ... plot of the approximation: see u_lambda.m" << endl;
//     
// //     double L2error = 0;
// //     for (int ii = 0; ii <= 2; ii++) {
// //       s[ii].add(-1.0, SampledMapping<2>(s[ii], *uexact));
// //       const double frob = frobenius_norm(s[ii].values());
// //       L2error += frob*frob;
// //     }
// //     L2error = sqrt(L2error/(N*N));
// //     cout << "  ... L_2 error " << L2error << endl;
// 
// //     cout << "- plot error into a file:" << endl;
// //     std::ofstream u_Lambda_error_stream("u_lambda_error.m");
// //     octave_output(u_Lambda_error_stream, s);
// //     u_Lambda_error_stream.close();
// //     cout << "  ... done, see file 'u_lambda_error.m'" << endl;
// //   }
// // #endif
// 
//   if (uexact) delete uexact;
//   if (f) delete f;

  return 0;
}

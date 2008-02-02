#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <ring/ring_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/ring_helmholtz.h>
#include <galerkin/generic_helmholtz.h>
#include "ring_functions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing discretizations of the Helmholtz equation with wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 3;
  const int dt = 3;
  const int s0 = 1;
  const int s1 = 1;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,s0,s1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  const int j0 = basis.j0();
  const int jmax = j0+1;

  cout << "j0=" << j0 << ", jmax=" << jmax << endl;

  Function<2> *uexact = 0, *f = 0;
  const int testcase = 1;
  switch(testcase) {
  case 1:
    uexact = new RingFunction4(r0, r1);
    f = new RingFunction4_RHS_Helmholtz(r0, r1);
    break;
  default:
    uexact = new RingFunction4(r0, r1);
    f = new RingFunction4_RHS_Helmholtz(r0, r1);
  }
  
  const int resi = 6;
  Point<2> sphi, y;
  Matrix<double> gridx((1<<resi)+1), gridy((1<<resi)+1);
  for (int i = 0; i <= 1<<resi; i++) {
    sphi[0] = i/(double)(1<<resi);
    for (int j = 0; j <= 1<<resi; j++) {
      sphi[1] = j/(double)(1<<resi);
      basis.chart().map_point(sphi, y);
      gridx.set_entry(j, i, y[0]);
      gridy.set_entry(j, i, y[1]);
    }
  }
  Grid<2> grid(gridx, gridy);

#if 1
  {
    // plot solution and right-hand side
    cout << "- plotting the exact solution..." << endl;
    SampledMapping<2> u_sm(grid, *uexact);
    std::ofstream u_stream("solution.m");
    u_sm.matlab_output(u_stream);
    u_stream.close();
    cout << "  ...done, see file solution.m!" << endl;
    cout << "- plotting the right-hand side..." << endl;
    SampledMapping<2> f_sm(grid, *f);
    std::ofstream f_stream("rhs.m");
    f_sm.matlab_output(f_stream);
    f_stream.close();
    cout << "  ...done, see file rhs.m!" << endl;
  }
#endif

  typedef GenericFullHelmholtz<Basis> Problem;
  ostringstream G_file, A_file;
  G_file << "ring_" << d << "_" << dt << "_" << s0 << "_" << s1 << "_" << jmax << "_G";
  A_file << "ring_" << d << "_" << dt << "_" << s0 << "_" << s1 << "_"
	 << r0 << "_" << r1 << "_" << jmax << "_A";

  Problem A(basis,
	    G_file.str().c_str(),
	    A_file.str().c_str(),
	    jmax,
	    1.0,
	    energy);

  cout << "- expand right-hand side..." << endl;
  Vector<double> fcoeffs;
  basis.expand(f, true, jmax, fcoeffs);
  cout << "  ... done!" << endl;
  cout << "fcoeffs=" << fcoeffs << endl;

  // setup preconditioned right-hand side
  Vector<double> rhs(fcoeffs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    rhs[i] = fcoeffs[i]/A.D(i);

  // solve Galerkin system
  Vector<double> ulambda(A.row_dimension()), residual(A.row_dimension(), false);
  unsigned int iterations;
  CG(A, rhs, ulambda, 1e-15, 500, iterations);

  //     cout << "  solution coefficients: " << ulambda;
  cout << "  Galerkin system solved with residual (infinity) norm ";
  A.apply(ulambda, residual);
  residual -= rhs;
  cout << linfty_norm(residual) << endl;

#if 1
  {
    InfiniteVector<double,Index> ucoeffs;
    unsigned int i = 0;
    for (Index lambda = basis.first_generator(basis.j0()); i < ulambda.size(); ++lambda, ++i)
      ucoeffs.set_coefficient(lambda, ulambda[i]/A.D(i));

    cout << "- compute L_infty error on a subgrid..." << endl;
    SampledMapping<2> u_sm(grid, *uexact);
    SampledMapping<2> u0_sm(basis.evaluate(ucoeffs, resi));
    cout << "  ... done, result: "
 	 << maximum_norm(u_sm.values()-u0_sm.values())
 	 << endl;
    
    cout << "- plot point values of the solution:" << endl;
    SampledMapping<2> u0_sm2(grid,u0_sm.values());
    std::ofstream u0_stream("u0.m");
    u0_sm2.matlab_output(u0_stream);
    u0_stream.close();
    cout << "  ...done, see file u0.m!" << endl;

    cout << "- plot the pointwise error..." << endl;
    SampledMapping<2> u0_err(grid,u_sm.values()-u0_sm.values());
    std::ofstream u0_err_stream("u0_error.m");
    u0_err.matlab_output(u0_err_stream);
    u0_err_stream.close();
    cout << "  ...done, see file u0_error.m!" << endl;
  }
#endif

  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}

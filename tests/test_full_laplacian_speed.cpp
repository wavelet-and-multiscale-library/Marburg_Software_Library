#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/function.h>
#include <interval/spline_basis.h>
#include <galerkin/full_laplacian.h>
#include <galerkin/full_gramian.h>
#include <Rd/cdf_utils.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main() {
  cout << "Testing speed of FullLaplacian ..." << endl;

  const unsigned int d = 3;
  const unsigned int dT = 3;

  SplineBasis<d,dT,P_construction,1,1,0,0> basis; // PBasis, complementary b.c.'s
  FullLaplacian<d,dT> delta(basis, dyadic);

  const int j0 = basis.j0();
  const int jmax = 10;
  
  clock_t tstart, tend;
  double time1, time2;

  cout << "speed measurements for dyadic preconditioning:" << endl;
  for (int j = j0; j <= jmax; j++) {
    cout << "j=" << j << ":" << endl;

    SparseMatrix<double> A;

    tstart = clock();
    delta.set_level(j);
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;

    cout << "- time for set_level(): " << time1 << "s" << endl;

    tstart = clock();
    delta.to_sparse(A);
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;

    cout << "- time for to_sparse(): " << time2 << "s" << endl;
  }

  FullLaplacian<d,dT> delta2(basis, energy);
  cout << "speed measurements for energy preconditioning:" << endl;
  for (int j = j0; j <= jmax; j++) {
    cout << "j=" << j << ":" << endl;

    SparseMatrix<double> A;

    tstart = clock();
    delta2.set_level(j);
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;

    cout << "- time for set_level(): " << time1 << "s" << endl;

    tstart = clock();
    delta2.to_sparse(A);
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;

    cout << "- time for to_sparse(): " << time2 << "s" << endl;
  }


  
  return 0;
}

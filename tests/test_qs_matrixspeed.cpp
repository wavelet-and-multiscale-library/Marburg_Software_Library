#include <iostream>
#include <map>
#include <cmath>
#include <time.h>
#include <algebra/matrix.h>
#include <algebra/vector.h>
#include <algebra/qs_matrix.h>
#include <algebra/sparse_matrix.h>
#include <utils/random.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Performance measurements for QuasiStationaryMatrix..." << endl;

  int j0 = 3;
  Matrix<double> M_l(5, 1, "0.625 -0.75 -0.25 0.25 0.125");
  Matrix<double> M_r; M_l.mirror(M_r);
  Vector<double> M_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
  QuasiStationaryMatrix<double> M(j0, 15, 8, M_l, M_r, M_band_lr, M_band_lr, 0, 0, M_SQRT1_2);

  cout << "* a QuasiStationaryMatrix M=" << endl << M << endl;
  
  clock_t tstart, tend;
  double time1, time2;

  for (int jtest = j0+4; jtest <= 10; jtest++) {
    M.set_level(jtest);
    Vector<double> x(M.column_dimension()), Mx(M.row_dimension());
    std::map<unsigned int,double> xsparse, Mxsparse;
    
    const unsigned int N = 2000;
    cout << "* applying M on level " << jtest << " " << N << " times to vectors of several density..." << endl;
    for (unsigned int density = 1; density <= 10; density++) {
      cout << "  density: " << density/100.0 << endl;
      
      x = 0;
      xsparse.clear();
      for (unsigned int i = 1; i <= x.size()*(density/100.0); i++) {
	// place entries at random points
	unsigned int k = random_integer(0, x.size()-1);
	x[k] = 1.0;
	xsparse[k] = 1.0;
      }
      
      tstart = clock();
      for (unsigned int n = 1; n <= N; n++)
	M.apply(x, Mx);
      tend = clock();
      time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
      
      tstart = clock();
      for (unsigned int n = 1; n <= N; n++)
	M.apply(xsparse, Mxsparse);
      tend = clock();
      time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;

      cout << "  time needed: " << time1 << "s for Vector<double>, "
	   << time2 << "s for std::map<unsigned int,double>" << endl;
    }
  }

  return 0;
}

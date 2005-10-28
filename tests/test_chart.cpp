#include <iostream>
#include <geometry/chart.h>
#include <geometry/point.h>
#include <algebra/matrix.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing some charts:" << endl;

  AffineLinearMapping<2> kappa_0;
  cout << "* the identity mapping: " << endl << kappa_0;
  Point<2> x(1.0, 2.0), y;
  kappa_0.map_point(x, y);
  cout << "  x=" << x << " is mapped to y=" << y << endl;
  kappa_0.map_point_inv(y, x);
  cout << "  y=" << y << " is mapped back to x=" << x << endl;

  Matrix<double> A(2, 2, "2 0 1 -4");
  Point<2> b(0, 1);
  AffineLinearMapping<2> kappa_1(A, b);
  cout << "* affine linear map: " << endl << kappa_1;
  kappa_1.map_point(x, y);
  cout << "  x=" << x << " is mapped to y=" << y << endl;
  kappa_1.map_point_inv(y, x);
  cout << "  y=" << y << " is mapped back to x=" << x << endl;
  
  Point<2> P(1.0, 2.0);
  kappa_1.map_point_inv(P, y);
  cout << "  test point: P=" << P
       << " (inverse: " << y << ")"
       << ", in patch: " << kappa_1.in_patch(P) << endl;

  // [0,1] -> [-1,1]
  Matrix<double> C(1, 1, "2");
  Point<1> d(-1.0);
  AffineLinearMapping<1> kappa_2(C, d);
  cout << "* affine linear map with C=" << endl << C
       << "  and d=" << d << ":" << endl;
  Point<1> x2(0.5), y2;
  kappa_2.map_point(x2, y2);
  cout << "  x=" << x2 << " is mapped to y=" << y2 << endl;
  kappa_2.map_point_inv(y2, x2);
  cout << "  y=" << y2 << " is mapped back to x=" << x2 << endl;

  cout << kappa_2 << endl;

  //############## tests Manuel ###############
  cout << "Testing class LinearBezierMapping..." << endl;
  cout << endl;
  LinearBezierMapping k_1(Point<2>(-0.2,-1),Point<2>(-0.5,0),
  			  Point<2>(1,-1.1), Point<2>(1.5,0));
  //LinearBezierMapping k_1(Point<2>(-0.2,-1),Point<2>(-0.2,0),
  //Point<2>(1,-1), Point<2>(1,0));


  clock_t tstart, tend;
  double time;
  
  const int J = 3;
  const double dx = 1.0 / (1<<J);
  Point<2> tmp;
  tstart = clock();
  Point<2> pc(0.5, 0.75);
  Matrix<double> R(2,2);
  R(0,0) = 1.;
  R(1,1) = 1.;
  for (int i = 0; i <= 1<<J; i++)
    {
      for (int j= 0; j<= 1<<J; j++)
	{ 
	  Point<2> pc(i*dx, j*dx);
	  //cout << k_1 << endl;
	  //cout << "###########################" << endl;
	  cout << "before = " << pc << endl;
	  k_1.map_point(pc,tmp);
	  cout << "in the meantime = " << tmp << endl;
	  k_1.map_point_inv(tmp,pc);
	  cout << "after = " << pc << endl;
	  cout << endl;
	  
	  // 	   Point<2> pc(i*dx, j*dx);
	  // 	   cout << k_2 << endl;
	  // 	   cout << "###########################" << endl;
	  // 	   cout << "before = " << pc << endl;
	  // 	   k_2.mapPoint(tmp,pc);
	  // 	   cout << "in the meantime = " << tmp << endl;
	  // 	   k_2.mapPointInv(pc, tmp);
	  // 	   cout << "after = " << pc << endl;
	  
	  // 	   if (j != 0)
	  // 	     continue;
	  // 	   Point<1> pc(i*dx);
	  // 	   cout << k_3 << endl;
	  // 	   cout << "###########################" << endl;
	  // 	   cout << "before = " << pc << endl;
	  // 	   k_3.mapPoint(tmp_1,pc);
	  // 	   cout << "in the meantime = " << tmp_1 << endl;
	  // 	   k_3.mapPointInv(pc, tmp_1);
	  //	   cout << "after = " << pc << endl;
	  
	}
    }
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "cpu time = " << time << "s" << endl;
  cout << "Gram_factor = " << k_1.Gram_factor(pc) << endl;
  cout << "Gram_D_factor = " << k_1.Gram_D_factor(1,pc) << endl;
  cout << "Dkappa_inv " << k_1.Dkappa_inv(1,0,pc) << endl;
  Point<2> p(1, 0);
  cout << "in_patch " << k_1.in_patch(p) << endl;

  cout << k_1 << endl;

  return 0;
}

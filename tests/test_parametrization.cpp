
#include <iostream>
#include "parametrization.h"

using std::cout;
using std::endl;

using FrameTL::LinearBezierMapping;
using FrameTL::AffinLinearMapping;

int main()
{
  
  cout << "Hello World, Frame software is alive again :-)!" << endl;
  LinearBezierMapping k_1(Point<2>(-0.5,-1),Point<2>(-0.5,0),
			Point<2>(1,-1), Point<2>(1,0));

  Matrix<double> A(2,2);
  A(0,0) = 2.3;
  A(1,1) = -0.21;
  Vector<double> b(2);
  b(0) = -1.9;
  b(1) = 5.78;
  AffinLinearMapping<2> k_2 (A,b);

  Matrix<double> AA(1,1);
  AA(0,0) = 2.3;
  Vector<double> bb(1);
  bb(0) = -1.9;
  AffinLinearMapping<1> k_3 (AA,bb);


  const int J = 3;
  const double dx = 1.0 / (1<<J);

   Point<2> tmp;
   Point<1> tmp_1;
//   Point<2> pc(0,0.5);
//   cout << "before = " << pc << endl;
//   b.mapPoint(tmp,pc);
//   cout << "in the meantime = " << tmp << endl;
//   b.mapPointInv(pc, tmp);
//   cout << "after = " << pc << endl;

   for (int i = 0; i <= 1<<J; i++)
     {
       for (int j= 0; j<= 1<<J; j++)
	 {
// 	   Point<2> pc(i*dx, j*dx);
// 	   cout << k_1 << endl;
// 	   cout << "###########################" << endl;
// 	   cout << "before = " << pc << endl;
// 	   k_1.mapPoint(tmp,pc);
// 	   cout << "in the meantime = " << tmp << endl;
// 	   k_1.mapPointInv(pc, tmp);
// 	   cout << "after = " << pc << endl;
	   
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
// 	   cout << "after = " << pc << endl;

	 }
     } 
   
  return 0;
}


#include <iostream>
#include "parametrization.h"
#include "atlas.h"
#include <misc.cpp>
#include <time.h> 

using std::cout;
using std::endl;


using FrameTL::LinearBezierMapping;
using FrameTL::AffineLinearMapping;
using FrameTL::Atlas;

using namespace FrameTL;

int main()
{
  
  cout << "Testing class Atlas..." << endl;
  //let us put together an  L-shaped domain by using both, LinearBezierMapping
  //and AffineLinearMapping
  Atlas<2,2> Lshaped;

  //##############################  
  Matrix<double> A(2,2);
  A(0,0) = 2.0;
  A(1,1) = 1.0;
  Vector<double> b(2);
  b[0] = -1.0;
  b[1] = -1.0;
  AffineLinearMapping<2> affineP(A,b);
  //##############################
  
  //##############################
  LinearBezierMapping bezierP(Point<2>(-1.,-1),Point<2>(-1.,1),
			      Point<2>(1.,-1.), Point<2>(1,0));
  //##############################

  Lshaped.add_chart(affineP);
  Lshaped.add_chart(bezierP);
  cout << Lshaped << endl;

  return 0;
}


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

using MathTL::SymmetricMatrix;
using MathTL::Array1D;



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

  //test copy constructor
  Atlas<2,2> Lshaped_copy(Lshaped);
  cout << Lshaped_copy << endl;

  //test constructor (does not really make sense)
  Atlas<3,3> dummy(4);
  //cout << dummy << endl;

  Array1D< Parametrization<2,2>* > ch(2);
  ch[0] = &bezierP;
  ch[1] = &affineP;
  
  Atlas<2,2> Lshaped2(ch);

  cout << Lshaped2 << endl;

  Lshaped.set_chart(1, bezierP);
  cout << Lshaped << endl;

  SymmetricMatrix<bool> adjaM(2);
  adjaM(0,0) = true;
  adjaM(1,1) = true;

  Lshaped.set_adjacency_matrix(adjaM);
  Lshaped.set_adjacent(1,2);
  
  cout << Lshaped << endl;

  Atlas<2,2> L3(3);
  L3.set_chart(1,affineP);
  L3.set_chart(2,affineP);
  L3.set_chart(3,bezierP);
  L3.set_adjacent(1,2);
  L3.set_adjacent(1,3);
  L3.set_adjacent(3,2);
  cout << L3 << endl;

  Point<2> p(10.,.0);
  Array1D<unsigned int> res(L3.get_circumjacent(p));
  cout << res << endl;

  cout << L3.adjacent(1,3) << endl;


  return 0;
}

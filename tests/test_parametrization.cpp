
#include <iostream>
#include "parametrization.h"

using std::cout;
using std::endl;

using FrameTL::LinearBezierMapping;

int main()
{
  
  cout << "Hello World, Frame software is alive again :-)!" << endl;
  LinearBezierMapping b(Point<2>(-0.5,-1),Point<2>(-0.5,0),
			Point<2>(1,-1), Point<2>(1,0));
  const int J = 3;
  const double dx = 1.0 / (1<<J);

   Point<2> tmp;
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
	    Point<2> pc(i*dx, j*dx);
	    cout << "###########################" << endl;
	    cout << "before = " << pc << endl;
	    b.mapPoint(tmp,pc);
	    cout << "in the meantime = " << tmp << endl;
	    b.mapPointInv(pc, tmp);
	    cout << "after = " << pc << endl;
	  }
    } 

  return 0;
}

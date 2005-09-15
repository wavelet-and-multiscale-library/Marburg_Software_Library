
#include <iostream>
#include "parametrization.h"
#include "atlas.h"
#include <misc.cpp>
#include <time.h> 

using std::cout;
using std::endl;


using FrameTL::LinearBezierMapping;
using FrameTL::AffinLinearMapping;
using FrameTL::Atlas;

using namespace FrameTL;

int main()
{
  
  cout << "Testing class Atlas..." << endl;
  Atlas<2,2> omega;
  return 0;
}

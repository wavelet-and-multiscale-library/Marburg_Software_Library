#include <fstream>
#include <iostream>
#include <biharmonic_equation.h>
#include <biharmonic_1d_testcase.h>
#include <biharmonic_1d_testcase.h>
#include <geometry/sampled_mapping.h>

using std::cout;
using std::endl;



using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing class Biharmonic1D_Solution..." << endl;
  
  const int resolution = 1 << 9;
  
  // We plot the exact solution.

  // grid resolution should be a multiple of three to get the re-entrant corners into the grid
  Grid<1> grid1D(0, 1, resolution);
  Biharmonic1D_Solution solution;
  
  SampledMapping<1> mapping_exact (grid1D, solution);
  std::ofstream ofexactsol("exact_solution_biharm_interval.m");
  mapping_exact.matlab_output(ofexactsol);
  ofexactsol.close();

   return 0;

}

#include <iostream>
#include <fstream>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/corner_singularity_biharmonic.h>
//#include <numerics/atez.h>

using namespace MathTL;
using namespace std;

int main()
{
  cout << "Testing CornerSingularity ..." << endl;

  Point<2> origin;

#if 1
   CornerSingularityBiharmonic s(origin, 0.5, 1.5);

  const bool matlab = true; // false -> Octave output

  if (matlab)
    cout << "- Matlab output of the corner singularity on a 2D grid..." << endl;
  else
    cout << "- Octave output of the corner singularity on a 2D grid..." << endl;
  
  Grid<2> grid(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), 1<<6, 1<<6);
  //   Grid<2> grid(Point<2>(-1.0, 0.0), Point<2>(0.0, 1.0), 1<<6, 1<<6); // only patch 0
  // Grid<2> grid(Point<2>(-1.0, -1.0), Point<2>(0.0, 1.0), 1<<6, 1<<6); // only patches 0+1
  //   Grid<1> grid(0.0,1.0,1<<10);
  SampledMapping<2> h(grid, s);
  std::ofstream fs("corner_biharmonic.m");

  if (matlab)
    h.matlab_output(fs);
  else
    h.octave_output(fs);

  if (matlab)
    fs << "surf(x,y,z)" << endl;
  else
    fs << "mesh(x,y,z)" << endl;

  fs.close();
  cout << "  ...done, see file corner_biharmonic.m!" << endl;

  CornerSingularityBiharmonicRHS rhs(origin, 0.5, 1.5);
  //CornerSingularityBiharmonic rhs(origin, 0.5, 1.5);
  if (matlab)
    cout << "- Matlab output of the corner singularity rhs on a 2D grid..." << endl;
  else
    cout << "- Octave output of the corner singularity rhs on a 2D grid..." << endl;
  
  SampledMapping<2> h_rhs(grid, rhs);
  std::ofstream fs_rhs("corner_biharmonic_rhs.m");
  if (matlab)
    h_rhs.matlab_output(fs_rhs);
  else
    h_rhs.octave_output(fs_rhs);
  
  if (matlab)
    fs_rhs << "surf(x,y,z)" << endl;
  else
    fs_rhs << "mesh(x,y,z)" << endl;
  
  fs_rhs.close();
  cout << "  ...done, see file corner_biharmonic_rhs.m!" << endl;
#endif
  
    //  Atez zeta;//(double r0 0.01 , double r1 = 0.99);
//      SampledMapping<1> h(grid, zeta);
//      std::ofstream fs_rhs("zeta.m");
//      h.matlab_output(fs_rhs);
//      fs_rhs.close();

  return 0;
}



  

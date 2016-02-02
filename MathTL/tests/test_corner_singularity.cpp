#include <iostream>
#include <fstream>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/corner_singularity.h>

using namespace MathTL;
using namespace std;

int main()
{
  cout << "Testing CornerSingularity ..." << endl;

  Point<2> origin;

#if 1
  CornerSingularity s(origin, 0.5, 1.5);

  const bool matlab = true; // false -> Octave output

  if (matlab)
    cout << "- Matlab output of the corner singularity on a 2D grid..." << endl;
  else
    cout << "- Octave output of the corner singularity on a 2D grid..." << endl;
  
   Grid<2> grid(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), 1<<6, 1<<6);
//  Grid<2> grid(Point<2>(-1.0, 0.0), Point<2>(0.0, 1.0), 1<<6, 1<<6); // only patch 0
//   Grid<2> grid(Point<2>(-1.0, -1.0), Point<2>(0.0, 1.0), 1<<6, 1<<6); // only patches 0+1
  SampledMapping<2> h(grid, s);
  std::ofstream fs("corner.m");

  if (matlab)
    h.matlab_output(fs);
  else
    h.octave_output(fs);

  if (matlab)
    fs << "surf(x,y,z)" << endl;
  else
    fs << "mesh(x,y,z)" << endl;

  fs.close();
  cout << "  ...done, see file corner.m!" << endl;

  CornerSingularityRHS rhs(origin, 0.5, 1.5);
  if (matlab)
    cout << "- Matlab output of the corner singularity rhs on a 2D grid..." << endl;
  else
    cout << "- Octave output of the corner singularity rhs on a 2D grid..." << endl;

  SampledMapping<2> h_rhs(grid, rhs);
  std::ofstream fs_rhs("corner_rhs.m");
  if (matlab)
    h_rhs.matlab_output(fs_rhs);
  else
    h_rhs.octave_output(fs_rhs);

  if (matlab)
    fs_rhs << "surf(x,y,z)" << endl;
  else
    fs_rhs << "mesh(x,y,z)" << endl;

  fs_rhs.close();
  cout << "  ...done, see file corner_rhs.m!" << endl;


  // we plot the first component of the gradient of the corner singularity
  CornerSingularityGradient corner_grad(origin, 0.5, 1.5);
  if (matlab)
    cout << "- Matlab output of the corner singularity rhs on a 2D grid..." << endl;
  else
    cout << "- Octave output of the corner singularity rhs on a 2D grid..." << endl;

  SampledMapping<2> h_rhs2(grid, corner_grad);
  std::ofstream fs_rhs2("corner_gradient.m");
  if (matlab)
    h_rhs2.matlab_output(fs_rhs2);
  else
    h_rhs2.octave_output(fs_rhs2);

  if (matlab)
    fs_rhs2 << "surf(x,y,z)" << endl;
  else
    fs_rhs2 << "mesh(x,y,z)" << endl;
  
  fs_rhs2.close();
  cout << "  ...done, see file corner_gradient.m!" << endl;


#endif

#if 0
  CornerTimeSingularity st(origin, 0.5, 1.5);
  st.set_time(1.0);

  cout << "- Matlab output of the time corner singularity on a 2D grid..." << endl;
  Grid<2> gridt(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), 1<<6, 1<<6);
  SampledMapping<2> ht(gridt, st);
  std::ofstream fst("corner_time.m");
  ht.matlab_output(fst);
  fst << "surf(x,y,z)" << endl;
  fst.close();
  cout << "  ...done, see file corner_time.m!" << endl;

  CornerTimeSingularityRHS st_rhs(origin, 0.5, 1.5);
  st_rhs.set_time(1.0);

  cout << "- Matlab output of the time corner singularity rhs on a 2D grid..." << endl;
  Grid<2> gridt_rhs(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), 1<<6, 1<<6);
  SampledMapping<2> ht_rhs(gridt_rhs, st_rhs);
  std::ofstream fst_rhs("corner_time_rhs.m");
  ht_rhs.matlab_output(fst_rhs);
  fst_rhs << "surf(x,y,z)" << endl;
  fst_rhs.close();
  cout << "  ...done, see file corner_time_rhs.m!" << endl;
#endif

  return 0;
}

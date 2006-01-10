#include <iostream>
#include <fstream>
#include <utils/plot_tools.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::get_color() ..." << endl;

  MatlabColorMap colormap = jet;

  double red, green, blue;
  for (double x = -1.0; x <= 1.0; x += 0.05) {
    get_color(x, colormap, red, green, blue);
    cout << "x=" << x << " yields RGB values "
	 << "[" << red << "," << green << "," << blue << "]" << endl;
  }

#if 1
  // create some colored rectangles (to visually check the colormap routine)
//   colormap = gray;
  colormap = jet;
  std::ofstream fs("colored_rectangles.m");
  int N = 10;
  double h = 2.0/N;
  for (int k = 0; k < N; k++) {
    double x = (k-N/2)*h;
    double y = 0.0;
    double width = h;
    double height = 1.0;
    get_color(x, colormap, red, green, blue);
    fs << "rectangle('position',["
       << x << "," << y<< "," << width << "," << height << "],"
       << "'LineWidth',0.125,"
       << "'FaceColor',[" << red << "," << green << "," << blue << "])" << endl;
  }
  fs.close();
#endif

  return 0;
}

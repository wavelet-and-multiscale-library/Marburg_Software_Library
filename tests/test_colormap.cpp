#include <iostream>
#include <utils/plot_tools.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::get_color() ..." << endl;

  MatlabColorMap colormap = gray;
  double red, green, blue;
  for (double x = -1.0; x <= 1.0; x += 0.2) {
    get_color(x, colormap, red, green, blue);
    cout << "x=" << x << " yields RGB values "
	 << "[" << red << "," << green << "," << blue << "]" << endl;
  }

  return 0;
}

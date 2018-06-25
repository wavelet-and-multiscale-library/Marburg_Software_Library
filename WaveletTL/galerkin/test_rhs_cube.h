
#include <utils/function.h>


/*
  Some test problems for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: u(x,y) = x(1-x)y(1-y), -Delta u(x,y) = 2(x(1-x)+y(1-y))
  2: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y)
  3: u(x,y) = x(1-x)^2y^2(1-y), -Delta u(x,y)= 4*(1-x)*y^2*(1-y)-2*x*y^2*(1-y)-2*x*(1-x)^2*(1-y)+4*x*(1-x)^2*y
*/
template <unsigned int N>
class TestRHS
  : public MathTL::Function<2,double>
{
public:
  virtual ~TestRHS() {}
  double value(const MathTL::Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
       return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
       break;     
    case 2:
      return
	(200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
	* exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
      break;
    case 3:
      return
	4*(1-p[0])*p[1]*p[1]*(1-p[1])
	- 2*p[0]*p[1]*p[1]*(1-p[1])
	- 2*p[0]*(1-p[0])*(1-p[0])*(1-p[1])
	+ 4*p[0]*(1-p[0])*(1-p[0])*p[1];
      break;
    default:
        return 1;
      break;
    }
  }
  void vector_value(const MathTL::Point<2>& p, MathTL::Vector<double>& values) const {
    values[0] = value(p);
  }
};

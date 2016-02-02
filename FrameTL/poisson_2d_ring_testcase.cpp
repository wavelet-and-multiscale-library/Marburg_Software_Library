// implementation for poisson_2d_ring_testcase.h

#include <cmath>

using namespace std;


namespace FrameTL
{

  double
  Poisson_Solution_Ring::value(const Point<2>& p, const unsigned int component) const
  {
    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    CornerSingularity sing0(origin, 0.5, 1.5);

    origin[0] = 1.0;
    origin[1] = 0.0;
    CornerSingularity sing1(origin, 1.0, 1.5);
    
    origin[0] = 1.0;
    origin[1] = 1.0;
    CornerSingularity sing2(origin, 1.5, 1.5);

    origin[0] = 0.0;
    origin[1] = 1.0;
    CornerSingularity sing3(origin, 0.0, 1.5);

    return sing0.value(p)+sing1.value(p)+sing2.value(p)+sing3.value(p);
    
    
  }

  double
  Poisson_RHS_Ring::value(const Point<2>& p, const unsigned int component) const
  {
    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    CornerSingularityRHS singRhs0(origin, 0.5, 1.5);

    origin[0] = 1.0;
    origin[1] = 0.0;
    CornerSingularityRHS singRhs1(origin, 1.0, 1.5);
    
    origin[0] = 1.0;
    origin[1] = 1.0;
    CornerSingularityRHS singRhs2(origin, 1.5, 1.5);

    origin[0] = 0.0;
    origin[1] = 1.0;
    CornerSingularityRHS singRhs3(origin, 0.0, 1.5);

    return singRhs0.value(p)+singRhs1.value(p)+singRhs2.value(p)+singRhs3.value(p);
  }
  
  void
  Poisson_SolutionGradient_Ring::vector_value(const Point<2> &p, Vector<double>& values) const
  {
    values[0] = 0.0;
    values[1] = 0.0;

    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    CornerSingularityGradient singGrad0(origin, 0.5, 1.5);

    origin[0] = 1.0;
    origin[1] = 0.0;
    CornerSingularityGradient singGrad1(origin, 1.0, 1.5);
    
    origin[0] = 1.0;
    origin[1] = 1.0;
    CornerSingularityGradient singGrad2(origin, 1.5, 1.5);

    origin[0] = 0.0;
    origin[1] = 1.0;
    CornerSingularityGradient singGrad3(origin, 0.0, 1.5);

    Vector<double> tmp_values(2);
    singGrad0.vector_value(p, tmp_values);
    values+=tmp_values;
    
    singGrad1.vector_value(p, tmp_values);
    values+=tmp_values;

    singGrad2.vector_value(p, tmp_values);
    values+=tmp_values;

    singGrad3.vector_value(p, tmp_values);
    values+=tmp_values; 
  }


}

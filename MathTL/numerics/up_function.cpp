// implementation for up_function.h

#include <iostream>
#include <algorithm>
#include <numerics/cardinal_splines.h>

namespace MathTL
{
  ApproximateUpFunction::ApproximateUpFunction(const unsigned int k)
    : k_(k), coeffs((1<<(k+1))-k)
  {
    // initialize spline expansion coefficients w.r.t. N_{k+1}(2^k.-m)
    coeffs[0] = coeffs[1] = 0.5;
    for (unsigned int i = 2; i < coeffs.size(); i++) coeffs[i] = 0;
    Array1D<double> help(coeffs.size());
    for (int j = 1; j <= (int)k; j++) {
      std::copy(coeffs.begin(), coeffs.end(), help.begin());
      
      for (unsigned int i = 0; i < coeffs.size(); i++)
	coeffs[i] = 0;

      for (int n = -(1<<j); n <= -1; n++) {
	for (int m = -(1<<(j-1));
	     m <= std::min(n+(1<<(j-1)),(1<<(j-1))-j); m++)
	  coeffs[(1<<j)+n] += help[(1<<(j-1))+m];
	coeffs[(1<<j)+n] /= (1<<(j-1));
      }

      for (int n = 0; n <= (1<<j)-j-1; n++) {
	for (int m = std::min(n-(1<<(j-1)),(1<<(j-1))-j)+1;
	     m <= std::min(n+(1<<(j-1)),(1<<(j-1))-j); m++)
	  coeffs[(1<<j)+n] += help[(1<<(j-1))+m];
	coeffs[(1<<j)+n] /= (1<<(j-1));
      }
      
    }
    
//     std::cout << "spline coefficients of u_" << k << ": " << coeffs << std::endl;
  }

  ApproximateUpFunction::~ApproximateUpFunction()
  {
  }

  double
  ApproximateUpFunction::value(const Point<1>& p,
			       const unsigned int component) const
  {
    double r = 0;
    for (unsigned int i = 0; i < coeffs.size(); i++) {
      r += coeffs[i]
	* EvaluateCardinalBSpline(k_+1, (-(1<<k_))+i, (1<<k_)*p[0]);
    }
    
    return r;
  }

  inline
  void
  ApproximateUpFunction::vector_value(const Point<1> &p,
				      Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
  
}

// implementation of quadrature.h inline functions

#include <algebra/polynomial.h>

namespace MathTL
{
  template <unsigned int DIM>
  double QuadratureRule<DIM>::integrate(const Function<DIM, double>& f) const
  {
    double r(0.0);
    
    for (unsigned int k(0); k < points_.size(); k++)
      r += weights_[k] * f.value(points_[k]);
    
    return r;
  }

  MidpointRule::MidpointRule()
  {
    points_.resize(1);
    points_[0] = 0.5;
    weights_.resize(1);
    weights_[0] = 1;
  }

  TrapezoidalRule::TrapezoidalRule()
  {
    points_.resize(2);
    points_[0] = 0.0;
    points_[1] = 1.0;
    weights_.resize(2);
    weights_[0] = weights_[1] = 0.5;
  }

  SimpsonRule::SimpsonRule()
  {
    points_.resize(3);
    points_[0] = 0.0;
    points_[1] = 0.5;
    points_[2] = 1.0;
    weights_.resize(3);
    weights_[0] = weights_[2] = 1.0/6.0;
    weights_[1] = 2.0/3.0;
  }

  template <unsigned int N>
  ClosedNewtonCotesRule<N>::ClosedNewtonCotesRule()
  {
    points_.resize(N+1);
    weights_.resize(N+1);
    for (unsigned int n(0); n <= N; n++)
      {
	points_[n]   = ((double) n)/ N;
	
	// construct n-th Lagrange polynomial
	Polynomial<double> p(1.0), q;
	q.set_coefficient(1, 1.0);
	for (int j(0); j <= (int) N; j++)
	  {
	    if (j != (int) n)
	      {
		q.set_coefficient(0, -j); // q(x)=x-j
		p *= 1./(((int) n)-j) * q;
	      }
	  }
	weights_[n] = p.integrate(0, N)/(double) N;
      }
  }
}

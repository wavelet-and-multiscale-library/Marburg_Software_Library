// implementation of quadrature.h inline functions

#include <algebra/polynomial.h>

namespace MathTL
{
  template <unsigned int DIM>
  QuadratureRule<DIM>::QuadratureRule()
    : points_(), weights_()
  {
  }

  template <unsigned int DIM>
  QuadratureRule<DIM>::QuadratureRule(const SubQuadratureRule& Q,
				      const QuadratureRule<1>& Q1)
    : points_(Q.get_N() * Q1.get_N()),
      weights_(Q.get_N() * Q1.get_N())
  {
    Array1D<Point<DIM-1> > Qpoints;
    Q.get_points(Qpoints);
    Array1D<Point<1> > Q1points;
    Q1.get_points(Q1points);
    Array1D<double> Qweights, Q1weights;
    Q.get_weights(Qweights);
    Q1.get_weights(Q1weights);

    // compute tensor product
    unsigned int current_index(0);
    for (unsigned int i(0); i < Q.get_N(); i++)
      for (unsigned int j(0); j < Q1.get_N(); j++, current_index++)
	{
	  for (unsigned int d(0); d < DIM-1; d++)
	    points_[current_index](d) = Qpoints[i](d);
	  points_[current_index](DIM-1) = Q1points[j](0);
	  weights_[current_index] = Qweights[i] * Q1weights[j];
	}
  }

  template <unsigned int DIM>
  inline
  unsigned int QuadratureRule<DIM>::get_N() const
  {
    return points_.size();
  }

  template <unsigned int DIM>
  inline
  void QuadratureRule<DIM>::get_points(Array1D<Point<DIM> >& points) const
  {
    points = points_;
  }

  template <unsigned int DIM>
  inline
  void QuadratureRule<DIM>::get_weights(Array1D<double>& weights) const
  {
    weights = weights_;
  }

  template <unsigned int DIM>
  double QuadratureRule<DIM>::integrate(const Function<DIM, double>& f) const
  {
    double r(0.0);
    
    for (unsigned int k(0); k < points_.size(); k++)
      r += weights_[k] * f.value(points_[k]);
    
    return r;
  }

  template <unsigned int DIM>
  double QuadratureRule<DIM>::integrate(const Function<DIM, double>& f,
					const Point<DIM>& a, const Point<DIM>& b) const
  {
    double r(0.0);
    
    Point<DIM> helppoint;
    double helpweight;
    for (unsigned int k(0); k < points_.size(); k++)
      {
	helpweight = weights_[k];
	for (unsigned int d(0); d < DIM; d++)
	  {
	    helppoint(d) = a(d) + points_[k](d)*(b(d)-a(d));
	    helpweight *= (b(d)-a(d));
	  }
	r += helpweight * f.value(helppoint);
      }
    
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

  template <unsigned int DIM, class QUADRATURE>
  CompositeRule<DIM, QUADRATURE>::CompositeRule(const unsigned int N)
    : Q_(), N_(N)
  {
//     points_.resize((N+1)*DIM);
//     cout << "after points.resize" << endl;
//     weights_.resize((N+1)*DIM);
  }

//   template <class Q>
//   CompositeRule<Q>::CompositeRule(const double a, const double b, const int N)
//     : a_(a), b_(b), N_(N) {}

//   template <class Q>
//   template <class FUNCTION>
//   double CompositeRule<Q>::integrate(FUNCTION f) const
//   {
//     double r(0);
//     for (int k(0); k < N_; k++)
//       {
// 	// TODO: use only one internal instance of Q
// 	// -> maybe we can use templates afterwards!
// 	r += Q(a_+k*(b_-a_)/N_, a_+(k+1)*(b_-a_)/N_).integrate(f);
//       }

//     return r;
//   }


}

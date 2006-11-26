// implementation of quadrature.h inline functions

#include <algebra/polynomial.h>

namespace MathTL
{
  template <unsigned int DIM>
  QuadratureRule<DIM>::QuadratureRule()
    : points_(), weights_()
  {
  }

  template <>
  QuadratureRule<0>::QuadratureRule()
  {
    // empty quadrature rule
  }

  template <unsigned int DIM>
  QuadratureRule<DIM>::QuadratureRule(const QuadratureRule<DIM>& Q)
    : points_(Q.points_), weights_(Q.weights_)
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
  template <>
  QuadratureRule<1>::QuadratureRule(const QuadratureRule<0>& Q,
				    const QuadratureRule<1>& Q1)
  {
    assert(false); // this constructor should never be called
  }

  template <unsigned int DIM>
  QuadratureRule<DIM>::~QuadratureRule()
  {
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
//     double helpweight;
    for (unsigned int k(0); k < points_.size(); k++)
      {
	double helpweight = weights_[k];
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

  ClosedNewtonCotesRule::ClosedNewtonCotesRule(const unsigned int N)
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

  template <unsigned int DIM>
  CompositeRule<DIM>::CompositeRule(const QuadratureRule<1>& Q,
				    const unsigned int N)
    : QuadratureRule<DIM>(CompositeRule<DIM-1>(Q, N),
			  CompositeRule<1>(Q, N)),
      Q_(Q), N_(N)
  {
  }

  template <>
  bool CompositeRule<1>::uses_both_endpoints(const QuadratureRule<1>& Q)
  {
    bool left(false), right(false);

    Array1D<Point<1> > points;
    Q.get_points(points);
    Array1D<double> weights;
    Q.get_weights(weights);

    for (unsigned int i(0); i < Q.get_N(); i++)
      {
	if (points[i][0] == 0.0) left = true;
	if (points[i][0] == 1.0) right = true;
      }

    return (left && right);
  }  

  template <>
  CompositeRule<1>::CompositeRule(const QuadratureRule<1>& Q,
				  const unsigned int N)
    : Q_(Q), N_(N)
  {
    Array1D<Point<1> > points;
    Q.get_points(points);
    Array1D<double> weights;
    Q.get_weights(weights);

    // check whether we have to take care of double points
    if (uses_both_endpoints(Q))
      {
	points_.resize(N*(Q.get_N()-1)+1);
	weights_.resize(N*(Q.get_N()-1)+1);

	// the double points have to be glued,
	// which changes their weight
	double double_point_weight(0.0);
	unsigned int n_end_points(0);
	for (unsigned int i(0); i < Q.get_N(); i++)
	  if (points[i][0] == 0.0 ||
	      points[i][0] == 1.0)
	    {
	      double_point_weight += weights[i];
	      n_end_points++;
	    }

	double_point_weight /= N;

	// the base quadrature formula should not have doubly used points:
	assert(n_end_points == 2);

	unsigned int next_point(0);
	for (unsigned int copy(0); copy < N; copy++)
	  for (unsigned int q_point(0); q_point < Q.get_N(); q_point++)
	    {
	      if (copy > 0 && points[q_point][0] == 0.0)
		continue;
	      
	      points_[next_point]
		= Point<1>(points[q_point](0)/N + (1.0*copy)/N);

	      if (copy != N-1 && points[q_point][0] == 1.0)
		weights_[next_point] = double_point_weight;
	      else
		weights_[next_point] = weights[q_point]/N;

	      next_point++;
	    }
      }
    else
      {
	points_.resize(N*Q.get_N());
	weights_.resize(N*Q.get_N());

	// we can just copy the points appropriately
	unsigned int next_point(0);
	for (unsigned int copy(0); copy < N; copy++)
	  for (unsigned int q_point(0); q_point < Q.get_N(); q_point++)
	    {
	      points_[next_point]
		= Point<1>(points[q_point](0)/N + (1.0*copy)/N);
	      weights_[next_point]
		= weights[q_point]/N;
	      
	      next_point++;
	    }
      }
  }
}

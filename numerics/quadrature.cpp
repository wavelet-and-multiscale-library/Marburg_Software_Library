// implementation of quadrature.h inline functions

#include <algebra/polynomial.h>

namespace MathTL
{
//   template <unsigned int DIM>
//   QuadratureRule<DIM>::QuadratureRule(const QuadratureRule<DIM-1>& Q,
// 				      const QuadratureRule<1>& Q1)
//     : 
//   {
//   }

// template <int dim>
// Quadrature<dim>::Quadrature (const SubQuadrature &q1,
// 			     const Quadrature<1> &q2)
// 		:
// 		n_quadrature_points (q1.n_quadrature_points *
// 				     q2.n_quadrature_points),
// 		quadrature_points (n_quadrature_points),
// 		weights (n_quadrature_points, 0)
// {
//   unsigned int present_index = 0;
//   for (unsigned int i=0; i<q1.n_quadrature_points; ++i)
//     for (unsigned int j=0; j<q2.n_quadrature_points; ++j)
//       {
// 					 // compose coordinates of
// 					 // new quadrature point by tensor
// 					 // product in the last component
// 	for (unsigned int d=0; d<dim-1; ++d)
// 	  quadrature_points[present_index](d)
// 	    = q1.point(i)(d);
// 	quadrature_points[present_index](dim-1)
// 	  = q2.point(j)(0);
					       
// 	weights[present_index] = q1.weight(i) * q2.weight(j);

// 	++present_index;
//       };

// #ifdef DEBUG
//   double sum = 0;
//   for (unsigned int i=0; i<n_quadrature_points; ++i)
//     sum += weights[i];
// 				   // we cant guarantee the sum of weights
// 				   // to be exactly one, but it should be
// 				   // near that. 
//   Assert ((sum>0.999999) && (sum<1.000001), ExcInternalError());
// #endif
// }


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
    cout << "integrate() called, with points " << endl;
    print_vector(points_, cout);
    cout << " and weights ";
    print_vector(weights_, cout);
    cout << endl;

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

    cout << "constructed closed NC rule, points: " << endl;
    print_vector(points_, cout);
    cout << ", weights: ";
    print_vector(weights_, cout);
    cout << endl;
  }

  template <unsigned int DIM>
  CompositeRule<DIM>::CompositeRule(const QuadratureRule<1>& Q,
				    const unsigned int N)
    : Q_(Q), N_(N)
  {
//     cout << "before points.resize" << endl;
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

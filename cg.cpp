// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>

using std::set;

namespace FrameTL{

/*!
*/
template<class VALUE = double>
class Singularity1D_RHS_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_RHS_2() {};
  virtual ~Singularity1D_RHS_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};

/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_2() {};
  virtual ~Singularity1D_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    if (0. <= p[0] && p[0] < 0.5)
      return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];

    if (0.5 <= p[0] && p[0] <= 1.0)
      return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);

    return 0.;

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};



  template <class PROBLEM>
  void cg_SOLVE(const PROBLEM& P,  const double epsilon,
		InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {
    //typedef DSBasis<2,2> Basis1D;
    typedef PBasis<2,2> Basis1D;

    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    
    CornerSingularity sing2D(origin, 0.5, 1.5);
    CornerSingularityRHS singRhs(origin, 0.5, 1.5);

//    Singularity1D_RHS_2<double> sing1D;
//    Singularity1D_2<double> exactSolution1D;

    unsigned int loops = 0;
    const int jmax = 4;
    typedef typename PROBLEM::Index Index;

    double a_inv     = P.norm_Ainv();
    double kappa     = P.norm_A()*a_inv;
    double omega_i   = a_inv*P.F_norm();
    cout << "a_inv = " << a_inv << endl;
    cout << "omega_i = " << omega_i << endl;

    InfiniteVector<double, Index> xk, rk, zk, pk, Apk, help, f;

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;

    double eta = 0.1;

    const double tol = 1.0e-6;
    const int maxiter = 100;

    // first (negative) residual
    P.RHS(0., f);
    APPLY_COARSE(P, xk, eta, rk, 1.0e-13, jmax, CDD1);
    rk -= f;

    const double normr0 = l2_norm_sqr(rk);
    double normrk = normr0, rhok = 0, oldrhok = 0;

    for (int iterations = 1; normrk/normr0 > tol*tol && iterations <= maxiter; iterations++) {
      zk = rk;
      rhok = rk * zk;

      if (iterations == 1) // TODO: shift this case upwards!
	pk = zk;
      else
	pk.sadd(rhok/oldrhok, zk);

      APPLY_COARSE(P, pk, eta, Apk, 1.0e-13, jmax, CDD1);
      //cout << pk << endl;
      const double alpha = rhok/(pk*Apk);
      xk.add(-alpha,  pk);
      rk.add(-alpha, Apk);
      normrk = l2_norm_sqr(rk);
      
      cout << "loop: " << iterations << ", " << "residual error = " << sqrt(normrk) << endl;
      
      oldrhok = rhok;     
      
      eta *= 0.8;
      cout << "eta = " << eta << endl;
      
      u_epsilon = xk;

    }
  }
}

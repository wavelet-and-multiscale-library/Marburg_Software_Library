// implementation for corner_singularity.h

#include <cmath>
#include <utils/tiny_tools.h>

namespace MathTL
{
  CornerSingularity::CornerSingularity(const Point<2>& x,
				       const double w0,
				       const double w,
				       const double t0,
				       const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w), r0(t0), r1(t1)
  {
  }
  
  double
  CornerSingularity::value(const Point<2>& p,
			   const unsigned int component) const
  {
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);

    if (r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;

    return zeta(r) * pow(r, 1.0/omega) * sin(theta/omega);
  }

  inline
  void
  CornerSingularity::vector_value(const Point<2> &p,
				  Vector<double>& values) const
  {
    values[0] = value(p);
  }

  double
  CornerSingularity::zeta(const double r) const {
    if (r <= r0)
      return 1.0;
    else {
      if (r >= r1)
	return 0.0;
      else {
	const double help1(r1-r);
	const double help0(r-r0);
	return exp(-1.0/(help1*help1))/(exp(-1.0/(help0*help0))+exp(-1.0/(help1*help1)));
      }
    }
  }

  CornerSingularityRHS::CornerSingularityRHS(const Point<2>& x,
					     const double w0,
					     const double w,
					     const double t0,
					     const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w), r0(t0), r1(t1)
  {
  }
  
  double
  CornerSingularityRHS::value(const Point<2>& p,
			      const unsigned int component) const
  {
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    if (r <= r0 || r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;

#if 1
    // Due to the representation of the Laplacian in polar coordinates,
    // we have to calculate
    //   -d^2/dr^2 s(r,phi) - 1/r * d/dr s(r,phi) - 1/r^2 * d^2/dphi^2 s(r,phi)
    // Since u(r,phi)=r^{1/omega}sin(phi/omega) is harmonic, this reduces to
    //   -zeta''(r)*u(r,phi)+2*zeta'(r)u_r(r,phi)+1/r*zeta'(r)u(r,phi)

    const double zeta_prime_r = zeta_prime(r);
    const double zeta_primeprime_r = zeta_primeprime(r);

    return
      -(zeta_primeprime_r+zeta_prime_r/r*(1+2/omega))
      * pow(r, 1/omega)
      * sin(theta/omega);
#else
    // code from Markus Juergens' diploma thesis
    const double t1 = 1/r;
    const double t2 = r-r1;
    const double t3 = t2*t2;
    const double t7 = exp(-1/t3);
    const double t8 = 1/t3/t2*t7;
    const double t9 = r0-r;
    const double t10 = t9*t9;
    const double t12 = exp(-1/t10);
    const double t13 = t7+t12;
    const double t14 = 1/t13;
    const double t15 = 1/omega;
    const double t16 = pow(r,t15);
    const double t19 = sin(theta*t15);
    const double t20 = t14*t16*t19;
    const double t23 = t13*t13;
    const double t24 = 1/t23;
    const double t25 = t7*t24;
    const double t26 = t16*t19;
    const double t30 = t8-1/t10/t9*t12;
    const double t31 = 2.0*t26*t30;
    const double t34 = t7*t14*t16;
    const double t35 = t15*t1;
    const double t38 = t3*t3;
    const double t40 = 1/t38*t7;
    const double t45 = 1/t38/t3*t7;
    const double t71 = t10*t10;
    const double t82 = omega*omega;
    const double t84 = r*r;
    const double t85 = 1/t84;
    const double t88 = t34/t82*t85*t19;
    return -t1*(2.0*t8*t20-t25*t31+t34*t35*t19+r*
		(-6.0*t40*t20+4.0*t45*t20-4.0*t8*t24*t31+4.0*t8*t14*t16*t15*t1*t19+8.0*t7/t23/t13*t26*t30*t30-4.0*t25*t16*
		 t35*t19*t30-t25*t26*(-6.0*t40+4.0*t45-6.0/t71*t12+4.0/t71/t10*t12)+t88-t34*t15*
		 t85*t19))+t88;
#endif
  }

  inline
  void
  CornerSingularityRHS::vector_value(const Point<2> &p,
				     Vector<double>& values) const
  {
    values[0] = value(p);
  }

  double
  CornerSingularityRHS::zeta(const double r) const {
    if (r <= r0)
      return 1.0;
    else {
      if (r >= r1)
	return 0.0;
      else {
	const double help1 = r1-r;
	const double help0 = r-r0;
	return exp(-1.0/(help1*help1))/(exp(-1.0/(help0*help0))+exp(-1.0/(help1*help1)));
      }
    }
  }

  double
  CornerSingularityRHS::zeta_prime(const double r) const {
    if (r <= r0 || r >= r1)
      return 0.0;
    else {
      const double help1 = r1-r;
      const double help0 = r-r0;
      return
	-2.0
	* exp(-(2*r*r-2*r*r0+r0*r0+r1*r1-2*r1*r)/(help1*help1*help0*help0))
	* (3*r*r*r0-3*r*r0*r0+r0*r0*r0-r1*r1*r1+3*r1*r1*r-3*r1*r*r)
	/ (exp(-1/(help0*help0))+exp(-1/(help1*help1)))
	/ (exp(-1/(help0*help0))+exp(-1/(help1*help1)))
	/ (-help0*help0*help0*help1*help1*help1);
    }
  }

  double
  CornerSingularityRHS::zeta_primeprime(const double r) const {
    if (r <= r0 || r >= r1)
      return 0.0;
    else {
      const double help0   = r-r0;
      const double help0_2 = help0*help0;
      const double help0_3 = help0_2*help0;
      const double help0_4 = help0_2*help0_2;
      const double help0_6 = help0_3*help0_3;
      const double help1   = r1-r;
      const double help1_2 = help1*help1;
      const double help1_3 = help1_2*help1;
      const double help1_4 = help1_2*help1_2;
      const double help1_6 = help1_3*help1_3;

      const double denom = exp(-1/help0_2)+exp(-1/help1_2);
      const double numer = 2*exp(-1/help0_2)/help0_3-2*exp(-1/help1_2)/help1_3;
      
      return
	-6 * exp(-1/help1_2)/(help1_4*denom)
	+ 4 * exp(-1/help1_2)/(help1_6*denom)
	+ 4 * exp(-1/help1_2)*numer/(help1_3*denom*denom)
	+ 2 * exp(-1/help1_2)*numer*numer/(denom*denom*denom)
	-exp(-1/help1_2)*(-6*exp(-1/help0_2)/help0_4+4*exp(-1/help0_2)/help0_6-6*exp(-1/help1_2)/help1_4+4*exp(-1/help1_2)/help1_6)/(denom*denom);
    }
  }

  CornerSingularityGradient::CornerSingularityGradient(const Point<2>& x,
						       const double w0,
						       const double w,
						       const double t0,
						       const double t1)
    : Function<2>(2), x0(x), theta0(w0), omega(w), r0(t0), r1(t1)
  {
  }

  inline
  void
  CornerSingularityGradient::vector_value(const Point<2> &p,
					  Vector<double>& values) const
  {
    values[0] = value(p,0);
    values[1] = value(p,1);
  }
  double
  CornerSingularityGradient::value(const Point<2>& p,
			      const unsigned int component) const
  {
    //cout << "com = " << component << endl;
    double res = 0.;
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);

    if (r >= r1) return 0.0;
    
    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;

    // note that the gradient in polar coordinates looks like
    // \nabla f(r,\theta) = \frac{\partial f}{\parial r} {\cos(\theta) \choose \sin(\theta)}
    // + \frac{1}{r} \frac{\partial f}{\parial \theta} {-\sin(\theta) \choose \cos(\theta)}

    // but due to the rotation about 90 degrees in this case it is
    // \frac{\partial f}{\parial r} {-\sin(\theta) \choose \cos(\theta)}
    // + \frac{1}{r} \frac{\partial f}{\parial \theta} {-\cos(\theta) \choose -\sin(\theta)}


    double power = pow(r, 1.0/omega);
    double Zeta = zeta(r);
    double Zeta_prime = zeta_prime(r);
    double theta_div_omega = theta/omega;
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    double dr = (Zeta_prime +  Zeta / omega / r) * sin(theta_div_omega);
    double dtheta = Zeta*cos(theta_div_omega)/omega;

//     double dr = (zeta_prime(r)*pow(r, 1.0/omega) +  zeta(r) / omega * pow(r, (1.0/omega) - 1.))*sin(theta/omega);
//     double dtheta = zeta(r)*pow(r,1.0/omega)*cos(theta/omega)/omega;

    // first derivative with respect to x
    if ( component == 1 ) {
      res = (dr * cos_theta - dtheta * sin_theta/r)*power;
      //res = dr * cos(theta) - dtheta * sin(theta)/r;
      
    }
    // first derivative with respect to y
    else if ( component == 0 ) {
      res = -power*(dr * sin_theta + dtheta * cos_theta/r);
      //res = -(dr * sin(theta) + dtheta * cos(theta)/r);

    }
    return res;
  }

  double
  CornerSingularityGradient::zeta(const double r) const {
    if (r <= r0)
      return 1.0;
    else {
      if (r >= r1)
	return 0.0;
      else {
	const double help1 = r1-r;
	const double help0 = r-r0;
	return exp(-1.0/(help1*help1))/(exp(-1.0/(help0*help0))+exp(-1.0/(help1*help1)));
      }
    }
  }

  double
  CornerSingularityGradient::zeta_prime(const double r) const {
    if (r <= r0 || r >= r1)
      return 0.0;
    else {
      const double help1 = r1-r;
      const double help0 = r-r0;
      return
	-2.0
	* exp(-(2*r*r-2*r*r0+r0*r0+r1*r1-2*r1*r)/(help1*help1*help0*help0))
	* (3*r*r*r0-3*r*r0*r0+r0*r0*r0-r1*r1*r1+3*r1*r1*r-3*r1*r*r)
	/ (exp(-1/(help0*help0))+exp(-1/(help1*help1)))
	/ (exp(-1/(help0*help0))+exp(-1/(help1*help1)))
	/ (-help0*help0*help0*help1*help1*help1);
    }
  }
  
  CornerTimeSingularity::CornerTimeSingularity(const Point<2>& x,
					       const double w0,
					       const double w)
    : Function<2>(), x0(x), theta0(w0), omega(w)
  {
  }
  
  double
  CornerTimeSingularity::value(const Point<2>& p,
			       const unsigned int component) const
  {
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;

    return pow(get_time(),0.75) * pow(r, 1.0/omega) * sin(theta/omega) * (1-x[0]*x[0]) * (1-x[1]*x[1]);
  }

  inline
  void
  CornerTimeSingularity::vector_value(const Point<2> &p,
				      Vector<double>& values) const
  {
    values[0] = value(p);
  }

  CornerTimeSingularityRHS::CornerTimeSingularityRHS(const Point<2>& x,
						     const double w0,
						     const double w)
    : Function<2>(), x0(x), theta0(w0), omega(w)
  {
  }
  
  double
  CornerTimeSingularityRHS::value(const Point<2>& p,
				  const unsigned int component) const
  {
    if (get_time() == 0) return 0;

    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;

    const double ct = cos(theta+theta0*M_PI);
    const double st = sin(theta+theta0*M_PI);

    return
      0.75*pow(get_time(),-0.25)*pow(r,1/omega)*sin(theta/omega)*(1-x[0]*x[0])*(1-x[1]*x[1]) // temporal derivative
      - pow(get_time(),0.75) * // spatial derivative:
      (2*1/omega*pow(r,1/omega-1)*sin(theta/omega)*(-2*r+4*r*r*r*ct*ct*st*st)           // 2*u_r*f_r(r,theta)
       + pow(r,1/omega)*sin(theta/omega)*(-2+12*r*r*ct*ct*st*st)                        // u*f_{rr}
       + pow(r,1/omega-1)*sin(theta/omega)*(-2*r+4*r*r*r*ct*ct*st*st)                   // 1/r*u*f_r
       + 2/(omega)*pow(r,1/omega+2)*cos(theta/omega)*(2*st*ct*ct*ct-2*st*st*st*ct)      // 2/r^2*u_theta*f_theta
       + pow(r,1/omega+2)*sin(theta/omega)*(2*ct*ct*ct*ct-12*ct*ct*st*st+2*st*st*st*st) // 1/r^2*u*f_{thetatheta}
       );
  }
  
  inline
  void
  CornerTimeSingularityRHS::vector_value(const Point<2> &p,
					 Vector<double>& values) const
  {
    values[0] = value(p);
  }

  CornerTimeSingularityRHSt::CornerTimeSingularityRHSt(const Point<2>& x,
						       const double w0,
						       const double w)
    : Function<2>(), x0(x), theta0(w0), omega(w)
  {
  }
  
  double
  CornerTimeSingularityRHSt::value(const Point<2>& p,
				   const unsigned int component) const
  {
    if (get_time() == 0) return 0;

    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    double theta = atan2(x[1],x[0]);
    
    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega * M_PI) return 0.0;
    
    const double ct = cos(theta+theta0*M_PI);
    const double st = sin(theta+theta0*M_PI);
    
    return
      -0.1875*pow(get_time(),-1.25)*pow(r,1/omega)*sin(theta/omega)*(1-x[0]*x[0])*(1-x[1]*x[1]) // temporal derivative
      -0.75*pow(get_time(),-0.25) * // spatial derivative:
      (2*1/omega*pow(r,1/omega-1)*sin(theta/omega)*(-2*r+4*r*r*r*ct*ct*st*st)           // 2*u_r*f_r(r,theta)
       + pow(r,1/omega)*sin(theta/omega)*(-2+12*r*r*ct*ct*st*st)                        // u*f_{rr}
       + pow(r,1/omega-1)*sin(theta/omega)*(-2*r+4*r*r*r*ct*ct*st*st)                   // 1/r*u*f_r
       + 2/(omega)*pow(r,1/omega+2)*cos(theta/omega)*(2*st*ct*ct*ct-2*st*st*st*ct)      // 2/r^2*u_theta*f_theta
       + pow(r,1/omega+2)*sin(theta/omega)*(2*ct*ct*ct*ct-12*ct*ct*st*st+2*st*st*st*st) // 1/r^2*u*f_{thetatheta}
       );
  }
  
  inline
  void
  CornerTimeSingularityRHSt::vector_value(const Point<2> &p,
					  Vector<double>& values) const
  {
    values[0] = value(p);
  }

}

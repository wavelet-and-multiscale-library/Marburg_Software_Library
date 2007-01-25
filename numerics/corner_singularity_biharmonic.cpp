// implementation for corner_singularity.h

#include <cmath>
#include <utils/tiny_tools.h>

namespace MathTL
{
  CornerSingularityBiharmonic::CornerSingularityBiharmonic(const Point<2>& x,
<<<<<<< corner_singularity_biharmonic.cpp
				       const double w0,
				       const double w,
				       const double t0,
				       const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w*M_PI), r0(t0), r1(t1)
=======
							   const double w0,
							   const double w,
							   const double t0,
							   const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(M_PI*w), r0(t0), r1(t1)
>>>>>>> 1.2
  {
  }
  
  double
  CornerSingularityBiharmonic::value(const Point<2>& p,
				     const unsigned int component) const
  {
    double res=0.0;
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    double rootz[2];
    rootz[0]=0.544483736782463929140876854601;
    rootz[1]=0.908529189846098818660368657869;
    if (r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
<<<<<<< corner_singularity_biharmonic.cpp
   if (theta >= omega) return 0.0;
 
    //   omega*=M_PI;
=======
    if (theta >= omega /** M_PI*/) return 0.0;

    

>>>>>>> 1.2
    for(int m=0;m<2;m++)    
<<<<<<< corner_singularity_biharmonic.cpp
	  res+= zeta(r) * pow(r, 1.0+rootz[m]) * ((1/(rootz[m]-1)*sin((rootz[m]-1)*omega)
					     -1/(rootz[m]+1)*sin((rootz[m]+1)*omega))
					    *(cos((rootz[m]-1)*theta)-cos((rootz[m]+1)*theta))
						  -(1/(rootz[m]-1)*sin((rootz[m]-1)*theta)
					      -1/(rootz[m]+1)*sin((rootz[m]+1)*theta))
					    *(cos((rootz[m]-1)*omega)-cos((rootz[m]+1)*omega)));
=======
	  res += zeta(r) * pow(r, 1.0+rootz[m]) * ((1./(rootz[m]-1)*sin((rootz[m]-1.)*omega)
					     -1./(rootz[m]+1)*sin((rootz[m]+1.)*omega))
					    *(cos((rootz[m]-1.)*theta)-cos((rootz[m]+1.)*theta))
						  -(1./(rootz[m]-1.)*sin((rootz[m]-1.)*theta)
					      -1./(rootz[m]+1.)*sin((rootz[m]+1.)*theta))
					    *(cos((rootz[m]-1.)*omega)-cos((rootz[m]+1.)*omega)));
>>>>>>> 1.2
    return res;
  }

  inline
  void
  CornerSingularityBiharmonic::vector_value(const Point<2> &p,
				  Vector<double>& values) const
  {
    values[0] = value(p);
  }

  double
  CornerSingularityBiharmonic::zeta(const double r) const {
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

  CornerSingularityBiharmonicRHS::CornerSingularityBiharmonicRHS(const Point<2>& x,
					     const double w0,
					     const double w,
					     const double t0,
					     const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w*M_PI), r0(t0), r1(t1)
  {
  }
  
  double
  CornerSingularityBiharmonicRHS::value(const Point<2>& p,
			      const unsigned int component) const
  {
    double res=0.0;
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    if (r <= r0 || r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega) return 0.0;

#if 1
    // Due to the representation of the Laplacian in polar coordinates,
    // we have to calculate
    //   -d^2/dr^2 s(r,phi) - 1/r * d/dr s(r,phi) - 1/r^2 * d^2/dphi^2 s(r,phi)
    // Since u(r,phi)=r^{1/omega}sin(phi/omega) is harmonic, this reduces to
    //   -zeta''(r)*u(r,phi)+2*zeta'(r)u_r(r,phi)+1/r*zeta'(r)u(r,phi)

    const double zeta_r=zeta(r);
    const double zeta_prime_r = zeta_prime(r);
    const double zeta_primeprime_r = zeta_primeprime(r);
    const double zeta_third_der_r=zeta_third_der(r);
    const double zeta_fourth_der_r = zeta_fourth_der(r);
    double root_z[2];
    root_z[0] = 0.544483736782463929140876854601;
    root_z[1] = 0.908529189846098818660368657869;
    const double r2 = r*r;
    const double r3 = r2*r;
    const double r4 = r3*r;
    double sin_th[2];
    double sin_th1[2];
    double sin_om[2];
    double sin_om1[2];
    double cos_th[2];
    double cos_th1[2];
    double cos_om[2];
    double cos_om1[2];
    double rest[2];
    double rest2[2];
    double rest4[2];
    double r_z[2];
    double root_z2[2];
    double root_z21[2];
    double root_z3[2];
    double root_z31[2];
    double root_z4[2];
    double root_z41[2];

      for(int i=0;i<2;i++)
	{
	  r_z[i] = pow(r,1+root_z[i]);
	  root_z2[i]=(root_z[i]-1)*(root_z[i]-1);
	  root_z21[i]=(root_z[i]+1)*(root_z[i]+1);
	  root_z3[i]=root_z2[i]*(root_z[i]-1);
	  root_z31[i]=root_z21[i]*(root_z[i]+1);
	  root_z4[2]=root_z2[i]*root_z2[i];
	  root_z41[2]=root_z21[i]*root_z21[i];
	  sin_th[i] = sin((root_z[i]-1)*theta);
	  sin_th1[i] = sin((root_z[i]+1)*theta);
	  sin_om[i] = sin((root_z[i]-1)*omega);
	  sin_om1[i] = sin((root_z[i]+1)*omega);
	  cos_th[i] = cos((root_z[i]-1)*theta);
	  cos_th1[i] = cos((root_z[i]+1)*theta);
	  cos_om[i] = cos((root_z[i]-1)*omega);
	  cos_om1[i] = cos((root_z[i]+1)*omega);
	  rest[i] = (sin_om[i]/(root_z[i]-1)-sin_om1[i]/(root_z[i]+1))*(cos_th[i]-cos_th1[i])
	    -(sin_th[i]/(root_z[i]-1)-sin_th1[i]/(root_z[i]+1))*(cos_om[i]-cos_om1[i]);
	  rest2[i] = (sin_om[i]/(root_z[i]-1)-sin_om1[i]/(root_z[i]+1))
	    *(-cos_th[i]*root_z2[i]+cos_th1[i]*root_z21[i])
	    -(-sin_th[i]*(root_z[i]-1)+sin_th1[i]*(root_z[i]+1))*(cos_om[i]-cos_om1[i]);
	  rest4[i]=(sin_om[i]/(root_z[i]-1)-sin_om1[i]/(root_z[i]+1))
	    *(cos_th[i]*root_z4[i]-cos_th1[i]*root_z41[i])
	    -(sin_th[i]*root_z3[i]-sin_th1[i]*root_z31[i])*(cos_om[i]-cos_om1[i]);
	}

      for(int i=0;i<2;i++)
	res += rest[i]*r_z[i]*(zeta_fourth_der_r+4*zeta_third_der_r*(1+root_z[i])/r
			       +(6*zeta_primeprime_r*root_z21[i]-6*zeta_primeprime_r*(1+root_z[i]))/r2
			       +(4*zeta_prime_r*root_z31[i]-12*zeta_prime_r*root_z21[i]
				 +8*zeta_prime_r*(1+root_z[i]))/r3
			       +(zeta_r*root_z41[i]-6*zeta_r*root_z31[i]
				 +11*zeta_r*root_z21[i]-6*zeta_r*(1+root_z[i]))/r4)
	  +2/r*rest[i]*r_z[i]*(zeta_third_der_r+3*zeta_primeprime_r*(1+root_z[i])/r
			       +(3*zeta_prime_r*root_z21[i]-3*zeta_prime_r*(1+root_z[i]))/r2
			       +(zeta_r*root_z31[i]-3*zeta_r*root_z21[i]
				 +2*zeta_r*(1+root_z[i]))/r3)
	  +1/r2*rest[i]*r_z[i]*(zeta_primeprime_r+2*zeta_prime_r*(1+root_z[i])/r
				+(zeta_r*root_z21[i]-zeta_r*(1+root_z[i]))/r2)
	  +2/r2*rest2[i]*r_z[i]*(zeta_primeprime_r+2*zeta_prime_r*(1+root_z[i])/r
				+(zeta_r*root_z21[i]-zeta_r*(1+root_z[i]))/r2)
	  +2/r3*rest2[i]*r_z[i]*(zeta_prime_r+zeta_r*(1+root_z[i])/r)
	  +1/r4*rest4[i]*zeta_r*r_z[i];

     
      return res;




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
  CornerSingularityBiharmonicRHS::vector_value(const Point<2> &p,
				     Vector<double>& values) const
  {
    values[0] = value(p);
  }


 double
  CornerSingularityBiharmonicRHS::zeta(const double r) const {
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


  double
  CornerSingularityBiharmonicRHS::zeta_prime(const double r) const {
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
  CornerSingularityBiharmonicRHS::zeta_primeprime(const double r) const {
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

   double
  CornerSingularityBiharmonicRHS::zeta_third_der(const double r) const {
    if (r <= r0 || r >= r1)
      return 0.0;
    else {
      const double help0 = r-r0;
      const double help0_2 = help0*help0;
      const double help0_3 = help0_2*help0;
      const double help0_4 = help0_2*help0_2;
      const double help0_5 = help0_4*help0;
      const double help0_6 = help0_3*help0_3;
      const double help0_7 = help0_6*help0;
      //  const double help0_8 = help0_7* help0;
      const double help0_9 = help0_7*help0_2;
      const double help1   = r1-r;
      const double help1_2 = help1*help1;
      const double help1_3 = help1_2*help1;
      const double help1_4 = help1_2*help1_2;
      const double help1_5 = help1_4*help1;
      const double help1_6 = help1_3*help1_3;
      const double help1_7 = help1_6*help1;
      //    const double help1_8 = help1_7*help1;
      const double help1_9 = help1_7*help1_2;
      const double h0 = -1/help0_2;
      const double h1 = -1/help1_2;
      

      return -24*exp(h1)/(help1_5*(exp(h0)+exp(h1)))
	+36*exp(h1)/(help1_7*(exp(h0)+exp(h1)))
	+18*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))/(help1_4*pow(exp(h0)+exp(h1),2))
	-8*exp(h1)/(help1_9*(exp(h0)+exp(h1)))
	-12*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))/(help1_6*pow(exp(h0)+exp(h1),2))
	-12*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),2)/(help1_3*pow(exp(h0)+exp(h1),3))
	+6*exp(h1)*(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)
	/(help1_3*pow(exp(h0)+exp(h1),2))
	-6*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),3)/pow(exp(h0)+exp(h1),4)
	+6*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))
	*(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)
	/pow(exp(h0)+exp(h1),3)
	-exp(h1)
	*(24*exp(h0)/help0_5-36*exp(h0)/help0_7+8*exp(h0)/help0_9
	  -24*exp(h1)/help1_5+36*exp(h1)/help1_7-8*exp(h1)/help1_9)
	/pow(exp(h0)+exp(h1),2);
    }
  }

   double
  CornerSingularityBiharmonicRHS::zeta_fourth_der(const double r) const {
    if (r <= r0 || r >= r1)
      return 0.0;
    else {
      const double help0 = r-r0;
      const double help0_2 = help0*help0;
      const double help0_3 = help0_2*help0;
      const double help0_4 = help0_2*help0_2;
      const double help0_5 = help0_4*help0;
      const double help0_6 = help0_3*help0_3;
      const double help0_7 = help0_6*help0;
      const double help0_8 = help0_7* help0;
      const double help0_9 = help0_8*help0;
      const double help0_10 = help0_9* help0;     
      const double help0_12 = help0_10* help0_2;
      const double help1   = r1-r;
      const double help1_2 = help1*help1;
      const double help1_3 = help1_2*help1;
      const double help1_4 = help1_2*help1_2;
      const double help1_5 = help1_4*help1;
      const double help1_6 = help1_3*help1_3;
      const double help1_7 = help1_6*help1;
      const double help1_8 = help1_7*help1;
      const double help1_9 = help1_8*help1;
      const double help1_10 = help1_9*help1;
      const double help1_12 = help1_10*help1_2;
      const double h0 = -1/help0_2;
      const double h1 = -1/help1_2;


     return 24*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),4)/pow(exp(h0)+exp(h1),5)
       -36*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),2)
       *(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)/pow(exp(h0)+exp(h1),4)
       +6*exp(h1)*pow(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6,2)/pow(exp(h0)+exp(h1),3)
       +16*exp(h1)*(exp(h0)/help0_3-exp(h1)/help1_3)
       *4*(6*exp(h0)/help0_5-9*exp(h0)/help0_7+2*exp(h0)/help0_9-6*exp(h1)/help1_5+9*exp(h1)/help1_7-2*exp(h1)/help1_9)
       /pow(exp(h0)+exp(h1),3)
       -72*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),2)
       /(help1_4*pow(exp(h0)+exp(h1),3))
       +48*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),2)
       /(help1_6*pow(exp(h0)+exp(h1),3))
       +48*exp(h1)*pow(2*(exp(h0)/help0_3-exp(h1)/help1_3),3)
       /(help1_3*pow(exp(h0)+exp(h1),4))
       -48*exp(h1)*2*(exp(h0)/help0_3-exp(h1)/help1_3)
       *(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)
       /(help1_3*pow(exp(h0)+exp(h1),3))
       -exp(h1)
       *(-120*exp(h0)/help0_6+300*exp(h0)/help0_8-144*exp(h0)/help0_10+16*exp(h0)/help0_12
	 -120*exp(h1)/help1_6+300*exp(h1)/help1_8-144*exp(h1)/help1_10+16*exp(h1)/help1_12)
       /pow(exp(h0)+exp(h1),2)
       -120*exp(h1)/(help1_6*(exp(h0)+exp(h1)))
       +300*exp(h1)/(help1_8*(exp(h0)+exp(h1)))
       -144*exp(h1)/(help1_10*(exp(h0)+exp(h1)))
       +16*exp(h1)/(help1_12*(exp(h0)+exp(h1)))
       +96*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))/(help1_5*pow(exp(h0)+exp(h1),2))
       -144*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))/(help1_7*pow(exp(h0)+exp(h1),2))
       +36*exp(h1)*(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)
       /(help1_4*pow(exp(h0)+exp(h1),2))
       +32*exp(h1)*(2*(exp(h0)/help0_3-exp(h1)/help1_3))
       /(help1_9*pow(exp(h0)+exp(h1),2))
       -24*exp(h1)*(-6*exp(h0)/help0_4+4*exp(h0)/help0_6-6*exp(h1)/help1_4+4*exp(h1)/help1_6)
       /(help1_6*pow(exp(h0)+exp(h1),2))
       +8*exp(h1)
       *(24*exp(h0)/help0_5-36*exp(h0)/help0_7+8*exp(h0)/help0_9
	 -24*exp(h1)/help1_5+36*exp(h1)/help1_7-8*exp(h1)/help1_9)
       /(help1_3*pow(exp(h0)+exp(h1),2));
    }
  }
}

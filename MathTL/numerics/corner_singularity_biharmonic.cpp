// implementation for corner_singularity.h

#include <cmath>
#include <utils/tiny_tools.h>

namespace MathTL
{
  CornerSingularityBiharmonic::CornerSingularityBiharmonic(const Point<2>& x,
				       const double w0,
				       const double w,
				       const double t0,
				       const double t1)
    : Function<2>(), x0(x), theta0(w0), omega(w*M_PI), r0(t0), r1(t1)
  {
  }
  
  double
  CornerSingularityBiharmonic::value(const Point<2>& p,
			   const unsigned int component) const
  {
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
    if (theta >= omega) return 0.0;
 
#if 0
   //seems to work
    double res=0.0;
    for(unsigned int m=0; m<2; m++)
      res+= zeta(r) * pow(r, 1.0+rootz[m]) * ( (1./(rootz[m]-1.)*sin((rootz[m]-1.)*omega)
					     - 1./(rootz[m]+1.)*sin((rootz[m]+1.)*omega))
					    *(cos((rootz[m]-1.)*theta)-cos((rootz[m]+1.)*theta))
						  - (1./(rootz[m]-1.)*sin((rootz[m]-1.)*theta)
					      - 1./(rootz[m]+1.)*sin((rootz[m]+1.)*theta))
					    *(cos((rootz[m]-1.)*omega)-cos((rootz[m]+1.)*omega)) );
    return res;
#else
    //Maple output
  double t56;
  double t58;
  double t21;
  double t32;
  double t61;
  double t62;
  double t35;
  double t36;
  double t16;
  double t24;
  double t25;
  double t26;
  double t42;
  double t27;
  double t43;
  double t44;
  double t17;
  double t1;
  double t46;
  double t3;
  double t47;
  double t5;
  double t18;
  double t6;
  double t7;
  double t50;
  double t51;
  double t9;
  double t52;
  double t11;
  double t53;
  double t12;
  double t15;
  double t20;
  double t30;
  t1 = 0.99e2 / 0.100e3 - r;
  t3 = t1 * t1;
  t5 = exp(-0.1e1 / t3);
  t6 = (0.0e0 < t1) ? t5 : 0;
  t7 = r - 0.1e1 / 0.100e3;
  t9 = t7 * t7;
  t11 = exp(-0.1e1 / t9);
  t12 = (0.0e0 < t7) ? t11 : 0;
  t15 = t6 / (t12 + t6);
  t16 = pow(r, 0.154448373678246e1);
  t17 = 0.683274394826304e0 * 0.3141592654e1;
  t18 = sin(t17);
  t20 = 0.231672560517369e1 * 0.3141592654e1;
  t21 = sin(t20);
  t24 = 0.455516263217536e0 * theta;
  t25 = cos(t24);
  t26 = 0.154448373678246e1 * theta;
  t27 = cos(t26);
  t30 = sin(t24);
  t32 = sin(t26);
  t35 = cos(t17);
  t36 = cos(t20);
  t42 = pow(r, 0.190852918984610e1);
  t43 = 0.137206215230852e0 * 0.3141592654e1;
  t44 = sin(t43);
  t46 = 0.286279378476915e1 * 0.3141592654e1;
  t47 = sin(t46);
  t50 = 0.91470810153901e-1 * theta;
  t51 = cos(t50);
  t52 = 0.190852918984610e1 * theta;
  t53 = cos(t52);
  t56 = sin(t50);
  t58 = sin(t52);
  t61 = cos(t43);
  t62 = cos(t46);
  return(t15 * t16 * ((0.219531130005438e1 * t18 - 0.6474655422e0 * t21) * (t25 - t27) - (0.219531130005438e1 * t30 - 0.6474655422e0 * t32) * (t35 - t36)) + t15 * t42 * ((0.109324493608123e2 * t44 - 0.5239636917e0 * t47) * (t51 - t53) - (0.109324493608123e2 * t56 - 0.5239636917e0 * t58) * (t61 - t62)));
#endif
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
    const Point<2> x(p-x0);
    const double r = hypot(x[0],x[1]);
    
    if (r <= r0 || r >= r1) return 0.0;

    double theta = atan2(x[1],x[0]);

    // shift theta to [0,2*pi]
    if (theta < 0) theta += 2.0 * M_PI;
    theta -= theta0 * M_PI;
    if (theta < 0) theta += 2.0 * M_PI;
    if (theta >= omega) return 0.0;

#if 0
    // Due to the representation of the Laplacian in polar coordinates,
    // we have to calculate
    //   -d^2/dr^2 s(r,phi) - 1/r * d/dr s(r,phi) - 1/r^2 * d^2/dphi^2 s(r,phi)
    // Since u(r,phi)=r^{1/omega}sin(phi/omega) is harmonic, this reduces to
    //   -zeta''(r)*u(r,phi)+2*zeta'(r)u_r(r,phi)+1/r*zeta'(r)u(r,phi)

    //error inside

    double res=0.0;

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

#endif


#if 0
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
#else
    //don't try to understand, it's just a Maple output
/* The options were    : operatorarrow */
  double t10;
  double t54;
  double t68;
  double t494;
  double t52;
  double t65;
  double t64;
  double t60;
  double t69;
  double t14;
  double t71;
  double t72;
  double t15;
  double t53;
  double t16;
  double t289;
  double t74;
  double t73;
  double t29;
  double t47;
  double t79;
  double t83;
  double t77;
  double t76;
  double t81;
  double t84;
  double t27;
  double t91;
  double t89;
  double t88;
  double t87;
  double t45;
  double t44;
  double t613;
  double t93;
  double t92;
  double t43;
  double t97;
  double t96;
  double t332;
  double t31;
  double t196;
  double t415;
  double t2;
  double t34;
  double t36;
  double t9;
  double t3;
  double t4;
  double t39;
  double t25;
  double t38;
  double t23;
  double t20;
  double t70;
  double t105;
  double t19;
  double t120;
  double t111;
  double t114;
  double t569;
  double t115;
  double t119;
  double t121;
  double t122;
  double t125;
  double t128;
  double t130;
  double t131;
  double t136;
  double t140;
  double t141;
  double t144;
  double t151;
  double t160;
  double t101;
  double t102;
  double t109;
  double t154;
  double t147;
  double t508;
  double t171;
  double t178;
  double t182;
  double t186;
  double t191;
  double t185;
  double t202;
  double t201;
  double t486;
  double t225;
  double t232;
  double t233;
  double t237;
  double t1;
  double t349;
  double t206;
  double t229;
  double t8;
  double t13;
  double t17;
  double t18;
  double t11;
  double t22;
  double t287;
  double t249;
  double t248;
  double t250;
  double t251;
  double t252;
  double t261;
  double t263;
  double t30;
  double t276;
  double t278;
  double t281;
  double t284;
  double t296;
  double t303;
  double t393;
  double t400;
  double t427;
  double t428;
  double t317;
  double t319;
  double t324;
  double t327;
  double t346;
  double t367;
  double t366;
  double t382;
  double t372;
  double t373;
  double t379;
  double t380;
  double t385;
  double t442;
  double t502;
  double t467;
  double t477;
  double t478;
  double t481;
  double t483;
  double t489;
  double t491;
  double t516;
  {
    t1 = 0.99-r;
    t2 = t1*t1;
    t3 = t2*t2;
    t4 = t3*t3;
    t8 = exp(-1/t2);
    t9 = t8/t4/t1;
    t10 = 0.1E-1-r;
    t11 = t10*t10;
    t13 = exp(-1/t11);
    t14 = t13+t8;
    t15 = t14*t14;
    t16 = 1/t15;
    t17 = t16*t9;
    t18 = pow(r,0.190852919E1);
    t19 = 0.914708102E-1*theta;
    t20 = cos(t19);
    t22 = 0.190852919E1*theta;
    t23 = cos(t22);
    t25 = sin(t19);
    t27 = sin(t22);
    t29 = 0.434888792E1*t20-0.434888792E1*t23-0.1986489872E2*t25+0.9520726166*
t27;
    t30 = t29*t18;
    t31 = t11*t10;
    t34 = t2*t1;
    t36 = t8/t34;
    t38 = -2.0*t13/t31-2.0*t36;
    t39 = t38*t30;
    t43 = t8/t3;
    t44 = 1/t14;
    t45 = pow(r,-0.9147081E-1);
    t47 = t29*t45*t44;
    t52 = t8/t3/t2;
    t53 = t16*t52;
    t54 = t11*t11;
    t60 = t13/t54/t11;
    t64 = -6.0*t13/t54+4.0*t60-6.0*t43+4.0*t52;
    t65 = t64*t30;
    t68 = 1/r;
    t69 = r*r;
    t70 = 1/t69;
    t71 = pow(r,0.1544483737E1);
    t72 = t71*t44;
    t73 = 0.4555162632*theta;
    t74 = cos(t73);
    t76 = 0.1544483737E1*theta;
    t77 = cos(t76);
    t79 = sin(t73);
    t81 = sin(t76);
    t83 = 0.1298288751E1*t74-0.1298288751E1*t77+0.2390622594E1*t79-0.7050689139
*t81;
    t84 = t83*t72;
    t87 = t16*t8;
    t88 = t83*t71;
    t89 = t38*t88;
    t91 = t44*t8;
    t92 = pow(r,0.544483737);
    t93 = t83*t92;
    t96 = t18*t44;
    t97 = t29*t96;
    t101 = pow(r,0.90852919);
    t102 = t29*t101;
    t105 = -2.0*t84*t36-t89*t87+0.1544483737E1*t93*t91-2.0*t97*t36-t39*t87+
0.190852919E1*t102*t91;
    t109 = t84*t52;
    t111 = t16*t36;
    t114 = t92*t44;
    t115 = t83*t114;
    t119 = 1/t15/t14;
    t120 = t119*t8;
    t121 = t38*t38;
    t122 = t121*t88;
    t125 = t38*t93;
    t128 = t64*t88;
    t130 = pow(r,-0.455516263);
    t131 = t83*t130;
    t136 = t97*t52;
    t140 = t101*t44;
    t141 = t29*t140;
    t144 = t121*t30;
    t147 = t38*t102;
    t151 = t29*t45;
    t154 = -6.0*t84*t43+4.0*t109+4.0*t89*t111-0.6177934948E1*t115*t36+2.0*t122*
t120-0.3088967474E1*t125*t87-t128*t87+0.8409462769*t131*t91-6.0*t97*t43+4.0*
t136+4.0*t39*t111-0.763411676E1*t141*t36+2.0*t144*t120-0.381705838E1*t147*t87-
t65*t87+0.1733954479E1*t151*t91;
    t160 = t16*t43;
    t171 = t54*t54;
    t178 = t8/t3/t1;
    t182 = t8/t3/t34;
    t185 = -24.0*t13/t54/t10+36.0*t13/t54/t31-8.0*t13/t171/t10-24.0*t178+36.0*
t182-8.0*t9;
    t186 = t185*t30;
    t191 = t8/t4/t3;
    t196 = t8/t4/t2;
    t201 = pow(r,-0.1455516263E1);
    t202 = t83*t201;
    t206 = pow(r,-0.109147081E1);
    t225 = t8/t4;
    t229 = -120.0*t60+300.0*t13/t171-144.0*t13/t171/t11+16.0*t13/t171/t54-120.0
*t52+300.0*t225-144.0*t196+16.0*t191;
    t232 = t119*t36;
    t233 = t38*t64;
    t237 = t64*t102;
    t248 = t15*t15;
    t249 = 1/t248;
    t250 = t249*t8;
    t251 = t121*t38;
    t252 = t251*t30;
    t261 = t71*t120;
    t263 = t64*t38*t83;
    t276 = t185*t88;
    t278 = t64*t93;
    t281 = t38*t131;
    t284 = t251*t88;
    t287 = -0.1040372687E2*t47*t36+6.0*t65*t111+0.2290235028E2*t147*t111+
0.2290235028E2*t141*t52-6.0*t252*t250-12.0*t144*t232-0.3830647054*t202*t91-8.0*
t97*t9+6.0*t263*t261+36.0*t97*t182-24.0*t97*t178-12.0*t39*t53-0.3435352542E2*
t141*t43+18.0*t39*t160-t276*t87-0.4633451211E1*t278*t87-0.2522838831E1*t281*t87
-6.0*t284*t250;
    t289 = t83*t130*t44;
    t296 = t121*t93;
    t303 = t29*t206;
    t317 = t18*t120;
    t319 = t64*t38*t29;
    t324 = t121*t102;
    t327 = t38*t151;
    t332 = -0.5045677661E1*t289*t36-12.0*t122*t232-8.0*t84*t9+0.9266902422E1*
t296*t120-24.0*t84*t178+0.1853380485E2*t125*t111-0.1586062207*t303*t91-12.0*t89
*t53+36.0*t84*t182+18.0*t89*t160-0.2780070726E2*t115*t43-t186*t87+6.0*t128*t111
+6.0*t319*t317-0.572558757E1*t237*t87+0.1145117514E2*t324*t120-0.5201863437E1*
t327*t87+0.1853380485E2*t115*t52;
    t346 = t16*t182;
    t349 = 32.0*t39*t17-0.6242236123E2*t47*t43-24.0*t65*t53+(-t105*t70+t154*t68
)*t68-2.0*t154*t70+36.0*t65*t160+8.0*t186*t111+16.0*t97*t191-144.0*t97*t196
-0.6107293408E2*t141*t9+0.1532258821E1*t38*t202*t87+0.1268849765E1*t29*t206*t44
*t36-t229*t30*t87-48.0*t233*t30*t232+0.4580470056E2*t237*t111+2.0*(t287+t332)*
t68+0.6344248828*t38*t303*t87-0.763411676E1*t185*t102*t87+8.0*t38*t185*t29*t317
-144.0*t89*t346;
    t366 = t69*t69;
    t367 = 1/t366;
    t372 = -0.2693885101*t74+0.3096976954E1*t77-0.4960423931*t79+0.1681892549E1
*t81;
    t373 = t372*t71;
    t379 = -0.3638674999E-1*t20+0.1584075323E2*t23+0.1662078022*t25
-0.3467908958E1*t27;
    t380 = t379*t18;
    t382 = t373*t91+t380*t91;
    t385 = t119*t43;
    t393 = t16*t178;
    t400 = t119*t52;
    t415 = -144.0*t39*t346+96.0*t89*t393-0.7413521938E2*t296*t232+300.0*t97*
t225+48.0*t122*t400+16.0*t84*t191+96.0*t39*t393-144.0*t84*t196-0.494234796E2*
t115*t9+32.0*t89*t17-0.1832188023E3*t141*t178;
    t427 = 1/t248/t14*t8;
    t428 = t121*t121;
    t442 = t249*t36;
    t467 = pow(r,-0.2455516263E1);
    t477 = 1/t69/r;
    t478 = t372*t72;
    t481 = t38*t373;
    t483 = t372*t92;
    t486 = t379*t96;
    t489 = t38*t380;
    t491 = t379*t101;
    t494 = -2.0*t478*t36-t481*t87+0.1544483737E1*t483*t91-2.0*t486*t36-t489*t87
+0.190852919E1*t491*t91;
    t502 = 0.3706760969E2*t278*t111-0.3706760969E2*t251*t93*t250+0.3706760968E2
*t263*t92*t120-0.4580470056E2*t251*t102*t250+8.0*t276*t111-0.6177934948E1*t185*
t93*t87+0.5575569085*t83*t467*t91-0.5045677662E1*t64*t131*t87+0.2748282034E3*
t141*t182-3.0*t494*t477-36.0*t64*t121*t29*t18*t250;
    t508 = pow(r,-0.209147081E1);
    t516 = t64*t64;
    t569 = -6.0*t478*t43+4.0*t478*t52+4.0*t481*t111-0.6177934948E1*t372*t114*
t36+2.0*t121*t373*t120-0.3088967474E1*t38*t483*t87-t64*t373*t87+0.8409462769*
t372*t130*t91-6.0*t486*t43+4.0*t486*t52+4.0*t489*t111-0.763411676E1*t379*t140*
t36+2.0*t121*t380*t120-0.381705838E1*t38*t491*t87-t64*t380*t87+0.1733954479E1*
t379*t45*t91;
    t613 = ((0.5589678671E-1*t74-0.7387621778E1*t77+0.1029263491*t79
-0.4012036966E1*t81)*t71*t91+(0.3044446303E-3*t20-0.5769968494E2*t23
-0.1390645576E-2*t25+0.1263180175E2*t27)*t18*t91)*t367+0.4580470056E2*t319*t101
*t120-0.9160940112E2*t147*t53+0.4161490749E2*t47*t52-120.0*t136-0.1040372687E2*
t64*t151*t87-0.9160940112E2*t324*t232+(-2.0*t382*t477+t494*t70)*t68+48.0*t252*
t442-72.0*t144*t385+48.0*t144*t400;
    return(6.0*t382*t367-0.7413521939E2*t125*t53+24.0*t428*t30*t427-36.0*t64*
t121*t83*t71*t250+2.0*t569*t70+8.0*t38*t185*t83*t261+0.3064517643E1*t83*t201*
t44*t36+0.2224056581E3*t115*t182+300.0*t84*t225+0.2018271065E2*t281*t111+48.0*
t284*t442+24.0*t428*t88*t427-0.3027406596E2*t289*t43+t502+0.2018271065E2*t289*
t52+0.1374141017E3*t147*t160+0.4161490749E2*t327*t111-24.0*t128*t53+2.0*t105*
t477+36.0*t128*t160-48.0*t233*t88*t232-t229*t88*t87+6.0*t516*t30*t120+
0.2080745374E2*t121*t151*t120-0.1482704387E3*t115*t178+0.1112028291E3*t125*t160
-72.0*t122*t385+0.1009135532E2*t121*t131*t120+0.1731140602*t29*t508*t91-120.0*
t109+t349+6.0*t516*t88*t120+t415+t613);
  }
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

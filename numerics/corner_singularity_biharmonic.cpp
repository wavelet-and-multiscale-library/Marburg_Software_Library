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
    for(int m=0;m<2;m++)    
	  res+= zeta(r) * pow(r, 1.0+rootz[m]) * ((1/(rootz[m]-1)*sin((rootz[m]-1)*omega)
					     -1/(rootz[m]+1)*sin((rootz[m]+1)*omega))
					    *(cos((rootz[m]-1)*theta)-cos((rootz[m]+1)*theta))
						  -(1/(rootz[m]-1)*sin((rootz[m]-1)*theta)
					      -1/(rootz[m]+1)*sin((rootz[m]+1)*theta))
					    *(cos((rootz[m]-1)*omega)-cos((rootz[m]+1)*omega)));
    return res;
#else
    //Maple output

  double t28;
  double t8;
  double t12;
  double t13;
  double t2;
  double t11;
  double t6;
  double t31;
  double t30;
  double t33;
  double t4;
  double t35;
  double t27;
  double t26;
  double t16;
  double t19;
  double t21;
  double t17;
  double t14;
  {
    t2 = pow(0.99-r,2.0);
    t4 = exp(-1/t2);
    t6 = pow(r-0.1E-1,2.0);
    t8 = exp(-1/t6);
    t11 = t4/(t8+t4);
    t12 = pow(r,0.1544483737E1);
    t13 = 0.4555162632*theta;
    t14 = cos(t13);
    t16 = 0.1544483737E1*theta;
    t17 = cos(t16);
    t19 = sin(t13);
    t21 = sin(t16);
    t26 = pow(r,0.190852919E1);
    t27 = 0.914708102E-1*theta;
    t28 = cos(t27);
    t30 = 0.190852919E1*theta;
    t31 = cos(t30);
    t33 = sin(t27);
    t35 = sin(t30);
    return(t11*t12*(0.1298288751E1*t14-0.1298288751E1*t17+0.2390622594E1*t19
-0.7050689139*t21)+t11*t26*(0.434888792E1*t28-0.434888792E1*t31-0.1986489872E2*
t33+0.9520726166*t35));
  }


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

    double t46;
  double t60;
  double t506;
  double t80;
  double t76;
  double t83;
  double t265;
  double t88;
  double t268;
  double t69;
  double t108;
  double t89;
  double t101;
  double t100;
  double t102;
  double t578;
  double t106;
  double t107;
  double t24;
  double t87;
  double t26;
  double t119;
  double t90;
  double t91;
  double t122;
  double t114;
  double t136;
  double t135;
  double t146;
  double t144;
  double t141;
  double t140;
  double t149;
  double t150;
  double t319;
  double t96;
  double t126;
  double t94;
  double t93;
  double t127;
  double t128;
  double t131;
  double t57;
  double t58;
  double t153;
  double t154;
  double t157;
  double t158;
  double t163;
  double t98;
  double t195;
  double t21;
  double t171;
  double t170;
  double t169;
  double t192;
  double t188;
  double t202;
  double t203;
  double t298;
  double t403;
  double t222;
  double t223;
  double t226;
  double t11;
  double t13;
  double t17;
  double t16;
  double t227;
  double t230;
  double t231;
  double t232;
  double t54;
  double t3;
  double t52;
  double t235;
  double t240;
  double t238;
  double t242;
  double t241;
  double t38;
  double t239;
  double t326;
  double t251;
  double t198;
  double t249;
  double t4;
  double t2;
  double t1;
  double t118;
  double t10;
  double t252;
  double t9;
  double t8;
  double t19;
  double t18;
  double t15;
  double t14;
  double t45;
  double t22;
  double t32;
  double t28;
  double t37;
  double t29;
  double t279;
  double t36;
  double t35;
  double t34;
  double t275;
  double t278;
  double t280;
  double t283;
  double t286;
  double t289;
  double t294;
  double t59;
  double t327;
  double t330;
  double t331;
  double t333;
  double t340;
  double t341;
  double t344;
  double t345;
  double t347;
  double t353;
  double t306;
  double t316;
  double t383;
  double t388;
  double t392;
  double t534;
  double t405;
  double t446;
  double t468;
  double t507;
  double t43;
  double t481;
  double t453;
  double t516;
  double t542;
  double t41;
  {
    t1 = 0.99-r;
    t2 = t1*t1;
    t3 = t2*t2;
    t4 = t3*t3;
    t8 = exp(-1/t2);
    t9 = 1/t4/t1*t8;
    t10 = r-0.1E-1;
    t11 = t10*t10;
    t13 = exp(-1/t11);
    t14 = t13+t8;
    t15 = 1/t14;
    t16 = pow(r,0.90852919);
    t17 = t15*t16;
    t18 = 0.914708102E-1*theta;
    t19 = cos(t18);
    t21 = 0.190852919E1*theta;
    t22 = cos(t21);
    t24 = sin(t18);
    t26 = sin(t21);
    t28 = 0.434888792E1*t19-0.434888792E1*t22-0.1986489872E2*t24+0.9520726166*
t26;
    t29 = t17*t28;
    t32 = t14*t14;
    t34 = 1/t32/t14;
    t35 = t8*t34;
    t36 = pow(r,-0.9147081E-1);
    t37 = t36*t28;
    t38 = t11*t10;
    t41 = t2*t1;
    t43 = 1/t41*t8;
    t45 = 2.0/t38*t13-2.0*t43;
    t46 = t45*t45;
    t52 = 1/t3/t2*t8;
    t54 = t15*t36*t28;
    t57 = pow(r,0.190852919E1);
    t58 = t35*t57;
    t59 = t28*t45;
    t60 = t11*t11;
    t69 = t60*t60;
    t76 = 1/t3/t1*t8;
    t80 = 1/t3/t41*t8;
    t83 = 24.0/t60/t10*t13-36.0/t60/t38*t13+8.0/t69/t10*t13-24.0*t76+36.0*t80
-8.0*t9;
    t87 = 1/t32;
    t88 = t52*t87;
    t89 = pow(r,0.544483737);
    t90 = 0.4555162632*theta;
    t91 = cos(t90);
    t93 = 0.1544483737E1*theta;
    t94 = cos(t93);
    t96 = sin(t90);
    t98 = sin(t93);
    t100 = 0.1298288751E1*t91-0.1298288751E1*t94+0.2390622594E1*t96
-0.7050689139*t98;
    t101 = t89*t100;
    t102 = t101*t45;
    t106 = 1/t3*t8;
    t107 = t106*t87;
    t108 = t57*t28;
    t114 = 1/t60/t11*t13;
    t118 = -6.0/t60*t13+4.0*t114-6.0*t106+4.0*t52;
    t119 = t108*t118;
    t122 = t118*t118;
    t126 = pow(r,0.1544483737E1);
    t127 = t126*t100;
    t128 = t127*t118;
    t131 = t8*t87;
    t135 = pow(r,-0.109147081E1);
    t136 = t135*t28;
    t140 = t43*t87;
    t141 = t37*t45;
    t144 = pow(r,-0.455516263);
    t146 = t15*t144*t100;
    t149 = t9*t87;
    t150 = t127*t45;
    t153 = t15*t89;
    t154 = t153*t100;
    t157 = t76*t87;
    t158 = t108*t45;
    t163 = t16*t28;
    t169 = 1/t4/t3*t8;
    t170 = t15*t57;
    t171 = t170*t28;
    t188 = 1/t4*t8;
    t192 = 1/t4/t2*t8;
    t195 = -120.0*t114+300.0/t69*t13-144.0/t69/t11*t13+16.0/t69/t60*t13-120.0*
t52+300.0*t188-144.0*t192+16.0*t169;
    t198 = -0.6107293408E2*t9*t29+0.2080745374E2*t35*t37*t46+0.4161490749E2*t52
*t54+8.0*t58*t59*t83-0.7413521939E2*t88*t102+36.0*t107*t119+6.0*t35*t108*t122
-24.0*t88*t128-0.1040372687E2*t131*t37*t118+0.6344248828*t131*t136*t45+
0.4161490749E2*t140*t141+0.2018271065E2*t52*t146+32.0*t149*t150-0.494234796E2*
t9*t154+96.0*t157*t158-0.6242236123E2*t106*t54-0.763411676E1*t131*t163*t83+16.0
*t169*t171-t131*t108*t195;
    t202 = t15*t126;
    t203 = t202*t100;
    t222 = t43*t34;
    t223 = t127*t46;
    t226 = t144*t100;
    t227 = t226*t45;
    t230 = t35*t126;
    t231 = t100*t45;
    t232 = t231*t118;
    t235 = t101*t46;
    t238 = t32*t32;
    t239 = 1/t238;
    t240 = t8*t239;
    t241 = t46*t45;
    t242 = t127*t241;
    t249 = t83*t127;
    t251 = 18.0*t107*t150-8.0*t9*t203+0.1853380485E2*t52*t154-12.0*t88*t150+
36.0*t80*t203-24.0*t76*t203-0.2780070726E2*t106*t154-0.5045677661E1*t43*t146+
6.0*t140*t128+0.1853380485E2*t140*t102-12.0*t222*t223-0.2522838831E1*t131*t227+
6.0*t230*t232+0.9266902422E1*t35*t235-6.0*t240*t242-0.3435352542E2*t106*t29+
18.0*t107*t158-t131*t249;
    t252 = t101*t118;
    t265 = t163*t45;
    t268 = t108*t46;
    t275 = t59*t118;
    t278 = t8*t15;
    t279 = pow(r,-0.1455516263E1);
    t280 = t279*t100;
    t283 = t163*t46;
    t286 = t108*t241;
    t289 = t163*t118;
    t294 = t108*t83;
    t298 = -0.4633451211E1*t131*t252-12.0*t88*t158+36.0*t80*t171-24.0*t76*t171
-8.0*t9*t171+0.2290235028E2*t52*t29+0.2290235028E2*t140*t265-12.0*t222*t268
-0.1040372687E2*t43*t54+6.0*t140*t119+6.0*t58*t275-0.3830647054*t278*t280+
0.1145117514E2*t35*t283-6.0*t240*t286-0.572558757E1*t131*t289-0.5201863437E1*
t131*t141-t131*t294-0.1586062207*t278*t136;
    t306 = pow(r,-0.2455516263E1);
    t316 = t80*t87;
    t319 = r*r;
    t326 = -0.2693885101*t91+0.3096976954E1*t94-0.4960423931*t96+0.1681892549E1
*t98;
    t327 = t202*t326;
    t330 = t126*t326;
    t331 = t330*t45;
    t333 = t89*t326;
    t340 = -0.3638674999E-1*t19+0.1584075323E2*t22+0.1662078022*t24
-0.3467908958E1*t26;
    t341 = t170*t340;
    t344 = t57*t340;
    t345 = t344*t45;
    t347 = t16*t340;
    t353 = t319*t319;
    t383 = t106*t34;
    t388 = t45*t118;
    t392 = t43*t239;
    t403 = 2.0/r*(t251+t298)+8.0*t140*t294+0.4580470056E2*t140*t289+
0.5575569085*t278*t306*t100+0.1268849765E1*t43*t15*t135*t28+0.1374141017E3*t107
*t265-144.0*t316*t158+2.0/t319/r*(-2.0*t43*t327-t131*t331+0.1544483737E1*t278*
t333-2.0*t43*t341-t131*t345+0.190852919E1*t278*t347)+1/t353*(t278*t126*(
0.5589678671E-1*t91-0.7387621778E1*t94+0.1029263491*t96-0.4012036966E1*t98)+
t278*t57*(0.3044446303E-3*t19-0.5769968494E2*t22-0.1390645576E-2*t24+
0.1263180175E2*t26))-36.0*t240*t57*t28*t46*t118-0.6177934948E1*t131*t101*t83-
t131*t127*t195+16.0*t169*t203-72.0*t383*t268-0.7413521938E2*t222*t235-48.0*t222
*t108*t388+48.0*t392*t242-48.0*t222*t127*t388-0.4580470056E2*t240*t163*t241
-0.1832188023E3*t76*t29;
    t405 = 1/t319;
    t446 = -6.0*t106*t327+4.0*t52*t327+4.0*t140*t331-0.6177934948E1*t43*t153*
t326+2.0*t35*t330*t46-0.3088967474E1*t131*t333*t45-t131*t330*t118+0.8409462769*
t278*t144*t326-6.0*t106*t341+4.0*t52*t341+4.0*t140*t345-0.763411676E1*t43*t17*
t340+2.0*t35*t344*t46-0.381705838E1*t131*t347*t45-t131*t344*t118+0.1733954479E1
*t278*t36*t340;
    t453 = t52*t203;
    t468 = t52*t171;
    t481 = -6.0*t106*t203+4.0*t453+4.0*t140*t150-0.6177934948E1*t154*t43+2.0*
t35*t223-0.3088967474E1*t131*t102-t131*t128+0.8409462769*t278*t226-6.0*t106*
t171+4.0*t468+4.0*t140*t158-0.763411676E1*t43*t29+2.0*t35*t268-0.381705838E1*
t131*t265-t131*t119+0.1733954479E1*t278*t37;
    t506 = t8/t238/t14;
    t507 = t46*t46;
    t516 = pow(r,-0.209147081E1);
    t534 = 2.0*t405*t446+0.2018271065E2*t140*t227+t405*t481+0.3064517643E1*t43*
t15*t279*t100-0.9160940112E2*t222*t283+8.0*t140*t249+0.3706760969E2*t140*t252
-36.0*t240*t126*t100*t46*t118+0.2748282033E3*t80*t29-0.3706760969E2*t240*t101*
t241-120.0*t453+24.0*t506*t127*t507+48.0*t392*t286+0.1009135532E2*t35*t226*t46+
0.1731140602*t278*t516*t28+24.0*t506*t108*t507+6.0*t35*t127*t122+0.3706760968E2
*t35*t89*t232+8.0*t230*t231*t83+300.0*t188*t171;
    t542 = t52*t34;
    t578 = -0.9160940112E2*t88*t265+0.1532258821E1*t131*t280*t45-144.0*t192*
t171+48.0*t542*t268-120.0*t468+32.0*t149*t158-24.0*t88*t119-72.0*t383*t223
-0.5045677662E1*t131*t226*t118+36.0*t107*t128+0.1112028291E3*t107*t102
-0.1482704387E3*t76*t154-0.3027406596E2*t106*t146-144.0*t316*t150+
0.4580470056E2*t35*t16*t275+96.0*t157*t150+0.2224056581E3*t80*t154+48.0*t542*
t223-144.0*t192*t203+300.0*t188*t203;
    return(t198+t403+t534+t578);
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

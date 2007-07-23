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
  double t314;
  double t8;
  double t7;
  double t87;
  double t334;
  double t6;
  double t24;
  double t131;
  double t461;
  double t424;
  double t258;
  double t279;
  double t287;
  double t207;
  double t574;
  double t420;
  double t1373;
  double t244;
  double t417;
  double t1312;
  double t412;
  double t413;
  double t468;
  double t2637;
  double t414;
  double t1434;
  double t469;
  double t472;
  double t195;
  double t1713;
  double t458;
  double t1484;
  double t88;
  double t331;
  double t409;
  double t324;
  double t228;
  double t477;
  double t457;
  double t866;
  double t867;
  double t456;
  double t1883;
  double t321;
  double t402;
  double t322;
  double t246;
  double t1266;
  double t191;
  double t319;
  double t165;
  double t205;
  double t72;
  double t89;
  double t1539;
  double t240;
  double t1404;
  double t277;
  double t396;
  double t397;
  double t226;
  double t256;
  double t395;
  double t2648;
  double t427;
  double t49;
  double t187;
  double t317;
  double t799;
  double t1441;
  double t2649;
  double t35;
  double t493;
  double t158;
  double t117;
  double t492;
  double t391;
  double t1407;
  double t911;
  double t914;
  double t390;
  double t185;
  double t519;
  double t104;
  double t224;
  double t73;
  double t448;
  double t184;
  double t164;
  double t1500;
  double t312;
  double t203;
  double t149;
  double t222;
  double t238;
  double t500;
  double t1129;
  double t765;
  double t171;
  double t385;
  double t311;
  double t77;
  double t2654;
  double t124;
  double t2212;
  double t384;
  double t308;
  double t954;
  double t1004;
  double t383;
  double t2655;
  double t744;
  double t51;
  double t275;
  double t145;
  double t254;
  double t90;
  double t182;
  double t453;
  double t2656;
  double t505;
  double t181;
  double t306;
  double t37;
  double t332;
  double t2657;
  double t694;
  double t510;
  double t1648;
  double t179;
  double t443;
  double t305;
  double t189;
  double t91;
  double t513;
  double t378;
  double t304;
  double t303;
  double t99;
  double t65;
  double t892;
  double t143;
  double t377;
  double t1797;
  double t516;
  double t221;
  double t236;
  double t273;
  double t160;
  double t376;
  double t301;
  double t309;
  double t440;
  double t524;
  double t520;
  double t122;
  double t1505;
  double t695;
  double t252;
  double t66;
  double t527;
  double t219;
  double t819;
  double t289;
  double t177;
  double t175;
  double t53;
  double t531;
  double t530;
  double t710;
  double t26;
  double t1337;
  double t92;
  double t539;
  double t265;
  double t1183;
  double t125;
  double t534;
  double t1546;
  double t373;
  double t113;
  double t84;
  double t93;
  double t1299;
  double t666;
  double t1142;
  double t430;
  double t310;
  double t623;
  double t722;
  double t100;
  double t234;
  double t39;
  double t1096;
  double t2046;
  double t1256;
  double t349;
  double t437;
  double t436;
  double t550;
  double t1186;
  double t403;
  double t919;
  double t1354;
  double t31;
  double t367;
  double t366;
  double t298;
  double t126;
  double t365;
  double t201;
  double t141;
  double t1848;
  double t79;
  double t140;
  double t156;
  double t297;
  double t329;
  double t2377;
  double t108;
  double t1114;
  double t872;
  double t875;
  double t764;
  double t1233;
  double t173;
  double t558;
  double t561;
  double t1307;
  double t2411;
  double t1317;
  double t18;
  double t362;
  double t904;
  double t361;
  double t654;
  double t392;
  double t484;
  double t55;
  double t295;
  double t1191;
  double t1151;
  double t1511;
  double t1079;
  double t824;
  double t85;
  double t266;
  double t599;
  double t119;
  double t1680;
  double t28;
  double t1815;
  double t2;
  double t359;
  double t464;
  double t1194;
  double t111;
  double t776;
  double t380;
  double t250;
  double t215;
  double t1047;
  double t358;
  double t2278;
  double t41;
  double t102;
  double t569;
  double t1035;
  double t1038;
  double t80;
  double t1017;
  double t609;
  double t20;
  double t1019;
  double t1072;
  double t169;
  double t81;
  double t657;
  double t354;
  double t984;
  double t989;
  double t353;
  double t291;
  double t232;
  double t691;
  double t348;
  double t110;
  double t154;
  double t1108;
  double t465;
  double t209;
  double t42;
  double t579;
  double t1026;
  double t346;
  double t650;
  double t653;
  double t634;
  double t2310;
  double t1200;
  double t57;
  double t1491;
  double t147;
  double t1363;
  double t151;
  double t95;
  double t585;
  double t217;
  double t2540;
  double t431;
  double t3;
  double t2013;
  double t1481;
  double t1056;
  double t262;
  double t136;
  double t199;
  double t1164;
  double t82;
  double t128;
  double t267;
  double t96;
  double t213;
  double t725;
  double t729;
  double t1702;
  double t1916;
  double t2080;
  double t44;
  double t19;
  double t1979;
  double t590;
  double t342;
  double t4;
  double t134;
  double t591;
  double t1059;
  double t1060;
  double t285;
  double t833;
  double t836;
  double t260;
  double t963;
  double t968;
  double t971;
  double t594;
  double t491;
  double t129;
  double t646;
  double t68;
  double t557;
  double t339;
  double t855;
  double t5;
  double t22;
  double t1998;
  double t338;
  double t1868;
  double t931;
  double t230;
  double t33;
  double t211;
  double t600;
  double t86;
  double t248;
  double t2245;
  double t2113;
  double t1844;
  double t268;
  double t197;
  double t336;
  double t283;
  double t152;
  double t335;
  double t2507;
  double t282;
  double t523;
  double t1368;
  double t1810;
  double t69;
  double t59;
  double t713;
  double t895;
  double t167;
  double t2148;
  double t193;
  double t618;
  double t120;
  double t162;
  double t1579;
  double t2474;
  double t1947;
  double t300;
  double t2181;
  double t749;
  double t2571;
  double t75;
  double t293;
  double t612;
  double t1391;
  double t97;
  double t582;
  double t10;
  double t11;
  double t2443;
  double t46;
  double t1320;
  double t302;
  double t61;
  double t132;
  double t2344;
  double t13;
  double t1067;
  double t70;
  double t734;
  double t737;
  double t803;
  double t806;
  double t14;
  double t685;
  double t688;
  double t1779;
  double t1615;
  double t848;
  double t841;
  double t269;
  double t1746;
  double t242;
  double t280;
  double t15;
  double t566;
  double t63;
  double t2605;
  double t47;
  double t327;
  double t480;
  double t71;
  double t115;
  double t30;
  double t635;
  double t326;
  double t106;
  double t1474;
  double t315;
  double t16;
  double t274;
  double t946;
  t2 = 0.154448373678246e1 * theta;
  t3 = sin(t2);
  t4 = r * r;
  t5 = r * t4;
  t6 = t4 * t4;
  t7 = t5 * t6;
  t8 = t3 * t7;
  t10 = 0.455516263217536e0 * theta;
  t11 = sin(t10);
  t13 = cos(t10);
  t14 = t6 * t6;
  t15 = r * t14;
  t16 = t13 * t15;
  t18 = r * t6;
  t19 = pow(r, 0.1e1 / 0.25000000000000e14);
  t20 = t19 * t19;
  t22 = t20 * t20;
  t24 = t22 * t22;
  t26 = t24 * t24;
  t28 = t26 * t26;
  t30 = t28 * t28;
  t31 = t30 * t30;
  t33 = t31 * t31;
  t35 = t33 * t33;
  t37 = t35 * t35;
  t39 = t37 * t37;
  t41 = t39 * t39;
  t42 = t41 * t41;
  t44 = t42 * t42;
  t46 = t44 * t44;
  t47 = t46 * t46;
  t49 = t47 * t47;
  t51 = t49 * t49;
  t53 = t51 * t51;
  t55 = t53 * t53;
  t57 = t55 * t55;
  t59 = t57 * t57;
  t61 = t59 * t59;
  t63 = t61 * t61;
  t65 = t63 * t63;
  t66 = t65 * t65;
  t68 = t66 * t66;
  t69 = t68 * t68;
  t70 = t69 * t69;
  t71 = t70 * t70;
  t72 = t71 * t71;
  t73 = t72 * t72;
  t75 = t73 * t73;
  t77 = t75 * t75;
  t79 = t77 * t77;
  t80 = t79 * t79;
  t81 = t80 * t80;
  t82 = t81 * t81;
  t84 = t82 * t82;
  t85 = t84 * t84;
  t86 = t85 * t85;
  t87 = t86 * t86;
  t88 = t87 * t87;
  t89 = t19 * t20 * t22 * t24 * t26 * t28 * t31 * t33 * t35 * t37 * t39 * t42 * t44 * t47 * t49 * t51 * t53 * t55 * t57 * t59 * t61 * t63 * t66 * t73 * t75 * t77 * t82 * t88;
  t90 = t89 * t18;
  t91 = 0.914708101539010e-1 * theta;
  t92 = sin(t91);
  t93 = t90 * t92;
  t95 = 0.190852918984610e1 * theta;
  t96 = cos(t95);
  t97 = t90 * t96;
  t99 = cos(t2);
  t100 = t99 * t15;
  t102 = t11 * r;
  t104 = t11 * t4;
  t106 = t11 * t6;
  t108 = t11 * t18;
  t110 = t4 * t6;
  t111 = t11 * t110;
  t113 = t11 * t7;
  t115 = t11 * t14;
  t117 = t11 * t15;
  t119 = sin(t95);
  t120 = t89 * t119;
  t122 = t18 * t99;
  t124 = t89 * t15;
  t125 = cos(t91);
  t126 = t124 * t125;
  t128 = t89 * r;
  t129 = t128 * t125;
  t131 = t89 * t110;
  t132 = t131 * t92;
  t134 = t5 * t3;
  t136 = 0.493709408956622e16 * t8 - 0.1732572320e16 * t11 + 0.734530333379934e15 * t16 - 0.9606796772e18 * t93 - 0.2931478419e19 * t97 - 0.6050226304e15 * t100 + 0.157506574530711e17 * t102 - 0.6363902001e17 * t104 + 0.340752905152216e18 * t106 - 0.5732014168e17 * t108 - 0.1545822022e18 * t111 + 0.669187022714911e17 * t113 - 0.1689866219e17 * t115 + 0.189659508302699e16 * t117 - 0.6533149021e15 * t120 + 0.825656773610300e17 * t122 + 0.292197523761556e15 * t126 + 0.242661343298357e16 * t129 + 0.111840894514604e18 * t132 + 0.939905220379247e17 * t134;
  t140 = t89 * t7;
  t141 = t140 * t125;
  t143 = t3 * t15;
  t145 = t13 * t14;
  t147 = t5 * t99;
  t149 = t140 * t92;
  t151 = t89 * t14;
  t152 = t151 * t92;
  t154 = t128 * t92;
  t156 = t140 * t119;
  t158 = t151 * t119;
  t160 = t124 * t119;
  t162 = t128 * t119;
  t164 = t89 * t6;
  t165 = t164 * t92;
  t167 = t164 * t96;
  t169 = t131 * t96;
  t171 = t140 * t96;
  t173 = t151 * t96;
  t175 = t124 * t96;
  t177 = -0.1278248422e15 * t3 - 0.6710061283e15 * t13 + 0.552698608019387e15 * t99 + 0.103097805493932e17 * t141 + 0.139926030494797e15 * t143 - 0.6544665270e16 * t145 + 0.104846525834210e18 * t147 - 0.4841597165e17 * t149 + 0.122262554675110e17 * t152 - 0.1139566906e17 * t154 + 0.252335702945362e17 * t156 - 0.6372113711e16 * t158 + 0.715164277299828e15 * t160 + 0.593922638270515e16 * t162 + 0.252437285919965e19 * t165 + 0.324659555733038e19 * t167 + 0.185691777908538e19 * t169 - 0.8038605104e18 * t171 + 0.202995078390548e18 * t173 - 0.2278283708e17 * t175;
  t179 = t128 * t96;
  t181 = t89 * t4;
  t182 = t181 * t119;
  t184 = t89 * t5;
  t185 = t184 * t119;
  t187 = t164 * t119;
  t189 = t90 * t119;
  t191 = t131 * t119;
  t193 = t184 * t92;
  t195 = t164 * t125;
  t197 = t90 * t125;
  t199 = t184 * t125;
  t201 = t131 * t125;
  t203 = t151 * t125;
  t205 = t6 * t99;
  t207 = t6 * t13;
  t209 = t13 * r;
  t211 = t13 * t4;
  t213 = t13 * t110;
  t215 = t13 * t7;
  t217 = t5 * t13;
  t219 = t5 * t11;
  t221 = -0.1892046784e18 * t179 - 0.6094867859e17 * t182 + 0.168533454120764e18 * t185 - 0.1988006049e18 * t187 + 0.124642882525077e18 * t189 - 0.5828954738e17 * t191 - 0.2444869577e19 * t193 - 0.5516603786e18 * t195 + 0.209321672343304e18 * t197 + 0.534589651703281e18 * t199 - 0.2381559307e17 * t201 - 0.2603479937e16 * t203 - 0.2359767244e18 * t205 + 0.220458496716970e18 * t207 + 0.610005571220435e16 * t209 - 0.2464668975e17 * t211 - 0.5986798003e17 * t213 + 0.259168744708442e17 * t215 - 0.9460449559e17 * t217 - 0.1311749490e18 * t219;
  t222 = t181 * t92;
  t224 = t181 * t125;
  t226 = t18 * t13;
  t228 = t99 * r;
  t230 = t89 * t125;
  t232 = t3 * r;
  t234 = t3 * t4;
  t236 = t89 * t92;
  t238 = t184 * t96;
  t240 = t181 * t96;
  t242 = t89 * t96;
  t244 = t3 * t6;
  t246 = t3 * t18;
  t248 = t3 * t110;
  t250 = t3 * t14;
  t252 = t99 * t4;
  t254 = t99 * t110;
  t256 = t99 * t7;
  t258 = t99 * t14;
  t260 = t124 * t92;
  t262 = 0.817038761884413e18 * t222 - 0.1785933626e18 * t224 - 0.6689069994e17 * t226 - 0.5024532800e16 * t228 - 0.2669274776e15 * t230 + 0.116204402026290e16 * t232 - 0.4695127355e16 * t234 + 0.125352359615179e16 * t236 - 0.2313246870e19 * t238 + 0.933252210759300e18 * t240 + 0.208125146199842e17 * t242 - 0.1842909755e18 * t244 + 0.101544230590683e18 * t246 - 0.1140468735e17 * t248 - 0.1246740932e16 * t250 + 0.203011426269747e17 * t252 + 0.493124396761490e17 * t254 - 0.2134737648e17 * t256 + 0.539075163716701e16 * t258 - 0.1372194777e16 * t260;
  t265 = 0.100e3 * r;
  t266 = -0.99e2 + t265;
  t267 = t266 * t266;
  t268 = t267 * t267;
  t269 = t268 * t268;
  t273 = pow(r, 0.1e1 / 0.50000000000000e14);
  t274 = t273 * t273;
  t275 = t274 * t274;
  t277 = t275 * t275;
  t279 = t277 * t277;
  t280 = t279 * t279;
  t282 = t280 * t280;
  t283 = t282 * t282;
  t285 = t283 * t283;
  t287 = t285 * t285;
  t289 = t287 * t287;
  t291 = t289 * t289;
  t293 = t291 * t291;
  t295 = t293 * t293;
  t297 = t295 * t295;
  t298 = t297 * t297;
  t300 = t298 * t298;
  t301 = t300 * t300;
  t302 = t301 * t301;
  t303 = t302 * t302;
  t304 = t303 * t303;
  t305 = t304 * t304;
  t306 = t305 * t305;
  t308 = t306 * t306;
  t309 = t308 * t308;
  t310 = t309 * t309;
  t311 = t310 * t310;
  t312 = t311 * t311;
  t314 = t312 * t312;
  t315 = t314 * t314;
  t317 = t315 * t315;
  t319 = t317 * t317;
  t321 = t319 * t319;
  t322 = t321 * t321;
  t324 = t322 * t322;
  t326 = t324 * t324;
  t327 = t326 * t326;
  t329 = t327 * t327;
  t331 = t329 * t329;
  t332 = t331 * t331;
  t334 = t332 * t332;
  t335 = t334 * t334;
  t336 = t335 * t335;
  t338 = t336 * t336;
  t339 = t338 * t338;
  t342 = 0.1e1 / t273 / t275 / t277 / t280 / t283 / t285 / t287 / t289 / t291 / t293 / t295 / t298 / t306 / t312 / t315 / t317 / t319 / t322 / t324 / t327 / t329 / t332 / t336 / t339 / t4;
  t346 = 0.1e1 / t267;
  t348 = exp(-0.10000e5 * t346);
  t349 = 0.10000e5 * t4;
  t353 = pow(t265 - 0.1e1, 0.2e1);
  t354 = 0.1e1 / t353;
  t358 = exp(-0.40000e5 * (t349 - 0.5100e4 * r + 0.2451e4) * t354 * t346);
  t359 = t13 * t358;
  t361 = t7 * t14;
  t362 = t361 * t99;
  t365 = t14 * t14;
  t366 = t18 * t365;
  t367 = t366 * t3;
  t373 = exp(-0.40000e5 * (0.7351e4 - 0.14900e5 * r + t349) * t354 * t346);
  t376 = t14 * t365;
  t377 = t89 * t376;
  t378 = t377 * t125;
  t380 = exp(-0.40000e5 * t346);
  t383 = t18 * t14;
  t384 = t89 * t383;
  t385 = t384 * t96;
  t390 = t4 * t365;
  t391 = t89 * t390;
  t392 = t391 * t119;
  t395 = t110 * t14;
  t396 = t89 * t395;
  t397 = t396 * t96;
  t402 = t5 * t14;
  t403 = t402 * t99;
  t409 = exp(-0.40000e5 * (t349 - 0.10000e5 * r + 0.4901e4) * t354 * t346);
  t412 = t5 * t365;
  t413 = t89 * t412;
  t414 = t413 * t125;
  t417 = t3 * t409;
  t420 = t395 * t11;
  t424 = exp(-0.40000e5 * t354);
  t427 = t391 * t125;
  t430 = 0.260430630129832e16 * t359 - 0.5039616543e58 * t362 * t358 - 0.1057303063e55 * t367 * t373 + 0.292197523761556e39 * t378 * t380 - 0.4398382562e59 * t385 * t373 - 0.3406458885e57 * t93 * t358 - 0.4776309285e57 * t392 * t358 - 0.7167729828e57 * t397 * t358 + 0.212136122393778e59 * t132 * t373 - 0.1533797525e59 * t403 * t409 + 0.199488608067453e57 * t414 * t373 + 0.839556182968782e39 * t417 * t376 - 0.2145445896e57 * t420 * t373 - 0.2153688088e59 * t152 * t424 - 0.1486608892e58 * t427 * t409;
  t431 = t383 * t13;
  t436 = t89 * t366;
  t437 = t436 * t125;
  t440 = t391 * t92;
  t443 = t412 * t3;
  t448 = t11 * t380;
  t453 = t390 * t99;
  t456 = t6 * t365;
  t457 = t89 * t456;
  t458 = t457 * t125;
  t461 = t361 * t3;
  t464 = t13 * t380;
  t465 = t7 * t365;
  t468 = t4 * t14;
  t469 = t468 * t13;
  t472 = t396 * t92;
  t477 = t456 * t3;
  t480 = 0.100965960838548e59 * t431 * t373 + 0.104618347277984e38 * t179 * t358 - 0.6466529895e41 * t437 * t380 - 0.7794278556e57 * t440 * t424 - 0.1100432030e57 * t443 * t424 - 0.7993905386e34 * t152 * t380 - 0.4197293887e42 * t448 * t366 + 0.213400309312559e43 * t234 * t358 + 0.161583356516153e58 * t453 * t373 + 0.221406402623160e56 * t458 * t358 + 0.211470323161505e58 * t461 * t424 - 0.8814364001e40 * t464 * t465 + 0.624990340163402e37 * t469 * t380 + 0.363331112326290e59 * t472 * t358 - 0.6871016147e37 * t162 * t373 + 0.174230482490933e56 * t477 * t424;
  t484 = t384 * t119;
  t491 = r * t365;
  t492 = t89 * t491;
  t493 = t492 * t92;
  t500 = t390 * t11;
  t505 = t384 * t125;
  t510 = t492 * t125;
  t513 = t456 * t11;
  t516 = t365 * t11;
  t519 = t89 * t361;
  t520 = t519 * t125;
  t523 = -0.2490141063e32 * t141 * t380 - 0.2674107925e58 * t484 * t358 - 0.1798980908e57 * t189 * t409 - 0.2736892671e58 * t461 * t358 - 0.1594574105e59 * t493 * t358 + 0.252448075863004e59 * t472 * t373 + 0.734530333379934e39 * t464 * t376 - 0.1156595550e43 * t500 * t380 - 0.1211661727e18 * t242 * t409 - 0.3242279453e41 * t505 * t380 - 0.1334706529e57 * t246 * t373 + 0.458954143064850e58 * t510 * t409 + 0.456912960209569e56 * t513 * t424 + 0.580811039895425e58 * t516 * t424 - 0.1383568893e59 * t520 * t424;
  t524 = t390 * t13;
  t527 = t436 * t96;
  t530 = t89 * t365;
  t531 = t530 * t92;
  t534 = t457 * t92;
  t539 = t391 * t96;
  t550 = t366 * t99;
  t557 = t89 * t402;
  t558 = t557 * t96;
  t561 = t384 * t92;
  t566 = -0.2962762139e57 * t524 * t409 + 0.504199676928859e43 * t527 * t380 + 0.424285272265945e59 * t531 * t409 - 0.1482807161e43 * t534 * t380 - 0.7891837277e18 * t464 * r + 0.140104172837588e58 * t539 * t424 - 0.1980908397e59 * t100 * t373 - 0.4474854304e58 * t169 * t409 + 0.105486029786671e29 * t93 * t380 + 0.772938249814366e58 * t149 * t424 - 0.1946880148e55 * t550 * t409 + 0.338912827796943e55 * t106 * t424 + 0.176435841379338e58 * t132 * t358 + 0.226011530693600e59 * t558 * t424 + 0.124122693242902e60 * t561 * t424 + 0.440269843357361e40 * t558 * t380;
  t569 = t413 * t119;
  t574 = t412 * t99;
  t579 = t557 * t92;
  t582 = t413 * t96;
  t585 = t3 * t380;
  t590 = t6 * t14;
  t591 = t590 * t11;
  t594 = t11 * t358;
  t599 = t89 * t590;
  t600 = t599 * t96;
  t609 = -0.8112574216e42 * t569 * t380 + 0.122666613827378e50 * t134 * t409 - 0.6078888866e57 * t574 * t373 + 0.631587900916018e59 * t126 * t373 - 0.3682466297e60 * t579 * t373 - 0.7417018430e57 * t582 * t424 + 0.154786836792634e42 * t585 * t456 + 0.392947717828996e43 * t252 * t358 + 0.202640251833905e59 * t591 * t424 - 0.9103656399e41 * t594 * t465 + 0.293812133351974e40 * t359 * t376 - 0.6091540287e41 * t600 * t380 + 0.180162093085962e44 * t493 * t380 + 0.664219207869541e56 * t458 * t373 + 0.467333348460534e42 * t385 * t380;
  t612 = t13 * t373;
  t618 = t395 * t13;
  t623 = t365 * t99;
  t634 = t110 * t365;
  t635 = t634 * t3;
  t646 = 0.236290953387731e36 * t102 * t424 + 0.260430630129832e16 * t612 - 0.1016717220e58 * t191 * t373 - 0.1421571812e58 * t117 * t358 - 0.4499347308e58 * t618 * t358 + 0.364144668788061e59 * t117 * t373 + 0.875133878339690e57 * t623 * t373 - 0.2756781471e59 * t126 * t409 + 0.688921669148688e59 * t579 * t358 + 0.162375767792857e57 * t569 * t358 - 0.2462640864e59 * t420 * t409 + 0.370069649003175e41 * t635 * t373 - 0.2290338716e37 * t162 * t358 - 0.8286218801e57 * t427 * t358 + 0.246772416175543e57 * t111 * t424 + 0.244960659639051e58 * t111 * t409;
  t650 = t599 * t125;
  t653 = t89 * t468;
  t654 = t653 * t96;
  t657 = t412 * t11;
  t666 = t99 * t380;
  t685 = -0.1841938336e55 * t207 * t358 - 0.6346971404e59 * t650 * t373 - 0.1697571904e59 * t654 * t424 - 0.1554259548e57 * t657 * t409 + 0.103599785846063e16 * t230 * t373 - 0.1235782229e56 * t527 * t358 + 0.321305175437331e44 * t224 * t373 + 0.457705719612654e28 * t666 * t18 + 0.245181512080547e57 * t122 * t409 - 0.8348937135e59 * t558 * t373 - 0.2026296289e57 * t574 * t424 + 0.161346668370079e57 * t500 * t373 - 0.8212747082e21 * t222 * t380 - 0.2290338716e37 * t162 * t424 + 0.320504399644355e42 * t458 * t380 + 0.749273246245596e56 * t97 * t424;
  t688 = t390 * t3;
  t691 = t653 * t119;
  t694 = t89 * t634;
  t695 = t694 * t119;
  t710 = t13 * t424;
  t713 = t519 * t119;
  t722 = t365 * t3;
  t725 = -0.2282004517e43 * t397 * t380 + 0.106476937322347e58 * t688 * t409 + 0.205162524255020e58 * t691 * t358 + 0.472858038179319e41 * t695 * t424 + 0.374675967341180e59 * t173 * t373 - 0.1041855666e43 * t585 * t491 + 0.752855626992498e49 * t217 * t424 - 0.2037712957e19 * t448 * r + 0.128323846284923e36 * t209 * t358 + 0.136320059687462e26 * t195 * t380 + 0.734530333379934e39 * t710 * t376 + 0.531881061179615e57 * t713 * t373 + 0.202557850248254e56 * t207 * t409 + 0.384971538854770e36 * t209 * t373 - 0.4159026287e50 * t219 * t373 + 0.346195625611024e58 * t722 * t409;
  t729 = t599 * t92;
  t734 = t456 * t99;
  t737 = t491 * t99;
  t744 = t519 * t92;
  t749 = t519 * t96;
  t764 = 0.289917667630530e60 * t729 * t373 - 0.7684097713e58 * t258 * t409 + 0.962465387584557e56 * t734 * t409 - 0.1762881758e56 * t737 * t424 + 0.116879009504622e40 * t378 * t373 - 0.1413424163e44 * t531 * t380 - 0.7922099434e59 * t744 * t409 + 0.125924571457733e60 * t260 * t409 - 0.2429532243e58 * t749 * t373 - 0.3099495764e57 * t165 * t409 - 0.1076851093e59 * t431 * t424 - 0.2339359517e58 * t654 * t409 - 0.6334307439e59 * t175 * t373 - 0.2318573903e59 * t600 * t424 - 0.1100042235e56 * t244 * t409;
  t765 = t530 * t125;
  t776 = t456 * t13;
  t799 = -0.8333042054e58 * t765 * t358 - 0.4426787714e57 * t169 * t424 - 0.4747009125e42 * t713 * t380 - 0.2988021632e42 * t666 * t412 + 0.725915968304717e43 * t749 * t380 + 0.744414811031304e56 * t776 * t409 + 0.320821795861477e56 * t734 * t358 + 0.726027156520810e40 * t666 * t465 + 0.169125029084278e58 * t171 * t424 + 0.925174122507938e40 * t635 * t424 + 0.386239402626190e57 * t169 * t358 - 0.3834297631e24 * t219 * t380 + 0.157174361503155e57 * t524 * t424 - 0.1876989510e60 * t579 * t409 + 0.115031645368518e39 * t666 * t402 + 0.179335671890435e59 * t558 * t358;
  t803 = t653 * t125;
  t806 = t11 * t373;
  t819 = t634 * t99;
  t824 = t412 * t13;
  t833 = t491 * t11;
  t836 = -0.4975236033e59 * t385 * t409 + 0.178881634412719e59 * t803 * t424 - 0.9103656399e41 * t806 * t465 + 0.589779626774411e37 * t228 * t373 + 0.576449749237333e58 * t539 * t409 + 0.472858038179319e41 * t695 * t380 + 0.339934303579609e26 * t207 * t380 - 0.1075782975e59 * t143 * t373 - 0.1600134812e42 * t819 * t358 + 0.248138270343821e56 * t776 * t358 - 0.8440804037e56 * t824 * t409 - 0.2093632890e44 * t240 * t424 + 0.108919477355032e59 * t484 * t409 + 0.489472445231188e40 * t591 * t380 + 0.241242774736184e58 * t833 * t358;
  t841 = t99 * t424;
  t848 = t468 * t3;
  t855 = t530 * t96;
  t866 = t89 * t465;
  t867 = t866 * t92;
  t872 = t468 * t11;
  t875 = 0.979649586039273e57 * t191 * t409 + 0.321305897552470e44 * t224 * t424 + 0.726027156520810e40 * t841 * t465 + 0.642194518291206e58 * t765 * t373 + 0.183182721230809e58 * t143 * t424 - 0.3275547296e58 * t848 * t424 + 0.450638563879818e39 * t579 * t380 - 0.3406213193e58 * t16 * t424 - 0.1548779283e44 * t855 * t380 - 0.6272079478e58 * t117 * t424 + 0.179833547188035e59 * t385 * t424 - 0.7527820685e49 * t147 * t358 + 0.727304204376191e58 * t100 * t409 + 0.164663373299812e41 * t867 * t380 - 0.5934600811e43 * t211 * t409 - 0.4765246983e59 * t872 * t373;
  t892 = t634 * t11;
  t895 = t492 * t96;
  t904 = t396 * t125;
  t911 = 0.242627340254370e57 * t362 * t373 - 0.2139179577e43 * t234 * t373 - 0.7240138926e58 * t16 * t409 + 0.864016295096288e58 * t403 * t424 - 0.1333173357e59 * t117 * t409 - 0.2067370658e58 * t737 * t358 + 0.925174122507938e40 * t585 * t634 + 0.501602363902581e42 * t892 * t373 + 0.230017094947168e44 * t895 * t380 + 0.593769478354827e43 * t211 * t424 - 0.2909816546e57 * t256 * t358 - 0.1713290425e42 * t585 * t361 + 0.325601166670929e42 * t904 * t380 - 0.5242918066e57 * t215 * t424 + 0.585773678119850e58 * t158 * t409;
  t914 = t413 * t92;
  t919 = t491 * t13;
  t931 = t530 * t119;
  t946 = -0.1092891782e59 * t397 * t424 - 0.9112262887e57 * t914 * t409 + 0.289415263046540e57 * t500 * t424 - 0.1052450202e58 * t919 * t424 + 0.744169899539267e15 * t417 + 0.650430785395217e49 * t185 * t424 - 0.3037420962e57 * t914 * t424 - 0.1541485003e58 * t623 * t424 - 0.9995591067e54 * t244 * t424 - 0.2650401381e58 * t931 * t409 - 0.1338054593e60 * t729 * t424 - 0.1274517569e58 * t713 * t424 + 0.146766796057449e45 * t222 * t358 - 0.2461040830e57 * t158 * t358 - 0.5936147797e43 * t211 * t358 - 0.3758381618e58 * t93 * t373;
  t954 = t395 * t99;
  t963 = t457 * t119;
  t968 = t361 * t11;
  t971 = t365 * t13;
  t984 = -0.1546412997e40 * t666 * t590 - 0.3825026892e59 * t558 * t409 - 0.9288590490e58 * t765 * t409 + 0.339917732102004e58 * t954 * t373 - 0.1679112366e40 * t585 * t465 - 0.3291353342e58 * t872 * t409 + 0.267322823438198e58 * t855 * t373 - 0.3514774245e56 * t963 * t358 - 0.6965128186e58 * t516 * t409 + 0.133455293521769e59 * t968 * t409 - 0.3782591021e58 * t971 * t409 + 0.134942346058522e55 * t187 * t424 - 0.3034027905e57 * t534 * t373 + 0.157013022682054e60 * t472 * t409 - 0.1022796340e59 * t126 * t424 - 0.6075870705e58 * t160 * t409;
  t989 = t396 * t119;
  t1004 = t3 * t358;
  t1017 = t13 * t409;
  t1019 = -0.5425372258e42 * t585 * t412 - 0.3936974214e58 * t484 * t424 + 0.239259405730010e58 * t989 * t424 + 0.113410046673884e22 * t448 * t4 + 0.113452744140412e57 * t254 * t358 + 0.175318514256933e40 * t378 * t409 + 0.139926030494797e39 * t585 * t376 - 0.3514774245e56 * t963 * t424 - 0.1206604427e56 * t246 * t358 + 0.559704121979188e39 * t1004 * t376 + 0.578034687729035e37 * t175 * t380 - 0.2971044995e50 * t238 * t424 - 0.3034027905e57 * t534 * t409 - 0.1951299080e39 * t654 * t380 + 0.116879009504622e40 * t378 * t358 + 0.390645945194749e16 * t1017;
  t1026 = t377 * t119;
  t1035 = t694 * t96;
  t1038 = t557 * t125;
  t1047 = t634 * t13;
  t1056 = t457 * t96;
  t1059 = -0.4865170767e16 * t236 * t373 + 0.429098566379896e40 * t1026 * t409 + 0.589779626774411e37 * t228 * t409 + 0.189143215271727e42 * t695 * t358 - 0.8077744849e17 * t242 * t373 - 0.9038243090e43 * t1035 * t409 - 0.9708394964e38 * t1038 * t380 - 0.2247240146e28 * t197 * t380 - 0.3067205893e57 * t392 * t424 - 0.1054484875e28 * t167 * t380 + 0.485662642066812e41 * t1047 * t424 - 0.1091072107e58 * t392 * t373 - 0.1490545945e43 * t472 * t380 + 0.493312609551463e57 * t931 * t424 + 0.160548249950962e57 * t1056 * t358;
  t1060 = t492 * t119;
  t1067 = t383 * t3;
  t1072 = t590 * t13;
  t1079 = t366 * t13;
  t1096 = 0.118742169328127e57 * t1060 * t424 - 0.2241355660e56 * t226 * t424 + 0.126138072166255e58 * t1060 * t373 + 0.510893945806757e57 * t1067 * t358 - 0.4135324523e58 * t215 * t409 + 0.380780480787768e58 * t1072 * t358 + 0.343539061308693e59 * t505 * t373 - 0.2258665608e50 * t217 * t373 - 0.5705689415e55 * t1079 * t424 - 0.6400352015e26 * t165 * t380 - 0.2221798355e56 * t122 * t358 - 0.4158935639e50 * t219 * t409 - 0.2275914100e41 * t448 * t465 + 0.309847563158691e57 * t165 * t373 + 0.258999464615157e15 * t230 * t424 - 0.9628253719e55 * t437 * t358;
  t1108 = t11 * t409;
  t1114 = t11 * t424;
  t1129 = -0.1692141431e58 * t141 * t424 + 0.877519748703892e57 * t688 * t373 + 0.114773734591781e59 * t258 * t373 + 0.100408140591232e35 * t115 * t380 - 0.2828858726e23 * t585 * t5 + 0.100866791361979e17 * t1108 - 0.4089241303e49 * t134 * t424 - 0.4000337030e41 * t819 * t424 + 0.189659508302699e40 * t1114 * t376 - 0.5180865160e56 * t657 * t424 - 0.1600490945e59 * t171 * t373 + 0.138628149806536e50 * t219 * t424 - 0.2933423282e59 * t531 * t373 - 0.8217402599e57 * t197 * t409 + 0.172394864267081e30 * t97 * t380;
  t1142 = t436 * t119;
  t1151 = t694 * t125;
  t1164 = 0.162375767792860e57 * t569 * t424 - 0.2067490504e32 * t169 * t380 + 0.110048972060320e59 * t1072 * t424 + 0.116767218929888e33 * t149 * t380 + 0.103539928954391e43 * t931 * t380 - 0.1291337427e43 * t414 * t380 + 0.811624334852520e55 * t1142 * t373 + 0.371637934927943e58 * t691 * t424 + 0.140446814223919e59 * t171 * t409 - 0.1404378832e59 * t141 * t409 + 0.772789985363420e41 * t1151 * t358 - 0.3674369864e58 * t989 * t373 + 0.962908144635886e58 * t484 * t373 - 0.1506373848e43 * t1035 * t380 + 0.141424611841618e59 * t115 * t409 - 0.3524343544e54 * t367 * t358;
  t1183 = t557 * t119;
  t1186 = t383 * t99;
  t1191 = t99 * t358;
  t1194 = t3 * t373;
  t1200 = 0.772789985363420e41 * t1151 * t373 - 0.1413941406e43 * t1060 * t380 + 0.158263946209984e59 * t561 * t358 + 0.439224587986642e21 * t464 * t4 + 0.337088936140765e58 * t954 * t358 + 0.616391329461637e55 * t195 * t424 + 0.392471128231606e43 * t252 * t409 + 0.475263837538731e57 * t722 * t373 + 0.182777633691709e59 * t1183 * t373 + 0.892448033258599e58 * t1186 * t424 - 0.4947917573e58 * t1183 * t424 + 0.290410862608327e41 * t1191 * t465 + 0.496113266359511e15 * t1194 + 0.100432354516343e59 * t175 * t424 - 0.2225105529e58 * t582 * t409;
  t1233 = 0.876234353183435e56 * t524 * t373 + 0.159706114225336e38 * t129 * t373 - 0.1165139272e57 * t618 * t373 + 0.501602363902581e42 * t892 * t358 + 0.650556841343834e30 * t191 * t380 + 0.629243479965558e37 * t691 * t380 + 0.196062834477694e58 * t453 * t409 - 0.1840561270e36 * t160 * t380 - 0.1134562589e57 * t213 * t358 - 0.1282769744e31 * t132 * t380 - 0.3934243074e43 * t252 * t424 - 0.5380498871e59 * t744 * t358 - 0.3074707109e58 * t156 * t409 - 0.2451814682e57 * t226 * t409 + 0.188522525499413e40 * t1072 * t380 - 0.4152370133e58 * t737 * t409;
  t1256 = t3 * t424;
  t1266 = 0.395810388460390e58 * t623 * t358 - 0.1729317224e58 * t848 * t358 - 0.4719830920e58 * t461 * t409 - 0.5852317763e57 * t931 * t373 + 0.658653493199249e41 * t867 * t373 + 0.496113266359511e15 * t1004 - 0.6025495394e43 * t1035 * t373 + 0.806177439099246e59 * t1038 * t373 + 0.811624334852488e55 * t1142 * t409 + 0.297161910976403e50 * t199 * t358 - 0.7295073696e38 * t154 * t409 + 0.124028316589878e15 * t1256 - 0.9654122318e57 * t113 * t424 - 0.1305383265e41 * t431 * t380 + 0.114099043999768e59 * t872 * t424 + 0.701155595328745e58 * t591 * t358;
  t1299 = 0.161467843073418e38 * t872 * t380 + 0.424913522219787e57 * t431 * t358 - 0.1467661364e45 * t222 * t373 + 0.180129450189575e57 * t189 * t373 + 0.768474374797446e57 * t453 * t358 + 0.417339665914342e57 * t688 * t358 - 0.1007467420e41 * t417 * t465 - 0.1011342635e57 * t534 * t358 + 0.389393910170803e58 * t362 * t424 - 0.5910045660e58 * t989 * t409 - 0.2652184542e60 * t561 * t409 + 0.103757717162143e43 * t666 * t390 - 0.5556817149e28 * t226 * t380 - 0.2112725552e59 * t115 * t373 - 0.6489600493e54 * t550 * t358 + 0.633910796190366e15 * t120 * t424;
  t1307 = t866 * t125;
  t1312 = t491 * t3;
  t1317 = t590 * t99;
  t1320 = t866 * t96;
  t1337 = 0.137073888062859e57 * t513 * t373 + 0.877729170298884e26 * t106 * t380 - 0.3506370285e40 * t1307 * t380 - 0.6716449464e40 * t1004 * t465 - 0.9573780313e55 * t1312 * t424 - 0.1340125593e57 * t254 * t424 + 0.218667212450653e59 * t1317 * t373 + 0.273394044970435e42 * t1320 * t380 + 0.483477664053694e58 * t215 * t373 + 0.244780146948558e20 * t179 * t380 - 0.1235782229e56 * t527 * t424 + 0.741608129613457e33 * t585 * t14 + 0.230090041394785e57 * t688 * t424 + 0.664962026891669e56 * t414 * t358 + 0.678551643339766e56 * t195 * t409;
  t1354 = t653 * t92;
  t1363 = t599 * t119;
  t1368 = t377 * t92;
  t1373 = 0.308978935778809e43 * t765 * t380 - 0.7310333781e59 * t149 * t373 - 0.4694399607e58 * t173 * t424 + 0.134016072755461e57 * t213 * t424 + 0.152270193746682e43 * t833 * t380 + 0.242177702148800e57 * t258 * t358 - 0.2400202218e42 * t819 * t409 + 0.407156336764161e51 * t193 * t373 - 0.1577170327e38 * t1354 * t380 + 0.320821795861496e56 * t734 * t424 + 0.821740505229881e57 * t97 * t409 + 0.109334593838060e44 * t104 * t424 - 0.1323774629e59 * t1363 * t409 + 0.758638033210798e40 * t806 * t376 - 0.1372194777e40 * t1368 * t380 - 0.8494825696e41 * t657 * t380;
  t1391 = t694 * t92;
  t1404 = t377 * t96;
  t1407 = -0.7297756151e16 * t236 * t409 + 0.768042529376457e58 * t145 * t409 + 0.447483463417546e58 * t201 * t409 + 0.423679595181998e57 * t453 * t424 - 0.1951268491e50 * t185 * t409 + 0.297035944209433e50 * t199 * t424 - 0.1840552486e55 * t205 * t424 + 0.278094844099418e58 * t115 * t424 - 0.5443672323e42 * t1391 * t409 + 0.293812133351974e40 * t612 * t376 - 0.3863324579e59 * t1354 * t358 + 0.822797066781517e57 * t197 * t373 - 0.1635267798e58 * t126 * t358 - 0.1625565576e42 * t464 * t366 - 0.2278283708e41 * t1404 * t424;
  t1434 = t361 * t13;
  t1441 = 0.273382714455458e30 * t201 * t380 + 0.122315812146063e24 * t666 * t5 + 0.606407968541319e42 * t666 * t365 + 0.302254481853269e58 * t469 * t358 + 0.555104473504763e41 * t635 * t409 + 0.580624822989402e59 * t505 * t409 + 0.487127303378530e57 * t569 * t373 + 0.189659508302699e40 * t448 * t376 - 0.8271534720e58 * t520 * t373 + 0.245768055748730e57 * t226 * t373 - 0.6852481113e42 * t516 * t380 + 0.394981151782001e58 * t143 * t409 - 0.1857700352e36 * t16 * t380 - 0.4072739814e58 * t1434 * t373 + 0.113795704981619e41 * t1108 * t376 - 0.4088175576e49 * t134 * t358;
  t1474 = -0.8786470507e35 * t126 * t380 - 0.1351053202e55 * t187 * t358 + 0.194265056826725e42 * t1047 * t373 + 0.744414811031288e56 * t776 * t373 - 0.5145767371e37 * t666 * t468 + 0.623308123486860e58 * t250 * t373 - 0.1711706825e56 * t1079 * t373 + 0.125400590975645e42 * t448 * t634 - 0.4572344206e43 * t182 * t409 + 0.650039886863762e18 * t666 * r + 0.385847106555705e39 * t585 * t590 + 0.222613255377080e58 * t440 * t373 + 0.456912960209592e56 * t513 * t358 + 0.105494900470259e41 * t666 * t383 + 0.321897070739432e58 * t713 * t409;
  t1481 = t402 * t3;
  t1484 = t366 * t11;
  t1491 = t383 * t11;
  t1500 = t468 * t99;
  t1505 = t590 * t3;
  t1511 = 0.303199096964296e60 * t729 * t409 - 0.2587889264e59 * t469 * t373 + 0.125400590975645e42 * t892 * t424 + 0.292951209135495e58 * t1481 * t358 - 0.3151875892e56 * t1484 * t409 - 0.1503373534e18 * t585 * r - 0.2888476116e56 * t437 * t409 + 0.185915118950513e59 * t1491 * t373 + 0.372982800445815e56 * t106 * t409 - 0.2019436212e17 * t242 * t424 + 0.196593208924805e37 * t228 * t424 - 0.6031475956e58 * t1500 * t424 - 0.1362353884e23 * t240 * t380 + 0.125838305844456e59 * t1505 * t409 - 0.6050226304e39 * t666 * t376 + 0.168111318936632e16 * t1114;
  t1539 = t436 * t92;
  t1546 = 0.225873927327254e50 * t147 * t409 - 0.2019436212e17 * t242 * t380 - 0.1392000957e58 * t213 * t373 - 0.9628253719e55 * t437 * t424 - 0.7614639075e58 * t113 * t409 + 0.167838276423527e59 * t397 * t373 + 0.650406192271234e49 * t185 * t358 + 0.160040192942715e59 * t141 * t373 + 0.193685075018479e42 * t464 * t412 - 0.6163911532e55 * t167 * t424 - 0.1216292692e16 * t236 * t424 - 0.8077744849e17 * t242 * t358 + 0.193197496340856e41 * t1151 * t424 + 0.131940134140556e57 * t1539 * t373 - 0.6506211994e58 * t1354 * t409 - 0.1337400113e59 * t618 * t409;
  t1579 = -0.6637732754e59 * t650 * t409 - 0.8284937646e58 * t420 * t358 - 0.1571856995e58 * t931 * t358 + 0.845771519518723e58 * t803 * t358 + 0.139926030494797e39 * t1256 * t376 + 0.439800447135278e56 * t1539 * t424 - 0.3213066197e44 * t224 * t358 - 0.1330325036e58 * t254 * t409 - 0.5288618400e41 * t1017 * t465 - 0.1356804099e51 * t193 * t424 + 0.987980239798872e41 * t867 * t409 - 0.5423918053e57 * t895 * t424 - 0.2957144884e41 * t666 * t395 - 0.2026296289e57 * t574 * t358 + 0.375355432339805e58 * t93 * t409 - 0.1506373848e43 * t1035 * t424;
  t1615 = 0.532353714084451e37 * t129 * t424 - 0.1141886077e58 * t1312 * t373 - 0.3301296090e57 * t443 * t373 - 0.1434806065e29 * t108 * t380 - 0.1011342635e57 * t534 * t424 - 0.2022162706e58 * t132 * t424 + 0.119971352970620e37 * t585 * t468 + 0.884616696744292e57 * t848 * t409 + 0.409114831334879e56 * t108 * t358 - 0.7295073696e38 * t154 * t373 - 0.1221270608e60 * t152 * t409 - 0.3213073418e44 * t224 * t409 - 0.3184303053e58 * t1500 * t358 - 0.5966248342e58 * t1434 * t424 - 0.9112262887e57 * t914 * t373;
  t1648 = 0.471492366139179e58 * t203 * t424 + 0.440718200027962e40 * t1017 * t376 + 0.270541444950892e55 * t1142 * t424 + 0.133895622268933e42 * t666 * t366 - 0.2431691232e38 * t154 * t424 + 0.283714822907591e42 * t695 * t409 - 0.5526674975e58 * t904 * t373 + 0.106764970732211e37 * t232 * t358 - 0.1937944262e58 * t833 * t424 - 0.1508209613e59 * t1038 * t358 + 0.273394044970435e42 * t1320 * t424 - 0.1058630944e28 * t585 * t18 + 0.590890608002908e58 * t493 * t373 - 0.2884974257e60 * t260 * t373 - 0.3729062622e56 * t106 * t373 + 0.481644749852981e57 * t1056 * t373;
  t1680 = -0.2420090522e40 * t1191 * t376 + 0.593924176986214e43 * t211 * t373 + 0.296408156822829e43 * t427 * t380 - 0.4865170767e16 * t236 * t358 - 0.8227973114e57 * t97 * t373 - 0.6078888866e57 * t574 * t409 + 0.258999464615157e15 * t230 * t380 - 0.7417018430e57 * t582 * t358 + 0.172036266169008e31 * t111 * t380 - 0.4873521430e57 * t427 * t373 + 0.507589702298482e58 * t1363 * t424 - 0.2888476116e56 * t437 * t373 + 0.236290953387731e36 * t102 * t358 - 0.8690928304e58 * t362 * t409 + 0.637472276607546e58 * t623 * t409;
  t1702 = t395 * t3;
  t1713 = 0.194265056826725e42 * t1047 * t358 - 0.5362829067e15 * t841 - 0.4577896542e43 * t182 * t358 - 0.1525367145e59 * t1363 * t373 - 0.3985364856e57 * t160 * t358 + 0.148826193945080e42 * t561 * t380 - 0.6489600493e54 * t550 * t424 + 0.104618347277984e38 * t179 * t424 + 0.278315545281073e24 * t193 * t380 + 0.617135588354615e55 * t167 * t358 + 0.487645727932519e32 * t666 * t7 + 0.184601019118166e58 * t1702 * t373 + 0.452548531596958e57 * t108 * t373 + 0.890258525506349e58 * t113 * t373 + 0.437202150408296e57 * t143 * t358 - 0.2699092995e59 * t403 * t373;
  t1746 = 0.349089323434265e58 * t510 * t358 - 0.1147369652e59 * t145 * t373 + 0.458344887866665e43 * t182 * t424 - 0.9044799266e41 * t666 * t361 - 0.8913598290e50 * t199 * t373 - 0.1093061082e44 * t104 * t358 + 0.103873073037893e58 * t141 * t358 + 0.225854303414122e50 * t147 * t373 + 0.121065329869654e59 * t855 * t409 + 0.163262124383637e56 * t189 * t358 + 0.678328330310744e56 * t167 * t373 + 0.464417225739364e58 * t169 * t373 + 0.164036426982262e43 * t1320 * t409 - 0.5347611132e58 * t1505 * t424 + 0.408369934806663e41 * t618 * t380;
  t1779 = 0.292931241857535e59 * t650 * t424 + 0.170634868574784e57 * t427 * t424 + 0.131013059530030e58 * t919 * t358 + 0.985889028817965e58 * t968 * t358 + 0.837387266041970e58 * t1183 * t409 - 0.3939008970e43 * t252 * t373 - 0.8455677999e56 * t191 * t358 + 0.717994591318360e58 * t855 * t358 - 0.8047609393e59 * t803 * t373 + 0.109357617988174e43 * t1320 * t373 - 0.1569221555e60 * t561 * t373 - 0.9543510386e59 * t472 * t424 + 0.715164277299828e39 * t1026 * t424 + 0.284006886456334e58 * t971 * t373 + 0.183791640082982e34 * t171 * t380 + 0.320294912196635e37 * t232 * t409;
  t1797 = t402 * t13;
  t1810 = t402 * t11;
  t1815 = 0.197758076782384e59 * t16 * t373 + 0.836710892658179e20 * t585 * t4 + 0.109273235633976e59 * t493 * t424 - 0.4817500207e58 * t1317 * t358 + 0.122655956559392e50 * t134 * t373 + 0.183064769114110e58 * t1702 * t358 + 0.100031155651345e55 * t244 * t358 + 0.410838690213511e36 * t260 * t380 - 0.9166021129e58 * t1797 * t424 - 0.2616266713e58 * t1505 * t358 - 0.3707346688e56 * t527 * t409 - 0.3744710266e59 * t203 * t373 - 0.3369899763e59 * t531 * t424 + 0.582175556326051e58 * t749 * t424 - 0.1687798437e59 * t1810 * t424 + 0.193197496340856e41 * t1151 * t380;
  t1844 = t99 * t373;
  t1848 = 0.155399678769094e16 * t230 * t409 - 0.2813601346e56 * t824 * t358 - 0.7499397824e58 * t968 * t373 - 0.3525745600e41 * t612 * t465 - 0.3702539870e57 * t156 * t424 + 0.737749361851325e58 * t765 * t424 + 0.131764983250504e57 * t461 * t373 + 0.532353714084451e37 * t129 * t358 - 0.1640333239e56 * t189 * t424 + 0.940741888819744e57 * t1186 * t358 - 0.1208013155e36 * t173 * t380 + 0.484666732243604e58 * t1067 * t424 + 0.137073888062860e57 * t513 * t409 - 0.2225105529e58 * t582 * t373 - 0.2145131627e16 * t1844 - 0.7224670392e57 * t248 * t409;
  t1868 = t99 * t409;
  t1883 = -0.1160047030e44 * t1056 * t380 + 0.592901623854378e43 * t914 * t380 + 0.969126692556203e56 * t191 * t424 + 0.171051293138232e60 * t152 * t373 - 0.6871016147e37 * t162 * t409 - 0.5761748928e58 * t895 * t373 - 0.3757634014e41 * t968 * t380 + 0.724763107781588e58 * t1434 * t409 + 0.435616293912487e41 * t1868 * t465 + 0.128323846284923e36 * t209 * t424 + 0.315424291685962e58 * t971 * t424 + 0.303676035184890e42 * t1539 * t380 + 0.696758905740781e59 * t600 * t373 - 0.1463460776e58 * t833 * t373 - 0.8233168665e40 * t1368 * t409;
  t1916 = -0.1092776225e44 * t104 * t409 + 0.342254987468301e57 * t93 * t424 + 0.173433164606601e59 * t520 * t409 - 0.1476486833e42 * t666 * t456 + 0.320294912196635e37 * t232 * t373 + 0.138672708396655e59 * t160 * t373 - 0.5920444518e32 * t215 * t380 + 0.103599785846063e16 * t230 * t358 + 0.127633008672214e43 * t392 * t380 - 0.3437379906e59 * t904 * t409 + 0.380346477714219e16 * t120 * t409 - 0.7529783077e49 * t147 * t424 - 0.5350166261e42 * t971 * t380 - 0.1787453622e58 * t469 * t409 + 0.587586257737852e42 * t585 * t365 + 0.664962026891508e56 * t414 * t424;
  t1947 = 0.227206987201318e57 * t156 * t358 + 0.253564318476145e16 * t120 * t373 - 0.2721871422e38 * t585 * t402 + 0.637013599139315e43 * t744 * t380 - 0.4644152644e58 * t201 * t373 + 0.631988141764236e59 * t744 * t424 - 0.3037420962e57 * t914 * t358 - 0.3549351438e35 * t585 * t15 + 0.752904855709666e49 * t217 * t358 + 0.458900121494707e43 * t182 * t373 - 0.1946880148e55 * t550 * t373 - 0.9582678043e40 * t729 * t380 - 0.1991743240e59 * t600 * t358 - 0.5362829067e15 * t666 + 0.212900205365897e59 * t431 * t409;
  t1979 = -0.2675703540e59 * t173 * t409 + 0.524132298043287e57 * t256 * t424 - 0.7720209345e57 * t16 * t358 + 0.121722314287707e56 * t246 * t424 - 0.7277894819e56 * t248 * t424 + 0.184055280131002e55 * t207 * t424 + 0.467194053010555e59 * t260 * t424 - 0.9072787205e41 * t1391 * t424 + 0.286065710919930e40 * t1026 * t358 - 0.1781756263e59 * t691 * t373 + 0.539429907275501e58 * t403 * t358 - 0.8912338623e50 * t199 * t409 + 0.168111318936632e16 * t448 + 0.616134146946717e56 * t248 * t358 - 0.2458516273e57 * t145 * t358 - 0.4514684102e57 * t108 * t409;
  t1998 = t866 * t119;
  t2013 = 0.170658887061555e34 * t203 * t380 + 0.267364513835535e59 * t203 * t409 - 0.9127971081e42 * t524 * t380 + 0.107739647934174e43 * t919 * t380 - 0.3202566699e34 * t666 * t14 - 0.1050625297e56 * t1484 * t358 - 0.5488779110e40 * t1368 * t373 + 0.141257207765667e58 * t713 * t358 - 0.8581971328e40 * t1998 * t380 + 0.248138270343800e56 * t776 * t424 + 0.891302597213288e50 * t238 * t409 + 0.313855041833953e38 * t179 * t409 + 0.164663373299812e41 * t867 * t424 - 0.4898030218e58 * t1797 * t358 - 0.1144711008e58 * t203 * t358;
  t2046 = 0.755963633080353e57 * t248 * t373 - 0.2403281077e59 * t1072 * t409 - 0.2625735732e58 * t8 * t373 + 0.758638033210798e40 * t594 * t376 - 0.3629114882e42 * t1391 * t373 + 0.224135127836936e56 * t122 * t424 - 0.5180865160e56 * t657 * t358 + 0.270541444950845e55 * t1142 * t358 + 0.133152091593212e57 * t246 * t409 + 0.522691447472816e56 * t477 * t373 + 0.666267970166994e30 * t213 * t380 - 0.6386619293e58 * t516 * t358 + 0.209109668796594e44 * t240 * t358 - 0.3558545755e58 * t1702 * t424 + 0.559704121979188e39 * t1194 * t376 + 0.195462403155893e44 * t582 * t380;
  t2080 = -0.3455751358e41 * t1491 * t380 + 0.708872860163194e36 * t102 * t409 + 0.208929617894522e59 * t904 * t424 + 0.370069649003175e41 * t635 * t358 - 0.2563180912e58 * t111 * t373 + 0.805048445273071e57 * t100 * t358 - 0.2096169093e44 * t240 * t373 - 0.6452362490e58 * t749 * t358 + 0.161269264084078e59 * t420 * t424 + 0.102389638300639e58 * t1060 * t358 + 0.460614853388067e25 * t238 * t380 - 0.1711706825e56 * t1079 * t409 - 0.1054432274e57 * t963 * t409 - 0.8814364001e40 * t710 * t465 + 0.752403545853867e42 * t892 * t409 - 0.2815559189e56 * t165 * t424;
  t2113 = -0.2546270232e59 * t1038 * t424 + 0.486773544644824e59 * t1810 * t373 + 0.782420603647732e57 * t1491 * t358 - 0.6744456628e58 * t1067 * t373 + 0.746959448481037e58 * t260 * t358 - 0.5455524524e57 * t500 * t409 + 0.213141485141087e43 * t234 * t409 + 0.131520795775812e57 * t250 * t358 - 0.2255050814e58 * t1312 * t409 - 0.3862590539e57 * t201 * t358 + 0.264354824520390e59 * t1797 * t373 - 0.7954162294e58 * t904 * t358 - 0.1392429766e43 * t520 * t380 - 0.2431691232e38 * t154 * t358 - 0.2813601346e56 * t824 * t424 - 0.6785516610e56 * t167 * t409;
  t2148 = 0.160548249950994e57 * t1056 * t424 - 0.1057303063e55 * t367 * t409 - 0.1216292692e16 * t236 * t380 + 0.745753003333204e56 * t197 * t358 + 0.286065710919930e40 * t1026 * t373 - 0.6171357648e55 * t195 * t358 + 0.641493729130501e59 * t149 * t409 + 0.498381133875070e58 * t539 * t373 - 0.2278283708e41 * t1404 * t380 + 0.512140156187915e57 * t691 * t409 + 0.377828085376019e59 * t744 * t373 + 0.439800447135287e56 * t1539 * t358 - 0.2096415853e59 * t493 * t409 + 0.350384657344789e58 * t156 * t373 + 0.556560621821184e58 * t872 * t358;
  t2181 = 0.122521924502821e42 * t420 * t380 - 0.2103822171e41 * t1307 * t409 - 0.1982875193e59 * t1491 * t424 - 0.1293596842e58 * t510 * t373 + 0.184193802022196e55 * t205 * t358 - 0.3151875892e56 * t1484 * t373 - 0.6783283127e56 * t195 * t373 + 0.522882370981901e58 * t152 * t358 + 0.133226459063894e59 * t954 * t409 + 0.407098797523860e51 * t193 * t409 + 0.290410862608327e41 * t1844 * t465 + 0.102771315419091e58 * t158 * t424 - 0.8440804037e56 * t824 * t373 - 0.3301296090e57 * t443 * t409 + 0.388834908578136e34 * t145 * t380 + 0.715164277299828e39 * t1026 * t380;
  t2212 = -0.1951243897e50 * t185 * t373 + 0.118752822974261e59 * t1505 * t373 + 0.139200446961659e58 * t254 * t373 - 0.1679112366e40 * t1256 * t465 - 0.2717334170e59 * t505 * t424 + 0.647578698838813e25 * t585 * t6 + 0.619644043648913e58 * t469 * t424 + 0.651076575324578e15 * t710 + 0.427680188051979e21 * t182 * t380 + 0.222180273837880e56 * t226 * t358 + 0.281895962872740e56 * t165 * t358 + 0.313855041833953e38 * t179 * t373 - 0.2198697607e58 * t160 * t424 + 0.522960150527700e58 * t516 * t373 - 0.1402548114e41 * t1307 * t373;
  t2245 = -0.3506370285e40 * t1307 * t424 + 0.253564318476145e16 * t120 * t358 - 0.1357379491e51 * t193 * t358 + 0.151026518430795e58 * t145 * t424 + 0.148550725588114e56 * t187 * t409 - 0.3617833104e21 * t666 * t4 + 0.413548409018292e58 * t256 * t409 - 0.3391679552e55 * t106 * t358 - 0.4174838577e57 * t500 * t358 - 0.2102628293e58 * t737 * t373 - 0.4834936120e58 * t256 * t373 + 0.708872860163194e36 * t102 * t373 + 0.813872612711812e59 * t654 * t373 + 0.487127303378521e57 * t569 * t409 - 0.1402548114e41 * t1307 * t358 + 0.633910796190366e15 * t120 * t380;
  t2278 = 0.182043940348002e58 * t175 * t358 - 0.1398246371e39 * t1797 * t380 - 0.3707346688e56 * t527 * t373 - 0.9113134832e41 * t1404 * t358 + 0.152991052254593e36 * t666 * t15 - 0.4127152366e56 * t108 * t424 + 0.723520363724951e58 * t1702 * t409 - 0.4000337030e41 * t666 * t634 - 0.3525745600e41 * t359 * t465 + 0.112415697513387e58 * t173 * t358 - 0.8202527531e58 * t158 * t373 - 0.4744725814e58 * t149 * t358 + 0.485662642066812e41 * t464 * t634 - 0.1122738517e58 * t1312 * t358 - 0.7947698079e57 * t919 * t373;
  t2310 = -0.3926074567e58 * t1183 * t358 + 0.101487945316179e43 * t585 * t390 + 0.106764970732211e37 * t232 * t424 + 0.159706114225336e38 * t129 * t409 - 0.7492756948e56 * t197 * t424 + 0.672445275746530e16 * t594 - 0.6552579224e58 * t954 * t424 - 0.2267252877e57 * t524 * t358 - 0.1969053175e41 * t484 * t380 + 0.284643451146915e57 * t8 * t424 + 0.246940049176821e41 * t585 * t395 - 0.2392245038e58 * t510 * t424 - 0.1100432030e57 * t443 * t358 + 0.410916617866632e59 * t1038 * t409 - 0.1485018333e56 * t187 * t373 - 0.1467664662e45 * t222 * t424;
  t2344 = 0.392026844259748e59 * t1491 * t409 - 0.5487900173e30 * t666 * t110 + 0.384971538854770e36 * t209 * t409 - 0.1600134812e42 * t819 * t373 - 0.1506589753e58 * t424 * t258 - 0.9113134832e41 * t1404 * t373 - 0.2025578534e56 * t205 * t409 - 0.1554259548e57 * t657 * t373 + 0.122148199995878e59 * t385 * t358 - 0.9846900395e58 * t1317 * t424 - 0.7683767227e18 * t162 * t380 + 0.138637214606924e50 * t219 * t358 + 0.117791622536915e59 * t520 * t358 - 0.3616225518e39 * t1810 * t380 - 0.4527024173e57 * t115 * t358 - 0.8371428603e57 * t722 * t424;
  t2377 = -0.5488779110e40 * t1368 * t358 - 0.1528744759e33 * t113 * t380 + 0.224588041604625e58 * t8 * t409 - 0.3139387450e18 * t129 * t380 + 0.260470633639946e59 * t1500 * t373 + 0.481644749852955e57 * t1056 * t409 + 0.962465387584567e56 * t734 * t373 + 0.218172788199872e58 * t539 * t358 + 0.212602037368247e42 * t464 * t456 - 0.1459189760e39 * t1183 * t380 - 0.3464765529e58 * t505 * t358 - 0.9371436663e58 * t654 * t358 - 0.1050625297e56 * t1484 * t424 + 0.291397585240088e42 * t1047 * t409 + 0.174230482490933e56 * t477 * t358 - 0.6716449464e40 * t1194 * t465;
  t2411 = 0.380637624260450e59 * t531 * t358 - 0.1054432274e57 * t963 * t373 - 0.2457680995e57 * t122 * t373 + 0.127927491786070e42 * t989 * t380 + 0.664344474786550e42 * t513 * t380 + 0.277534299650751e59 * t175 * t409 - 0.1128378719e32 * t585 * t7 + 0.367599599338031e60 * t1354 * t373 + 0.156781098899844e58 * t919 * t409 - 0.2258616380e50 * t217 * t409 + 0.126946510420674e30 * t585 * t110 - 0.5927038928e23 * t199 * t380 + 0.221406402623176e56 * t458 * t424 - 0.8170975246e59 * t1354 * t424 + 0.536097006996510e57 * t113 * t358;
  t2443 = 0.189143215271727e42 * t695 * t373 - 0.2044019955e59 * t132 * t409 + 0.382188846961351e34 * t158 * t380 - 0.4425316054e59 * t591 * t409 + 0.331212849606813e26 * t187 * t380 + 0.337306128961239e58 * t100 * t424 - 0.2799989408e26 * t666 * t6 + 0.891291363538679e50 * t238 * t373 - 0.2145131627e16 * t1191 - 0.2089143084e57 * t111 * t358 - 0.1261982446e58 * t392 * t409 - 0.1465811491e59 * t1481 * t373 + 0.199488608067435e57 * t414 * t409 - 0.4173045813e58 * t250 * t409 - 0.5149182797e41 * t1998 * t409 + 0.109363079492251e44 * t104 * t373;
  t2474 = 0.124028316589878e15 * t585 + 0.141455240133989e59 * t848 * t373 - 0.9072787205e41 * t1391 * t380 - 0.6025495394e43 * t1035 * t358 - 0.1365548460e42 * t1108 * t465 - 0.1580250305e57 * t8 * t358 - 0.2025162842e56 * t207 * t373 - 0.1366970225e42 * t1404 * t409 + 0.396766731471060e42 * t963 * t380 - 0.2970932658e50 * t238 * t358 + 0.156918261056378e57 * t989 * t358 - 0.6050226304e39 * t841 * t376 - 0.3629114882e42 * t1391 * t358 + 0.535412254994095e58 * t1434 * t358 + 0.174888128293040e21 * t224 * t380;
  t2507 = -0.3468416969e58 * t971 * t358 + 0.147429417190605e19 * t154 * t380 + 0.109357617988174e43 * t1320 * t358 + 0.291141192387267e57 * t215 * t358 - 0.3096658205e41 * t585 * t366 + 0.604675251254038e59 * t600 * t409 - 0.2040274031e59 * t1186 * t409 - 0.1484980709e24 * t217 * t380 - 0.2253358839e58 * t855 * t424 - 0.2275914100e41 * t1114 * t465 - 0.3937209394e43 * t510 * t380 + 0.658653493199249e41 * t867 * t358 + 0.158046547802360e59 * t1797 * t409 + 0.378498826500953e58 * t440 * t358 + 0.133032152025704e58 * t213 * t409 - 0.3785922959e59 * t591 * t373;
  t2540 = 0.338281559387530e37 * t803 * t380 + 0.146767125906010e45 * t222 * t409 - 0.7457505525e56 * t97 * t358 + 0.291021275787978e59 * t1810 * t409 + 0.142436098544504e58 * t803 * t409 - 0.1098603689e59 * t968 * t424 - 0.1357548789e44 * t440 * t380 + 0.139088386166631e59 * t650 * t358 - 0.2056042302e59 * t1072 * t373 - 0.8181921021e57 * t250 * t424 - 0.4797736259e36 * t117 * t380 - 0.3432788531e41 * t1998 * t373 - 0.2420090522e40 * t1844 * t376 - 0.8581971328e40 * t1998 * t424 + 0.566656400770185e41 * t1434 * t380;
  t2571 = -0.1187175466e43 * t666 * t491 + 0.231714167876961e59 * t1317 * t409 - 0.5705689415e55 * t1079 * t358 - 0.1470365030e59 * t749 * t409 - 0.1241900187e59 * t1186 * t373 + 0.202516287353545e56 * t205 * t373 - 0.1108023001e59 * t1067 * t409 - 0.3330384295e40 * t585 * t383 + 0.116308815617617e60 * t579 * t424 + 0.207639801598430e40 * t650 * t380 + 0.651076575324578e15 * t464 - 0.3217697440e16 * t1868 - 0.3432788531e41 * t1998 * t358 - 0.1446262728e24 * t185 * t380 + 0.875813989410341e58 * t618 * t424 + 0.215820017109399e40 * t1363 * t380;
  t2605 = -0.3630135783e40 * t1868 * t376 + 0.109981649978393e56 * t244 * t373 - 0.2136591335e43 * t234 * t424 + 0.436038875878512e58 * t1363 * t358 + 0.679054867440064e58 * t440 * t409 + 0.442698403732952e57 * t201 * t424 - 0.9896322439e58 * t895 * t409 - 0.1582707177e42 * t1142 * t380 + 0.288691123315963e58 * t833 * t409 + 0.208856048560765e44 * t240 * t409 - 0.1037838610e58 * t171 * t358 + 0.115918497804513e42 * t1151 * t409 + 0.292197523761556e39 * t378 * t424 + 0.664219207869572e56 * t458 * t409 - 0.6353294815e59 * t729 * t358 - 0.1372194777e40 * t1368 * t424;
  t2637 = -0.4676965317e58 * t895 * t358 + 0.216653493534073e58 * t1060 * t409 + 0.214954955822726e58 * t722 * t358 - 0.9019058136e58 * t1810 * t358 - 0.3524343544e54 * t367 * t424 + 0.131940134140557e57 * t1539 * t409 + 0.469226149583779e58 * t1481 * t424 - 0.5794565209e32 * t156 * t380 + 0.269959724770053e59 * t397 * t409 + 0.162890163214325e58 * t1500 * t409 - 0.8329679786e58 * t1481 * t409 + 0.522691447472823e56 * t477 * t409 + 0.196593208924805e37 * t228 * t358 - 0.5418438531e28 * t189 * t380 - 0.2463370715e44 * t539 * t380 + 0.672445275746530e16 * t806;
  t2648 = t353 * t353;
  t2649 = t2648 * t2648;
  t2654 = exp(-0.10000e5 * t354);
  t2655 = t2654 + t348;
  t2656 = t2655 * t2655;
  t2657 = t2656 * t2656;
  if (r <= 0.100000000000000e-1)
    return 0.100000000000000e-10 * (t136 + t177 + t221 + t262) / t266 / t269 * t342;
  else if (r < 0.990000000000000e0)
    return 0.100000000000000e-4 * t348 * (t725 + t1373 + t2080 + t836 + t2637 + t1713 + t2113 + t875 + t911 + t2571 + t946 + t566 + t2605 + t2148 + t984 + t1019 + t1059 + t2181 + t1096 + t1129 + t1164 + t2212 + t2245 + t609 + t2278 + t1200 + t1233 + t1266 + t2310 + t1299 + t1337 + t1407 + t2344 + t646 + t1441 + t1474 + t2377 + t1511 + t1546 + t1579 + t1615 + t2411 + t1648 + t685 + t1680 + t1746 + t1779 + t2443 + t1815 + t1848 + t2474 + t1883 + t1916 + t430 + t1947 + t2507 + t764 + t480 + t1979 + t2013 + t2046 + t523 + t799 + t2540) * t342 / t268 / t269 / t2648 / t2649 / t2655 / t2657;
  else
    return 0.0;
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

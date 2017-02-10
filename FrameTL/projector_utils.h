#ifndef _FRAMETL_PROJECTOR_UTILS_H
#define _FRAMETL_PROJECTOR_UTILS_H

#include <cmath>
#include <set>

namespace FrameTL
{

  // smooth partition of unity
  double zeta(double r, double r0, double r1) {
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


  // linear partition of unity - C0
  double linear(double r, double r0, double r1){
      if(r <= r0)
          return 1.0;
      else{
          if(r>=r1)
              return 0.0;
          else{
              return (r1-r)/(r1-r0);
          }
      }
                 
  }

  // C1- partition of unity
  double C1(double r, double r0, double r1){
      if (r<= r0) return 1.0;
      if (r>= r1) return 0.0;
      if (r<= (r0+r1)/2) return 1.0 - 2.0*(1.0/(r1-r0)*(r-r0))*(1.0/(r1-r0)*(r-r0));
      else return 2.0*(1.0/(r1-r0)*(r-r0)-1.0)*(1.0/(r1-r0)*(r-r0)-1.0);
  }

  // C2- function (needs degree 4 - degree 3 not sufficient)
  double C2(double r, double r0, double r1){
      if (r<= r0) return 1.0;
      if (r>= r1) return 0.0;
      if (r<= (r0+r1)/2) return 1.0 - 8.0*pow((r-r0)/(r1-r0),4);
      else return 8.0*pow((r-r0)/(r1-r0)-1.0,4);
  }

  // C^\infty - function, possibly better than zeta
  double Cinfty(double r, double r0, double r1){
      if (r<= r0) return 1.0;
      if (r>= r1) return 0.0;
      else return pow(2*(zeta(r,r0,r1)-0.5),100);
  }

  // Partition of unity on L-Domain, cf Manuels thesis
  double sigmaL(Point<2> x, char funktion)
  {
      const double PI = 3.14159265;
      if (x[1] >= 0) return 0;
      else if (x[0] >= 0) return 1;
      else{
          double alpha = atan(x[0]/x[1]);
          if (funktion == '1') return C1(alpha, 0, 0.5 * PI);
          else if (funktion == '2') return C2(alpha, 0*PI, 0.5 * PI);
          else if (funktion == 'z') return zeta(alpha,0,0.5*PI);
          else return linear(alpha, 0, 0.5 * PI);
      }
  }
  
}

#endif
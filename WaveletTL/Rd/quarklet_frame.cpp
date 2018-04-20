// implementation of quarklet_frame.h

#include <cassert>
#include <numerics/cardinal_splines.h>
#include <utils/tiny_tools.h>

namespace WaveletTL
{
  template <int d, int dt>
  double
  QuarkletFrame<d,dt>::evaluate(const unsigned int derivative,
			   const RQIndex& lambda,
			   const double x)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;
    const int j = lambda.j();
    const int p = lambda.p();
    const int k = lambda.k();
    const int e = lambda.e();
//    const double pfkt = (d+1)/2;
    const double y = 2*((1<<j)*x-k-((d % 2) == 0 ? 0 : 0.5))/d;
    
    if (e == 0) // generator
      {
	r = (derivative == 0
	     ? MathTL::EvaluateCardinalBSpline_td<d>(j, k, x)*pow(y,p)
	     : MathTL::EvaluateCardinalBSpline_td_x<d>(j, k, x)*pow(y,p)+MathTL::EvaluateCardinalBSpline_td<d>(j, k, x)*pow(y,p-1)*(1<<j)*p*2./d);
      }
    else // wavelet
      {
	// old generic version: readable and slow
 	InfiniteVector<double, RQIndex> gcoeffs;
        QuarkletFrame<d, dt> tempframe;
	tempframe.reconstruct_1(lambda, j+1, gcoeffs);
	//QuarkletFrame<d,dt>::reconstruct_1(lambda, j+1, gcoeffs); Hier gab es Fehlermeldung
 	for (typename InfiniteVector<double, RQIndex>::const_iterator it(gcoeffs.begin());
 	     it != gcoeffs.end(); ++it)
 	  r += *it * evaluate(derivative, it.index(), x);
	  }
    
      return r;
  }

  template <int d, int dt>
  void
  QuarkletFrame<d,dt>::evaluate
  (const unsigned int derivative,
   const RQIndex& lambda,
   const Array1D<double>& points, Array1D<double>& values)
  {
    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++) {
      values[i] = evaluate(derivative, lambda, points[i]);
    }
  }

  // 
  // some fast template specializations without reconstruct_1() calls
  

  template <>
  double
  QuarkletFrame<1,1>::evaluate(const unsigned int derivative,
 			  const RQIndex& lambda,
 			  const double x)
  {
    if (derivative == 0) {
      const double y = (1<<lambda.j())*x-lambda.k();
      const double factor = twotothejhalf(lambda.j());
      if (lambda.e() == 0) // generator
 	return (y >= 0 && y < 1 ? factor * pow(y, lambda.p()) : 0.);
      else // wavelet
 	{
 	  if (y >= 0 && y < 1)
 	    return (y < 0.5 ? factor * pow(2 * y, lambda.p()) : -factor * pow(2 * y-1, lambda.p()));
 	}
    }
    return 0.;
  }
  
  /*template <>
  double
  QuarkletFrame<1,3>::evaluate(const unsigned int derivative,
			  const RQIndex& lambda,
			  const double x)
  {
    if (derivative == 0) {
      const double y = (1<<lambda.j())*x-lambda.k();
      const double factor = twotothejhalf(lambda.j());
      if (lambda.e() == 0) // generator
	return (y >= 0 && y < 1 ? factor : 0.);
      else // wavelet
	{
	  if (y >= -1 && y < 2)
	    {
	      switch ((int) floor(y))
		{
		case -1:
		  return -0.125 * factor;
		case 0:
		  return (y < 0.5 ? factor : -factor);
		case 1:
		  return 0.125 * factor;
		}
	    }
	}
    }
    return 0.;
  }
  
  template <>
  double
  QuarkletFrame<1,5>::evaluate(const unsigned int derivative,
			  const RQIndex& lambda,
			  const double x)
  {
    if (derivative == 0) {
      const double y = (1<<lambda.j())*x-lambda.k();
      const double factor = twotothejhalf(lambda.j());     
      if (lambda.e() == 0) // generator
	return (y >= 0 && y < 1 ? factor : 0.);
      else // wavelet
	{
	  if (y >= -2 && y < 3)
	    {
	      switch((int) floor(y))
		{
		case -2:
		  return 0.0234375 * factor;  // 3/128
		case -1:
		  return -0.171875 * factor;  // -11/64
		case 0:
		  return (y < 0.5 ? factor : -factor);
		case 1:
		  return 0.171875 * factor;   // 11/64
		case 2:
		  return -0.0234375 * factor; // -3/128
		}
	    }
	}
    }
    return 0.;
  }*/
  
//  template <>
//  double
//  QuarkletFrame<2,2>::evaluate(const unsigned int derivative,
// 			  const RQIndex& lambda,
// 			  const double x)
//  {
//    const int j = lambda.j();
//    const int p = lambda.p();
//    const double y = (1<<j)*x-lambda.k();
//    if (lambda.e() == 0) // generator
//      {
//	if (y >= -1 && y < 1) {
//	  switch (derivative) {
//	  case 0: {
//	    const double factor = twotothejhalf(j);
//	    return (y < 0 ? pow(y,p)*factor*(y+1) : 
//		    pow(y,p)*factor*(-y+1));
//	  }
//	  case 1: {
//	    const double factor = twotothejhalf(3*j);
//	    return (y < 0 ? factor*(p*pow(y,p-1)+(p+1)*pow(y,p)) : 
//                    factor*(p*pow(y,p-1)-(p+1)*pow(y,p)));
//	  }
//	  }
//	}
//      }
//    else // wavelet noch auf neues quarklet-setting anpassen @PHK
//      {
//	if (y >= -1 && y < 2) {
//	  switch (derivative) {
//	  case 0: {
//	    const double factor = twotothejhalf(j);
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor*(0.5*y+0.5);
//	    case 0:
//	      return (y < 0.5 ? factor*(-4*y+0.5) : factor*(4*y-3.5));
//	    case 1:
//	      return factor*(-0.5*y+1);
//	    }
//	  }
//	  case 1: {
//	    const double factor = twotothejhalf(3*j);
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor*0.5;
//	    case 0:
//	      return (y < 0.5 ? -factor*4 : factor*4);
//	    case 1:
//	      return -factor*0.5;
//	    }
//	  }
//	  }
//	}
//      }
//      return 0.;
//      }
  
  /*template <>
  double
  QuarkletFrame<2,4>::evaluate(const unsigned int derivative,
 			  const RQIndex& lambda,
 			  const double x)
  {
    const int j = lambda.j();
    const double y = (1<<j)*x-lambda.k();
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y < 1) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    return (y < 0 ? factor*(y+1) : factor*(-y+1));
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    return (y < 0 ? factor : -factor);
	  }
	  }
	}
      }
    else // wavelet
      {
	if (y >= -2 && y < 3) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);	    
	    switch((int) floor(y)) {
	    case -2:
	      return factor*(-0.09375*y-0.1875); // -(3/32)*y-3/16
	    case -1:
	      return factor*(0.6875*y+0.59375);  // (11/16)*y+19/32
	    case 0:
	      return (y < 0.5
		      ? factor*(-4*y+0.59375)    // -4*y+19/32
		      : factor*(4*y-3.40625));   // 4*y-109/32
	    case 1:
	      return factor*(-0.6875*y+1.28125); // -(11/16)*y+41/32
	    case 2:
	      return factor*(0.09375*y-0.28125); // (3/32)*y-9/32
	    }
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);	    
	    switch((int) floor(y)) {
	    case -2:
	      return -(1<<j)*factor*0.09375;
	    case -1:
	      return (1<<j)*factor*0.6874;
	    case 0:
	      return (y < 0.5 ? -(1<<j)*factor*4 : (1<<j)*factor*4);
	    case 1:
	      return -(1<<j)*factor*0.6874;
	    case 2:
	      return (1<<j)*factor*0.09375;
	    }
	  }
	  }
	}
      }
    return 0.;
  }*/

  
//  template <>
//  double
//  QuarkletFrame<3,3>::evaluate(const unsigned int derivative,
// 			  const RQIndex& lambda,
// 			  const double x)
//  {
//    const int j = lambda.j();
//    const int p = lambda.p();
//    const double y = (1<<j)*x-lambda.k();
//    if (lambda.e() == 0) // generator
//      {
//	if (y >= -1 && y < 2) {
//	  switch (derivative) {
//	  case 0: {
//	    const double factor = twotothejhalf(j) * pow(0.5 * y, p);
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor*(0.5*y*y+y+0.5);
//	    case 0:
//	      return factor*(-y*y+y+0.5);
//	    case 1:
//	      return factor*(0.5*y*y-2*y+2);
//	    }
//	  }
//	  case 1: {
//           if(p==0){
//	    const double factor = twotothejhalf(3*j) ;
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor*(y+1);
//	    case 0:
//	      return factor*(-2*y+1);
//	    case 1:
//	      return factor*(y-2);
//            }
//           }
//           else{
//            const double factor = twotothejhalf(3*j) * pow(0.5 * y, p-1) * 0.5;
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor*(y*y*(0.5*p+1)+y*(1+p)+0.5*p);
//	    case 0:
//	      return factor*(-y*y*(p+2)+y*(1+p)+0.5*p);
//	    case 1:
//	      return factor*(y*y*(0.5*p+1)-2*y*(1+p)+2*p);
//            }
//           }
//	  }
//	  /*case 2: {
//	    const double factor = twotothejhalf(5*j);
//	    switch((int) floor(y)) {
//	    case -1:
//	      return factor;
//	    case 0:
//	      return -2*factor;
//	    case 1:
//	      return factor;
//	    }
//	  }*/
//	  }
//	}
//      }
//    else // wavelet
//      {
//	if (y >= -2 && y < 3) {
//	  switch (derivative) {
//	  case 0: {
//	    const double factor = twotothejhalf(j);
//	    switch((int) floor(2*y)) {
//            case -4:
//              return factor*(-0.09375* pow(y+1.5, p)*((y+1.5)*(2*y+3)+2*y+3.5));
//            case -3:
//              return factor*(-0.28125* pow(y+1, p)*((y+1)*(2*y+2)+2*y+2.5)
//                      -0.09375* pow(y+1.5, p)*((-2*y-3)*(2*y+3)+2*y+3.5));
//	    case -2:
//	      return factor*(0.21875* pow(y+0.5, p)*((y+0.5)*(2*y+1)+2*y+1.5)
//                      -0.28125* pow(y+1, p)*((-2*y-2)*(2*y+2)+2*y+2.5)
//                      -0.09375* pow(y+1.5, p)*((y+1.5)*(2*y+3)-4*y-4));
//	    case -1:
//	      return factor*(1.40625* pow(y, p)*(y*2*y+2*y+0.5)
//                      +0.21875* pow(y+0.5, p)*((-2*y-1)*(2*y+1)+2*y+1.5)
//                      -0.28125* pow(y+1, p)*((y+1)*(2*y+2)-4*y-2));
//	    case 0:
//	      return factor*(-1.40625* pow(y-0.5, p)*((y-0.5)*(2*y-1)+2*y-0.5)
//                      +1.40625* pow(y, p)*(-4*y*y+2*y+0.5)
//                      +0.21875* pow(y+0.5, p)*((y+0.5)*(2*y+1)-4*y));
//	    case 1:
//	      return factor*(-0.21875* pow(y-1, p)*((y-1)*(2*y-2)+2*y-1.5)
//                      -1.40625* pow(y-0.5, p)*((-2*y+1)*(2*y-1)+2*y-0.5)
//                      +1.40625* pow(y, p)*(y*2*y-4*y+2));
//	    case 2:
//	      return factor*(0.28125* pow(y-1.5, p)*((y-1.5)*(2*y-3)+2*y-2.5)
//                      -0.21875* pow(y-1, p)*((-2*y+2)*(2*y-2)+2*y-1.5)
//                      -1.40625* pow(y-0.5, p)*((y-0.5)*(2*y-1)-4*y+4));
//            case 3:
//              return factor*(0.09375* pow(y-2, p)*((y-2)*(2*y-4)+2*y-3.5)
//                      +0.28125* pow(y-1.5, p)*((-2*y+3)*(2*y-3)+2*y-2.5)
//                      -0.21875* pow(y-1, p)*((y-1)*(2*y-2)-4*y+6));
//            case 4: 
//              return factor*(0.09375* pow(y-2, p)*((-2*y+4)*(2*y-4)+2*y-3.5)
//                      +0.28125* pow(y-1.5, p)*((y-1.5)*(2*y-3)-4*y+8));
//            case 5:
//              return factor*(0.09375* pow(y-2, p)*((y-2)*(2*y-4)-4*y+10));
//                            
//	    }
//	  }
//	  case 1: {
//           if(p == 0){
//              const double factor = twotothejhalf(3*j);
//	    switch((int) floor(y)) {
//	    case -2:
//	      return factor*(-0.375*y-0.75);
//	    case -1:
//	      return factor*(2.75*y+2.375);
//	    case 0:
//	      return (y < 0.5
//		      ? factor*(-16*y+2.375)
//		      : factor*(16*y-13.625));
//	    case 1:
//	      return factor*(-2.75*y+5.125);
//	    case 2:
//	      return factor*(0.375*y-1.125);
//	    }    
//           }
//           else{
//	    const double factor = twotothejhalf(3*j);
//	    switch((int) floor(2*y)) {
//            case -4:
//              return factor*(-0.09375* pow(y+1.5, p-1)*(p*((y+1.5)*(2*y+3)+2*y+3.5)+(2*y+3)*(2*y+4)));
//            case -3:
//              return factor*(-0.28125* pow(y+1, p-1)*(p*((y+1)*(2*y+2)+2*y+2.5)+(2*y+2)*(2*y+3))
//                      -0.09375* pow(y+1.5, p-1)*(p*((-2*y-3)*(2*y+3)+2*y+3.5)+(2*y+3)*(-4*y-5)));
//	    case -2:
//	      return factor*(0.21875* pow(y+0.5, p-1)*(p*((y+0.5)*(2*y+1)+2*y+1.5)+(2*y+1)*(2*y+2))
//                      -0.28125* pow(y+1, p-1)*(p*((-2*y-2)*(2*y+2)+2*y+2.5)+(2*y+2)*(-4*y-3))
//                      -0.09375* pow(y+1.5, p-1)*(p*((y+1.5)*(2*y+3)-4*y-4)+(2*y+3)*(2*y+1)));
//	    case -1:
//	      return factor*(1.40625* pow(y, p-1)*(p*(y*2*y+2*y+0.5)+2*y*(2*y+1))
//                      +0.21875* pow(y+0.5, p-1)*(p*((-2*y-1)*(2*y+1)+2*y+1.5)+(2*y+1)*(-4*y-1))
//                      -0.28125* pow(y+1, p-1)*(p*((y+1)*(2*y+2)-4*y-2)+(2*y+2)*(2*y)));
//            case 0:
//	      return factor*(-1.40625* pow(y-0.5, p-1)*(p*((y-0.5)*(2*y-1)+2*y-0.5)+(2*y-1)*2*y)
//                      +1.40625* pow(y, p-1)*(p*(-4*y*y+2*y+0.5)+(2*y)*(-4*y+1))
//                      +0.21875* pow(y+0.5, p-1)*(p*((y+0.5)*(2*y+1)-4*y)+(2*y+1)*(2*y-1)));
//            case 1:
//	      return factor*(-0.21875* pow(y-1, p-1)*(p*((y-1)*(2*y-2)+2*y-1.5)+(2*y-2)*(2*y-1))
//                      -1.40625* pow(y-0.5, p-1)*(p*((-2*y+1)*(2*y-1)+2*y-0.5)+(2*y-1)*(-4*y+3))
//                      +1.40625* pow(y, p-1)*(p*(y*2*y-4*y+2)+2*y*(2*y-2)));
//	    case 2:
//	      return factor*(0.28125* pow(y-1.5, p-1)*(p*((y-1.5)*(2*y-3)+2*y-2.5)+(2*y-3)*(2*y-2))
//                      -0.21875* pow(y-1, p-1)*(p*((-2*y+2)*(2*y-2)+2*y-1.5)+(2*y-2)*(-4*y+5))
//                      -1.40625* pow(y-0.5, p-1)*(p*((y-0.5)*(2*y-1)-4*y+4)+(2*y-1)*(2*y-3)));
//            case 3:
//              return factor*(0.09375* pow(y-2, p-1)*(p*((y-2)*(2*y-4)+2*y-3.5)+(2*y-4)*(2*y-3))
//                      +0.28125* pow(y-1.5, p-1)*(p*((-2*y+3)*(2*y-3)+2*y-2.5)+(2*y-3)*(-4*y+7))
//                      -0.21875* pow(y-1, p-1)*(p*((y-1)*(2*y-2)-4*y+6)+(2*y-2)*(2*y-4)));
//            case 4: 
//              return factor*(0.09375* pow(y-2, p-1)*(p*((-2*y+4)*(2*y-4)+2*y-3.5)+(2*y-4)*(-4*y+9))
//                      +0.28125* pow(y-1.5, p-1)*(p*((y-1.5)*(2*y-3)-4*y+8)+(2*y-3)*(2*y-5)));
//            case 5:
//              return factor*(0.09375* pow(y-2, p-1)*(p*((y-2)*(2*y-4)-4*y+10)+(2*y-4)*(2*y-6)));
//                            
//	    }
//           }
//	  }
//	  /*case 2: {
//	    const double factor = twotothejhalf(5*j);
//	    switch((int) floor(y)) {
//	    case -2:
//	      return -factor*0.375;
//	    case -1:
//	      return factor*2.75;
//	    case 0:
//	      return (y < 0.5 ? -factor*16 : factor*16);
//	    case 1:
//	      return -factor*2.75;
//	    case 2:
//	      return factor*0.375;
//	    }
//	  }*/
//	  }
//	}
//      }
//    return 0.;
//  }

  /*template <>
  double
  QuarkletFrame<3,5>::evaluate(const unsigned int derivative,
 			  const RQIndex& lambda,
 			  const double x)
  {
    const int j = lambda.j();
    const double y = (1<<j)*x-lambda.k();
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y < 2) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    switch((int) floor(y)) {
	    case -1:
	      return factor*(0.5*y*y+y+0.5);
	    case 0:
	      return factor*(-y*y+y+0.5);
	    case 1:
	      return factor*(0.5*y*y-2*y+2);
	    }
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    switch((int) floor(y)) {
	    case -1:
	      return factor*(y+1);
	    case 0:
	      return factor*(-2*y+1);
	    case 1:
	      return factor*(y-2);
	    }
	  }
	  case 2: {
	    const double factor = twotothejhalf(5*j);
	    switch((int) floor(y)) {
	    case -1:
	      return factor;
	    case 0:
	      return -2*factor;
	    case 1:
	      return factor;
	    }
	  }
	  }
	}
      }
    else // wavelet
      {
	if (y >= -3 && y < 4) {
 	  switch (derivative) {
 	  case 0: {
 	    const double factor = twotothejhalf(j);
 	    switch((int) floor(y)) {
 	    case -3:
 	      return factor*(0.0390625*y*y+0.234375*y+0.3515625);  // (5/128)*y^2+(15/64)*y+(45/128)
 	    case -2:
 	      return factor*(-0.34375*y*y-1.296875*y-1.1796875);   // -(11/32)*y^2-(83/64)*y-(151/128)
	    case -1:
 	      return factor*(1.5703125*y*y+2.53125*y+0.734375);    // (201/128)*y^2+(81/32)*y+(47/64)
 	    case 0:
 	      return (y < 0.5
 		      ? factor*(-8*y*y+2.53125*y+0.734375)         // -8*y^2+(81/32)*y+(47/64)
 		      : factor*(8*y*y-13.46875*y+4.734375));       // 8*y^2-(431/32)*y+(303/64)
 	    case 1:
 	      return factor*(-1.5703125*y*y+5.671875*y-4.8359375); // -(201/128)*y^2+(363/64)*y-(619/128)
 	    case 2:
 	      return factor*(0.34375*y*y-1.984375*y+2.8203125);    // (11/32)*y^2-(127/64)*y+(361/128)
 	    case 3:
 	      return factor*(-0.0390625*y*y+0.3125*y-0.625);       // (-5/128)*y^2+(5/16)*y-(5/8)
 	    }
	  }
 	  case 1: {
 	    const double factor = twotothejhalf(3*j);
 	    switch((int) floor(y)) {
 	    case -3:
 	      return factor*(0.078125*y+0.234375);
 	    case -2:
 	      return factor*(-0.6875*y-1.296875);
	    case -1:
 	      return factor*(3.140625*y+2.53125);
 	    case 0:
 	      return (y < 0.5
 		      ? factor*(-16*y+2.53125)
 		      : factor*(16*y-13.46875));
 	    case 1:
 	      return factor*(-3.140625*y+5.671875);
 	    case 2:
 	      return factor*(0.6875*y-1.984375);
 	    case 3:
 	      return factor*(-0.078125*y+0.3125);
 	    }
 	  }
 	  case 2: {
 	    const double factor = twotothejhalf(5*j);
 	    switch((int) floor(y)) {
 	    case -3:
 	      return factor*0.078125;
 	    case -2:
 	      return -factor*0.6875;
	    case -1:
 	      return factor*3.140625;
 	    case 0:
 	      return (y < 0.5 ? -factor*16 : factor*16);
 	    case 1:
 	      return -factor*3.140625;
 	    case 2:
 	      return factor*0.6875;
 	    case 3:
 	      return -factor*0.078125;
 	    }
 	  }
 	  }
	}
      }
    return 0.;
    }*/
  
  template <int d, int dt>
  double
  QuarkletFrame<d, dt>:: integrate(const RQIndex& lambda){
      
    double r = 0;

    // first we compute the support of psi_lambda
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(lambda, k1, k2);
    const int length = k2 > k1 ? k2-k1 : k2+(1<<j)-k1;
    //cout << k1 << ", " << k2 << endl;
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = 1.0/(1<<j);

    Array1D<double> gauss_points (N_Gauss*length);
    int k = k1;
    for (int patch = 0; patch < length; patch++, k++){ // work on 2^{-j}[k,k+1]
        //cout << "k: " << k << endl; 
        for (unsigned int n = 0; n < N_Gauss; n++)
 	gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
    }
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
 	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
 	for (int patch = k1; patch < k1+length; patch++)
 	  {
 	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
            //cout << t << endl;
	    
 	    r += evaluate(0, lambda, t)
 		 * gauss_weight;
 	  }
      }
    
    return r;
      
  }
  
}

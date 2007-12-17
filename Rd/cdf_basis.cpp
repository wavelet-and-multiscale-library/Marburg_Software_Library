// implementation of cdf_basis.h

#include <cassert>
#include <numerics/cardinal_splines.h>
#include <utils/tiny_tools.h>

namespace WaveletTL
{
  template <int d, int dt>
  double
  CDFBasis<d,dt>::evaluate(const unsigned int derivative,
			   const RIndex& lambda,
			   const double x) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) // generator
      {
	r = (derivative == 0
	     ? MathTL::EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x)
	     : MathTL::EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), x));
      }
    else // wavelet
      {
	// old generic version: readable and slow
 	InfiniteVector<double, RIndex> gcoeffs;
	RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::reconstruct_1
	  (lambda, lambda.j()+1, gcoeffs);
 	for (typename InfiniteVector<double, RIndex>::const_iterator it(gcoeffs.begin());
 	     it != gcoeffs.end(); ++it)
 	  r += *it * evaluate(derivative, it.index(), x);
      }
    
    return r;
  }

  //
  // some fast template specializations without reconstruct_1() calls

  template <>
  double
  CDFBasis<1,1>::evaluate(const unsigned int derivative,
 			  const RIndex& lambda,
 			  const double x) const
  {
    if (derivative == 0) {
      const double y = (1<<lambda.j())*x-lambda.k();
      const double factor = twotothejhalf(lambda.j());
      if (lambda.e() == 0) // generator
 	return (y >= 0 && y < 1 ? factor : 0.);
      else // wavelet
 	{
 	  if (y >= 0 && y < 1)
 	    return (y < 0.5 ? factor : -factor);
 	}
    }
    return 0.;
  }
  
  template <>
  double
  CDFBasis<1,3>::evaluate(const unsigned int derivative,
			  const RIndex& lambda,
			  const double x) const
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
  CDFBasis<1,5>::evaluate(const unsigned int derivative,
			  const RIndex& lambda,
			  const double x) const
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
  }

  template <>
  double
  CDFBasis<2,2>::evaluate(const unsigned int derivative,
 			  const RIndex& lambda,
 			  const double x) const
  {
    const int j = lambda.j();
    const double y = (1<<j)*x-lambda.k();
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y < 1) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    return (y <= 0 ? factor*(y+1) : factor*(-y+1));
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    return (y <= 0 ? factor : -factor);
	  }
	  }
	}
      }
    else // wavelet
      {
	if (y >= -1 && y < 2) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    switch((int) floor(y)) {
	    case -1:
	      return factor*(0.5*y+0.5);
	    case 0:
	      return (y < 0.5 ? factor*(-4*y+0.5) : factor*(4*y-3.5));
	    case 1:
	      return factor*(-0.5*y+1);
	    }
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    switch((int) floor(y)) {
	    case -1:
	      return factor*0.5;
	    case 0:
	      return (y < 0.5 ? -factor*4 : factor*4);
	    case 1:
	      return -factor*0.5;
	    }
	  }
	  }
	}
      }
    return 0.;
  }
  
  template <>
  double
  CDFBasis<2,4>::evaluate(const unsigned int derivative,
 			  const RIndex& lambda,
 			  const double x) const
  {
    const int j = lambda.j();
    const double y = (1<<j)*x-lambda.k();
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y < 1) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    return (y <= 0 ? factor*(y+1) : factor*(-y+1));
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    return (y <= 0 ? factor : -factor);
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
	      return factor*(0.6874*y+0.59375);  // (11/16)*y+19/32
	    case 0:
	      return (y < 0.5
		      ? factor*(-4*y+0.59375)    // -4*y+19/32
		      : factor*(4*y-3.40625));   // 4*y-109/32
	    case 1:
	      return factor*(-0.6874*y+1.28125); // -(11/16)*y+41/32
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
  }
  
  template <>
  double
  CDFBasis<3,3>::evaluate(const unsigned int derivative,
 			  const RIndex& lambda,
 			  const double x) const
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
	if (y >= -2 && y < 3) {
	  switch (derivative) {
	  case 0: {
	    const double factor = twotothejhalf(j);
	    switch((int) floor(y)) {
	    case -2:
	      return factor*(-0.1875*y*y-0.75*y-0.75);   // -(3/16)*y^2-(3/4)*y-(3/4)
	    case -1:
	      return factor*(1.375*y*y+2.375*y+0.8125);  // (11/8)*y^2+(19/8)*y+(13/16)
	    case 0:
	      return (y < 0.5
		      ? factor*(-8*y*y+2.375*y+0.8125)   // -8*y^2+(19/8)*y+(13/16)
		      : factor*(8*y*y-13.625*y+4.8125)); // 8*y^2-(109/8)*y+(77/16)
	    case 1:
	      return factor*(-1.375*y*y+5.125*y-4.5625); // -(11/8)*y^2+(41/8)*y-(73/16)
	    case 2:
	      return factor*(0.1875*y*y-1.125*y+1.6875); // (3/16)*y^2-(9/8)*y+(27/16)
	    }
	  }
	  case 1: {
	    const double factor = twotothejhalf(3*j);
	    switch((int) floor(y)) {
	    case -2:
	      return factor*(-0.375*y-0.75);
	    case -1:
	      return factor*(2.75*y+2.375);
	    case 0:
	      return (y < 0.5
		      ? factor*(-16*y+2.375)
		      : factor*(16*y-13.625));
	    case 1:
	      return factor*(-2.75*y+5.125);
	    case 2:
	      return factor*(0.375*y-1.125);
	    }
	  }
	  case 2: {
	    const double factor = twotothejhalf(5*j);
	    switch((int) floor(y)) {
	    case -2:
	      return -factor*0.375;
	    case -1:
	      return factor*2.75;
	    case 0:
	      return (y < 0.5 ? -factor*16 : factor*16);
	    case 1:
	      return -factor*2.75;
	    case 2:
	      return factor*0.375;
	    }
	  }
	  }
	}
      }
    return 0.;
  }
  
  template <>
  double
  CDFBasis<3,5>::evaluate(const unsigned int derivative,
 			  const RIndex& lambda,
 			  const double x) const
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
  }
  
}

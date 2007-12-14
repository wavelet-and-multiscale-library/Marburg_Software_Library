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
    const double factor = twotothejhalf(j);
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y <= 1) {
	  switch (derivative) {
	  case 0:
	    return (y <= 0 ? factor*(y+1) : factor*(-y+1));
	  case 1:
	    return (y <= 0 ? (1<<j)*factor : -(1<<j)*factor);
	  }
	}
      }
    else // wavelet
      {
	if (y >= -1 && y <= 2) {
	  switch (derivative) {
	  case 0:
	    switch((int) floor(y)) {
	    case -1:
	      return factor*(0.5*y+0.5);
	    case 0:
	      return (y < 0.5 ? factor*(-4*y+0.5) : factor*(4*y-3.5));
	    case 1:
	      return factor*(-0.5*y+1);
	    }
	  case 1:
	    switch((int) floor(y)) {
	    case -1:
	      return (1<<j)*factor*0.5;
	    case 0:
	      return (y < 0.5 ? -(1<<j)*factor*4 : (1<<j)*factor*4);
	    case 1:
	      return -(1<<j)*factor*0.5;
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
    const double factor = twotothejhalf(j);
    if (lambda.e() == 0) // generator
      {
	if (y >= -1 && y <= 1) {
	  switch (derivative) {
	  case 0:
	    return (y <= 0 ? factor*(y+1) : factor*(-y+1));
	  case 1:
	    return (y <= 0 ? (1<<j)*factor : -(1<<j)*factor);
	  }
	}
      }
    else // wavelet
      {
	if (y >= -2 && y <= 3) {
	  switch (derivative) {
	  case 0:
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
	  case 1:
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
    return 0.;
  }
  
}

// implementation for frame_support.h


#include <utils/fixed_array1d.h>
#include <typeinfo>
#include <cube/cube_support.h>
#include <map>
#include <time.h>
#include <math.h>
#include <iostream>

using namespace WaveletTL;

using std::type_info;
using std::map;


namespace FrameTL
{

//        Attention the calculation of the support between two Patches only works if the cube is mapped on a rectangle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  // #####################################################################################
  // First some helpers.
  // #####################################################################################

  // check two doubles for equality
  inline bool eq (const double x, const double y)
  {
    const double eps = 1.0e-12;
    return fabs(x-y) < eps;
  }

  // relation "less than"
  inline bool lt (const double x, const double y)
  {
    return x < y;
//     const double eps = 1.0e-12;
//     return x-y < eps;
  }

  // relation "greater than"
  inline bool gt (const double x, const double y)
  {
    return ! (lt(x,y) || eq(x,y));
  }

  // relation "less than or equal"
  inline bool leq (const double x, const double y)
  {
    return lt(x,y) || eq(x,y);
  }
  // #####################################################################################


  // This routine determines the location of the point
  // p relative to the line given by the points p1 and p2.
  // We compute the determinant of the matrix (p2-p1, p-p1).

  // If the result is negative, then p lies on the right of the line,
  // when looking from point p1 into the direction of p2. We return 0.

  // If the result is positive, then p lies on the left of the line,
  // when looking from point p1 into the direction of p2. We return 1.

  // If the result is 0, then p lies on the line. We return 2.
  template <unsigned int DIM>
  inline
  unsigned short int pos_wrt_line (const Point<DIM>& p,
				   const Point<DIM>& p1, const Point<DIM>&  p2)
  {
    assert ( DIM == 2 );
    double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    if (fabs(d) < 1.0e-15)
      return 2;
    if( d > 0.0 )
      return 1;
    else if (d < 0.0)
      return 0;
    // dummy
    return 3;
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const Point<DIM_m>& p)
  {

    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
	
    typename CUBEBASIS::Support supp_lambda;

    // create a CubeIndex
    typename CUBEBASIS::Index lambda_c(lambda.j(),
				       lambda.e(),
				       lambda.k(),frame.bases()[lambda.p()]);

    // get the support of the reference wavelet corresponding to the index lambda
    WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[lambda.p()], lambda_c, supp_lambda);

    // granularity of the dyadic grid with respect to which the current wavelet is
    // a piecewise polynomial
    const double dx = 1.0 / (1 << supp_lambda.j);
    
    // Determine at runtime which type of parametric mappings are used.
    // Start with the case of a one-dimensional AffineLinearMapping.
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>))
	 )
      {
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx), b);
	//cout << a[0] << " " << p[0] << " " <<  b[0] << " " << (leq(a[0],p[0]) && leq(p[0],b[0])) << endl;
	//return (a[0] <= p[0]) && (p[0] <= b[0]);
	//return leq(a[0],p[0]) && leq(p[0],b[0]);

	// At the moment, if the point lies on the boundary of the support,
	// we return false. This should not be the case when free boundary
	// conditions are chosen!!!
	return (a[0] < p[0]) && (p[0] < b[0]);
	
      }

    // The case of a two-dimensional AffineLinearMapping.
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>))
	 )
      {
	const double dx = 1.0 / (1 << supp_lambda.j);
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx,supp_lambda.a[1]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx,supp_lambda.b[1]*dx), b);
	return (a[0] <= p[0] && p[0] <= b[0]) && (a[1] <= p[1] && p[1] <= b[1]);

      }

    // The case of a two-dimensional LinearBezierMapping.
    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(LinearBezierMapping))
      {
	
	// this array represents the quadrangle representing the support
	// of the frame element
	FixedArray1D<Point<DIM_d>,4 > poly;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx,supp_lambda.a[1]*dx), poly[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx,supp_lambda.a[1]*dx), poly[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx,supp_lambda.b[1]*dx), poly[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx,supp_lambda.b[1]*dx), poly[3]);


	//make sure to walk through the vertices counter clockwise!!!
	unsigned short int res = 
	  FrameTL::pos_wrt_line(p, poly[0], poly[1]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[1], poly[2]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[2], poly[3]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[3], poly[0]);
	if (res == 0)
	  return false;
	
	return true;

      }
    return true;
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
		  const Point<DIM_m>& p)
  {

    // Determine at runtime which type of parametric mappings are used.
    // Start with the case of a one-dimensional AffineLinearMapping.
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>))
	 )
      {
	const double dx = 1.0 / (1 << supp_lambda->j);
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx), b);
	//return (a[0] < p[0]) && (p[0] < b[0]);
	return (leq(a[0],p[0]) && leq(p[0],b[0]));
	
      }

    // The case of a two-dimensional AffineLinearMapping.
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>))
	 )
      {
	const double dx = 1.0 / (1 << supp_lambda->j);
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx,supp_lambda->a[1]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx,supp_lambda->b[1]*dx), b);
	return (a[0] <= p[0] && p[0] <= b[0]) && (a[1] <= p[1] && p[1] <= b[1]);

      }
    
    // The case of a two-dimensional LinearBezierMapping.
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(LinearBezierMapping))
	 )
      {
	assert ( DIM_d == 2 && DIM_m == 2 );
	
	const double dx = 1.0 / (1 << supp_lambda->j);

	FixedArray1D<Point<DIM_m>,4 > poly;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx,supp_lambda->a[1]*dx), poly[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx,supp_lambda->a[1]*dx), poly[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx,supp_lambda->b[1]*dx), poly[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx,supp_lambda->b[1]*dx), poly[3]);

	//make sure to walk through the vertices counter clockwise!!!
	unsigned short int res = 
	  FrameTL::pos_wrt_line(p, poly[0], poly[1]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[1], poly[2]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[2], poly[3]);
	if (res == 0)
	  return false;
	res = 
	  FrameTL::pos_wrt_line(p, poly[3], poly[0]);
	if (res == 0)
	  return false;
	
	return true;

      }
    return true;
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool
  intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
		     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
		     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu)
  {    
    const double dx1 = 1.0 / (1 << supp_lambda->j);
    const double dx2 = 1.0 / (1 << supp_mu->j);

    // Determine at runtime which type of parametric mappings are used.
    // Start with the case of a one-dimensional AffineLinearMapping.
    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(AffineLinearMapping<1>) &&
	 typeid(*frame.atlas()->charts()[mu.p()])     == 
	 typeid(AffineLinearMapping<1>) )
      {
	assert ( DIM_d == 1 );
	
	Point<DIM_d> a;
	Point<DIM_d> b;

	Point<DIM_d> c;
	Point<DIM_d> d;

	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1), b);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->a[0]*dx2), c);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->b[0]*dx2), d);

	return (a[0] < d[0]) && (b[0] > c[0]);

      }

    // The case of a two-dimensional AffineLinearMapping or SimpleAffineLinearMapping.
    if ( ( 
	  (typeid(*frame.atlas()->charts()[lambda.p()])  == typeid(AffineLinearMapping<2>))
	  &&
	  (typeid(*frame.atlas()->charts()[mu.p()])      == typeid(AffineLinearMapping<2>))
	  )
	 ||
	 (
	  (typeid(*frame.atlas()->charts()[lambda.p()])  == typeid(SimpleAffineLinearMapping<2>))
	  &&
	  (typeid(*frame.atlas()->charts()[mu.p()])      == typeid(SimpleAffineLinearMapping<2>))
	  )
	 )
      {
	const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda.p()];
	const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu.p()];
	
	Point<DIM_m> a_la, b_la, a_mu, b_mu;
	
	chart_la->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->a[1]*dx1), a_la);
	chart_la->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->b[1]*dx1), b_la);
	
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), a_mu);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), b_mu);
	
	return (a_la[0] < b_mu[0] && b_la[0] > a_mu[0]) && (a_la[1] < b_mu[1] && b_la[1] > a_mu[1]);
      }
    
    // The case of a two-dimensional LinearBezierMapping.
    if ( ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(LinearBezierMapping) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(LinearBezierMapping) )
	 )
      {

	assert ( DIM_d == 2 && DIM_m == 2);

	// these arrays represent the quadrangles representing the supports
	// of the frame elements 
	FixedArray1D<Point<DIM_m>,4 > poly1;
	FixedArray1D<Point<DIM_m>,4 > poly2;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->a[1]*dx1), poly1[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->a[1]*dx1), poly1[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->b[1]*dx1), poly1[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->b[1]*dx1), poly1[3]);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), poly2[0]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->a[1]*dx2), poly2[1]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), poly2[2]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->b[1]*dx2), poly2[3]);

	// now check if the two quadrangles have a non-trivial intersection
	return quadrangles_intersect (poly1,poly2);
      }
    else
      return false;
  }


  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void
  precompute_supports_simple(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
			     Array1D<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support>& all_patch_supports)
  {
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;

    // loop over all degrees of freedom
    for (int i = 0; i < frame->degrees_of_freedom(); i++) {

      // get the Index corresponding to number i
      const Index* ind = frame->get_wavelet(i);

      // get the support of the reference wavelet on the unit cube
      const typename CubeBasis<IBASIS,DIM_d>::Support* supp_ind = &((frame->all_supports)[ind->number()]);

      const double dx1 = 1.0 / (1 << supp_ind->j);
      const Chart<DIM_d,DIM_m>* chart_la = frame->atlas()->charts()[ind->p()];

      // map the reference support to the patch and store the result
      Point<DIM_m> a, b;
      if (DIM_m == 2) {
	chart_la->map_point(Point<DIM_d>(supp_ind->a[0]*dx1,supp_ind->a[1]*dx1), a);
	chart_la->map_point(Point<DIM_d>(supp_ind->b[0]*dx1,supp_ind->b[1]*dx1), b);
      }
      else if (DIM_m == 1) {
	chart_la->map_point(Point<DIM_d>(supp_ind->a[0]*dx1), a);
	chart_la->map_point(Point<DIM_d>(supp_ind->b[0]*dx1), b);
      }
      all_patch_supports[i].j = supp_ind->j;
      all_patch_supports[i].a = a;
      all_patch_supports[i].b = b;
    }
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool
  intersect_supports_simple(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			    const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			    const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu)
  {
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support SuppType;

    // We fetch the suppports of the two wavelets from the collection of
    // precomputed supports.
    const SuppType* supp_la = &(frame.all_patch_supports[lambda.number()]);
    const SuppType* supp_mu = &(frame.all_patch_supports[mu.number()]);

    // The 2D case.
    if (DIM_m == 2)
      return (supp_la->a[0] < supp_mu->b[0] && supp_la->b[0] > supp_mu->a[0]) && 
      (supp_la->a[1] < supp_mu->b[1] && supp_la->b[1] > supp_mu->a[1]);
    // The 1D case.
    else if (DIM_m == 1)
      return (supp_la->a[0] < supp_mu->b[0] && supp_la->b[0] > supp_mu->a[0]);

//       const double dx1 = 1.0 / (1 << supp_lambda->j);
//       const double dx2 = 1.0 / (1 << supp_mu->j);
//       const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda.p()];
//       const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu.p()];
      
//       Point<DIM_m> a_la, b_la, a_mu, b_mu;
      
//       chart_la->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->a[1]*dx1), a_la);
//       chart_la->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->b[1]*dx1), b_la);
      
//       chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), a_mu);
//       chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), b_mu);
      
//       return (a_la[0] < b_mu[0] && b_la[0] > a_mu[0]) && (a_la[1] < b_mu[1] && b_la[1] > a_mu[1]);
  }


 template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void
  support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				 typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support& SuppType)
  {
    SuppType = frame.all_patch_supports[lambda.number()];
  }


 template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void
  support_on_cube(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				 typename CubeBasis<IBASIS,DIM_d>::Support& SuppType)
  {
    SuppType = (frame.all_supports)[lambda.number()];
  }



  /*!
    This function checks whether the convex qudrangles given by the vertices in poly1 and poly2
    have a non-trivial intersection.
   */
  template <unsigned int DIM>
  bool quadrangles_intersect (FixedArray1D<Point<DIM>, 4> poly1, FixedArray1D<Point<DIM>, 4> poly2)
  {
    bool result;
    Point<DIM> p11;
    Point<DIM> p12;
	
    unsigned short int tmp = 0;
    int countspez = 0;
	
    //check, if two arbitrary edges of the quadrangle intersect
    for (unsigned int i = 1; i <= 3; i++) {	
      result = true;
      p11 = poly1[i-1];
      p12 = poly1[i];
      for (unsigned int j = 1; j <= 3; j++) {
	tmp = edgesIntersect(p11, p12, poly2[j-1], poly2[j]);
	if (tmp == 1)
	  countspez++;
	if (tmp == 3) {
	  return true;
	}
      }
      tmp = edgesIntersect(p11, p12, poly2[3], poly2[0]);
      if (tmp == 1)
	countspez++;
      if (tmp == 3)
	return true;
	
    }//end outer for
    p11 = poly1[3];
    p12 = poly1[0];
    for (unsigned int j = 1; j <= 3; j++) {
      tmp = edgesIntersect(p11, p12, poly2[j-1], poly2[j]);
      if (tmp == 1)
	countspez++;
      if (tmp == 3)
	return true;
    }
    tmp = edgesIntersect(p11, p12, poly2[3], poly2[0]);

    if (tmp == 1)
      countspez++;
    if (tmp == 3)
      return true;
    // if control reaches this point, then no 'ordinary intersection' exists
    // ATTENTION!!!!
    // this is only correct in case of convex polygons!!!!!!!
    if (countspez >= 2) {
      return true;
    }

    // remaining possible cases: one quadrangle contains the other
    // or they do not intersect
    countspez = 0;
    tmp = 0;
	
    //loop over all knots of the second polygon
    for (unsigned int i = 0; i <= 3; i++) {
      result = true;
      tmp = FrameTL::pos_wrt_line(poly2[i], poly1[3], poly1[0]);
      if (tmp == 1)
	countspez++;
      if (tmp == 0 || tmp ==2) {
	result = false;
	continue;
      }
      for (unsigned int j = 1; j<= 3 ; j++) {
	tmp = FrameTL::pos_wrt_line(poly2[i], poly1[j-1], poly1[j]);
	if (tmp == 1)
	  countspez++;
	if (tmp == 0 || tmp == 2) {
	  result = false;
	  break;
	}
      }//end inner for
      
      if (result)
	return result;
    }//end outer for
	
    //the same with reversed roles
    //loop over all knots of the first polygon
    for (unsigned int i = 0; i <= 3; i++) {
      result = true;
      tmp = FrameTL::pos_wrt_line(poly1[i], poly2[3], poly2[0]);
      if (tmp == 1)
	countspez++;
      if (tmp == 0 || tmp ==2) {
	result = false;
	continue;
      }
		
      for (unsigned int j = 1; j<= 3 ; j++) {
	tmp = FrameTL::pos_wrt_line(poly1[i], poly2[j-1], poly2[j]);
	if (tmp == 1)
	  countspez++;
	if (tmp == 0 || tmp ==2) {
	  result = false;
	  break;
	}
      }//end inner for
		
      if (result)
	return result;

    }//end outer for

    return false;

  }
  
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports_1D(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			     const Index1D<IBASIS>& lambda,
 			     const Index1D<IBASIS>& mu,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			     const int dir,
			     Array1D<double>& supp_intersect)
  {
    const double dx1 = 1.0 / (1 << supp_lambda->j);
    const double dx2 = 1.0 / (1 << supp_mu->j);

    double x;
    double x_patch;
    double y0;
    double y1;

    const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda.p()];
    const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu.p()];
    
    x = supp_mu->a[dir] * dx2;
      
    x_patch = chart_mu->map_point(x, dir);
    y0 = chart_la->map_point_inv(x_patch, dir); 

    x = supp_mu->b[dir] * dx2;

    x_patch = chart_mu->map_point(x, dir);
    y1 = chart_la->map_point_inv(x_patch, dir);
    
    // compute lower left and right end of  intersection interval
    FixedArray1D<double,2> hyperCube_intersect;

    hyperCube_intersect[0] = std::max(supp_lambda->a[dir]*dx1,y0);
    hyperCube_intersect[1] = std::min(supp_lambda->b[dir]*dx1,y1);

    if ( hyperCube_intersect[0] >= hyperCube_intersect[1] )
      return false;
    
    // we just need some kind of a sorted list
    std::list<double> irregular_grid;
    std::list<double> irregular_grid_tmp;

    // collect gridpoints for all directions,
    // we have to do this for both patches  
    for (int k = supp_lambda->a[dir]; k <= supp_lambda->b[dir]; k++) {
      double d = k*dx1;
      //if ( hyperCube_intersect[0] <= d && d <= hyperCube_intersect[1]) {
      if ( leq(hyperCube_intersect[0],d) && leq(d,hyperCube_intersect[1])) {
	irregular_grid.push_back(d);
      }
    }
    
    // left end of support interval
    x = supp_mu->a[dir]*dx2;  
    
    // second patch
    //for (unsigned int i = 0; i < DIM_d; i++) {
    for (int k = supp_mu->a[dir]; k <= supp_mu->b[dir]; k++) {
      x = k*dx2;
      x_patch = chart_mu->map_point(x, dir);
      y0 = chart_la->map_point_inv(x_patch, dir);

      //if ( hyperCube_intersect[0] <= y0 &&  y0 <= hyperCube_intersect[1]){
      if ( leq(hyperCube_intersect[0],y0) &&  leq(y0,hyperCube_intersect[1])){
	irregular_grid.push_back(y0);
      }
    }
    x = supp_mu->a[dir]*dx2;
    irregular_grid.sort();
    irregular_grid.unique();    
    
    supp_intersect.resize(irregular_grid.size());
    unsigned int j = 0;
    for (list<double>::const_iterator it = irregular_grid.begin();
	 it != irregular_grid.end(); ++it) {
      supp_intersect[j] = *it;
      j++;
    }

//     double old = -1.;
//     for (list<double>::const_iterator it = irregular_grid.begin();
// 	 it != irregular_grid.end(); ++it) {
//       if (! (fabs(old-*it) < 1.0e-13))
// 	irregular_grid_tmp.push_back(*it);
//       old = *it;
//     }
    
//     supp_intersect.resize(irregular_grid_tmp.size());
//     unsigned int j = 0;
//     for (list<double>::const_iterator it = irregular_grid_tmp.begin();
// 	 it != irregular_grid_tmp.end(); ++it) {
//       supp_intersect[j] = *it;
//       j++;
//     }
    return true;

  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			  FixedArray1D<Array1D<double>,DIM_d >& supp_intersect)
  {    
    const double dx1 = 1.0 / (1 << supp_lambda->j);
    const double dx2 = 1.0 / (1 << supp_mu->j);

    Point<DIM_d> x;
    Point<DIM_d> x_patch;
    Point<DIM_d> y0;
    Point<DIM_d> y1;

    const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda.p()];
    const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu.p()];
    
    for (unsigned int i = 0; i < DIM_d; i++)
      x[i] = supp_mu->a[i] * dx2;
      
    chart_mu->map_point(x,x_patch);
    chart_la->map_point_inv(x_patch,y0);
  
    for (unsigned int i = 0; i < DIM_d; i++) {
      x[i] = supp_mu->b[i] * dx2;
    }

    chart_mu->map_point(x,x_patch);
    chart_la->map_point_inv(x_patch,y1);

    // compute lower left and upper right corner of intersecting cube
    FixedArray1D<FixedArray1D<double,2>,DIM_d > hyperCube_intersect;
    for (unsigned int i = 0; i < DIM_d; i++) {
      assert ( y0[i] <= y1[i] );
      hyperCube_intersect[i][0] = std::max(supp_lambda->a[i]*dx1,y0[i]);
      hyperCube_intersect[i][1] = std::min(supp_lambda->b[i]*dx1,y1[i]);
      

      if ( hyperCube_intersect[i][0] >= hyperCube_intersect[i][1]  )
	return false;      
    }

#if 0
    for (unsigned int i = 0; i < DIM_d; i++) {
      cout << "i = " << i << " intersect_cube[i] = " << hyperCube_intersect[i] << endl;
    }
#endif    

    // we just need some kind of a sorted list
    // FixedArray1D<map<double,bool>, DIM_d> irregular_grid;
    FixedArray1D<list<double>, DIM_d> irregular_grid;
    FixedArray1D<list<double>, DIM_d> irregular_grid_tmp;

    // collect gridpoints for all directions,
    // we have to do this for both patches
    
    for (unsigned int i = 0; i < DIM_d; i++) {
      for (int k = supp_lambda->a[i]; k <= supp_lambda->b[i]; k++) {
	double d = k*dx1;
	if ( hyperCube_intersect[i][0] <= d && d <= hyperCube_intersect[i][1]) {
	  irregular_grid[i].push_back(d);
	}
      }
    }
    
    // lower left corner of support cube
    for (unsigned int j = 0; j < DIM_d; j++) {
      x[j] = supp_mu->a[j]*dx2;
    }
    
    // second patch
    for (unsigned int i = 0; i < DIM_d; i++) {
      for (int k = supp_mu->a[i]; k <= supp_mu->b[i]; k++) {
	x[i] = k*dx2;
	chart_mu->map_point(x,x_patch);
	chart_la->map_point_inv(x_patch,y0);
	if ( hyperCube_intersect[i][0] <= y0[i] &&  y0[i] <= hyperCube_intersect[i][1]){
	  irregular_grid[i].push_back(y0[i]);
	}
      }
      x[i] = supp_mu->a[i]*dx2;
      irregular_grid[i].sort();    
      
      double old = -1.;
      for (list<double>::const_iterator it = irregular_grid[i].begin();
	   it != irregular_grid[i].end(); ++it) {
	if (! (fabs(old-*it) < 1.0e-13))
	  irregular_grid_tmp[i].push_back(*it);
	old = *it;
      }
    }

    for (unsigned int i = 0; i < DIM_d; i++) {
      supp_intersect[i].resize(irregular_grid_tmp[i].size());
      unsigned int j = 0;
      for (list<double>::const_iterator it = irregular_grid_tmp[i].begin();
	   it != irregular_grid_tmp[i].end(); ++it) {
	supp_intersect[i][j] = *it;
	j++;
      }
    }
    return true;
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  void intersecting_wavelets (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			      const int j, const bool generators,
			      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting)
  {
    //intersecting.erase(intersecting.begin(),intersecting.end());

#if 1

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }
    std::list<typename Frame::Index> intersect_diff;


    if (! generators) {
      unsigned int p = lambda.p();
      FixedArray1D<int,DIM_d>
      minkwavelet, maxkwavelet, minkgen, maxkgen;
      typedef typename IBASIS::Index Index1D;
      int minkkkk;
      int maxkkkk;

      // prepare all intersecting wavelets and generators in the i-th coordinate direction
      for (unsigned int i = 0; i < DIM_d; i++) {
        get_translation_wavelets(*frame.bases()[p]->bases()[i],
	    	Index1D(lambda.j(),
		    lambda.e()[i],
		    lambda.k()[i],
		    frame.bases()[p]->bases()[i]),
		    j, true, minkkkk,maxkkkk);
        minkgen[i]=minkkkk;
        maxkgen[i] = maxkkkk;
        if (!(generators))
	  get_translation_wavelets(*frame.bases()[p]->bases()[i],
		      Index1D(lambda.j(),
			      lambda.e()[i],
			      lambda.k()[i],
			      frame.bases()[p]->bases()[i]),
		      j, false, minkkkk,maxkkkk);
        minkwavelet[i] = minkkkk;
        maxkwavelet[i] = maxkkkk;
       } // end for

      unsigned int result = 0;
      int deltaresult = 0;
      const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();

      MultiIndex<int,DIM_d> type;
      type[DIM_d-1] = 1;
      unsigned int tmp = 1;
      bool exit = 0;

      while(!exit){
      FixedArray1D<int,DIM_d> help1, help2;
      
      for(unsigned int i = 0; i<DIM_d; i++)
         help1[i]=0;

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches vor p liegen 
      for (unsigned int pp = 0; pp < p; pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }
            // fügt die Indizes aus diesen Patches ein falls sie sich überlappen
            for (unsigned int i = result; i < result + tmp; i++) {
              const Index* ind = &((*full_collection_levelwise)[k][i]);
              if ( intersect_supports_simple(frame, lambda, *ind) ){
	        intersecting.push_back(*ind);
              }
            }
            result += tmp;
	  }



      // berechnet wie viele indices mit einem zu keinem translationstyp es gibt, so dass sich die Wavelets nicht schneiden
      unsigned int result2 = 0;
      for (unsigned int i = 0; i < DIM_d; i++) {  // begin for1
        int tmp = 1;

        for (unsigned int l = i+1; l < DIM_d; l++) {
	  if (type[l] == 0)
	    tmp *= ((frame.bases()[p])->bases())[l]->Deltasize(j);
	  else
	    tmp *= ((frame.bases()[p])->bases())[l]->Nablasize(j);
        }

        help2[i] = tmp;

        if (type[i] == 0) {
	  if (minkgen[i] == ((frame.bases()[p])->bases())[i]->DeltaLmin())
	    continue;
        }
        else
	  if (minkwavelet[i] == ((frame.bases()[p])->bases())[i]->Nablamin())
	    continue;
      
        
        if (type[i] == 0) {
	tmp *= minkgen[i]-((frame.bases()[p])->bases())[i]->DeltaLmin();
        }
        else
	  tmp *= minkwavelet[i]-((frame.bases()[p])->bases())[i]->Nablamin();

        result2 += tmp;
      }  // end for1

      int tmp = 0;

      if (type[DIM_d-1] == 0) {
	tmp = maxkgen[DIM_d-1] - minkgen[DIM_d-1]+1;
      }
      else{
	tmp = maxkwavelet[DIM_d-1] - minkwavelet[DIM_d-1]+1; 
      }
     
      bool exit2 = 0;
      
      while(!exit2){

      // fügt die Indizes aus diesem Patch ein falls sie sich überlappen
      for (unsigned int i = result + result2; i < result + result2 + tmp; i++) {
        const Index* ind = &((*full_collection_levelwise)[k][i]);
	intersecting.push_back(*ind);
      }

      for (unsigned int i = DIM_d-2; i >= 0; i--) {
            if(type[i]==0){
	      if ( help1[i] < maxkgen[i]-minkgen[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                for (unsigned int j = i+1; j<=DIM_d-2;j++){
                    if(type[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    } 
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
            }
            else {
              if ( help1[i] < maxkwavelet[i] - minkwavelet[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                for (unsigned int j = i+1; j<=DIM_d-2;j++){
                    if(type[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    }
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
	    }
	  } //end for
      } //end while 2


      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches p liegen 
      tmp = 1;
      for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[p])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[p])->bases())[i]->Nablasize(j);
	    }
         
      result += tmp;

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches nach p liegen     
      for (int pp = p+1; pp < frame.n_p(); pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }

            // fügt die Indizes aus diesen Patches ein falls sie sich überlappen
            for (unsigned int i = result; i < result + tmp; i++) {
               const Index* ind = &((*full_collection_levelwise)[k][i]);
               if ( intersect_supports_simple(frame, lambda, *ind) ){
	          intersecting.push_back(*ind);
                }
             }
	    result += tmp;
	  }

     // berechnet den nächsten Typ 
     for (unsigned int i = DIM_d-1; i >= 0; i--) {
	    if ( type[i] == 1 ) {
	      type[i] = 0;
              exit = (i == 0);
	      if(exit)
               break;
	    }
	    else {
	      type[i]++;
	      break;
	    }
	  } //end for
       } // end while 1
      } // end if
    // } // end if
    
//*/
    else{  //if generator

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }
    std::list<typename Frame::Index> intersect_diff;

    const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();

    for (unsigned int i = 0; i < (*full_collection_levelwise)[k].size(); i++) {
      const Index* ind = &((*full_collection_levelwise)[k][i]);
      if ( intersect_supports_simple(frame, lambda, *ind) ){
	intersecting.push_back(*ind);
      }
    }
    } //end else
#else

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    std::list<typename Frame::Index> intersect_diff;

    const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();

    int k = -1;
    if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }

    for (unsigned int i = 0; i < (*full_collection_levelwise)[k].size(); i++) {
      const Index* ind = &((*full_collection_levelwise)[k][i]);
      if (frame.atlas()->get_adjacency_matrix().get_entry(lambda.p(), ind->p()) && 
	  intersect_supports_simple(frame, lambda, *ind) ){
	//intersect_diff.push_back(*ind);
	intersecting.push_back(*ind);
      }
    }
#endif
  }
  
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  void intersecting_wavelets (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			      const int p,
			      const std::set<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& Lambda,	      
			      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting)
  {
    //intersecting.erase(intersecting.begin(),intersecting.end());

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;
    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;
    const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda = &((frame.all_supports)[lambda.number()]);
    // loop over the set Lambda
    for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	 it1 != itend; ++it1) {
      const typename CubeBasis<IBASIS,DIM_d>::Support* supp_ind = &((frame.all_supports)[(*it1).number()]);
      if ( (*it1).p()==p && intersect_supports_simple(frame, lambda, *it1) ){
	intersecting.push_back(*it1);
      }
    }
  }




  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  void intersecting_wavelets_on_patch (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				       const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				       const int p,
				       const int j, const bool generators,
				       std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting)
  {
    //intersecting.erase(intersecting.begin(),intersecting.end());

#if 1

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }
    std::list<typename Frame::Index> intersect_diff;


    if (! generators) {
      if(p == lambda.p()){

        FixedArray1D<int,DIM_d>
            minkwavelet, maxkwavelet, minkgen, maxkgen;
        typedef typename IBASIS::Index Index1D;
        int minkkkk;
        int maxkkkk;

        // prepare all intersecting wavelets and generators in the i-th coordinate direction
        for (unsigned int i = 0; i < DIM_d; i++) {
          get_translation_wavelets(*frame.bases()[p]->bases()[i],
			    Index1D(lambda.j(),
				    lambda.e()[i],
				    lambda.k()[i],
				    frame.bases()[p]->bases()[i]),
			    j, true, minkkkk,maxkkkk);
          minkgen[i]=minkkkk;
          maxkgen[i] = maxkkkk;
          if (!(generators))
	     get_translation_wavelets(*frame.bases()[p]->bases()[i],
			      Index1D(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      frame.bases()[p]->bases()[i]),
			      j, false, minkkkk,maxkkkk);
          minkwavelet[i] = minkkkk;
          maxkwavelet[i] = maxkkkk;
         } // end for

        unsigned int result = 0;
        int deltaresult = 0;
        const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();

        MultiIndex<int,DIM_d> type;
        type[DIM_d-1] = 1;
        unsigned int tmp = 1;
        bool exit = 0;

      while(!exit){
      FixedArray1D<int,DIM_d> help1, help2;
      
      for(unsigned int i = 0; i<DIM_d; i++)
         help1[i]=0;

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches vor p liegen 
      for (int pp = 0; pp < p; pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }
	    result += tmp;
	  }

      // berechnet wie viele indices mit einem zu keinem translationstyp es gibt, so dass sich die Wavelets nicht schneiden
      unsigned int result2 = 0;
      for (unsigned int i = 0; i < DIM_d; i++) {  // begin for1
        int tmp = 1;

        for (unsigned int l = i+1; l < DIM_d; l++) {
	  if (type[l] == 0)
	    tmp *= ((frame.bases()[p])->bases())[l]->Deltasize(j);
	  else
	    tmp *= ((frame.bases()[p])->bases())[l]->Nablasize(j);
        }

        help2[i] = tmp;

        if (type[i] == 0) {
	  if (minkgen[i] == ((frame.bases()[p])->bases())[i]->DeltaLmin())
	    continue;
        }
        else
	  if (minkwavelet[i] == ((frame.bases()[p])->bases())[i]->Nablamin())
	    continue;
      
        
        if (type[i] == 0) {
	tmp *= minkgen[i]-((frame.bases()[p])->bases())[i]->DeltaLmin();
        }
        else
	  tmp *= minkwavelet[i]-((frame.bases()[p])->bases())[i]->Nablamin();

        result2 += tmp;
      }  // end for1

      int tmp = 0;

      if (type[DIM_d-1] == 0) {
	tmp = maxkgen[DIM_d-1] - minkgen[DIM_d-1]+1;
      }
      else{
	tmp = maxkwavelet[DIM_d-1] - minkwavelet[DIM_d-1]+1; 
      }
     
      bool exit2 = 0;
      
      while(!exit2){

      // fügt die Indizes aus diesem Patch ein falls sie sich überlappen
      for (unsigned int i = result + result2; i < result + result2 + tmp; i++) {
        const Index* ind = &((*full_collection_levelwise)[k][i]);
	intersecting.push_back(*ind);
      }

      for (unsigned int i = DIM_d-2; i >= 0; i--) {
            if(type[i]==0){
	      if ( help1[i] < maxkgen[i]-minkgen[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                for (unsigned int j = i+1; j<=DIM_d-2;j++){
                    if(type[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    } 
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
            }
            else {
              if ( help1[i] < maxkwavelet[i] - minkwavelet[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                for (unsigned int j = i+1; j<=DIM_d-2;j++){
                    if(type[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    }
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
	    }
	  } //end for
      } //end while 2


      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches p liegen 
      tmp = 1;
      for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[p])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[p])->bases())[i]->Nablasize(j);
	    }
         

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches nach p liegen     
      result += tmp;
      for (int pp = p+1; pp < frame.n_p(); pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }
	    result += tmp;
	  }
     

     // berechnet den nächsten Typ 
     for (unsigned int i = DIM_d-1; i >= 0; i--) {
	    if ( type[i] == 1 ) {
	      type[i] = 0;
              exit = (i == 0);
	      if(exit)
               break;
	    }
	    else {
	      type[i]++;
	      break;
	    }
	  } //end for
       } // end while 1
    //  } // end if
    // } // end if
    
//*/
    }  // end if von if(lamda.p() == p) //////////////////////////////////////////////////////////////////////////////////////////

     else{  
        unsigned int result = 0;
        int deltaresult = 0;
        const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();

        MultiIndex<int,DIM_d> type;
        type[DIM_d-1] = 1;
        unsigned int tmp = 1;
        bool exit = 0;

      while(!exit){

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches vor p liegen 
      for (int pp = 0; pp < p; pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }
	    result += tmp;
	  }
      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches p liegen 
      tmp = 1;
      for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[p])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[p])->bases())[i]->Nablasize(j);
	    }
      
      // fügt die Indizes aus diesem Patch ein falls sie sich überlappen
      for (unsigned int i = result; i < result + tmp; i++) {
        const Index* ind = &((*full_collection_levelwise)[k][i]);
        if ( intersect_supports_simple(frame, lambda, *ind) ){
	  intersecting.push_back(*ind);
        }
      }

      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches nach p liegen     
      result += tmp;
      for (int pp = p+1; pp < frame.n_p(); pp++) {
	    tmp = 1;
	    for (unsigned int i = 0; i < DIM_d; i++) {
	      if (type[i] == 0)
	        tmp *= ((frame.bases()[pp])->bases())[i]->Deltasize(j);
	      else
	        tmp *= ((frame.bases()[pp])->bases())[i]->Nablasize(j);
	    }
	    result += tmp;
	  }

     // berechnet den nächsten Typ 
     for (unsigned int i = DIM_d-1; i >= 0; i--) {
	    if ( type[i] == 1 ) {
	      type[i] = 0;
              exit = (i == 0);
	      if(exit)
               break;
	    }
	    else {
	      type[i]++;
	      break;
	    }
	  } //end for
       } // end while
      } // end else  // lambda.p = p
     } // end if

    else{  //if generator

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }
    std::list<typename Frame::Index> intersect_diff;

    const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();
    for (unsigned int i = 0; i < (*full_collection_levelwise)[k].size(); i++) {
      const Index* ind = &((*full_collection_levelwise)[k][i]);
      if ( (ind->p() == p) && 
	   intersect_supports_simple(frame, lambda, *ind) ){
	intersecting.push_back(*ind);
      }
    }
    } //end else
#else

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-frame.j0()+1;
    }
    std::list<typename Frame::Index> intersect_diff;

    const Array1D<Array1D<Index> >* full_collection_levelwise = frame.get_full_collection_levelwise();
    for (unsigned int i = 0; i < (*full_collection_levelwise)[k].size(); i++) {
      const Index* ind = &((*full_collection_levelwise)[k][i]);
      if ( (ind->p() == p) && 
	   intersect_supports_simple(frame, lambda, *ind) ){
	intersecting.push_back(*ind);
      }
    }
#endif
  }
  
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool intersect_singular_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu)
  {
   
    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
    const typename CUBEBASIS::Support* supp_lambda = &((frame.all_supports)[lambda.number()]);
    const typename CUBEBASIS::Support* supp_mu = &((frame.all_supports)[mu.number()]);
    
    if (! intersect_supports<IBASIS,DIM_d,DIM_m>(frame, lambda, mu, supp_lambda, supp_mu) )
      return false;

    typedef typename IBASIS::Index Index_1D;
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;

    if (lambda.p() == mu.p()) {
      for (unsigned int i = 0; i < DIM_d; i++) {
	if ( intersect_singular_support(*frame.bases()[lambda.p()]->bases()[i],
					Index_1D(lambda.j(),lambda.e()[i],lambda.k()[i],
						 frame.bases()[lambda.p()]->bases()[i]),
					Index_1D(mu.j(),mu.e()[i],mu.k()[i],
						 frame.bases()[mu.p()]->bases()[i]))
	     )
	  return true;
      }
      return false;
    }
    else {
      // The following appraoach will only work for those cases where the support of a wavelet 
      // is some convex quadrangle.
      // So far, all kinds of parametrizations we use satisfy this assumption.
      // Thus, we do not check the typeid here.

      const typename CUBEBASIS::Support* tmp_supp;

      Index lambda_c = lambda;
      Index mu_c = mu;
      
      // Assume: mu is the index with the higher resolution!

      // swap indices and supports if necessary
      if (supp_mu->j < supp_lambda->j) {
	lambda_c = mu;
	mu_c = lambda;
	
	tmp_supp = supp_lambda;
	supp_lambda = supp_mu;
	supp_mu = tmp_supp;
      }
      
      const double dx1 = 1.0 / (1 << supp_lambda->j);
      const double dx2 = 1.0 / (1 << supp_mu->j);
      
      if (DIM_d == 1 && DIM_m == 1) {
	
	const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda_c.p()];
	const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu_c.p()];
	
	FixedArray1D<Point<DIM_m>,2> poly_mu;
	FixedArray1D<Point<DIM_m>,2> patch_lambda;   

	// map the knots of the unit cube to patch
	// get the corners of the patch of index lambda
	chart_la->map_point(Point<DIM_m>(0.), patch_lambda[0]);
	chart_la->map_point(Point<DIM_m>(1.), patch_lambda[1]);
	
	chart_mu->map_point(Point<DIM_m>(supp_mu->a[0]*dx2), poly_mu[0]);
	chart_mu->map_point(Point<DIM_m>(supp_mu->b[0]*dx2), poly_mu[1]);
	
	bool fully_contained = 1;
	if (! ((patch_lambda[0][0] <= poly_mu[0][0]) && (poly_mu[1][0] <= patch_lambda[1][0])))
	  fully_contained = false;
	
	if (! fully_contained)
	  // the supports intersect, but mu does not correspond to a wavelet
	  // that is fully contained in the patch of lambda
	  return true;
		
	FixedArray1D<Point<DIM_m>,2> poly_mu_pulled_back;
	

	chart_la->map_point_inv(poly_mu[0], poly_mu_pulled_back[0]);
	chart_la->map_point_inv(poly_mu[1], poly_mu_pulled_back[1]);

	FixedArray1D<double, 2> k1;
	FixedArray1D<double, 2> k2;

	k1[0] = floor(poly_mu_pulled_back[0][0] / dx1);
	k2[0] = floor(poly_mu_pulled_back[0][0] / dx1);

	k1[1] = floor(poly_mu_pulled_back[1][0] / dx1);
	k2[1] = floor(poly_mu_pulled_back[1][0] / dx1);


	if ((k1[1] != k1[0]) || (k2[1] != k2[0]))
	  return true;
	
	return false;
      }

      if (DIM_d == 2 && DIM_m == 2) {
	const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda_c.p()];
	const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu_c.p()];

	FixedArray1D<Point<DIM_m>,4> poly_mu;
	FixedArray1D<Point<DIM_m>,4> patch_lambda;   

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01     
	// get the corners of the patch of index lambda
	chart_la->map_point(Point<DIM_d>(0.,0.), patch_lambda[0]);
	chart_la->map_point(Point<DIM_d>(1.,0.), patch_lambda[1]);
	chart_la->map_point(Point<DIM_d>(1.,1.), patch_lambda[2]);
	chart_la->map_point(Point<DIM_d>(0.,1.), patch_lambda[3]);
	
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), poly_mu[0]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->a[1]*dx2), poly_mu[1]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), poly_mu[2]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->b[1]*dx2), poly_mu[3]);  
	
	bool fully_contained = 1;
	
	// first check whether all of the four corners of the support of
	// mu are contained in the patch of lambda
	for (int i = 0; i < 4; i++) {
	  
	  //make sure to walk through the vertices counter clockwise!!!
	  unsigned short int res =  FrameTL::pos_wrt_line(poly_mu[i], patch_lambda[0], patch_lambda[1]);
	  if (res == 0) {
	    fully_contained = false;
	 break;
	  }
	  res =  FrameTL::pos_wrt_line(poly_mu[i], patch_lambda[1], patch_lambda[2]);
	  if (res == 0) {
	    fully_contained = false;
	    break;
	  }
	  res = FrameTL::pos_wrt_line(poly_mu[i], patch_lambda[2], patch_lambda[3]);
	  if (res == 0) {
	    fully_contained = false;
	    break;
	  }
	  res = FrameTL::pos_wrt_line(poly_mu[i], patch_lambda[3], patch_lambda[0]);
	  if (res == 0) {
	    fully_contained = false;
	    break;
	  }
	}
	
	if (! fully_contained)
	  // the supports intersect, but mu does not correspond to a wavelet
	  // that is fully contained in the patch of lambda
	  return true;
		
	FixedArray1D<Point<DIM_m>,4> poly_mu_pulled_back;
	
	chart_la->map_point_inv(poly_mu[0], poly_mu_pulled_back[0]);
	chart_la->map_point_inv(poly_mu[1], poly_mu_pulled_back[1]);
	chart_la->map_point_inv(poly_mu[2], poly_mu_pulled_back[2]);
	chart_la->map_point_inv(poly_mu[3], poly_mu_pulled_back[3]);
	
	FixedArray1D<double, 4> k1;
	FixedArray1D<double, 4> k2;
	for (int i = 0; i < 4; i++) {
	  // determine that square [k1*2^-j1,(k1+1)*2^-j]x[k2*2^-j2,(k2+1)*2^-j2]
	  // where psi_lambda is smooth and that contains poly_mu_pulled_back[i]
	  // if these squares are NOT all the same, the singular supports intersect
	  
	  k1[i] = floor(poly_mu_pulled_back[i][0] / dx1);
	  k2[i] = floor(poly_mu_pulled_back[i][1] / dx1);
	  
	  if (i > 0) {
	    //compare with the preceeding point
	    if ((k1[i] != k1[i-1]) || (k2[i] != k2[i-1]))
	      return true;
	  }
	}
	return false;
      }
    }// end two dimensional case

    return true;
  }
 
  template <unsigned int DIM>
  inline
  int edgesIntersect (const Point<DIM>& A, const Point<DIM>& B,
		      const Point<DIM>& C, const Point<DIM>& D) {
    
    // so far only the 2D case is implemented
    assert ( DIM == 2 );
    // The points A and B define the line f(s) = A + s*(B-A).
    // The points C and D define the line g(t) = C + t*(D-C).
    // ==>[ f(s) = g(t) <==> C-A = t*(C-D) + s*(B-A) ]

    double s = 0;
    double t = 0;
        
    double tmp = 0.;
    double u1 = B[0] - A[0];
    double u2 = B[1] - A[1];
 
    double v1 = 0.;
    double v2 = 0.;
	
    double w1 = 0.;
    double w2 = 0.;

    double x1 = 0.;
    double x2 = 0.;

    double y1 = 0.;

    if (u1 == 0.) { // ==> u2 != 0.
      // reverse role of lines
      tmp = u1;
      u1 = u2;
      u2 = tmp;
      
      v1 = C[1] - D[1];
      v2 = C[0] - D[0];
      
      w1 = C[1] - A[1];
      w2 = C[0] - A[0];
      
      x1 = C[1] - B[1];
      x2 = C[0] - B[0];
      
      y1 = D[1] - A[1];

    }
    else {
      v1 = C[0] - D[0];
      v2 = C[1] - D[1];
      
      w1 = C[0] - A[0];
      w2 = C[1] - A[1];
      

      x1 = C[0] - B[0];
      x2 = C[1] - B[1];
      
      y1 = D[0] - A[0];
    }

    // we have to solve the 2 by 2 linear system
    //   s    t
    // (u1   v1 | w1)
    // (u2   v2 | w2)

    // u1 != 0. is now guaranteed

    double v2_bak = v2;
    double w2_bak = w2;
//    if (!(u2 == 0.)) {
    if (!eq(u2,0.0)) {
      // gauss elimination
      v2 += -(u2 / u1)*v1;
      w2 += -(u2 / u1)*w1;
      
    }
    
    //      if (v2 == 0. && !(w2 == 0.)) {
    if (eq(v2,0.0) && !eq(w2,0.0)) {
      // no solutions at all
      return 0;
    }
    //    if (!(v2 == 0.)) {
    if (!eq(v2,0.0)) {
      // unique solution
      t = w2 / v2;
      s = (w1 - t*v1) / u1;
      
      // check whether the intersecting point is a knot
      // for any of the two edges
      if (eq(s,0.) || eq(s,1.) || eq(t,0.) || eq(t,1.))
	return 2;
      else if (lt(0.0,s) && lt(s,1.0) && lt(0.0,t) && lt(t,1.0)) {
	return 3;
      }
      else return 0;
    }
    if (eq(v2,0.0) && eq(w2,0.0)) {
      // infinitely many solutions --> most complicated case
      // (the lines given by the edges coincide)
      // test now if the edges also intersect
      
      int n_ex_case = 0;
      // first: check if there exists a knots of 'f' between the knots of 'g'
      // f(0) = A, g(t) = C + t(D-C) = A <==> A-C = t(D-C) <==> [-w1 == -t*v1]
      if (! eq(v1,0.0) ) {
	t = w1 / v1;
	if (lt(0.0,t) && lt(t,1.0)) {
	  return 1;
	}
	if (eq(t,0.0) || eq(t,1.0))
	  n_ex_case++;
      }
      else  { // v1 == 0 (==> w1 == 0)
	// then use the other component
	// observe that v2 == 0 is impossible, otherwise C==D, which is never true
	t = w2_bak / v2_bak;
	if (lt(0.0,t) && lt(t, 1.0))
	  return 1;
	if (eq(t,0.0) || eq(t,1.0))
	  n_ex_case++;
      }
      // f(1) = B, g(t) = C + t(D-C) = B <==> B-C = t(D-C) <==> [-x1 == -t*v1]
      if (! eq(v1,0.0) ) {
	t = x1 / v1;
	if (lt(0.0,t) && lt(t,1.0))
	  return 1;
	if (eq(t,0.0) || eq(t,1.0))
	  n_ex_case++;
      }
      else { // v1 == 0
	t = x2 / v2_bak;
	if (lt(0.0,t) && lt(t,1.0))
	  return 1;
	if (eq(t,0.0) || eq(t,1.0))
	  n_ex_case++;
      }
      if (n_ex_case == 2)
	// the two edges coincide
	return 1;

      n_ex_case = 0;      

      // second: check if there exists a knots of 'g' between the knots of 'f'
      // g(0) = C, f(s) = A + s*(B-A) = C  <==> C-A = s*(B-A)[w1 == t*u1]
      // keep in mind: u1 != 0. is guaranteed
      s = w1 / u1;
      if (lt(0.0,s) && lt(s,1.0))
	return 1;
      if (eq(s,0.0) || eq(s,1.0))
	n_ex_case++;
      
      // g(1) = D, f(s) = A + s*(B-A) = D  <==> D-A = s*(B-A)[y1 == t*u1]
      // keep in mind: u1 != 0. is guaranteed
      s = y1 / u1;
      if (lt(0.0,s) && lt(s,1.0))
	return 1;
      if (eq(s,0.0) || eq(s,1.0))
	n_ex_case++;

      if (n_ex_case == 1)
	return 2;
      
    }   

    // dummy
    return -1;
													 
  }

}

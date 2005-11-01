// implementation for frame_support.h


#include <utils/fixed_array1d.h>
#include <typeinfo>
#include <cube/cube_support.h>

using namespace WaveletTL;

using std::type_info;

namespace FrameTL
{

  template <unsigned int DIM>
  inline
  unsigned short int pos_wrt_line (const Point<DIM>& p,
				   const Point<DIM>& p1, const Point<DIM>&  p2)
  {
   
    assert ( DIM == 2 );

    double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    
    if( d > 0.0 )
      return 1;
    else if (d < 0.0)
      return 0;
    else
      return 2;
  }


  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const Point<DIM_m>& p)
  {
    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(LinearBezierMapping))
      {
	typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
	
	typename CUBEBASIS::Support supp_lambda;
	
	typename CubeBasis<IBASIS,DIM_d>::Index lambda_c(lambda.j(),
							 lambda.e(),
							 lambda.k(),frame.bases()[lambda.p()]);

	WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[lambda.p()], lambda_c, supp_lambda);
	  
	const double dx = 1.0 / (1 << supp_lambda.j);
	
	FixedArray1D<Point<2>,4 > poly;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx,supp_lambda.a[1]*dx), poly[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx,supp_lambda.a[1]*dx), poly[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx,supp_lambda.b[1]*dx), poly[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx,supp_lambda.b[1]*dx), poly[3]);

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
	
	return true;;

      }
    return true;
  }

  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
		  const Point<DIM_m>& p)
  {
    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(LinearBezierMapping))
      {
	const double dx = 1.0 / (1 << supp_lambda.j);

	FixedArray1D<Point<2>,4 > poly;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx,supp_lambda.a[1]*dx), poly[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx,supp_lambda.a[1]*dx), poly[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx,supp_lambda.b[1]*dx), poly[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx,supp_lambda.b[1]*dx), poly[3]);

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
	
	return true;;

      }
    return true;
  }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool
  intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
		     typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
		     typename CubeBasis<IBASIS,DIM_d>::Support& supp_mu)
  {
    // both charts are LinearBezierMappings
    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(LinearBezierMapping) &&
	 typeid(*frame.atlas()->charts()[mu.p()])     == 
	 typeid(LinearBezierMapping))
      {


	typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
	
	typename CubeBasis<IBASIS,DIM_d>::Index lambda_c(lambda.j(),
							 lambda.e(),
							 lambda.k(),frame.bases()[lambda.p()]);
	
	typename CubeBasis<IBASIS,DIM_d>::Index mu_c(mu.j(),
						     mu.e(),
						     mu.k(),frame.bases()[mu.p()]);

	WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[lambda.p()], lambda_c, supp_lambda);
	WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[mu.p()], mu_c, supp_mu);
  
	const double dx1 = 1.0 / (1 << supp_lambda.j);
	const double dx2 = 1.0 / (1 << supp_mu.j);

	FixedArray1D<Point<2>,4 > poly1;
	FixedArray1D<Point<2>,4 > poly2;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx1,supp_lambda.a[1]*dx1), poly1[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx1,supp_lambda.a[1]*dx1), poly1[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.b[0]*dx1,supp_lambda.b[1]*dx1), poly1[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<2>(supp_lambda.a[0]*dx1,supp_lambda.b[1]*dx1), poly1[3]);

	frame.atlas()->charts()[mu.p()]->map_point(Point<2>(supp_mu.a[0]*dx2,supp_mu.a[1]*dx2), poly2[0]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<2>(supp_mu.b[0]*dx2,supp_mu.a[1]*dx2), poly2[1]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<2>(supp_mu.b[0]*dx2,supp_mu.b[1]*dx2), poly2[2]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<2>(supp_mu.a[0]*dx2,supp_mu.b[1]*dx2), poly2[3]);

//  	for (int i = 0; i < 4; i++)
//  	  cout << "i=" << i << " " << poly1[i] << endl;

//  	for (int i = 0; i < 4; i++)
//  	  cout << "i=" << i << " " << poly2[i] << endl;

// 	double a1 = poly1[0][0], b1 = poly1[2][0], a2 = poly1[0][1], b2 = poly1[2][1];

// 	double v1 = poly2[0][0], w1 = poly2[2][0], v2 = poly2[0][1], w2 = poly2[2][1];

// 	cout << "res must be " << ((a1 < w1 && b1 > v1) && (a2 < w2 && b2 > v2)) << endl;

 	bool result;
 	Point<DIM_m> p11;
 	Point<DIM_m> p12;
	
 	unsigned short int tmp = 0;
 	int countspez = 0;
	
 	//check, if two arbritrary edges of the quadrangle intersect
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
		
	//if control reaches this point, then no 'ordinary intersection' exists
	//ATTENTION!!!! But sufficient so far.
	//this is only correct in case of convex polygons!!!!!!!
	if (countspez >= 2) {
	  return true;
	}
	
	
	//remaining case: qudrangles contain each other
	countspez = 0;
	tmp = 0;
	
	//loop over all points knots of the second polygon
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
	    if (tmp == 0 || tmp ==2) {
	      result = false;
	      break;
	    }
	  }//end inner for
		
	  if (result)
	    return result;

	}//end outer for
	
	//the same with swapped order of the 
	//loop over all points knots of the second polygon
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
    else return false;
  }
 


  template <unsigned int DIM>
  inline int edgesIntersect  (const Point<DIM>& A, const Point<DIM>& B,
			      const Point<DIM>& C, const Point<DIM>& D) {
    
    // so far only the 2D case is implemented
    assert ( DIM == 2 );

    double s = 0;
    double t = 0;
        
    double u1 = B[0] - A[0];
    double u2 = B[1] - A[1];
	
    double v1 = C[0] - D[0];
    double v2 = C[1] - D[1];
	
    double w1 = C[0] - A[0];
    double w2 = C[1] - A[1];

    double x1 = D[0]-A[0];
    double x2 = D[1]-A[1];
	
    double t1, t2;
	
    //u1 = u2 = 0 cannot happen!!
    if (u1 == 0.0 && ! (v1 == 0.0)) {
      s = w1 / v1;
      t = (w2 - v2*s) / u2;
      if (   (((s == 0.0) || (s == 1.0)) && (t < 0.0 || t > 1.0))
	     || (((t == 0.0) || (t == 1.0)) && (s < 0.0 || s > 1.0))) {

	return 0;// to return 0 seems to be incorrect
	//return 2;
      }
      else 	if (   (((s == 0.0) || (s == 1.0)) && ((0.0 <= t) || (t <= 1.0)))
		       || ((( t== 0.0) || (t == 1.0)) && ((0.0 <= s) || (s <= 1.0))))

	return 2;
      else if (s < 0.0 || s > 1.0 || t < 0.0 || t > 1.0) {
	return 0;
      }
      else
	return 3;

    }//u1 = 0 ==> u2 != 0
    else if ((u1 == 0.0) && (v1 == 0.0)) {
      //if (eq(w1,0.0)) {
      if ((w1 == 0.0)) {
	//determine if one knot of second edge lies on first edge
	double d = w2 / u2;
	double e = x2 / u2;
	if ( (0 < d) && (d < 1) ) {
	  return 1;
	}
	else if ( (0 < e) && (e < 1)) {
	  return 1;
	}
	else if ( ((d == 0.0) &&  (e == 1.0)) || ((d == 1.0) && (e == 0.0)))
	  return 1;
	else if ((d < 0 && (e == 1.0)) || ((e < 0) && (d == 1.0))) {
	  return 1;
	}
	else if ((d > 1 && (e == 0.0)) || ((e > 1) && (d == 0.0))) {
	  return 1;
	}
	else {
	  return 0;
	}
      }
      else {
	return 0;
      }
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
    t1 = ((-u2*v1) / u1) + v2;
    t2 = w2 - (u2*w1/u1);

    if (! (t1 == 0.0)) {
      s = t2 / t1;
      t = (w1 - v1*s) / u1;
      //intersecting point is a knot for each edge:
      if (   (((s == 0.0) || (s == 1.0)) && (t < 0.0 || t > 1.0))
	     || (((t == 0.0) || (t == 1.0)) && (s < 0.0 || s > 1.0))) {

	return 0;
      }
      else 	if (   (((s == 0.0) || (s == 1.0)) && ((0.0 <= t) || (t <= 1.0)))
		       || (((t == 0.0) || (t == 1.0)) && ((0.0 <= s) || (s <= 1.0))))

	return 2;
      else if (s < 0.0 || s > 1.0 || t < 0.0 || t > 1.0) {
	return 0;
      }
      else
	return 3;
    }

    //straight lines are the same
    //--> find out, wether the edges overlap or not
    //	if (eq(t2,0.0) && eq(t1,0.0)) {
    if ((t2 == 0.0) && (t1 == 0.0)) {
      //determine if one knot of second edge lies on first edge
      double d = w1 / u1;
      double e = x1 / u1;
      if ( (0 < d) && (d < 1) ) {
	return 1;
      }
      else if ( (0 < e) && (e < 1)) {
	return 1;
      }
      else if ( ((d == 0.0) &&  (e == 1.0)) || ((d == 1.0) && (e == 0.0)))
	return 1;
      else if ((d < 0 && (e == 1.0)) || ((e < 0) && (d == 1.0))) {
	return 1;
      }
      else if ((d > 1 && (e == 0.0)) || ((e > 1) && (d == 0.0))) {
	return 1;
      }
      else {
	return 0;
      }
    }
	
    //parallel straight lines
    if ((t1 == 0.0) && ! (t2 == 0.0)) {
      return 0;
    }

    //dummy
    return -1;
														 
  }

}

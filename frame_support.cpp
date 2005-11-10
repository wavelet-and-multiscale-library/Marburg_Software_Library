// implementation for frame_support.h


#include <utils/fixed_array1d.h>
#include <typeinfo>
#include <cube/cube_support.h>
#include <map>

using namespace WaveletTL;

using std::type_info;
using std::map;

namespace FrameTL
{

  template <unsigned int DIM>
  inline
  unsigned short int pos_wrt_line (const Point<DIM>& p,
				   const Point<DIM>& p1, const Point<DIM>&  p2)
  {
   
    assert ( DIM == 2 );

    double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    
    //cout << "d = " << d << endl;    
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

    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
	
    typename CUBEBASIS::Support supp_lambda;
    
    typename CUBEBASIS::Index lambda_c(lambda.j(),
				       lambda.e(),
				       lambda.k(),frame.bases()[lambda.p()]);
    
    WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[lambda.p()], lambda_c, supp_lambda);
    
    const double dx = 1.0 / (1 << supp_lambda.j);
    
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>))
	 )
      {
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx), b);
	
	return (a[0] <= p[0]) && (p[0] <= b[0]);
	
      }

    if ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	 typeid(LinearBezierMapping))
      {
	
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

  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
		  const Point<DIM_m>& p)
  {

    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>))
	 )
      {
	const double dx = 1.0 / (1 << supp_lambda.j);
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx), b);

	return (a[0] <= p[0]) && (p[0] <= b[0]);
	
      }

    // just to get a better performance in the 
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
    
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(LinearBezierMapping))
	 ||
	 (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>))
	 )
      {
	assert ( DIM_d == 2 && DIM_m == 2 );
	
	const double dx = 1.0 / (1 << supp_lambda.j);

	FixedArray1D<Point<DIM_m>,4 > poly;

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
    const double dx1 = 1.0 / (1 << supp_lambda.j);
    const double dx2 = 1.0 / (1 << supp_mu.j);


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

	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1), b);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2), c);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2), d);

	return (a[0] < d[0]) && (b[0] > c[0]);

      }

    // both charts are LinearBezierMappings
    if ( ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(LinearBezierMapping) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(LinearBezierMapping) )
	 ||
	 ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(AffineLinearMapping<2>) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(AffineLinearMapping<2>)
	   )
	 )
      {

	assert ( DIM_d == 2 && DIM_m == 2);

	FixedArray1D<Point<DIM_m>,4 > poly1;
	FixedArray1D<Point<DIM_m>,4 > poly2;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1,supp_lambda.a[1]*dx1), poly1[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1,supp_lambda.a[1]*dx1), poly1[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1,supp_lambda.b[1]*dx1), poly1[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1,supp_lambda.b[1]*dx1), poly1[3]);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2,supp_mu.a[1]*dx2), poly2[0]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2,supp_mu.a[1]*dx2), poly2[1]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2,supp_mu.b[1]*dx2), poly2[2]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2,supp_mu.b[1]*dx2), poly2[3]);

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
 
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda)
  {
    
    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
    typedef typename CUBEBASIS::Index CubeIndex;
    typename CUBEBASIS::Support supp_mu;

    WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[mu.p()], 
				     CubeIndex(mu.j(),
					       mu.e(),
					       mu.k(),
					       frame.bases()[mu.p()]),
				     supp_mu);
    
    const double dx1 = 1.0 / (1 << supp_lambda.j);
    const double dx2 = 1.0 / (1 << supp_mu.j);


    //cout << supp_mu.a[0] << " " << supp_mu.b[0] << " " << supp_mu.j << endl;

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

	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1), b);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2), c);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2), d);

	return (a[0] < d[0]) && (b[0] > c[0]);

      }

    // both charts are LinearBezierMappings
    if ( ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(LinearBezierMapping) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(LinearBezierMapping) )
	 ||
	 ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(AffineLinearMapping<2>) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(AffineLinearMapping<2>)
	   )
	 )
      {

	assert ( DIM_d == 2 && DIM_m == 2);

	FixedArray1D<Point<DIM_m>,4 > poly1;
	FixedArray1D<Point<DIM_m>,4 > poly2;

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1,supp_lambda.a[1]*dx1), poly1[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1,supp_lambda.a[1]*dx1), poly1[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.b[0]*dx1,supp_lambda.b[1]*dx1), poly1[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda.a[0]*dx1,supp_lambda.b[1]*dx1), poly1[3]);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2,supp_mu.a[1]*dx2), poly2[0]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2,supp_mu.a[1]*dx2), poly2[1]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.b[0]*dx2,supp_mu.b[1]*dx2), poly2[2]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu.a[0]*dx2,supp_mu.b[1]*dx2), poly2[3]);

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

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
			  typename CubeBasis<IBASIS,DIM_d>::Support& supp_mu,
			  FixedArray1D<Array1D<double>,DIM_d >& supp_intersect)
  {    

    // template parameter DIM_d has to be specified when
    // the typeid shall really contain the information that an
    // affine mapping is used
//     assert ( typeid(frame.atlas()->charts()[lambda.p()])
// 	     ==  typeid(AffineLinearMapping<DIM_d>)
// 	     &&
// 	     typeid(frame.atlas()->charts()[mu.p()])
// 	     ==  typeid(AffineLinearMapping<DIM_d>));

    // WE ALSO ASSUME THAT THE MATRIX 'A' OF THE CHART IS A DIAGONAL ONE
    // WITH POSITIVE ENTRIES ONLY

//     for (unsigned int i = 0; i < DIM_d; i++) 
//       cout << supp_lambda.a[i] << " " << supp_lambda.b[i] << " " << supp_lambda.j << endl;
//     for (unsigned int i = 0; i < DIM_d; i++)
//       cout << supp_mu.a[i] << " " << supp_mu.b[i] << " " << supp_mu.j << endl;

    const unsigned int n_points_la = (1 << supp_lambda.j)+1;
    const unsigned int n_points_mu = (1 << supp_mu.j)+1;

    const double dx1 = 1.0 / (n_points_la-1);
    const double dx2 = 1.0 / (n_points_mu-1);

    Point<DIM_d> x;
    Point<DIM_d> x_patch;
    Point<DIM_d> y0;
    Point<DIM_d> y1;
    
    for (unsigned int i = 0; i < DIM_d; i++)
      x[i] = supp_mu.a[i] * dx2;
      
    frame.atlas()->charts()[mu.p()]->map_point(x,x_patch);
    frame.atlas()->charts()[lambda.p()]->map_point_inv(x_patch,y0);
  
    for (unsigned int i = 0; i < DIM_d; i++) {
      x[i] = supp_mu.b[i] * dx2;
    }
    
    frame.atlas()->charts()[mu.p()]->map_point(x,x_patch);
    frame.atlas()->charts()[lambda.p()]->map_point_inv(x_patch,y1);
    
 //    cout << y0 << " " << y1 << endl;


    FixedArray1D<FixedArray1D<double,2>,DIM_d > hyperCube_intersect;
    for (unsigned int i = 0; i < DIM_d; i++) {
      assert ( y0[i] <= y1[i] );
    
      hyperCube_intersect[i][0] = std::max(supp_lambda.a[i]*dx1,y0[i]);
      hyperCube_intersect[i][1] = std::min(supp_lambda.b[i]*dx1,y1[i]);
      
 
      if ( hyperCube_intersect[i][0] >= hyperCube_intersect[i][1]  )
	return false;      
    }

//     for (unsigned int i = 0; i < DIM_d; i++) {
//       cout << "i = " << i << " intersect_cube[i] = " << hyperCube_intersect[i] << endl;
//     }
    
    // the bool is alway just a dummy
    // we just need some kind of a sorted list
    FixedArray1D<map<double,bool>, DIM_d> irregular_grid;

    // collect gridpoints for all directions,
    // we have to do this for both patches
    
    for (unsigned int i = 0; i < DIM_d; i++) {
      for (int k = supp_lambda.a[i]; k <= supp_lambda.b[i]; k++) {
	double d = k*dx1;
	if ( hyperCube_intersect[i][0] <= d && d <= hyperCube_intersect[i][1])
	  irregular_grid[i][d] = 1;
      }
    }
    
    // lower left corner of support cube
    for (unsigned int j = 0; j < DIM_d; j++) {
      x[j] = supp_mu.a[j]*dx2;
    }
    
    // second patch
    for (unsigned int i = 0; i < DIM_d; i++) {
      for (int k = supp_mu.a[i]; k <= supp_mu.b[i]; k++) {
	x[i] = k*dx2;
	frame.atlas()->charts()[mu.p()]->map_point(x,x_patch);
	frame.atlas()->charts()[lambda.p()]->map_point_inv(x_patch,y0);
	if ( hyperCube_intersect[i][0] <= y0[i] &&
	     y0[i] <= hyperCube_intersect[i][1])
	  irregular_grid[i][y0[i]] = 1;
      }
      x[i] = supp_mu.a[i]*dx2;
    }

    for (unsigned int i = 0; i < DIM_d; i++) {
      supp_intersect[i].resize(irregular_grid[i].size());
      unsigned int j = 0;
      for (map<double,bool>::const_iterator it = irregular_grid[i].begin();
	   it != irregular_grid[i].end(); it++) {
	supp_intersect[i][j] = it->first;
	j++;
      }
    }
    return true;
  }


  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void intersecting_wavelets (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			      const int j, const bool generators,
			      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting)
  {
    intersecting.erase(intersecting.begin(),intersecting.end());

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;
    
    typename CubeBasis<IBASIS,DIM_d>::Support supp_lambda;
    
    WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[lambda.p()], 
				     CubeIndex(lambda.j(),
					       lambda.e(),
					       lambda.k(),
					       frame.bases()[lambda.p()]),
				     supp_lambda);
    
#if 0
    // test
    std::list<Index> intersect_du;
    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;
    if ( generators ) {
      for (Index ind = frame.first_generator(j);
	 ind <= frame.last_generator(j); ++ind)
	{
	  //intersect_du.push_back(ind);
	}
      
    }
    else {
      for (Index ind = frame.first_wavelet(j);
	   ind <= frame.last_wavelet(j); ++ind)
	{
	  //intersect_du.push_back(ind);
	}
    }
    //end test
#endif
    //cout << supp_lambda.a[0] << " " << supp_lambda.b[0] << " " << supp_lambda.j << endl;

    std::list<CubeIndex> intersect_same_cube;
    WaveletTL::intersecting_wavelets<IBASIS,DIM_d>(*(frame.bases()[lambda.p()]),
						   CubeIndex(lambda.j(),
							     lambda.e(),
							     lambda.k(),
							     frame.bases()[lambda.p()]),
						   j, generators,
						   intersect_same_cube);

    // ################ brute force approach ##################
    std::list<typename Frame::Index> intersect_same;

    // create list of FrameIndices
    for (typename std::list<CubeIndex>::const_iterator  it = intersect_same_cube.begin();
	 it != intersect_same_cube.end(); ++it) {
      intersecting.push_back( FrameIndex<IBASIS,DIM_d,DIM_m>(&frame,*it,lambda.p()) );
    }

    std::list<typename Frame::Index> intersect_diff;

    if ( generators ) {
      for (Index ind = FrameTL::first_generator<IBASIS,DIM_d,DIM_m,Frame>(&frame, j);
	 ind <= FrameTL::last_generator<IBASIS,DIM_d,DIM_m,Frame>(&frame, j); ++ind)
	{
// 	  if ( (lambda.p() != ind.p()) &&
// 	       frame.atlas()->get_adjacency_matrix().get_entry(lambda.p(),ind.p()) && 
// 	       FrameTL::intersect_supports<IBASIS,DIM_d,DIM_m>(frame,lambda,ind,supp_lambda) )
// 	    intersect_diff.push_back(ind);
	  intersecting.push_back(ind);
	}
      
    }
    else {
      for (Index ind = FrameTL::first_wavelet<IBASIS,DIM_d,DIM_m,Frame>(&frame, j);
	   ind <= FrameTL::last_wavelet<IBASIS,DIM_d,DIM_m,Frame>(&frame, j); ++ind)
	{
// 	  if ( (lambda.p() != ind.p()) &&
// 	       frame.atlas()->get_adjacency_matrix().get_entry(lambda.p(),ind.p()) && 
// 	       FrameTL::intersect_supports<IBASIS,DIM_d,DIM_m>(frame,lambda,ind,supp_lambda) )
// 	    intersect_diff.push_back(ind);
	  intersecting.push_back(ind);
	}
    }
    intersecting.unique();
    //intersecting.merge(intersect_diff);
//     intersecting.unique();
     intersecting.sort();
     
    //    cout << lambda << endl;
//     for (typename std::list<Index>::const_iterator  it = intersecting.begin();
// 	 it != intersecting.end(); ++it) {
//       cout << *it << endl;
//     }
//     cout << "###############" << endl;


  }
  
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_singular_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu)
  {

    return true;

    switch (lambda.p() == mu.p()) {
    case 0:
      // different patches
      // TODO
      return true;
    case 1:
      // same patches
      //      return true;
      return WaveletTL::intersect_singular_support<IBASIS,DIM_d>
	(
	 *frame.bases()[lambda.p()],
	 typename CubeBasis<IBASIS,DIM_d>::Index(lambda.j(),
						 lambda.e(),
						 lambda.k(),
						 frame.bases()[lambda.p()]),
	 typename CubeBasis<IBASIS,DIM_d>::Index(mu.j(),
						 mu.e(),
						 mu.k(),
						 frame.bases()[mu.p()])
	 );
    }
    // dummy
    return true;
  }

  template <unsigned int DIM>
  inline int edgesIntersect (const Point<DIM>& A, const Point<DIM>& B,
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

// implementation for frame_support.h


#include <utils/fixed_array1d.h>
#include <typeinfo>
#include <cube/cube_support.h>
#include <map>
#include <time.h>

using namespace WaveletTL;

using std::type_info;
using std::map;

namespace FrameTL
{

  inline bool eq (const double x, const double y)
  {
    const double eps = 1.0e-10;
    return fabs(x-y) < eps;
  }

  inline bool lt (const double x, const double y)
  {
    return x < y;
//     const double eps = 1.0e-15;
//     return x-y < eps;
  }

  inline bool gt (const double x, const double y)
  {
    return ! (lt(x,y) || eq(x,y));
  }



  template <unsigned int DIM>
  inline
  unsigned short int pos_wrt_line (const Point<DIM>& p,
				   const Point<DIM>& p1, const Point<DIM>&  p2)
  {
   
    assert ( DIM == 2 );

    double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    
    //    cout << "d = " << d << endl;

    if (fabs(d) < 1.0e-15)
      return 2;

    if( d > 0.0 )
      return 1;
    else if (d < 0.0)
      return 0;

    // dummy
    return 3;
  }


  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
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
  inline
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
		  const Point<DIM_m>& p)
  {

    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>))
	 )
      {
	const double dx = 1.0 / (1 << supp_lambda->j);
	Point<DIM_d> a;
	Point<DIM_d> b;
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx), a);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx), b);

	return (a[0] <= p[0]) && (p[0] <= b[0]);
	
      }

    // just to get a better performance in the 
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
    
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(LinearBezierMapping))
	 ||
	 (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>))
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
	
	return true;;

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

    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>) &&
	  typeid(*frame.atlas()->charts()[mu.p()])     == 
	  typeid(AffineLinearMapping<2>)
	  ||
	  (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(SimpleAffineLinearMapping<2>) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(SimpleAffineLinearMapping<2>)))
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
    
    // both charts are LinearBezierMappings
    if ( ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(LinearBezierMapping) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(LinearBezierMapping) )
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
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->a[1]*dx1), poly1[0]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->a[1]*dx1), poly1[1]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->b[1]*dx1), poly1[2]);
	frame.atlas()->charts()[lambda.p()]->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->b[1]*dx1), poly1[3]);

	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), poly2[0]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->a[1]*dx2), poly2[1]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), poly2[2]);
	frame.atlas()->charts()[mu.p()]->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->b[1]*dx2), poly2[3]);

// 	cout << quadrangles_intersect2 (poly1,poly2) << " " << quadrangles_intersect (poly1,poly2) << endl;
//  	if ((quadrangles_intersect2 (poly1,poly2) != quadrangles_intersect (poly1,poly2))/* && (mu.p()!=lambda.p())*/) {
// 	  cout << quadrangles_intersect2 (poly1,poly2) << " " << quadrangles_intersect (poly1,poly2) << endl;
// 	  cout << lambda << " " << mu << endl;
//  	  cout << poly1[0] << " " << poly1[1] << " " << poly1[2] << " " << poly1[3] << endl;
//  	  cout << poly2[0] << " " << poly2[1] << " " << poly2[2] << " " << poly2[3] << endl;
// 	  abort();
// 	}
// 	//cout << quadrangles_intersect (poly1,poly2) << endl;;
// 	//	abort();
	return quadrangles_intersect (poly1,poly2);
      }
    else return false;
  }
 
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda)
  {
    
    typedef WaveletTL::CubeBasis<IBASIS,DIM_d> CUBEBASIS;
    typedef typename CUBEBASIS::Index CubeIndex;
    const typename CUBEBASIS::Support* supp_mu = &((frame.all_supports)[mu.number()]);;

//     WaveletTL::support<IBASIS,DIM_d>(*frame.bases()[mu.p()], 
// 				     CubeIndex(mu.j(),
// 					       mu.e(),
// 					       mu.k(),
// 					       frame.bases()[mu.p()]),
// 				     supp_mu);


    
    const double dx1 = 1.0 / (1 << supp_lambda->j);
    const double dx2 = 1.0 / (1 << supp_mu->j);


    //cout << supp_mu->a[0] << " " << supp_mu->b[0] << " " << supp_mu->j << endl;

    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<1>) &&
	  typeid(*frame.atlas()->charts()[mu.p()])     == 
	  typeid(AffineLinearMapping<1>)) ||
	 (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(SimpleAffineLinearMapping<1>) &&
	  typeid(*frame.atlas()->charts()[mu.p()])     == 
	  typeid(SimpleAffineLinearMapping<1>)) )
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
#if 1
    if ( (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	  typeid(AffineLinearMapping<2>) &&
	  typeid(*frame.atlas()->charts()[mu.p()])     == 
	  typeid(AffineLinearMapping<2>)
	  ||
	  (typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(SimpleAffineLinearMapping<2>) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(SimpleAffineLinearMapping<2>)))
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
#endif
    if ( ( typeid(*frame.atlas()->charts()[lambda.p()]) ==
	   typeid(LinearBezierMapping) &&
	   typeid(*frame.atlas()->charts()[mu.p()])     == 
	   typeid(LinearBezierMapping) )
	 )
      {

	assert ( DIM_d == 2 && DIM_m == 2);

	FixedArray1D<Point<DIM_m>,4 > poly1;
	FixedArray1D<Point<DIM_m>,4 > poly2;

	const Chart<DIM_d,DIM_m>* chart_la = frame.atlas()->charts()[lambda.p()];
	const Chart<DIM_d,DIM_m>* chart_mu = frame.atlas()->charts()[mu.p()];

	// map the knots of the unit cube to patch
	// 0 -- 00
	// 1 -- 10
	// 2 -- 11
	// 3 -- 01
	chart_la->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->a[1]*dx1), poly1[0]);
	chart_la->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->a[1]*dx1), poly1[1]);
	chart_la->map_point(Point<DIM_d>(supp_lambda->b[0]*dx1,supp_lambda->b[1]*dx1), poly1[2]);
	chart_la->map_point(Point<DIM_d>(supp_lambda->a[0]*dx1,supp_lambda->b[1]*dx1), poly1[3]);

	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), poly2[0]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->a[1]*dx2), poly2[1]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), poly2[2]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->b[1]*dx2), poly2[3]);

//  	if (!((mu.number() == 0 && lambda.number()==2) || (mu.number() == 2 && lambda.number() == 0)))
//  	  return false;

// 	bool b1 = quadrangles_intersect2 (poly1,poly2);
// 	cout << "11111111111111111111" << endl;
// 	bool b2 = quadrangles_intersect (poly1,poly2);

//  	if ((b1 != b2)/* && (mu.p()!=lambda.p())*/) {
// 	  cout << "222222222222222" << b1 << " " << b2 << endl;
// 	  cout << lambda << " " << mu << endl;
//  	  cout << poly1[0] << " " << poly1[1] << " " << poly1[2] << " " << poly1[3] << endl;
//  	  cout << poly2[0] << " " << poly2[1] << " " << poly2[2] << " " << poly2[3] << endl;
// 	  //cout << "33333333333" << quadrangles_intersect2 (poly1,poly2) << " " << quadrangles_intersect (poly1,poly2) << endl;  
// 	  abort();
// 	}
	
	
	return quadrangles_intersect (poly1,poly2);

      }
    else return false;
  }
  
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
// 	cout << p11 << " " << p12 << " " << poly2[j-1] << " " << poly2[j] << endl;
// 	cout << "tmp = " << tmp << endl;
	//abort();
	if (tmp == 1)
	  countspez++;
	if (tmp == 3) {
	  return true;
	}
      }
      tmp = edgesIntersect(p11, p12, poly2[3], poly2[0]);
//       cout << p11 << " " << p12 << " " << poly2[3] << " " << poly2[0] << endl;
//       cout << "tmp = " << tmp << endl;
      if (tmp == 1)
	countspez++;
      if (tmp == 3)
	return true;
	
    }//end outer for
    p11 = poly1[3];
    p12 = poly1[0];
    for (unsigned int j = 1; j <= 3; j++) {
      tmp = edgesIntersect(p11, p12, poly2[j-1], poly2[j]);
//       cout << p11 << " " << p12 << " " << poly2[j-1] << " " << poly2[j] << endl;
//       cout << "tmp = " << tmp << endl;
      
      if (tmp == 1)
	countspez++;
      if (tmp == 3)
	return true;
    }
    tmp = edgesIntersect(p11, p12, poly2[3], poly2[0]);
//     cout << p11 << " " << p12 << " " << poly2[3] << " " << poly2[0] << endl;
//     cout << "tmp = " << tmp << endl;

    if (tmp == 1)
      countspez++;
    if (tmp == 3)
      return true;
    //if control reaches this point, then no 'ordinary intersection' exists
    //ATTENTION!!!!
    //this is only correct in case of convex polygons!!!!!!!
    if (countspez >= 2) {
      return true;
    }

    // remaining possible cases: one quadrangle contains the other
    // or they do not intersect
    countspez = 0;
    tmp = 0;
    //     cout << "passing the point " << endl;
	
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
	
    //    cout << "6.5 ----" << endl;
    //the same with reversed roles
    //loop over all knots of the first polygon
    for (unsigned int i = 0; i <= 3; i++) {
      //cout << "i = " << i << endl;
      //cout << poly1[i] << endl;
      //cout << poly2[3] << " " << poly2[0] << endl;
      result = true;
      tmp = FrameTL::pos_wrt_line(poly1[i], poly2[3], poly2[0]);
      //cout << "tmp1=" << tmp << endl;
      if (tmp == 1)
	countspez++;
      if (tmp == 0 || tmp ==2) {
	result = false;
	continue;
      }
		
      for (unsigned int j = 1; j<= 3 ; j++) {
	//cout << poly2[j-1] << " " << poly2[j] << endl;
	tmp = FrameTL::pos_wrt_line(poly1[i], poly2[j-1], poly2[j]);
	//cout << "tmp2=" << tmp << endl;
	if (tmp == 1)
	  countspez++;
	if (tmp == 0 || tmp ==2) {
	  result = false;
	  break;
	}
      }//end inner for
		
      if (result)
	return result;
      //cout << "7777777777777" << endl;
    }//end outer for

    return false;

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
//       cout << supp_lambda->a[i] << " " << supp_lambda->b[i] << " " << supp_lambda->j << endl;
//     for (unsigned int i = 0; i < DIM_d; i++)
//       cout << supp_mu->a[i] << " " << supp_mu->b[i] << " " << supp_mu->j << endl;

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
    
 //    cout << y0 << " " << y1 << endl;

    // compute lower left and upper right corner of intersecting cube
    FixedArray1D<FixedArray1D<double,2>,DIM_d > hyperCube_intersect;
    for (unsigned int i = 0; i < DIM_d; i++) {
      assert ( y0[i] <= y1[i] );
//       cout << supp_lambda->a[i]*dx1 << " " << y0[i] << endl;
//       cout << supp_lambda->b[i]*dx1 << " " << y1[i] << endl;
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
	  //cout << "inserting " << d << endl;
	  //cout << "in " << irregular_grid[i][d] << endl;
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
	//	cout << "x[i] = " << x[i] << endl;
	chart_mu->map_point(x,x_patch);
	chart_la->map_point_inv(x_patch,y0);
	//	cout << "y0[i] = " << y0[i] << endl;
	if ( hyperCube_intersect[i][0] <= y0[i] &&  y0[i] <= hyperCube_intersect[i][1]){
	  //cout << "in2 " << irregular_grid[i][y0[i]] << endl;
	  irregular_grid[i].push_back(y0[i]);
	  //cout << "inserting " << y0[i] << endl;
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

  //  static double time = 0;

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  void intersecting_wavelets (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			      const int j, const bool generators,
			      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting)
  {

    intersecting.erase(intersecting.begin(),intersecting.end());

    typedef AggregatedFrame<IBASIS,DIM_d,DIM_m> Frame;
    typedef typename Frame::Index Index;

    typedef typename CubeBasis<IBASIS,DIM_d>::Index CubeIndex;
    
    //cout << "LEVEL = " << frame.all_supports[lambda.number()].j << endl;


    const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda = &((frame.all_supports)[lambda.number()]);
    //cout << lambda.number() << endl;
        
    //    cout << supp_lambda->a[0] << " " << supp_lambda->b[0] << " " << supp_lambda->j << endl;
#if 1

    std::list<typename Frame::Index> intersect_diff;

    const Array1D<Array1D<Index> >* indices_levelwise = frame.indices();

    if ( generators ) {
      for (unsigned int i = 0; i < (*indices_levelwise)[0].size(); i++) {
	Index ind = (*indices_levelwise)[0][i];
	if (frame.atlas()->get_adjacency_matrix().get_entry(lambda.p(), ind.p()) && 
	    intersect_supports(frame, lambda, ind, supp_lambda) ){
	  intersect_diff.push_back(ind);
	}
      }
    }
    else {
      for (unsigned int i = 0; i < (*indices_levelwise)[j-frame.j0()+1].size(); i++) {
	Index ind = (*indices_levelwise)[j-frame.j0()+1][i];
	if (frame.atlas()->get_adjacency_matrix().get_entry(lambda.p(), ind.p()) && 
	    intersect_supports(frame,lambda, ind, supp_lambda) ){
	  intersect_diff.push_back(ind);
	}
      }
    }
    intersecting.merge(intersect_diff);
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
      // Thus, we do not check the typeid here. BE CAREFUL ANYWAY.

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
	
	//       cout << patch_lambda[0] << endl;
	//       cout << patch_lambda[1] << endl;
	
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
	
	//       cout << patch_lambda[0] << endl;
	//       cout << patch_lambda[1] << endl;
	//       cout << patch_lambda[2] << endl;
	//       cout << patch_lambda[3] << endl;
	
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->a[1]*dx2), poly_mu[0]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->a[1]*dx2), poly_mu[1]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->b[0]*dx2,supp_mu->b[1]*dx2), poly_mu[2]);
	chart_mu->map_point(Point<DIM_d>(supp_mu->a[0]*dx2,supp_mu->b[1]*dx2), poly_mu[3]);  
	
	bool fully_contained = 1;
	
	// first check whether all of the four corners of the support of
	// mu are contained in the patch of lambda
	for (int i = 0; i < 4; i++) {
	  //cout << poly_mu[i] << endl;
	  
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
	
	//       cout << poly_mu_pulled_back[0] << endl;
	//       cout << poly_mu_pulled_back[1] << endl;
	//       cout << poly_mu_pulled_back[2] << endl;
	//       cout << poly_mu_pulled_back[3] << endl;
	
	FixedArray1D<double, 4> k1;
	FixedArray1D<double, 4> k2;
	for (int i = 0; i < 4; i++) {
	  // determine that square [k1*2^-j1,(k1+1)*2^-j]x[k2*2^-j2,(k2+1)*2^-j2]
	  // where psi_lambda is smooth and that contains poly_mu_pulled_back[i]
	  // if these squares are NOT all the same, the singular supports intersect
	  
	  k1[i] = floor(poly_mu_pulled_back[i][0] / dx1);
	  k2[i] = floor(poly_mu_pulled_back[i][1] / dx1);
	  
	  // 	cout << "k1 = " << k1[i] << endl;
	  // 	cout << "k2 = " << k2[i] << endl;
	  
	  if (i > 0) {
	    //compare with the preceeding point
	    if ((k1[i] != k1[i-1]) || (k2[i] != k2[i-1]))
	      return true;
	  }
	}
	//cout << "#################### " << lambda << " " << mu << endl;
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
//       cout << u1 << " " << v1 << " " << w1 << endl;
//       cout << u2 << " " << v2 << " " << w2 << endl;
//    if (!(u2 == 0.)) {
    if (!eq(u2,0.0)) {
      // gauss elimination
      v2 += -(u2 / u1)*v1;
      w2 += -(u2 / u1)*w1;
      //       cout << u1 << " " << v1 << " " << w1 << endl;
      //       cout << u2 << " " << v2 << " " << w2 << endl;
      
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
	//	cout << "s= " << s << " t=  " << t << endl;
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

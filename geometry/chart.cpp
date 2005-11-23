// implementation for chart.h

#include <iostream>
#include <cassert>
#include <sstream>

using std::cout;
using std::endl;

namespace MathTL
{
  template <unsigned int DIM>
  const string AffineLinearMapping<DIM>::className = "AffineLinearMapping";

  template <unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream& s, const Chart<DIM_d,DIM_m>& c)
  {
    s << c.to_string();
    return s;
  }

  template <unsigned int DIM>
  AffineLinearMapping<DIM>::AffineLinearMapping()
    : A_(DIM, DIM), A_inv(DIM, DIM), det_A(1.0), b_()
  {
    A_.diagonal(DIM, 1.0);
    A_inv.diagonal(DIM, 1.0);
  }

  template <unsigned int DIM>
  AffineLinearMapping<DIM>::AffineLinearMapping(const Matrix<double>& A, const Point<DIM>& b)
    : A_(A), A_inv(DIM, DIM), b_(b)
  {
    assert(A.row_dimension() == DIM && A.column_dimension() == DIM);
    assert(DIM == 1 || DIM == 2);
    switch (DIM) {
    case 1:
      det_A = A_(0,0);
      A_inv(0,0) = 1./det_A;
      break;
    case 2:
      det_A = A_(0,0)*A_(1,1)-A_(1,0)*A_(0,1);
      A_inv(0,0) = A_(1,1)/det_A;
      A_inv(0,1) = -A_(0,1)/det_A;
      A_inv(1,0) = -A_(1,0)/det_A;
      A_inv(1,1) = A_(0,0)/det_A;
      break;
    }
  }

  template <unsigned int DIM>
  inline
  void
  AffineLinearMapping<DIM>::map_point(const Point<DIM>& x, Point<DIM>& y) const {
    A_.apply(x, y);
    y += b_;
  }
    
  template <unsigned int DIM>
  inline
  void
  AffineLinearMapping<DIM>::map_point_inv(const Point<DIM>& x, Point<DIM>& y) const {
    A_inv.apply(Point<DIM>(x-b_), y);
  }

  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Gram_factor(const Point<DIM>& x) const {
    return sqrt(fabs(det_A));
  }

  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Gram_D_factor(const unsigned int i, const Point<DIM>& x) const {
    return 0.0;
  }
  
  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Dkappa_inv(const unsigned int i, const unsigned int j,
				       const Point<DIM>& x) const
  {
    return A_inv(i, j);
  }

  template <unsigned int DIM>
  inline
  const bool
  AffineLinearMapping<DIM>::in_patch(const Point<DIM>& x) const {
    Point<DIM> y;
    map_point_inv(x, y);
    for (unsigned int i = 0; i < DIM; i++)
      if (y[i] < 0 || y[i] > 1) return false;
    return true;
  }

  template <unsigned int DIM>
  const string
  AffineLinearMapping<DIM>::to_string() const {
    std::stringstream strs;
 
    strs << className << ":" << endl << "A =" << endl
	 << A()
	 << "b = "
	 << b() << endl;

    return strs.str();
  }
//################### SimpleAffineLinear ###################
  template <unsigned int DIM>
  const string SimpleAffineLinearMapping<DIM>::className = "SimpleAffineLinearMapping";
  
  template <unsigned int DIM>
  SimpleAffineLinearMapping<DIM>::SimpleAffineLinearMapping()
    : det_A(1.0), b_()
  {
  }

  template <unsigned int DIM>
  SimpleAffineLinearMapping<DIM>::SimpleAffineLinearMapping(
							    const FixedArray1D<double,DIM>& A,
							    const Point<DIM>& b
							    )
    : A_(A), b_(b)
  {
    det_A = 1.;
    for (unsigned int i = 0; i < DIM; i++) {
      det_A *= A_[i];
      A_inv[i] = 1./A_[i];
    }
  }

  template <unsigned int DIM>
  inline
  void
  SimpleAffineLinearMapping<DIM>::map_point(const Point<DIM>& x, Point<DIM>& y) const {
    for (unsigned int i = 0; i< DIM; i++)
      y[i] = A_[i] * x[i] + b_[i];
  }
    
  template <unsigned int DIM>
  inline
  void
  SimpleAffineLinearMapping<DIM>::map_point_inv(const Point<DIM>& x, Point<DIM>& y) const {
    for (unsigned int i = 0; i< DIM; i++)
      y[i] = A_inv[i] * (x[i]-b_[i]);
  }

  template <unsigned int DIM>
  inline
  const double
  SimpleAffineLinearMapping<DIM>::Gram_factor(const Point<DIM>& x) const {
    return sqrt(fabs(det_A));
  }

  template <unsigned int DIM>
  inline
  const double
  SimpleAffineLinearMapping<DIM>::Gram_D_factor(const unsigned int i, const Point<DIM>& x) const {
    return 0.0;
  }
  
  template <unsigned int DIM>
  inline
  const double
  SimpleAffineLinearMapping<DIM>::Dkappa_inv(const unsigned int i, const unsigned int j,
				       const Point<DIM>& x) const
  {
    if (i != j)
      return 0.;
    return A_inv[i];
  }

  template <unsigned int DIM>
  inline
  const bool
  SimpleAffineLinearMapping<DIM>::in_patch(const Point<DIM>& x) const {
    Point<DIM> y;
    map_point_inv(x, y);
    for (unsigned int i = 0; i < DIM; i++)
      if (y[i] < 0 || y[i] > 1) return false;
    return true;
  }

  template <unsigned int DIM>
  const string
  SimpleAffineLinearMapping<DIM>::to_string() const {
    std::stringstream strs;
 
    strs << className << ":" << endl << "A =" << endl
	 << A()
	 << "b = "
	 << b() << endl;

    return strs.str();
  }
//################### LinearBezierMapping ###################
  const string LinearBezierMapping::className = "LinearBezierMapping";

  inline
  LinearBezierMapping::LinearBezierMapping ()
  {
  }

  inline
  LinearBezierMapping:: LinearBezierMapping (const LinearBezierMapping& kappa)
    : b_00(kappa.get_b_00()),  b_01(kappa.get_b_01()),
      b_10(kappa.get_b_10()),  b_11(kappa.get_b_11())
  {
    setup();
  }

  inline
  LinearBezierMapping:: LinearBezierMapping (const Point<2> & _b_00, const Point<2> & _b_01,
					     const Point<2> & _b_10, const Point<2> & _b_11)
    : b_00(_b_00),  b_01(_b_01), b_10(_b_10),  b_11(_b_11)
  {
    setup();
  }

  inline
  void LinearBezierMapping::setup ()
  {

    d_ds_d_dt_kappa_r.resize(2);
    d_dt_d_ds_kappa_r.resize(2);

    d_ds_d_dt_kappa_r[0] = b_00[0] - b_01[0] - b_10[0] + b_11[0];
    d_ds_d_dt_kappa_r[1] = b_00[1] - b_01[1] - b_10[1] + b_11[1];
    d_dt_d_ds_kappa_r = d_ds_d_dt_kappa_r;
	
    min_b00_plus_b10.resize(2);
    min_b00_plus_b01.resize(2);

    min_b00_plus_b10[0] = b_10[0] - b_00[0];
    min_b00_plus_b10[1] = b_10[1] - b_00[1];
	
    min_b00_plus_b01[0] = b_01[0] - b_00[0];
    min_b00_plus_b01[1] = b_01[1] - b_00[1];
	
    rot_angle = 0;
	
    //ATTENTION with this initial value!!!!!!
    cos_rot_angle = 1.0;
    sin_rot_angle = 0.0;
	
    shearing_param = 0.0;
    scaleX = 0.0;
    scaleY = 0.0;
	
    //get the signum of det (D kappa),
    //this is the same for each point since det (D kappa) is
    //continuous and never equals zero!!
    sgn_det_D = (int) (det_DKappa(Point<2>(0,0)) / det_DKappa(Point<2>(0,0)));

    //setup the generic quadrangle belonging to the parametrized patch,
    //we do this once and for all to be able to invert the bezier mapping
	
    b_gen_00[0] = 0;
    b_gen_00[0] = 0;

    b_gen_01[0]  = b_01[0] - b_00[0];
    b_gen_01[1]  = b_01[1] - b_00[1];

    b_gen_10[0]  = b_10[0] - b_00[0];
    b_gen_10[1]  = b_10[1] - b_00[1];
	
    b_gen_11[0]  = b_11[0] - b_00[0];
    b_gen_11[1]  = b_11[1] - b_00[1];

    //rotate the coordinate system in such a way that b_01 lies ion the second axis.
    Matrix<double> R(2,2);

    //rotate the whole quadrangle
    if (! (b_gen_01[0] == 0.0) )
      {
	//rotation angle
	rot_angle = atan(b_gen_01[0] / b_gen_01[1]);
	//cout << "rot_angle = " << rot_angle << endl;

	cos_rot_angle = cos(rot_angle);
	sin_rot_angle = sin(rot_angle);
	R(0,0) = cos_rot_angle;
	R(0,1) = -sin_rot_angle; // = -sin(rot_angle)
	R(1,0) = -R(0,1); // = sin(rot_angle)
	R(1,1) =  R(0,0);

	cout << "R = " << R << endl;

	Point<2> tmp;

	R.apply(b_gen_01,tmp);
	b_gen_01 = tmp; 

	R.apply(b_gen_10,tmp);
	b_gen_10 = tmp;

	R.apply(b_gen_11,tmp);
	b_gen_11 = tmp;
      }

    //rescale the second coordinate by 1/b_01(2) ==> b_01 = (0,1)^T
    if ( b_gen_01[1] != 0)
      {
	scaleY = 1.0 / b_gen_01[1];
	b_gen_00[1] *= scaleY;
	b_gen_01[1] *= scaleY;
	b_gen_10[1] *= scaleY;
	b_gen_11[1] *= scaleY;	
      }
    //rescale the first coordinate by 1/b_10(1) ==> b_10(1) = 1
    if ( b_gen_10[0] != 0)
      {
	scaleX = 1.0 / b_gen_10[0];
	b_gen_00[0] *= scaleX;
	b_gen_01[0] *= scaleX;
	b_gen_10[0] *= scaleX;
	b_gen_11[0] *= scaleX;
      }
  
    //finally apply a shearing so that b_10 = (1,0)^T
    if ( b_gen_10[1] != 0)
      {
	shearing_param = -b_gen_10[1];
	R(0,0) = 1;
	R(0,1) = 0;
	R(1,0) = shearing_param;
	R(1,1) = 1;
	
	Point<2> tmp;

	R.apply(b_gen_00,tmp);
	b_gen_00 = tmp; 

	R.apply(b_gen_01,tmp);
	b_gen_01 = tmp; 

	R.apply(b_gen_10,tmp);
	b_gen_10 = tmp;

	R.apply(b_gen_11,tmp);
	b_gen_11 = tmp; 

      }	
    cout << b_gen_00[0] << " " <<  b_gen_00[1] << endl
	 << b_gen_01[0] << " " <<  b_gen_01[1] << endl
	 << b_gen_10[0] << " " <<  b_gen_10[1] << endl
	 << b_gen_11[0] << " " <<  b_gen_11[1] << endl;

  }
  inline
  const Point<2>& LinearBezierMapping::get_b_00() const
  {
    return b_00;
  }
  
  inline
  const Point<2>& LinearBezierMapping::get_b_10() const
  {
    return b_10;
  }

  inline
  const Point<2>& LinearBezierMapping::get_b_01() const
  {
    return b_01;
  }

  inline
  const Point<2>& LinearBezierMapping::get_b_11() const
  {
    return b_11;
  }

//   inline
//   const unsigned short int LinearBezierMapping::get_sgn_det_D() const
//   {
//     return sgn_det_D;
//   }

  inline
  void LinearBezierMapping::map_point(const Point<2>& pc, Point<2>& pp) const
  {    
    double t1 = 1-pc[0];
    double t2 = 1-pc[1];
    double t3 = t1*t2;
    double t4 = pc[0]*t2;
    double t5 = t1*pc[1];
    t1 = pc[0]*pc[1];
	
    pp[0]= t3*b_00[0] + t4*b_10[0] +  t5*b_01[0] + t1*b_11[0];
    pp[1]= t3*b_00[1] + t4*b_10[1] +  t5*b_01[1] + t1*b_11[1];
  }

  inline
  void LinearBezierMapping::map_point_inv(const Point<2>& pp, Point<2>& pc) const
  {
    pc[0] = pp[0];
    pc[1] = pp[1];
    
    //shift the quadrangle so that b_00 lies in the origin.
    pc[0] -= b_00[0];
    pc[1] -= b_00[1];

    //rotate the coordinate system in such a way that b_01 lies ion the second axis
    Matrix<double> R(2,2);
    R(0,0) = cos_rot_angle ;
    R(0,1) = -sin_rot_angle ;
    R(1,0) = -R(0,1) ;
    R(1,1) = cos_rot_angle;

    //apply the rotation
    Point<2> tmpP;
    R.apply(pc,tmpP);
    pc = tmpP;

    //rescale the second coordinate by 1/b_01(2) ==> b_01 = (0,1)^T
    pc[1] *= scaleY;
    
    //rescale the first coordinate by 1/b_10(1) ==> b_10(1) = 1
    pc[0] *= scaleX;
	
    //finally apply a shearing so that b_10 = (1,0)^T
    R(0,0) = 1;
    R(0,1) = 0;
    R(1,0) = shearing_param;
    R(1,1) = 1;

    R.apply(pc,tmpP);
    pc = tmpP;

    if (0 == pc[0] && 0 == pc[1]) {

      pc[0] = 0;
      pc[1] = 0;
      return;
    }
    double s = 0, t = 0;
    double tmp = 0;
    short int i,j;
    //we have to take care of the correct role of the variables
    //in our calculations we assumed coordPatch(0) != 0
    //if this is not the case we obtain coordPatch(1) != 0 and
    //the role of the variables has to reversed
   
    bool swapped_var = 0;
    if (! (0 == pp[0]) )
      {
	i=0;
	j=1;
	//do nothing else		
      }
    else
      {
	swapped_var = 1;
	i = 1;
	j = 0;
	//cout << "swapping coordinates in inverse linear bezier mapping" << endl;
	//swap coordinates
	tmp = pc[0];
	pc[0] = pc[1];
	pc[1]= tmp;
      }
    double h1 = 1 + pc[0]*(b_gen_11[j]-1);
    double h2 = b_gen_11[i] - 1;
    double h3 = h1 / h2 - pc[1];
    double h4  = h3 * h3;
		
    if (b_gen_11[i] == 1)
      {
	t = pc[1] / h1;
	//s = xcube / (1 - t + t*bgen_11[i]);
	s = pc[0] / (1 + t * h2);
	
	if (swapped_var) {
	  pc[0] = t;
	  pc[1] = s;
	}
	else {
	  pc[0] = s;
	  pc[1] = t;
	}
	return;
      }
	
    if (b_gen_11[i] > 1)
      t = -0.5 * h3 + sqrt(0.25 * h4 + pc[1]/h2);
    else
      t = -0.5 * h3 - sqrt(0.25 * h4 + pc[1]/h2);

    s = pc[0] / (1 + t * h2);

    if (swapped_var) {
      pc[0] = t;
      pc[1] = s;
    }
    else {
      pc[0] = s;
      pc[1] = t;
    }
  }


  inline
  const double LinearBezierMapping::det_DKappa(const Point<2>& p) const
  {
    return (p[1] * d_ds_d_dt_kappa_r[0] + min_b00_plus_b10[0]) * 
      (p[0] * d_ds_d_dt_kappa_r[1] + min_b00_plus_b01[1]) -
      (p[1] * d_ds_d_dt_kappa_r[1] + min_b00_plus_b10[1]) *
      (p[0] * d_ds_d_dt_kappa_r[0] + min_b00_plus_b01[0]);
    
  }

  inline
  const double LinearBezierMapping::Gram_factor(const Point<2>& p) const
  {
    //anyway: ((det ((D Kappa)^T) * D Kappa))^1/2)^1/2 = | det D Kappa |^1/2
    //in the square case
    return sqrt(fabs(det_DKappa(p)));
    
  }

  inline
  const double LinearBezierMapping::Gram_D_factor(const unsigned int i,
						  const Point<2>& p) const
  {
    double partial_i_det_DKappa = 0.0;

    switch(i) {
    case 0:
        partial_i_det_DKappa = (p[1] * d_ds_d_dt_kappa_r[0] + min_b00_plus_b10[0]) * d_ds_d_dt_kappa_r[1]
	 -
	 (p[1] * d_ds_d_dt_kappa_r[1] + min_b00_plus_b10[1]) *  d_ds_d_dt_kappa_r[0];
      break;
    case 1:
       partial_i_det_DKappa = d_dt_d_ds_kappa_r[0] * (p[0] * d_ds_d_dt_kappa_r[1] + min_b00_plus_b01[1])
      -
      d_dt_d_ds_kappa_r[1] * (p[0] * d_ds_d_dt_kappa_r[0] + min_b00_plus_b01[0]);
      break;
    };
 
    return 0.5 * (1.0 / Gram_factor(p)) * sgn_det_D * partial_i_det_DKappa;
   
  }
  //(dim = i, j = direc)
  inline
  const double LinearBezierMapping::partial_i_Kappa_j(const unsigned int i, const unsigned int j,
						      const Point<2>& p) const
  {
    return (i==1) ? p[0] * d_ds_d_dt_kappa_r[j] + min_b00_plus_b01[j] : 
      p[1] * d_ds_d_dt_kappa_r[j] + min_b00_plus_b10[j];
  }

  inline
  const double LinearBezierMapping::Dkappa_inv(const unsigned int i, const unsigned int j,
					       const Point<2>& x) const
  {

    switch ((i << 1) + j) {
    case 0: // i == 0 && j == 0
      return partial_i_Kappa_j(1,1,x)/det_DKappa(x);
    case 1: // i == 0 && j == 1
      return -partial_i_Kappa_j(1,0,x)/det_DKappa(x);
    case 2: // i == 1 && j == 0
      return -partial_i_Kappa_j(0,1,x)/det_DKappa(x);
    case 3: // i == 1 && j == 1
      return partial_i_Kappa_j(0,0,x)/det_DKappa(x);
    }

//     if (i==0 && j==0)
//       return partial_i_Kappa_j(1,1,x)/det_DKappa(x);
//     if (i==0 && j==1)
//      return -partial_i_Kappa_j(1,0,x)/det_DKappa(x);
//     if (i==1 && j==0)
//      return -partial_i_Kappa_j(0,1,x)/det_DKappa(x);
//     if (i==1 && j==1)
//       return partial_i_Kappa_j(0,0,x)/det_DKappa(x);

    //dummy
    return 0;

  }

  /*!
    Checks whether p lies left of right of or on the line specified by the
    Points p1 and p2. This line is oriented by the vector starting in p1 and
    ending in p2.
    returning 0 means RIGHT OF LINE
    returning 1 means LEFT OF LINE
    returning 2 means ON LINE
   */
  template <unsigned int DIM>
  inline
  unsigned short int pos_wrt_line (const Point<DIM>& p,
					 const Point<DIM>& p1, const Point<DIM>&  p2)
  {
    double d = (p(1)-p1(1)) * (p2(0)-p1(0)) - (p(0)-p1(0)) * (p2(1)-p1(1));
    
    if( d > 0.0 )
      return 1;
    else if (d < 0.0)
      return 0;
    else
      return 2;
  }

  inline
  const bool LinearBezierMapping::in_patch(const Point<2>& x) const
  {
    //make sure to walk through the vertices counter clockwise!!!
    unsigned short int res = 
      pos_wrt_line(x, b_00, b_10);
    if (res == 0)
      return false;
    res = 
      pos_wrt_line(x, b_10, b_11);
    if (res == 0)
      return false;
    res = 
      pos_wrt_line(x, b_11, b_01);
    if (res == 0)
      return false;
    res = 
      pos_wrt_line(x, b_01, b_00);
    if (res == 0)
      return false;

    return true;;
  }
  
  const string LinearBezierMapping::to_string() const
  {   
    std::stringstream strs;
 
    strs << className << ":" << endl
	 << "b_00=" << b_00 << endl
	 << "b_01=" << b_01 << endl
	 << "b_10=" << b_10 << endl
	 << "b_11=" << b_11;
 
    return strs.str();
  }
}

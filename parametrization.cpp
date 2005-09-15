// implementation for parametrization.h
#include <math.h>
#include <cmath>

#include <algebra/matrix.h>
#include <misc.cpp>

using MathTL::Matrix;
using std::swap;

namespace FrameTL
{

//################### LinearBezierMapping ###################
  
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

    d_ds_d_dt_kappa_r(0) = b_00(0) - b_01(0) - b_10(0) + b_11(0);
    d_ds_d_dt_kappa_r(1) = b_00(1) - b_01(1) - b_10(1) + b_11(1);
    d_dt_d_ds_kappa_r = d_ds_d_dt_kappa_r;
	
    min_b00_plus_b10.resize(2);
    min_b00_plus_b01.resize(2);

    min_b00_plus_b10(0) = b_10(0) - b_00(0);
    min_b00_plus_b10(1) = b_10(1) - b_00(1);
	
    min_b00_plus_b01(0) = b_01(0) - b_00(0);
    min_b00_plus_b01(1) = b_01(1) - b_00(1);
	
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
    sgn_det_D = (int) (det_D(Point<2>(0,0)) / det_D(Point<2>(0,0)));

    //setup the generic quadrangle belonging to the parametrized patch,
    //we do this once and for all to be able to invert the bezier mapping
	
    b_gen_00(0) = 0;
    b_gen_00(0) = 1;

    b_gen_01(0)  = b_01(0) - b_00(0);
    b_gen_01(1)  = b_01(1) - b_00(1);

    b_gen_10(0)  = b_10(0) - b_00(0);
    b_gen_10(1)  = b_10(1) - b_00(1);
	
    b_gen_11(0)  = b_11(0) - b_00(0);
    b_gen_11(1)  = b_11(1) - b_00(1);

    //rotate the coordinate system in such a way that b_01 lies ion the second axis.
    Matrix<double> R(2,2);

    //rotate the whole quadrangle
    if (! (b_gen_01(0) == 0.0) )
      {
	//rotation angle
	rot_angle = atan(b_gen_01(0) / b_gen_01(1));

	cos_rot_angle = cos(rot_angle);
	sin_rot_angle = sqrt(1 - cos(rot_angle)*cos(rot_angle));
	R(0,0) = cos_rot_angle;
	R(0,1) = -sin_rot_angle; // = -sin(rot_angle)
	R(1,0) = -R(0,1); // = sin(rot_angle)
	R(1,1) =  R(0,0);

	b_gen_00 = apply(R, b_gen_00);
	b_gen_01 = apply(R, b_gen_01);
	b_gen_10 = apply(R, b_gen_10);
	b_gen_11 = apply(R, b_gen_11);

      }
 
    //rescale the second coordinate by 1/b_01(2) ==> b_01 = (0,1)^T
    if ( b_gen_01(1) != 0)
      {
	scaleY = 1.0 / b_gen_01(1);
	b_gen_00(1) *= scaleY;
	b_gen_01(1) *= scaleY;
	b_gen_10(1) *= scaleY;
	b_gen_11(1) *= scaleY;	
      }
    //rescale the first coordinate by 1/b_10(1) ==> b_10(1) = 1
    if ( b_gen_10(0) != 0)
      {
	scaleX = 1.0 / b_gen_10(0);
	b_gen_00(0) *= scaleX;
	b_gen_01(0) *= scaleX;
	b_gen_10(0) *= scaleX;
	b_gen_11(0) *= scaleX;
      }
  
    //finally apply a shearing so that b_10 = (1,0)^T
    if ( b_gen_10(1) != 0)
      {
	shearing_param = -b_gen_10(1);
	R(0,0) = 1;
	R(0,1) = 0;
	R(1,0) = shearing_param;
	R(1,1) = 1;
	b_gen_00 = apply(R, b_gen_00);
 	b_gen_01 = apply(R, b_gen_01);
 	b_gen_10 = apply(R, b_gen_10);
 	b_gen_11 = apply(R, b_gen_11);
      }	
//     cout << b_gen_00[0] << " " <<  b_gen_00[1] << endl
// 	 << b_gen_01[0] << " " <<  b_gen_01[1] << endl
// 	 << b_gen_10[0] << " " <<  b_gen_10[1] << endl
// 	 << b_gen_11[0] << " " <<  b_gen_11[1] << endl;

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
  
  inline
  const unsigned short int LinearBezierMapping::get_sgn_det_D() const
  {
    return sgn_det_D;
  }
  
  inline
  void LinearBezierMapping::mapPoint(Point<2>& pp, const Point<2>& pc) const
  {    
    double t1 = 1-pc(0);
    double t2 = 1-pc(1);
    double t3 = t1*t2;
    double t4 = pc(0)*t2;
    double t5 = t1*pc(1);
    t1 = pc(0)*pc(1);
	
    pp(0)= t3*b_00(0) + t4*b_10(0) +  t5*b_01(0) + t1*b_11(0);
    pp(1)= t3*b_00(1) + t4*b_10(1) +  t5*b_01(1) + t1*b_11(1);
  }

  inline
  void LinearBezierMapping::mapPointInv(Point<2>& pc, const Point<2>& pp) const
  {
    pc(0) = pp(0);
    pc(1) = pp(1);
    
    //shift the quadrangle so that b_00 lies in the origin.
    pc(0) -= b_00(0);
    pc(1) -= b_00(1);

    //rotate the coordinate system in such a way that b_01 lies ion the second axis
    Matrix<double> R(2,2);
    R(0,0) = cos_rot_angle ;
    R(0,1) = -sin_rot_angle ;
    R(1,0) = -R(0,1) ;
    R(1,1) = cos_rot_angle;
    
    //apply the rotation
    pc = apply(R, pc);

    //rescale the second coordinate by 1/b_01(2) ==> b_01 = (0,1)^T
    pc(1) *= scaleY;
    
    //rescale the first coordinate by 1/b_10(1) ==> b_10(1) = 1
    pc(0) *= scaleX;
	
    //finally apply a shearing so that b_10 = (1,0)^T
    R(0,0) = 1;
    R(0,1) = 0;
    R(1,0) = shearing_param;
    R(1,1) = 1;

    pc = apply(R, pc);

    if (0 == pc(0) && 0 == pc(1)) {

      pc(0) = 0;
      pc(1) = 0;
      return;
    }
    double s = 0, t = 0;
    double tmp = 0;
    short int i,j;
    //we have to take care of the correct role of the variables
    //in our calculations we assumed coordPatch(1) != 0
    //if this is not the case we obtain coordPatch(1) != 0 and
    //the role of the variables has to reversed
   

    if (! (0 == pp(0)) )
      {
	i=0;
	j=1;
	//do nothing else		
      }
    else
      {
	i = 1;
	j = 0;
	cout << "swapping coordinates in inverse linear bezier mapping" << endl;
	//swap coordinates
	tmp = pc(0);
	pc(0) = pc(1);
	pc(1)= tmp;
      }
    double h1 = 1 + pc(0)*(b_gen_11(j)-1);
    double h2 = b_gen_11(i) - 1;
    double h3 = h1 / h2 - pc(1);
    double h4  = h3 * h3;
		
    if (b_gen_11(i) == 1)
      {
	t = pc(1) / h1;
	//s = xcube / (1 - t + t*bgen_11[i]);
	s = pc(0) / (1 + t * h2);
	pc(0) = s;
	pc(1) = t;
	return;
      }
	
    if (b_gen_11(i) > 1)
      t = -0.5 * h3 + sqrt(0.25 * h4 + pc(1)/h2);
    else
      t = -0.5 * h3 - sqrt(0.25 * h4 + pc(1)/h2);

    s = pc(0) / (1 + t * h2);

    pc(0) = s;
    pc(1) = t;
  }

  inline
  const double LinearBezierMapping::det_D(const Point<2>& p) const
  {
    return (p(1) * d_ds_d_dt_kappa_r(0) + min_b00_plus_b10(0)) * 
      (p(0) * d_ds_d_dt_kappa_r(1) + min_b00_plus_b01(1)) -
      (p(1) * d_ds_d_dt_kappa_r(1) + min_b00_plus_b10(1)) *
      (p(0) * d_ds_d_dt_kappa_r(0) + min_b00_plus_b01(0));
    
  }
  inline    
  const double LinearBezierMapping::abs_Det_D(const Point<2>& p) const
  {
    return fabs(det_D(p));
  }

  inline
  const double LinearBezierMapping::d_x_det_D(const Point<2,double>& p) const
  {
    return (p(1) * d_ds_d_dt_kappa_r(0) + min_b00_plus_b10(0)) * d_ds_d_dt_kappa_r(1)
      -
      (p(1) * d_ds_d_dt_kappa_r(1) + min_b00_plus_b10(1)) *  d_ds_d_dt_kappa_r(0);
  }

  inline
  const double LinearBezierMapping::d_y_det_D(const Point<2,double>& p) const
  {
    return d_dt_d_ds_kappa_r(0) * (p(0) * d_ds_d_dt_kappa_r(1) + min_b00_plus_b01(1))
      -
      d_dt_d_ds_kappa_r(1) * (p(0) * d_ds_d_dt_kappa_r(0) + min_b00_plus_b01(0));
  }

  inline
  const double LinearBezierMapping::d_dim_kappa_direc( const unsigned short int& dim,
						       const unsigned short int& direc,
						       const Point<2>& p) const
  {
    return (dim==1) ? p(0) * d_ds_d_dt_kappa_r((direc==1) ? 1 : 0) + min_b00_plus_b01((direc==1) ? 1 : 0) : 
      p(1) * d_ds_d_dt_kappa_r((direc==1) ? 1 : 0) + min_b00_plus_b10((direc==1) ? 1 : 0);
  }

//   inline
//   const double LinearBezierMapping::d_dim_kappa_direc(const bool& dim, const bool& direc,
// 						      const Point<2>& p) const
//   {
//     return dim ? p(0) * d_ds_d_dt_kappa_r(direc ? 1 : 0) + min_b00_plus_b01(direc ? 1 : 0) : 
//       p(1) * d_ds_d_dt_kappa_r(direc ? 1 : 0) + min_b00_plus_b10(direc ? 1 : 0);
//   }

   LinearBezierMapping&
   LinearBezierMapping::operator = (const LinearBezierMapping& x)
   {
     b_00 = x.b_00;
     b_01 = x.b_01;
     b_10 = x.b_10;
     b_11 = x.b_11;
     
     setup();

     return *this;
   }

  inline
  std::ostream& operator << (std::ostream& os, const LinearBezierMapping& p)
  {
    
    os << "b_00=" << p.get_b_00() << endl
    << "b_01=" << p.get_b_01() << endl
    << "b_10=" << p.get_b_10() << endl
    << "b_11=" << p.get_b_11();
    
    return os;
  }
  
  //################### AffinLinearMapping ###################
 
 template <unsigned int DIM>
  inline
  AffinLinearMapping<DIM>::AffinLinearMapping ()
  {
  }

  template<unsigned int DIM>
  inline
  AffinLinearMapping<DIM>::AffinLinearMapping (const AffinLinearMapping<DIM>& kappa)
    : A(kappa.get_A()), b(kappa.get_b())
  {
    Point<DIM> p;//p equals origin is assumed!!
    sgn_det_D = (int) ( det_D(p) / det_D(p) );
  }

  template <unsigned int DIM>
  inline
  AffinLinearMapping<DIM>::AffinLinearMapping (const Matrix<double>& _A,
					       const Vector<double>& _b)
    : A(_A), b(_b)
  {
    Point<DIM> p;//p equals origin is assumed!!
    sgn_det_D = (int) ( det_D(p) / det_D(p) );
  }
  
  template <unsigned int DIM>
  inline
  AffinLinearMapping<DIM>& AffinLinearMapping<DIM>::operator = (const AffinLinearMapping& kappa)
  {
    A = kappa.get_A();
    b = kappa.get_b();
    
    Point<DIM> p;//p equals origin is assumed!!
    sgn_det_D = (int) ( det_D(p) / det_D(p) );

    return *this;
  }

  template <unsigned int DIM>
  inline
  const Matrix<double>& AffinLinearMapping<DIM>::get_A() const
  {
    return A;
  }
  
  template <unsigned int DIM>
  inline
  const Vector<double>&  AffinLinearMapping<DIM>::get_b() const
  {
    return b;
  }
  
  template <unsigned int DIM>
  inline
  const unsigned short int AffinLinearMapping<DIM>::get_sgn_det_D() const
  {
    return sgn_det_D;
  }

  template <unsigned int DIM>
  inline
  void  AffinLinearMapping<DIM>::mapPoint(Point<DIM>& pp,
					  const Point<DIM>& pc) const
  {
    pp = add(apply(A,pc),b);
  }
  
  template <unsigned int DIM>
  inline
  void  AffinLinearMapping<DIM>::mapPointInv(Point<DIM>& pc,
					     const Point<DIM>& pp) const
  {
    assert(DIM == 1 || DIM == 2);
    switch (DIM)
      {
      case 1:
	pc  =  sub(pp,b);
	pc *= 1.0/A(0,0);
	return;
      case 2:
	Matrix<double> A_Inv_m_detA(2,2);
	A_Inv_m_detA(0,0) =  A(1,1);
	A_Inv_m_detA(1,0) = -A(0,1);
	A_Inv_m_detA(0,1) = -A(1,0);
	A_Inv_m_detA(1,1) =  A(0,0);
	pc  =  apply(A_Inv_m_detA,sub(pp,b));
	pc /=  A(0,0)*A(1,1) - A(1,0)*A(0,1);
	return;
      };
  }
  
  template <unsigned int DIM>
  inline
  const double  AffinLinearMapping<DIM>::det_D(const Point<DIM>& p) const
  {
    assert(DIM == 1 || DIM == 2);
    switch (DIM)
      {
      case 1:
	return A(0,0);
      case 2:
	return  A(0,0)*A(1,1) - A(1,0)*A(0,1);
      };
  }
  
  template <unsigned int DIM>
  inline
  const double AffinLinearMapping<DIM>::abs_Det_D(const Point<DIM>& p) const
  {
    return  fabs(det_D(p));
  }
  
  template <unsigned int DIM>
  inline
  const double  AffinLinearMapping<DIM>::d_x_det_D(const Point<DIM>&) const
  {
    return 0;
  }
  
  template <unsigned int DIM>
  inline
  const double  AffinLinearMapping<DIM>::d_y_det_D(const Point<DIM>&) const
  {
    return 0;
  }
  
  template <unsigned int DIM>
  inline
  const double AffinLinearMapping<DIM>::d_dim_kappa_direc(const unsigned short int& dim,
							  const unsigned short int& direc,
							  const Point<DIM>& p) const
  {
    return 0;
  }

  template <unsigned int DIM>
  inline
  std::ostream& operator << (std::ostream& os, const AffinLinearMapping<DIM>& kappa)
  {
    os << "A = " << endl
       << kappa.get_A() << endl
       << "b = " << endl
       << kappa.get_b();
    return os;
}

}

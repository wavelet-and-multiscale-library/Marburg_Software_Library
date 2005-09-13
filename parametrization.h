// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_PARAMETRIZATION_H
#define _FRAMETL_PARAMETRIZATION_H

#include <iostream>
using std::cout;
using std::endl;

namespace Frame_TL
{

  class LinearBezierMapping
  {
  public:
    /*!
     */
    LinearBezierMapping () {};
    /*
    igpm::vector b_00;
    igpm::vector b_01;
    igpm::vector b_10;
	igpm::vector b_11;

	igpm::vector bgen_00;
	igpm::vector bgen_01;
	igpm::vector bgen_10;
	igpm::vector bgen_11;
	
	//equal(!)
	igpm::vector d_ds_d_dt_kappa_r;
	igpm::vector d_dt_d_ds_kappa_r;
	
	igpm::vector min_b00_plus_b10;
	igpm::vector min_b00_plus_b01;
	
	double cos_rot_angle;
	double sin_rot_angle;
	double rot_angle;
	double shearing_param;
	double scaleX;
	double scaleY;
	
	bool signumOfdetDKappa;
	
public:
	
	igpm::vector get_b_00() const { return b_00; };
	igpm::vector get_b_10() const { return b_10; };
	igpm::vector get_b_01() const { return b_01; };
	igpm::vector get_b_11() const { return b_11; };
	
    */

  }
}
// include implementation
#include <parametrization.cpp>

#endif

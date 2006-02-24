// implementation for p_basis.h

#include <cassert>
#include <cmath>
#include <numerics/schoenberg_splines.h>

namespace WaveletTL
{
  template <int d, int dT>
  PBasis<d,dT>::PBasis(const int s0, const int s1, const int sT0, const int sT1) {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
        
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;

    setup();
  }


  template <int d, int dT>
  void
  PBasis<d,dT>::setup() {

    // setup the refinement matrix block for all "real" primal boundary B-splines,
    // obeying the block structure (3.15)
    // (ignoring the current values of s0 and s1
    MathTL::SchoenbergKnotSequence<d> sknots;
    Matrix<double> ML_0;
    MathTL::compute_Bspline_refinement_matrix<d>(&sknots, ML_0);

    ML_.resize(3*dT+2*d-5, d+dT-2);
    ML_.set_block(0, 0, ML_0);
    for (int k = d; k <= d+dT-2; k++)
      for (int n = 2*k-d; n <= 2*k; n++)
	ML_.set_entry(n-1,k-1,cdf.a().get_coefficient(MultiIndex<int,1>(-(d/2)+n+d-2*k)));
    
    // compute the Gramian of all primal boundary B-splines
    // TODO!


   



    // for simplicity, we do not implement a generic setup here
    // but only setup the parameters for several important special cases from [P]
    switch(d) {
    case 2:
      switch(dT) {
      case 2:
	j0_ = 2;
	if (s0 == 0) {
	  MLT_.resize(5,2);
	  MLT_(0,0) =  5./4.;
	  MLT_(0,1) = -1./8.;
	  MLT_(1,0) =  3./2.;
	  MLT_(1,1) =  1./4.;
	  MLT_(2,0) = -3./4.;
	  MLT_(2,1) = 13./8.;
	  MLT_(3,1) =  1./2.;
	  MLT_(4,1) = -1./4.;
	}
	break;
      case 4:
	j0_ = 3;
	if (s0 == 0) {
	  MLT_.resize(11,4);
	  MLT_( 0,0) =   93./ 64 ;
	  MLT_( 0,1) = -241./768.;
	  MLT_( 0,2) =   41./384.;
	  MLT_( 0,3) = -  5./256.;
	  MLT_( 1,0) =   35./ 32.;
	  MLT_( 1,1) =  241./384.;
	  MLT_( 1,2) = - 41./192.;
	  MLT_( 1,3) =    5./128.;
	  MLT_( 2,0) = -  5./ 16.;
	  MLT_( 2,1) =  245./192.;
	  MLT_( 2,2) = - 13./ 96.;
	  MLT_( 2,3) =    1./ 64.;
	  MLT_( 3,0) = - 15./ 32.;
	  MLT_( 3,1) =  105./128.;
	  MLT_( 3,2) =   31./ 64.;
	  MLT_( 3,3) = -  9./128.;
	  MLT_( 4,0) =   15./ 64.;
	  MLT_( 4,1) = - 93./256.;
	  MLT_( 4,2) =  187./128.;
	  MLT_( 4,3) = - 67./256.;
	  MLT_( 5,1) = -  3./ 32.;
	  MLT_( 5,2) =   19./ 32.;
	  MLT_( 5,3) =   19./ 32.;
	  MLT_( 6,1) =    3./ 64.;
	  MLT_( 6,2) = -  1./  4.;
	  MLT_( 6,3) =   45./ 32.;
	  MLT_( 7,2) = -  3./ 32.;
	  MLT_( 7,3) =   19./ 32.;
	  MLT_( 8,2) =    3./ 64.;
	  MLT_( 8,3) = -  1./  4.;
	  MLT_( 9,3) = -  3./ 32.;
	  MLT_(10,3) =    3./ 64.;
	}
	break;
      case 6:
	j0_ = 3;
      default:
	break;
      }
      break;
    case 3:
      switch(dT) {
      case 3:
	j0_ = 3;
	if (s0 == 0) {
	  MLT_.resize(10,4);
	  MLT_(0,0) =    3./  2.;
	  MLT_(0,1) = - 35./ 96.;
	  MLT_(0,2) =   53./576.;
	  MLT_(0,3) = -  1./ 64.;
	  MLT_(1,0) =    1.;
	  MLT_(1,1) =   35./ 48.;
	  MLT_(1,2) = - 53./288.;
	  MLT_(1,3) =    1./ 32.;
	  MLT_(2,0) = -  3./  4.;
	  MLT_(2,1) =   71./ 32.;
	  MLT_(2,2) = - 53./192.;
	  MLT_(2,3) =    3./ 64.;
	  MLT_(3,0) =    1./  4.;
	  MLT_(3,1) = - 11./ 96.;
	  MLT_(3,2) =  689./576.;
	  MLT_(3,3) = - 13./ 64.;
	  MLT_(4,1) = - 45./ 64.;
	  MLT_(4,2) =  213./128.;
	  MLT_(4,3) = - 37./128.;
	  MLT_(5,1) =   15./ 64.;
	  MLT_(5,2) = - 39./128.;
	  MLT_(5,3) =  183./128.;
	  MLT_(6,2) = -  9./ 32.;
	  MLT_(6,3) =   45./ 32.;
	  MLT_(7,2) =    3./ 32.;
	  MLT_(7,3) = -  7./ 32.;
	  MLT_(8,3) = -  9./ 32.;
	  MLT_(9,3) =    3./ 32.;
	}
	break;
      case 5:
	j0_ = 4;
	if (s0 == 0) {
// 	  ML_.resize(16,6);
//  	  ML_( 0,0) = 1.0;
//  	  ML_( 1,0) = ML_( 1,1) = 0.5;
// 	  ML_( 2,1) = 0.75;
//  	  ML_( 2,2) = ML_( 3,1) = 0.25;
//  	  ML_( 3,2) = ML_( 4,2) = 0.75;
//  	  ML_( 4,3) = ML_( 5,2) = 0.25;
//  	  ML_( 5,3) = ML_( 6,3) = 0.75;
// 	  ML_( 6,4) = ML_( 7,3) = 0.25;
// 	  ML_( 7,4) = ML_( 8,4) = 0.75;
// 	  ML_( 8,5) = ML_( 9,4) = 0.25;
// 	  ML_( 9,5) = ML_(10,5) = 0.75;
// 	  ML_(11,5) = 0.25;

	  MLT_.resize(16,6);
	  MLT_( 0,0) =       5./     3.;
	  MLT_( 0,1) = -  2359./  3840.;
	  MLT_( 0,2) =  101909./345600.;
	  MLT_( 0,3) = - 17611./115200.;
	  MLT_( 0,4) =    3119./ 57600.;
	  MLT_( 0,5) = -    61./  6912.;
	  MLT_( 1,0) =       2./     3.;
	  MLT_( 1,1) =    2359./  1920.;
	  MLT_( 1,2) = -101909./172800.;
	  MLT_( 1,3) =   17611./ 57600.;
	  MLT_( 1,4) = -  3119./ 28800.;
	  MLT_( 1,5) =      61./  3456.;
	  MLT_( 2,0) = -     1./     2.;
	  MLT_( 2,1) =    2503./  1280.;
	  MLT_( 2,2) = - 18413./115200.;
	  MLT_( 2,3) =    1027./ 38400.;
	  MLT_( 2,4) = -    83./ 19200.;
	  MLT_( 2,5) =       1./  2304.;
	  MLT_( 3,0) =       1./     6.;
	  MLT_( 3,1) = -  1243./  3840.;
	  MLT_( 3,2) =  573353./345600.;
	  MLT_( 3,3) = - 79687./115200.;
	  MLT_( 3,4) =   13223./ 57600.;
	  MLT_( 3,5) = -   253./  6912.;
	  MLT_( 4,1) = -    77./   256.;
	  MLT_( 4,2) =    2583./  2560.;
	  MLT_( 4,3) =     789./  2560.;
	  MLT_( 4,4) = -   181./  1280.;
	  MLT_( 4,5) =      19./   768.;
	  MLT_( 5,1) = -    21./   256.; 
	  MLT_( 5,2) =     399./  2560.;
	  MLT_( 5,3) =    2877./  2560.;
	  MLT_( 5,4) = -   333./  1280.;
	  MLT_( 5,5) =       9./   256.;
	  MLT_( 6,1) =     105./   512.;
	  MLT_( 6,2) = -   547./  1024.;
	  MLT_( 6,3) =    1523./  1024.;
	  MLT_( 6,4) = -    79./   512.;
	  MLT_( 6,5) =      43./   512.;
	  MLT_( 7,1) = -    35./   512.;
	  MLT_( 7,2) =     129./  1024.;
	  MLT_( 7,3) = -   145./  1024.; 
	  MLT_( 7,4) =     709./   512.;
	  MLT_( 7,5) = -   587./  1536.;
	  MLT_( 8,2) =      15./   256.;
	  MLT_( 8,3) = -    97./   256.;
	  MLT_( 8,4) =     175./   128.;
	  MLT_( 8,5) = -    13./   128.;
	  MLT_( 9,2) = -     5./   256.;
	  MLT_( 9,3) =      19./   256.;
	  MLT_( 9,4) = -    13./   128.;
	  MLT_( 9,5) =     175./   128.;
	  MLT_(10,3) =      15./   256.;
	  MLT_(10,4) = -    97./   256.; 
	  MLT_(10,5) =     175./   128.;
	  MLT_(11,3) = -     5./   256.;
	  MLT_(11,4) =      19./   256.; 
	  MLT_(11,5) = -    13./   128.;
	  MLT_(12,4) =      15./   256.;
	  MLT_(12,5) = -    97./   256.;
	  MLT_(13,4) = -     5./   256.;
	  MLT_(13,5) =      19./   256.;
	  MLT_(14,5) =      15./   256.;
	  MLT_(15,5) = -     5./   256.;
	}
	break;
      case 7:
	j0_ = 4;
      default:
	break;
      }
      break;
    case 4:
      switch(dT) {
      case 6:
	j0_ = 4;
	break;
      case 8:
	j0_ = 5;
	break;
      default:
	break;
      }
      break;
    default:
      break;
    }

    cout << "ML=" << endl << ML_;
    cout << "MLT=" << endl << MLT_;
  }

}

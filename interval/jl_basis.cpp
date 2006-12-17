// implementation for jl_basis.h

#include <cmath>

namespace WaveletTL
{
  JLBasis::JLBasis()
  {
    setup();
  }
  
  void
  JLBasis::setup() {
    j0_ = 1;
  }

  inline
  JLBasis::Index
  JLBasis::first_generator(const int j) const
  {
    assert(j >= j0());
    return JLIndex(j, 0, 0, 1);
  }
  
  inline
  JLBasis::Index
  JLBasis::last_generator(const int j) const
  {
    assert(j >= j0());
    return JLIndex(j, 0, 1, 1<<j);
  }
  
  inline
  JLBasis::Index
  JLBasis::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return JLIndex(j, 1, 0, 1);
  }
  
  inline
  JLBasis::Index
  JLBasis::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return JLIndex(j, 1, 1, 1<<j);
  }
  
  void
  JLBasis::reconstruct(const InfiniteVector<double, Index>& c,
		       const int j,
		       InfiniteVector<double, Index>& v) const {
    v.clear();
    for (InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      v.add(*it, help);
    }
  }
  
  void
  JLBasis::reconstruct_1(const Index& lambda,
			 const int j,
			 InfiniteVector<double, Index>& c) const {
    c.clear();
    
    if (lambda.j() >= j)
      c.set_coefficient(lambda, 1.0); 
//     else {
//       const double phi_0_factor = sqrt(5./24.); // [JL], p.10
//       const double phi_1_factor = sqrt(15./8.); // dito
//       if (lambda.e() == 0) {
// 	// generator
	
// 	const int first_half = (1<<lambda.j())-get_s1();
// 	if (lambda.k() <= first_half) {
// 	  // type phi_0
// 	  // phi_0(x) = 1/2*phi_0(2*x+1)+phi_0(2*x)+1/2*phi_0(2*x-1)+3/4*phi_1(2*x+1)-3/4*phi_1(2*x-1)
// 	  const int offset = (1<<(lambda.j()+1))-get_s1()+1;

// 	  for (int m = 2*lambda.k()-1; m <= 2*lambda.k()+1; m++) {
// 	    switch (m-2*lambda.k()) {
// 	    case -1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(0.5*M_SQRT1_2, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(0.75*M_SQRT1_2*phi_0_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    case 0:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(M_SQRT1_2, dhelp);
// 	      }
// 	      break;
// 	    case 1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(0.5*M_SQRT1_2, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(-0.75*M_SQRT1_2*phi_0_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    default:
// 	      break;
// 	    }
// 	  }
// 	} else {
// 	  // type phi_1
// 	  // phi_1(x) = -1/8*phi_0(2*x+1)+1/8*phi_0(2*x-1)-1/8*phi_1(2*x+1)+1/2*phi_1(2*x)-1/8*phi_1(2*x-1)
	  
// 	  const int offset = (1<<(lambda.j()+1))-get_s1()+1;
// 	  const int lambdak = lambda.k()-first_half-1;
// 	  for (int m = 2*lambdak-1; m <= 2*lambdak+1; m++) {
// 	    switch (m-2*lambdak) {
// 	    case -1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(-0.125*M_SQRT1_2*phi_1_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(-0.125*M_SQRT1_2, dhelp);
// 	      }
// 	      break;
// 	    case 0:
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(0.5*M_SQRT1_2, dhelp);
// 	      }
// 	      break;
// 	    case 1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(0.125*M_SQRT1_2*phi_1_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(-0.125*M_SQRT1_2, dhelp);
// 	      }
// 	      break;
// 	    default:
// 	      break;
// 	    }
// 	  }
// 	}
//       } else {
// 	// wavelet
// 	const int first_half = (1<<lambda.j())-1;
// 	if (lambda.k() <= first_half) {
// 	  // type psi_0
// 	  // psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
// 	  const int offset = (1<<(lambda.j()+1))-get_s1()+1;
// 	  const double psi_0_factor = 1./sqrt(729.6);
	  
// 	  for (int m = 2*lambda.k()-1; m <= 2*lambda.k()+1; m++) {
// 	    switch (m-2*lambda.k()) {
// 	    case -1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(-2.0*M_SQRT1_2*psi_0_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(-21.0*M_SQRT1_2*psi_0_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    case 0:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(4.0*M_SQRT1_2*psi_0_factor/phi_0_factor, dhelp);
// 	      }
// 	      break;
// 	    case 1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(-2.0*M_SQRT1_2*psi_0_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(21.0*M_SQRT1_2*psi_0_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    default:
// 	      break;
// 	    }
// 	  }
// 	} else {
// 	  // type psi_1
// 	  // psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
// 	  const int offset = (1<<(lambda.j()+1))-get_s1()+1;
// 	  const int lambdak = lambda.k()-first_half-1;
// 	  const double psi_1_factor = (lambdak==0||lambdak==1<<lambda.j()) ? 1./sqrt(76.8) : 1./sqrt(153.6);
	  
// 	  for (int m = 2*lambdak-1; m <= 2*lambdak+1; m++) {
// 	    switch (m-2*lambdak) {
// 	    case -1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(M_SQRT1_2*psi_1_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(9.0*M_SQRT1_2*psi_1_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    case 0:
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(12.0*M_SQRT1_2*psi_1_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    case 1:
// 	      if (m >= get_s0() && m < offset) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, m, this), j, dhelp);
// 		c.add(-M_SQRT1_2*psi_1_factor/phi_0_factor, dhelp);
// 	      }
// 	      if (m >= 0 && m+offset <= DeltaRmax(lambda.j()+1)) {
// 		InfiniteVector<double, Index> dhelp;
// 		reconstruct_1(Index(lambda.j()+1, 0, offset+m, this), j, dhelp);
// 		c.add(9.0*M_SQRT1_2*psi_1_factor/phi_1_factor, dhelp);
// 	      }
// 	      break;
// 	    default:
// 	      break;
// 	    }
// 	  }
// 	}
//       }
//     }
}
  
}

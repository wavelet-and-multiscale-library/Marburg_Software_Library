// implementation for jl_basis.h

#include <cmath>
#include <interval/jl_utils.h>

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
    else {
#if _JL_PRECOND==0
      const double phi0factor = sqrt(35./26.);
      const double phi1factor = sqrt(105./2.);
#else
      const double phi0factor = sqrt(5./12.);
      const double phi1factor = sqrt(15./4.);
#endif
      if (lambda.e() == 0) {
 	// generator
	if (lambda.c() == 0) {
 	  // type phi_0
 	  // phi_0(x) = 1/2*phi_0(2*x+1)+phi_0(2*x)+1/2*phi_0(2*x-1)+3/4*phi_1(2*x+1)-3/4*phi_1(2*x-1)

	  int m = 2*lambda.k()-1; // m-2k=-1
	  { // phi_0(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(0.5*M_SQRT1_2, dhelp);
	  }
	  if (m >= 0) { // phi_1(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(0.75*M_SQRT1_2 * phi0factor/phi1factor, dhelp);
	  }

	  // m=2k <-> m-2k=0
	  m++;
	  { // phi_0(2x)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(M_SQRT1_2, dhelp);
	  }
	  
	  // m=2k+1 <-> m-2k=1
	  m++;
	  { // phi_0(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(0.5*M_SQRT1_2, dhelp);
	  }
	  if (m <= (1<<(lambda.j()+1))-1) { // phi_1(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(-0.75*M_SQRT1_2 * phi0factor/phi1factor, dhelp);
	  }
	} else { // lambda.c() == 1
  	  // type phi_1
 	  // phi_1(x) = -1/8*phi_0(2*x+1)+1/8*phi_0(2*x-1)-1/8*phi_1(2*x+1)+1/2*phi_1(2*x)-1/8*phi_1(2*x-1)
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  { // phi_0(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(-0.125*M_SQRT1_2 * phi1factor/phi0factor, dhelp);
	  }
	  if (m >= 0) { // phi_1(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(-0.125*M_SQRT1_2, dhelp);
	  }

	  // m=2k <-> m-2k=0
	  m++;
	  { // phi_1(2x)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(0.5*M_SQRT1_2, dhelp);
	  }

	  // m=2k+1 <-> m-2k=1
	  m++;
	  { // phi_0(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(0.125*M_SQRT1_2 * phi1factor/phi0factor, dhelp);
	  }
	  if (m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(-0.125*M_SQRT1_2, dhelp);
	  }
	} // end if (lambda.c() == 0)
      } else { // lambda.e() == 1
	// wavelet
	if (lambda.c() == 0) {
 	  // type psi_0
 	  // psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
#if _JL_PRECOND==0
	  const double factor = sqrt(35./352.);
#else
	  const double factor = sqrt(5./3648.);
#endif

	  int m = 2*lambda.k()-1; // m-2k=-1
	  { // phi_0(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(-2.0*M_SQRT1_2 * factor/phi0factor, dhelp);
	  }
	  if (m >= 0) { // phi_1(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(-21.0*M_SQRT1_2 * factor/phi1factor, dhelp);
	  }

	  // m=2k <-> m-2k=0
	  m++;
	  { // phi_0(2x)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(4.0*M_SQRT1_2 * factor/phi0factor, dhelp);
	  }

	  // m=2k+1 <-> m-2k=1
	  m++;
	  { // phi_0(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(-2.0*M_SQRT1_2 * factor/phi0factor, dhelp);
	  }
	  if (m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(21.0*M_SQRT1_2 * factor/phi1factor, dhelp);
	  }
	} else { // lambda.c() == 1
 	  // type psi_1
 	  // psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
#if _JL_PRECOND==0
	  const double factor = sqrt(35./48.);
#else
	  const double factor = sqrt(5./768.);
#endif
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  { // phi_0(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(M_SQRT1_2 * factor/phi0factor, dhelp);
	  }
	  if (m >= 0) { // phi_1(2x+1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(9.0*M_SQRT1_2 * factor/phi1factor, dhelp);
	  }

	  // m=2k <-> m-2k=0
	  m++;
	  { // phi_1(2x)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(12.0*M_SQRT1_2 * factor/phi1factor, dhelp);
	  }

	  // m=2k+1 <-> m-2k=1
	  m++;
	  { // phi_0(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
	    c.add(-M_SQRT1_2 * factor/phi0factor, dhelp);
	  }
	  if (m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
	    InfiniteVector<double, Index> dhelp;
	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
	    c.add(9.0*M_SQRT1_2 * factor/phi1factor, dhelp);
	  }
	}
      }
    }
  }
 
  void
  JLBasis::setup_full_collection()
  {
    if (jmax_ == -1 || jmax_ < j0_) {
      cout << "JLBasis::setup_full_collection(): specify a mximal level of resolution first!" << endl;
      abort();
    }   

    int degrees_of_freedom = Deltasize(jmax_+1);
    cout << "total degrees of freedom between j0_ and jmax_ is " << degrees_of_freedom << endl;

    cout << "setting up collection of wavelet indices..." << endl;
    full_collection.resize(degrees_of_freedom);
    int k = 0;
    for (Index ind = first_generator(j0_); ind <= last_wavelet(jmax_); ++ind) {
      full_collection[k] = ind;
      k++;
    }
    cout << "done setting up collection of wavelet indices..." << endl;

  }

 
}

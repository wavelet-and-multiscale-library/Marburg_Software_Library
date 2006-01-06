// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_BASIS_H
#define _WAVELETTL_P_BASIS_H

#include <interval/i_index.h>

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval as introduced in [P].

    The primal generators are exactly those B-splines associated with the
    Schoenberg knot sequence

      t^j_{-d+1} = ... = t_0 = 0          (knot with multiplicity d at x=0)
      t^j_k = k * 2^{-j}, 1 <= k <= 2^j-1
      t^j_{2^j} = ... = t_{2^j+d-1} = 1   (knot with multiplicity d at x=1)

    i.e.

      phi_{j,k}(x) = (t^j_{k+d}-t^j_k)[t^j_k,...,t^j_{k+d}](t-x)^{d-1}_+

    with supp(phi_{j,k}(x) = [t^j_k, t^j_{k+d}].
    In other words, we have exactly

     d-1     left boundary splines  (k=-d+1,...,-1),
     2^j-d+1 inner splines          (k=0,...,2^j-d),
     d-1     right boundary splines (k=2^j-d+1,...,2^j-1)

    References:
    [P] Primbs:
        Biorthogonale Wavelet-Basen auf dem Intervall
	Dissertation, Univ. Duisburg-Essen, 2006
  */
  template <int d, int dT>
  class PBasis
  {
  public:
    /*!
      constructor
      
      You can specify the order of either the primal (s) or the dual (sT) boundary conditions at
      the left and right end of the interval [0,1]. Several combinations are possible:
      
      si=sTi=0  : no b.c.
      si=s>0=sTi: primal basis with b.c., dual basis without b.c.
      si=0<s=sTi: primal basis without b.c., dual basis with b.c.
      si=sTi>0  : b.c. for primal and dual basis, like [CTU] (not recommended, loss of approximation order)
    */
    PBasis(const int s0 = 0, const int s1 = 0, const int sT0 = 0, const int sT1 = 0);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef IntervalIndex<PBasis<d,dT> > Index;
    
    //! extremal generator indices
    inline const int DeltaLmin() const { return 1-d; }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1; }

    //! boundary indices in \nabla_j
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }


  protected:
    //! coarsest possible level
    int j0_;

    //! boundary condition orders at 0 and 1
    int s0, s1, sT0, sT1;

    //! general setup routine which is shared by the different constructors
    void setup();
  };
}

#include <interval/p_basis.cpp>

#endif

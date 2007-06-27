// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_FUNCTIONAL_H
#define _FRAMETL_FUNCTIONAL_H

#include <utils/function.h>
#include <aggregated_frame.h>

using MathTL::Function;


namespace FrameTL
{
  /*!
    basis class for representation of a functional with respect to a given
    wavelet frame basis
    The function to which the functional is to be applied is given as a
    frame (index).
    This class models functionals of the form v -> \int g(x) v(x) dx.
    Other (more complicated) functionals need to be implemented in an
    inherited class.
   */
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class Functional
  {
  public:
    /*!
      Constructor with a given function f and an instance of a frame
     */
    Functional(const Function<DIM_d>* g, const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame)
      : g_(g), frame_(frame)
    { }

    /*!
      Evaluate the functional with the wavelet with given index lambda
      Here, the functional \psi_\lambda -> \int f(x) \psi_\lambda(x) dx
      is evalueted.
     */
    virtual double evaluate(const typename AggregatedFrame<IBASIS, DIM_d, DIM_m>::Index& lambda) const;

    //! reading access to the frame
    const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame() const { return frame_; }

    //! reading access to the function g
    const Function<DIM_d>* g() const { return g_; }


  protected:
    //! the function representing the simple functional
    const Function<DIM_d>* g_;

    //! an instance of the wavelet basis
    const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame_;
  };
}

#include <functional.cpp>

#endif // _FRAMETL_FUNCTIONAL_H

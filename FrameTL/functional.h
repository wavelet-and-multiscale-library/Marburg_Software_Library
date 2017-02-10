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
    Base class for representation of a functional with respect to a given
    wavelet (frame) basis. The function to which the functional is to be applied is given as a
    (frame) index.
    This class models functionals of the form \f$v \mapsto \int g(x) v(x) dx\f$.
    Other (more complicated) functionals need to be implemented in an
    inherited class.
   */
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class Functional
  {
  public:
    /*!
      Constructor with a given function f and an instance of a frame.
     */
    Functional(const Function<DIM_d>* g, const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame)
      : g_(g), frame_(frame)
    { }

    /*
      Virtual destructor.
     */
    virtual ~Functional()
    { }

    /*!
      Evaluate the functional with the wavelet with given index lambda
      Here, the functional \f$\psi_\lambda \mapsto \int f(x) \psi_\lambda(x) dx\f$
      is evaluated.
     */
    virtual double evaluate(const typename AggregatedFrame<IBASIS, DIM_d, DIM_m>::Index& lambda) const;

    //! Reading access to the frame.
    const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame() const { return frame_; }

    //! Reading access to the function g.
    const Function<DIM_d>* g() const { return g_; }


  protected:
    //! The function representing the simple functional.
    const Function<DIM_d>* g_;

    //! An instance of the wavelet basis.
    const AggregatedFrame<IBASIS, DIM_d, DIM_m>* frame_;
  };
}

#include <functional.cpp>

#endif // _FRAMETL_FUNCTIONAL_H

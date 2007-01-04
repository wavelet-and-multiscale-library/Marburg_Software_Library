// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GOERTZEL_REINSCH_H
#define _MATHTL_GOERTZEL_REINSCH_H

#include <cmath>
#include <algebra/vector.h>
#include <geometry/point.h>
#include <utils/function.h>

namespace MathTL
{
  // point evaluation of sine (and cosine) series according to Goertsch
  // (adjoint summmation for p_{2k}=cos(k.), p_{2k+1}=sin(k.))
  //
  // Note that this is potentially unstable for x\in l\mathbb Z!
  template <class C>
  class Goertzel : public Function<1,C>
  {
  public:
    Goertzel(const Vector<C>& coeffs, const bool sine = true)
      : coeffs_(coeffs), sine_(sine)
    {
    }
    
    C value(const Point<1,C>&p, const unsigned int component = 0) const;
    
    void vector_value(const Point<1,C> &p, Vector<C>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }

  protected:
    Vector<C> coeffs_;
    bool sine_;
  };

  // stable point evaluation of sine (or cosine) series
  //
  // references:
  // * Stoer, Numerische Mathematik 1
  // * Deuflhard/Bornemann, Numerische Mathematik 1

  template <class C>
  class GoertzelReinsch : public Function<1,C>
  {
  public:
    GoertzelReinsch(const Vector<C>& coeffs, const bool sine = true)
      : coeffs_(coeffs), sine_(sine)
    {
    }
    
    C value(const Point<1,C>&p, const unsigned int component = 0) const;
    
    void vector_value(const Point<1,C> &p, Vector<C>& values) const {
      values.resize(1, false);
      values[0] = value(p);
    }

    void set_coefficients(const Vector<C>& coeffs)
    {
      coeffs_ = coeffs;
    }

  protected:
    Vector<C> coeffs_;
    bool sine_;
  };
}

// include implementation of inline functions
#include <numerics/goertzel_reinsch.cpp>

#endif

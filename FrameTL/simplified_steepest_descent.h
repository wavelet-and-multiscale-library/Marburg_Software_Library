// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_SIMPLIFIED_STEEPEST_DESCENT_H
#define _FRAMETL_SIMPLIFIED_STEEPEST_DESCENT_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{
  /*!
    Implementation of a simplified steepest descent solver as described in
    [DRWFS]
    Simplified means that the appearing constants are "guessed" and not
    estimated as in above publication. That is, we have a standard steepest
    descent algorithm with APPLY and RHS instead of the normal matrix vector
    operations. The main iteration step is
     u := u + alpha * (RHS[eta] - APPLY[eta, u])
    with the acceleration parameter
     alpha = res^T res / (res^T APPLY[eta, res])
    and the residual
     res = (RHS[eta] - APPLY[eta, u])

    [DRWFS] Stephan Dahlke, Thorsten Raasch, Manuel Werner, Massimo Fornasier and Rob Stevenson
            Adaptive frame methods for elliptic operator equations: the steepest descent approach
  */
  template <class PROBLEM>
  void simplified_steepest_descent_SOLVE(const PROBLEM& problem, const double epsilon,
                                         InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
                                         const int jmax);
}

#include <simplified_steepest_descent.cpp>

#endif // _FRAMETL_SIMPLIFIED_STEEPEST_DESCENT_H

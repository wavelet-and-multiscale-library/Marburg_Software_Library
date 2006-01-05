// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ATLAS_H
#define _MATHTL_ATLAS_H

#include <iostream>
#include <geometry/chart.h>
#include <utils/array1d.h>
#include <algebra/symmetric_matrix.h>

namespace MathTL
{
  /*!
    Base class for an overlapping or non-overlapping atlas
      A = (kappa_i)_{i=1,...,n}
    of a d-dimensional manifold Omega in R^m.
    Essentially, an atlas is a collecition of mappings, together
    with the topological information of the manifold given by a symmetric
    adjacency matrix.
   */
  template <unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class Atlas {
  public:
    //! default constructor, yields empty atlas
    Atlas() {}

    //! constructor from a number of charts and adjacency relations
    Atlas (const Array1D<Chart<DIM_d,DIM_m>* >&,
	   const SymmetricMatrix<bool>&);
    
    //! read access to the charts
    const Array1D<Chart<DIM_d,DIM_m>* >& charts() const { return charts_; }

    //! read access to the adjacency matrix
    const SymmetricMatrix<bool>& get_adjacency_matrix() const { return adjacency_matrix; }
    
  protected:
    //! pointers to the parametrizations
    Array1D<Chart<DIM_d,DIM_m>* > charts_;

    //! adjacency relations
    SymmetricMatrix<bool> adjacency_matrix;
  };

  /*!
    stream output of an atlas
  */
  template <unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream&, const Atlas<DIM_d, DIM_m>&);
}

#include "geometry/atlas.cpp"

#endif

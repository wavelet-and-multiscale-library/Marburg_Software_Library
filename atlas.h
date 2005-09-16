// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_ATLAS_H
#define _FRAMETL_ATLAS_H

#include <iostream>
#include <list>

using std::cout;
using std::endl;
using std::list;

using MathTL::SymmetricMatrix;


namespace FrameTL
{
  /*!
    This class provides a realization of an arbitrary bounded domain
    in \mathbb R^{DIM_d}, that is covered with (overlapping) patches, each of
    which is represented by a parametric mapping from \mathbb R^{DIM_m} onto
    \mathbb R^{DIM_d}.
  */
  template <unsigned int DIM_m, unsigned int DIM_d>
  class Atlas
  {

  public:
    /*!
      default constructor
     */
    Atlas ();

    /*!
      constructor setting number of patches
     */
    Atlas (const int);

    /*!
      copy constructor
     */
    Atlas (const Atlas &);

    /*!
      constructor
     */
    Atlas (const Array1D<Parametrization<DIM_m, DIM_d> *>& charts);

    /*!
      access to the collection of parametrizations
     */
    const Array1D<Parametrization<DIM_m, DIM_d> *>& get_charts() const;

    /*!
      access to parametrization of patch with number p.
     */
    const Parametrization<DIM_m, DIM_d>& get_chart(const unsigned int p) const;

    /*!
      add a new patch
     */
    void add_chart(Parametrization<DIM_m, DIM_d>& kappa);

    /*!
      sets parametrization of patch p by to kappa, overwrites a possible
      existing one for patch p
     */
    void set_chart(const int p, const Parametrization<DIM_m, DIM_d>& kappa);

    /*!
      access to total number of patches
     */
    const unsigned int get_num_patches() const;

    /*!
      declares patches p1 and p2 to be (overlapping) neighbors
     */
    void set_adjacent(const int p1, const int p2);

    /*!
      checks wether patches p1 and p2 are (overlapping) neighbors
    */
    const bool adjacent(const int p1, const int p2);

    /*!
      returns adjacency matrix
     */
    const SymmetricMatrix<bool>& get_adjacency_matrix() const;

    /*!
      collects those patches, p is contained in
     */
    const Array1D<unsigned int>& get_circumjacent (const Point<DIM_d>& p);

  protected:
    /*!
      number of patches
     */
    unsigned int num_patches;

    /*!
      the parametrizations of the patches
     */
    Array1D<Parametrization<DIM_m, DIM_d> *> charts;

  private:
    /*!
      a matrix representing the adjacency relation of the patches,
      (patch collection can be viewed as an undirected graph)
     */
    SymmetricMatrix<bool> adjacency_matrix;
  };

  /*!
    stream output for Atlas
   */
  template <unsigned int DIM_m, unsigned int DIM_d>
  std::ostream& operator << (std::ostream&, const Atlas<DIM_m, DIM_d>&);

}

// include implementation
#include "atlas.cpp"

#endif


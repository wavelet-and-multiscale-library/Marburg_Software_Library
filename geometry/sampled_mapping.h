// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SAMPLED_MAPPING_H
#define _MATHTL_SAMPLED_MAPPING_H

#include <iostream>

#include <geometry/point.h>
#include <geometry/grid.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>
#include <algebra/matrix.h>
#include <algebra/infinite_vector.h>

namespace MathTL
{
  /*!
    Abstract base class for a mapping from R^n to R,
    represented by finite many samples on a rectangular grid (Matlab style).
    Meaningful values for n are n=1 and n=2 so far (see specializations).
  */
  template <unsigned int DIM>
  class SampledMapping
    : public Grid<DIM>
  {
  public:
    /*!
      default constructor, yields empty mapping
    */
    SampledMapping() {}

    /*!
      copy constructor
    */
    SampledMapping(const SampledMapping<DIM>& sm);

    /*!
      constructor from a given grid and given values
    */
    SampledMapping(const Grid<DIM>& grid, const Array1D<double>& values);

    /*!
      constructor from a fixed grid and a Function object
    */
    SampledMapping(const Grid<DIM>& grid, const Function<DIM>& f);

    /*!
      constructor from given values on 2^{-resolution}\mathbb Z^d,
      clipped to the cuboid <a,b> in \mathbb Z^d
    */
    SampledMapping(const MultiIndex<int, DIM>& a,
		   const MultiIndex<int, DIM>& b,
		   const InfiniteVector<double, MultiIndex<int, DIM> >& values,
		   const int resolution);

    /*!
      constructor from given tensor product values on a uniform subgrid
      of the cuboid <a,b>,
      the number of grid points in the i-th direction is taken from the length
      of the i-th component array of 'values'
    */
    SampledMapping(const Point<DIM>& a,
		   const Point<DIM>& b,
		   const FixedArray1D<Array1D<double>,DIM>& values);

    /*!
      assignment operator
    */
    SampledMapping<DIM>& operator = (const SampledMapping<DIM>& sm);

    /*!
      pointwise in-place summation *this += s
      of two sampled mappings over the same grid
    */
    void add(const SampledMapping<DIM>& s);
    
    /*!
      pointwise in-place summation *this += alpha * s
      of two sampled mappings over the same grid
    */
    void add(const double alpha, const SampledMapping<DIM>& s);

    /*!
      Matlab output of the sampled mapping onto a stream
    */
    void matlab_output(std::ostream& os) const;
  };

  //
  //
  // template specializations to one and two space dimensions

  template <>
  class SampledMapping<1>
    : public Grid<1>
  {
  public:
    /*!
      default constructor, yields empty mapping
    */
    SampledMapping();

    /*!
      copy constructor
    */
    SampledMapping(const SampledMapping<1>& sm);

    /*!
      constructor from a given grid and given values
    */
    SampledMapping(const Grid<1>& grid, const Array1D<double>& values);

    /*!
      constructor from a fixed grid and a Function object
    */
    SampledMapping(const Grid<1>& grid, const Function<1>& f);

    /*!
      constructor from given values on 2^{-resolution}\mathbb Z, clipped to [a,b]
    */
    SampledMapping(const MultiIndex<int,1>& a,
		   const MultiIndex<int,1>& b,
		   const InfiniteVector<double, MultiIndex<int, 1> >& values,
		   const int resolution);

    /*!
      constructor from given tensor product values on a uniform subgrid
      of the cuboid <a,b>, a,b integer vectors,
      the number of grid points in the i-th direction is taken from the length
      of the i-th component array of 'values'
      (this constructor does not make much sense in 1 space dimension...)
    */
    SampledMapping(const Point<1>& a,
		   const Point<1>& b,
		   const FixedArray1D<Array1D<double>,1>& values);

    /*!
      assignment operator
    */
    SampledMapping<1>& operator = (const SampledMapping<1>& sm);

    /*!
      pointwise in-place summation *this += s
      of two sampled mappings over the same grid
    */
    void add(const SampledMapping<1>& s);
    
    /*!
      pointwise in-place summation *this += alpha * s
      of two sampled mappings over the same grid
    */
    void add(const double alpha, const SampledMapping<1>& s);

    /*!
      reading access to the function values
    */
    inline const Array1D<double>& values() const { return values_; }

    /*!
      Matlab output of the sampled mapping onto a stream
    */
    void matlab_output(std::ostream& os) const;

  protected:
    /*!
      internal storage for the function values
    */
    Array1D<double> values_;
  };

  template <>
  class SampledMapping<2>
    : public Grid<2>
  {
  public:
    /*!
      default constructor, yields empty mapping
    */
    SampledMapping();

    /*!
      copy constructor
    */
    SampledMapping(const SampledMapping<2>& sm);

    /*!
      constructor from a given grid and given values
    */
    SampledMapping(const Grid<2>& grid, const Matrix<double>& values);

    /*!
      constructor from a fixed grid and a Function object
    */
    SampledMapping(const Grid<2>& grid, const Function<2>& f);

    /*!
      constructor from given values on 2^{-resolution}\mathbb Z^2,
      clipped to [a,b]^2
      (CURRENTLY NOT IMPLEMENTED!)
    */
    SampledMapping(const MultiIndex<int, 2>& a,
		   const MultiIndex<int, 2>& b,
		   const InfiniteVector<double, MultiIndex<int, 2> >& values,
		   const int resolution);

    /*!
      constructor from given tensor product values on a uniform subgrid
      of the cuboid <a,b>, a,b integer vectors,
      the number of grid points in the i-th direction is taken from the length
      of the i-th component array of 'values'
    */
    SampledMapping(const Point<2>& a,
		   const Point<2>& b,
		   const FixedArray1D<Array1D<double>,2>& values);

    /*!
      assignment operator
    */
    SampledMapping<2>& operator = (const SampledMapping<2>& sm);

    /*!
      pointwise in-place summation *this += s
      of two sampled mappings over the same grid
    */
    void add(const SampledMapping<2>& s);
    
    /*!
      pointwise in-place summation *this += alpha * s
      of two sampled mappings over the same grid
    */
    void add(const double alpha, const SampledMapping<2>& s);

    /*!
      reading access to the function values
    */
    inline const Matrix<double>& values() const { return values_; }

    /*!
      Matlab output of the sampled mapping onto a stream
    */
    void matlab_output(std::ostream& os) const;

  protected:
    /*!
      internal storage for the function values
    */
    Matrix<double> values_;
  };
}

#include "geometry/sampled_mapping.cpp"

#endif

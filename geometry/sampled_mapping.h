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

#include <geometry/grid.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <utils/multiindex.h>
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
      assignment operator
    */
    SampledMapping<DIM>& operator = (const SampledMapping<DIM>& sm);

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
    SampledMapping() : Grid<1>(), values_() {}

    /*!
      copy constructor
    */
    SampledMapping(const SampledMapping<1>& sm)
      : Grid<1>(sm), values_(sm.values_)
    {
    }

    /*!
      constructor from a given grid and given values
    */
    SampledMapping(const Grid<1>& grid, const Array1D<double>& values)
      : Grid<1>(grid), values_(values)
    {
    }

    /*!
      constructor from a fixed grid and a Function object
    */
    SampledMapping(const Grid<1>& grid, const Function<1>& f)
      : Grid<1>(grid), values_(grid.size())
    {
      for (unsigned int n(0); n < grid.size(); n++)
	values_[n] = f.value(Point<1>(grid_[n]));
    }

    /*!
      constructor from given values on 2^{-resolution}\mathbb Z, clipped to [a,b]
    */
    SampledMapping(const MultiIndex<int, 1>& a,
		   const MultiIndex<int, 1>& b,
		   const InfiniteVector<double, MultiIndex<int, 1> >& values,
		   const int resolution)
      : Grid<1>(a[0], b[0], (1<<resolution)*(b[0]-a[0]))
    {
      values_.resize(Grid<1>::size());
      for (int k(a[0]<<resolution), n(0); k <= (b[0]<<resolution); k++, n++)
	values_[n] = values.get_coefficient(MultiIndex<int, 1>(k));
    }

    /*!
      assignment operator
    */
    SampledMapping<1>& operator = (const SampledMapping<1>& sm)
    {
      Grid<1>::operator = (sm);
      values_ = sm.values_;
      return *this;
    }

    /*!
      reading access to the function values
    */
    inline const Array1D<double>& values() const { return values_; }

    /*!
      Matlab output of the sampled mapping onto a stream
    */
    void matlab_output(std::ostream& os) const
    {
      Grid<1>::matlab_output(os);
      os << "y = " // here we can take y
	 << values_
	 << ";"
	 << std::endl;
    }

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
    SampledMapping() : Grid<2>(), values_() {}

    /*!
      copy constructor
    */
    SampledMapping(const SampledMapping<2>& sm)
      : Grid<2>(sm), values_(sm.values_)
    {
    }

    /*!
      constructor from a given grid and given values
    */
    SampledMapping(const Grid<2>& grid, const Matrix<double>& values)
      : Grid<2>(grid), values_(values)
    {
    }

    /*!
      constructor from a fixed grid and a Function object
    */
    SampledMapping(const Grid<2>& grid, const Function<2>& f)
      : Grid<2>(grid), values_()
    {
      values_.resize(gridx_.row_dimension(), gridx_.column_dimension());
      for (unsigned int m(0); m < values_.row_dimension(); m++)
	for (unsigned int n(0); n < values_.column_dimension(); n++)
	  values_(m,n) = f.value(Point<2>(gridx_(m,n), gridy_(m,n)));
    }

    /*!
      constructor from given values on 2^{-resolution}\mathbb Z^2,
      clipped to [a,b]^2
    */
    SampledMapping(const MultiIndex<int, 2>& a,
		   const MultiIndex<int, 2>& b,
		   const InfiniteVector<double, MultiIndex<int, 2> >& values,
		   const int resolution)
    {
    }

    /*!
      assignment operator
    */
    SampledMapping<2>& operator = (const SampledMapping<2>& sm)
    {
      Grid<2>::operator = (sm);
      values_ = sm.values_;
      return *this;
    }

    /*!
      Matlab output of the sampled mapping onto a stream
    */
    void matlab_output(std::ostream& os) const
    {
      Grid<2>::matlab_output(os);
      os << "z = ";
      print_matrix(values_, os);
      os << ";"
	 << std::endl;
    }

  protected:
    /*!
      internal storage for the function values
    */
    Matrix<double> values_;
  };
}
#endif

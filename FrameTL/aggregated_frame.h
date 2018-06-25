// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_AGGREGATED_FRAME_H
#define _FRAMETL_AGGREGATED_FRAME_H

#include <list>

#include <geometry/chart.h>
#include <cube/mapped_cube_basis.h>
#include <geometry/atlas.h>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <frame_index.h>
#include <frame_support.h>

using std::list;

using WaveletTL::MappedCubeBasis;
using MathTL::Atlas;
using MathTL::Array1D;
using MathTL::FixedArray1D;

//using namespace WaveletTL;

namespace FrameTL
{
  //forward declaration of class AggregatedFrame
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class FrameIndex;

  /*!
    Class for frames on \f$d\f$-dimensional manifolds in \f$R^m\f$.
    The construction principle is the following.
    We write the domain or manifold as an overlapping union
    of subdomains (patches), each of them being the smooth
    parametric image of a reference domain, i.e., the
    \f$d\f$-dimensional hypercube. By lifting a wavelet basis
    on the reference domain to the subdomains and taking
    the union of these lifted bases, a frame is obtained.

    The manifold is given by an appropriate Atlas. The
    corresponding reference bases, or their
    lifted versions, are then internally constructed.
    For each lifted cube, the user may specify
    \f$2d\f$ boundary conditions.

    @tparam IBASIS The underlying one-dimensional wavelet basis on the interval.
    @tparam DIM_d Dimension of the domain or manifold.
    @tparam DIM_m Dimension of the surrounding space.

  */
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class AggregatedFrame
  {
  public:

    /*!
      Destructor.
     */
    virtual ~AggregatedFrame();

    /*!
      Constructor, creates a frame from a given Atlas, sets primal and dual boundary conditions,
      and fixes a maximal level of resolution that will be regarded.
      The fixation of a maximal level is included only for performance reasons.
      A lot of information about the wavelets between the minimal and maximal level
      is preprocessed. The maximal level should be chosen large enough.
      In a future iteration of FrameTL, we should finally get rid of this.

      @param Atlas The underlying Atlas of parametric mappings representing the domain covering.
      @param bc One-dimensional array of one-dimensional arrays of fixed length 2*DIM_d.
      The outer array runs over the subdomains. The inner one runs over the wavelet bases in each spatial
      direction. For each of the latter, the boundary conditions at the left and at the right end of the interval
      are specified. A boundary condition of the value r means that all derivatives up to degree r-1 vanish.
      @param bcT The same as <em>bc</em> for the dual side.
      @param jmax The maximal level of resolution.
    */
    AggregatedFrame(const Atlas<DIM_d, DIM_m>* a,
		    const Array1D<FixedArray1D<int,2*DIM_d> >& bc,
		    const Array1D<FixedArray1D<int,2*DIM_d> >& bcT,
		    const int jmax);

    /*!
      Constructor. Does the same as the first constructor, but no dual boundary
      conditions are specified.
     */
    AggregatedFrame(const Atlas<DIM_d, DIM_m>*,
		    const Array1D<FixedArray1D<int,2*DIM_d> >&,
		    const int jmax);

    /*!
      Frame index class.
    */
    typedef FrameIndex<IBASIS, DIM_d, DIM_m> Index;


    /*!
      Interval basis.
     */
    typedef IBASIS IntervalBasis;


    /*!
      Read only access to the local bases.
     */
    const Array1D<MappedCubeBasis<IBASIS, DIM_d, DIM_m>* > bases() const
    {
      return lifted_bases;
    }

//     /*!
//       geometric type of the support sets.
//       support is only implicitly given by the
//       support of the corresponding wavelet on the
//       hypercube and the respective patch chart.
//     */
//     typedef struct {
//       int j;
//       int a[DIM_d];
//       int b[DIM_d];
//       const Chart<DIM_d,DIM_m>* ch;
//     } Support;

    /*!
      This structure represents a rectangular support of a wavelet
      given by the lower left and upper right corner of the rectangle.
      This is tailored to the special case of rectangular subdomains.
      In the general case, full support information is given by the support of
      the reference wavelet and the respective parametric mapping (see the commented out
      code above).
    */
    typedef struct {
      int j;
      Point<DIM_d> a;
      Point<DIM_d> b;
    } Support;


    /*!
      Coarsest level.
     */
    const int j0() const { return j0_; }

    /*!
      Finest level.
     */
    const int jmax() const { return jmax_; }


    /*!
      Number of patches.
     */
    const int n_p() const { return lifted_bases.size(); }

    /*!
      Read only access to underlying atlas.
     */
    const Atlas<DIM_d, DIM_m>* atlas() const { return atlas_;  }

    /*!
      Read only access to collection of wavelet frame indices.
      In the constructors, between minimal and maximal level,
      on each level an array of all wavelet indices is preprocessed.
    */
    const Array1D<Array1D<Index> >* get_full_collection_levelwise() const { return &full_collection_levelwise; }

    /*!
      Read only access to collection of wavelet frame indices.
      In the constructors, between minimal and maximal level,
      an array of all wavelet indices is preprocessed.
     */
    const Array1D<Index>* get_full_collection() const { return &full_collection; }

    /*!
      Critical Sobolev regularity for the primal generators/wavelets.
    */
    static double primal_regularity() { return IBASIS::primal_regularity(); }

    /*!
      Number of vanishing moments of the primal wavelets.
    */
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    /*!
      Number of vanishing moments of the dual wavelets.
    */
    static unsigned int dual_vanishing_moments() { return IBASIS::dual_vanishing_moments(); }


    /*!
      Given the number of a frame element, according to the canonical lexicographical ordering
      (level, type, patch, shift parameter), the corresponding wavelet index is returned.
      This index is taken from the preprocessed collection of all wavelet indices between
      minimal and maximal level.
    */
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

    /*!
      Set finest possible level.
    */
    void set_jmax(const int jmax) { jmax_ = jmax; }

    //! Total number of wavelets between coarsest and finest level.
    const int degrees_of_freedom() const { return full_collection.size(); };

    /*!
      Point evaluation of a single frame element.
    */
    double evaluate(const Index& lambda, const Point<DIM_m>& x) const;

    /*!
      Point evaluation of the gradient of a single frame element.
    */
    void evaluate_gradient(const Index& lambda, const Point<DIM_d>& x, Vector<double>& values) const;


    //! Index of first generator on level \f$j \geq j_0\f$.
    Index first_generator(const int j) const;

    //! Index of last generator on level \f$j \geq j_0\f$.
    Index last_generator(const int j) const;

    //! Index of first wavelet on level \f$j \geq j_0\f$.
    Index first_wavelet(const int j) const;

    //! Index of last wavelet on level \f$j \geq j_0\f$.
    Index last_wavelet(const int j) const;


    /*!
      All supports of the reference wavelets on the unit cube,
      preprocessed between minimal and maximal level.
    */

    Array1D<typename WaveletTL::CubeBasis<IBASIS,DIM_d>::Support> all_supports;

    /*!
      All supports of the frame elements,
      preprocessed between minimal and maximal level.
    */
    Array1D<Support> all_patch_supports;


  protected:
    //! Pointer to the underlying atlas.
    const Atlas<DIM_d, DIM_m>* atlas_;

    //! Primal boundary conditions.
    Array1D<FixedArray1D<int,2*DIM_d> > bc_;

    //! Dual boundary conditions.
    Array1D<FixedArray1D<int,2*DIM_d> > bcT_;

    /*!
      Collection of all wavelet indices between coarsest and finest level
      stored level by level.
    */

    Array1D<Array1D<Index> > full_collection_levelwise;

    //! Collection of all wavelet indices between coarsest and finest level.
    Array1D<Index> full_collection;

  private:
    /*!
      Collection of mapped cube bases forming the aggregated frame.
     */
    Array1D<MappedCubeBasis<IBASIS, DIM_d, DIM_m>* > lifted_bases;


    //! Coarsest possible level \f$j_0\f$.
    int j0_;

    //! Finest possible level.
    int jmax_;


  };

  /*!
    stream output of an AggregatedFrame
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream&,
			     const AggregatedFrame<IBASIS, DIM_d, DIM_m>&);

}

#include "aggregated_frame.cpp"

#endif

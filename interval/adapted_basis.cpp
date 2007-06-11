// implementation for adapted_basis.h

#include <cassert>
#include <algebra/infinite_vector.h>

using namespace MathTL;

namespace WaveletTL
{
  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::encode_infinite_vector(const InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex>& c_multi, InfiniteVector<double, typename AdaptedBasis<IBASIS>::Index>& c) const
  {
    c.clear();
    for (typename InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex>::const_iterator it(c_multi.begin()), itend(c_multi.end()); it != itend; ++it)
      c.set_coefficient(typename AdaptedBasis<IBASIS>::Index(it.index(), this), *it);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::decode_infinite_vector(const InfiniteVector<double, Index>& c_adapt, InfiniteVector<double, MultiIndex>& c_multi) const
  {
    c_multi.clear();
    for (typename InfiniteVector<double, typename AdaptedBasis<IBASIS>::Index>::const_iterator it(c_adapt.begin()), itend(c_adapt.end()); it != itend; ++it)
      c_multi.set_coefficient(typename AdaptedBasis<IBASIS>::MultiIndex(*it.index().multi_index()), *it);
  }


  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::decompose_1(const Index& lambda, const int jmin, InfiniteVector<double, Index>& c) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi;
    basis_->decompose_1(typename AdaptedBasis<IBASIS>::Index(lambda, this), jmin, c_multi);
    encode_infinite_vector(c_multi, c);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::decompose(const InfiniteVector<double, Index>& c, const int jmin, InfiniteVector<double, Index>& v) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi, v_multi;
    decode_infinite_vector(c, c_multi);
    basis_->decompose(c_multi, jmin, v_multi);
    encode_infinite_vector(v_multi, v);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::decompose_t_1(const Index& lambda, const int jmin, InfiniteVector<double, Index>& c) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi;
    basis_->decompose_t_1(typename AdaptedBasis<IBASIS>::Index(lambda, this), jmin, c_multi);
    encode_infinite_vector(c_multi, c);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::decompose_t(const InfiniteVector<double, Index>& c, const int jmin, InfiniteVector<double, Index>& v) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi, v_multi;
    decode_infinite_vector(c, c_multi);
    basis_->decompose_t(c_multi, jmin, v_multi);
    encode_infinite_vector(v_multi, v);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::reconstruct_1(const Index& lambda, const int j, InfiniteVector<double, Index>& c) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi;
    basis_->reconstruct_1(typename AdaptedBasis<IBASIS>::Index(lambda, this), j, c_multi);
    encode_infinite_vector(c_multi, c);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::reconstruct(const InfiniteVector<double, Index>& c, const int j, InfiniteVector<double, Index>& v) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi, v_multi;
    decode_infinite_vector(c, c_multi);
    basis_->reconstruct(c_multi, j, v_multi);
    encode_infinite_vector(v_multi, v);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::reconstruct_t_1(const Index& lambda, const int j, InfiniteVector<double, Index>& c) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi;
    basis_->reconstruct_t_1(typename AdaptedBasis<IBASIS>::Index(lambda, this), j, c_multi);
    encode_infinite_vector(c_multi, c);
  }

  template <class IBASIS>
  void
  AdaptedBasis<IBASIS>::reconstruct_t(const InfiniteVector<double, Index>& c, const int j, InfiniteVector<double, Index>& v) const
  {
    InfiniteVector<double, typename AdaptedBasis<IBASIS>::MultiIndex> c_multi, v_multi;
    decode_infinite_vector(c, c_multi);
    basis_->reconstruct_t(c_multi, j, v_multi);
    encode_infinite_vector(v_multi, v);
  }
}

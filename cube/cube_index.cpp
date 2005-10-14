// implementation for cube_index.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>::CubeIndex(const CubeBasis<IBASIS,DIM>* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = 0; // invalid (e and k are initialized by zero automatically)
    } else {
      j_ = basis_->j0(); // coarsest level;
      // e_ is zero by default: generator
      for (unsigned int i = 0; i < DIM; i++)
	k_[i] = basis_->bases()[i]->DeltaLmin();
    }
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>::CubeIndex(const int j,
				   const type_type& e,
				   const translation_type& k,
				   const CubeBasis<IBASIS,DIM>* basis)
    : basis_(basis), j_(j), e_(e), k_(k)
  {
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>::CubeIndex(const CubeIndex& lambda)
    : basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_)
  {
  }

  template <class IBASIS, unsigned int DIM>
  bool
  CubeIndex<IBASIS,DIM>::operator == (const CubeIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IBASIS, unsigned int DIM>
  bool
  CubeIndex<IBASIS,DIM>::operator < (const CubeIndex& lambda) const
  {
    // standard lexicographic order on (j,e,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>
  first_generator(const CubeBasis<IBASIS,DIM>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM>::type_type e;
    typename CubeIndex<IBASIS,DIM>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(first_generator(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>
  last_generator(const CubeBasis<IBASIS,DIM>* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM>::type_type e;
    typename CubeIndex<IBASIS,DIM>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_generator(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>
  first_wavelet(const CubeBasis<IBASIS,DIM>* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM>::type_type e;
    typename CubeIndex<IBASIS,DIM>::translation_type k;
    for (unsigned int i = 0; i < DIM-1; i++) {
      typename IBASIS::Index lambda(first_generator(basis->bases()[i], j));
      k[i] = lambda.k();
    }
    typename IBASIS::Index lambda(first_wavelet(basis->bases()[DIM-1], j));
    k[DIM-1] = lambda.k();
    e[DIM-1] = 1;

    return CubeIndex<IBASIS,DIM>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM>
  CubeIndex<IBASIS,DIM>
  last_wavelet(const CubeBasis<IBASIS,DIM>* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM>::type_type e;
    typename CubeIndex<IBASIS,DIM>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_wavelet(basis->bases()[i], j));
      k[i] = lambda.k();
      e[i] = 1;
    }

    return CubeIndex<IBASIS,DIM>(j, e, k, basis);
  }
}

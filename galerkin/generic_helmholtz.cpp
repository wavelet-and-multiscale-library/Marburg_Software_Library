// implementation for generic_helmholtz.h

namespace WaveletTL
{
  template <class WBASIS>
  GenericFullHelmholtz<WBASIS>::GenericFullHelmholtz
  (const WBASIS& basis,
   const char* filenameG,
   const char* filenameA,
   const int jmax,
   const double alpha,
   const PreconditioningType precond)
    : basis_(basis), jmax_(jmax), precond_(precond)
  {
    G_.matlab_input(filenameG);
    A_.matlab_input(filenameA);
    set_alpha(alpha);
  }
  
  template <class WBASIS>
  inline
  const
  typename GenericFullHelmholtz<WBASIS>::size_type
  GenericFullHelmholtz<WBASIS>::row_dimension() const
  {
    return WBASIS::Deltasize(jmax_+1);
  }
  
  template <class WBASIS>
  inline const
  typename GenericFullHelmholtz<WBASIS>::size_type
  GenericFullHelmholtz<WBASIS>::column_dimension() const
  {
    return WBASIS::Deltasize(jmax_+1);
  }

  template <class WBASIS>
  void
  GenericFullHelmholtz<WBASIS>::set_alpha(const double alpha) const
  {
    assert(alpha >= 0);
    alpha_ = alpha;
    setup_D();
  }

  template <class WBASIS>
  void
  GenericFullHelmholtz<WBASIS>::setup_D() const
  {
    D_.resize(WBASIS::Deltasize(jmax_+1));
    if (precond_ == no_precond) {
      D_ = 1.0;
    } else {
      if (precond_ == dyadic) {
	for (int k(0); k < WBASIS::Deltasize(WBASIS::j0()); k++)
	  D_[k] = alpha_ + (1<<WBASIS::j0());
	for (int j = WBASIS::j0(); j < jmax_; j++) {
	  for (int k(WBASIS::Deltasize(j)); k < WBASIS::Deltasize(j+1); k++)
	    D_[k] = alpha_ + (1<<j);
	}
      } else {
 	for (size_type k(0); k < D_.size(); k++)
 	  D_[k] = sqrt(diagonal(k));
      }
    }
  }
  
  template <class WBASIS>
  inline
  const double
  GenericFullHelmholtz<WBASIS>::diagonal(const size_type row) const
  {
    return (alpha_*G_.get_entry(row,row)+A_.get_entry(row,row));
  }
  
  template <class WBASIS>
  inline
  double
  GenericFullHelmholtz<WBASIS>::D(const size_type k) const {
    return D_[k];
  }

  template <class WBASIS>
  template <class VECTOR>
  void
  GenericFullHelmholtz<WBASIS>::apply(const VECTOR& x, VECTOR& Mx,
				      const bool preconditioning) const
  {
    assert(Mx.size() == row_dimension());
    
    VECTOR y(x), yhelp(x.size(), false);
    
    if (preconditioning) {
      // apply diagonal preconditioner D^{-1}
      for (size_type k(0); k < y.size(); k++)
	y[k] /= D(k);
    }
    
    // apply Gramian
    G_.apply(y, yhelp);
    
    // apply unpreconditioned Laplacian
    A_.apply(y, Mx);
    
    // add Gramian part, y+=alpha*yhelp
    Mx.add(alpha_, yhelp);
    
    if (preconditioning) {
      // apply diagonal preconditioner D^{-1}
      for (size_type k(0); k < y.size(); k++)
	Mx[k] /= D(k);
    }
  }
  
}

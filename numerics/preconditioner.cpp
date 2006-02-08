// implementation for preconditioner.h

namespace MathTL
{
  template <class VECTOR>
  Preconditioner<VECTOR>::~Preconditioner()
  {
  }

  template <class MATRIX, class VECTOR>
  IdentityPreconditioner<MATRIX,VECTOR>::IdentityPreconditioner(const MATRIX& M)
    : A(M)
  {
  }

  template <class MATRIX, class VECTOR>
  inline
  void
  IdentityPreconditioner<MATRIX,VECTOR>::apply(const VECTOR& x, VECTOR& Px) const
  {
    Px = x;
  }
  
  template <class MATRIX, class VECTOR>
  void
  IdentityPreconditioner<MATRIX,VECTOR>::apply_preconditioner(const VECTOR& Px, VECTOR& x) const
  {
    x = Px;
  }
  
  template <class MATRIX, class VECTOR>
  JacobiPreconditioner<MATRIX,VECTOR>::JacobiPreconditioner(const MATRIX& M)
    : A(M)
  {
  }

  template <class MATRIX, class VECTOR>
  void
  JacobiPreconditioner<MATRIX,VECTOR>::apply(const VECTOR& x, VECTOR& Px) const
  {
    Px = x;
    unsigned int i(0);
    for (typename VECTOR::iterator it(Px.begin()), itend(Px.end());
	 it != itend; ++it, ++i)
      *it *= A.get_entry(i, i);
  }
  
  template <class MATRIX, class VECTOR>
  void
  JacobiPreconditioner<MATRIX,VECTOR>::apply_preconditioner(const VECTOR& Px, VECTOR& x) const
  {
    x = Px;
    unsigned int i(0);
    for (typename VECTOR::iterator it(x.begin()), itend(x.end());
	 it != itend; ++it, ++i)
      *it /= A.get_entry(i, i);
  }
}

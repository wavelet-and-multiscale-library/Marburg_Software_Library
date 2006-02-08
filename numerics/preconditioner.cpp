// implementation for preconditioner.h

namespace MathTL
{
  template <class C>
  Preconditioner<C>::~Preconditioner()
  {
  }

  template <class C, class MATRIX>
  JacobiPreconditioner<C,MATRIX>::JacobiPreconditioner(const MATRIX& M)
    : A(M)
  {
  }

  template <class C, class MATRIX>
  void
  JacobiPreconditioner<C,MATRIX>::apply(const Vector<C>& x, Vector<C>& Px) const
  {
    Px = x;
    unsigned int i(0);
    for (typename Vector<C>::iterator it(Px.begin()), itend(Px.end());
	 it != itend; ++it, ++i)
      *it *= A.get_entry(i, i);
  }
  
  template <class C, class MATRIX>
  void
  JacobiPreconditioner<C,MATRIX>::apply_preconditioner(const Vector<C>& Px, Vector<C>& x) const
  {
    x = Px;
    unsigned int i(0);
    for (typename Vector<C>::iterator it(x.begin()), itend(x.end());
	 it != itend; ++it, ++i)
      *it /= A.get_entry(i, i);
  }
}

namespace MathTL
{
  template <class VECTOR, class IVP>
  Rosenbrock<VECTOR, IVP>::Rosenbrock()
    : alpha_(1), gamma_(1), b_(1)
  {
    gamma_(0, 0) = 1;
    b_[0] = 1;
  }

  template <class VECTOR, class IVP>
  Rosenbrock<VECTOR, IVP>::Rosenbrock(const LowerTriangularMatrix<double>& alpha,
				 const LowerTriangularMatrix<double>& gamma,
				 const Vector<double>& b)
    : alpha_(alpha), gamma_(gamma), b_(b)
  {
  }

  template <class VECTOR, class IVP>
  void Rosenbrock<VECTOR, IVP>::increment(const IVP& ivp,
					  const double t, const VECTOR& um,
					  const double tau,
					  VECTOR& umplus1) const
  {
  }
}

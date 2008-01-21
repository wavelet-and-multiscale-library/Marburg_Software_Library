// implementation for ring_laplacian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingLaplacian<d,dt,s0,s1>::RingLaplacian(const RingBasis<d,dt,s0,s1>& basis,
					   const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  inline
  double
  RingLaplacian<d,dt,s0,s1>::a(const Index& lambda,
			       const Index& mu) const
  {
    // Some notes concerning the strategy are in order.
    // We have to compute the integral
    //   I = int_R nabla psi_lambda(x) nabla psi_mu(x) dx
    //     = int_0^1 int_0^1 nabla(psi^0_lambda(phi,s)) M(phi,s) (nabla(psi^0_mu(phi,s)))^T |det Dkappa(phi,s)| dphi ds
    // with the matrix
    //   M(phi,s) = D(kappa^{-1})(kappa(phi,s)) (D(kappa^{-1})(kappa(phi,s)))^T.
    // In the case of the ring-shaped domain, we have the chart
    //   kappa(phi,s) = r(s)*(cos(2*pi*phi),sin(2*pi*phi)), r(s)=r0+s*(r1-r0),
    // so that
    //   M(phi,s) = diag(1/(4*pi^2*r(s)^2), 1/(r1-r0)^2).
    // Using the fact that the wavelets psi_lambda and psi_mu are tensor products,
    //   psi_lambda(phi,s) = f1(phi)*f2(s)
    //   psi_mu    (phi,s) = g1(phi)*g2(s),
    // some calculation (cf. handritten notes+Maple) yields the following
    // representation as sum of products of 1D integrals:
    //   I = 1/(4*pi^2) * I1 * I3
    //       +1/(r1-r0)^2 * I2 * I4
    //       +1/4 * I2 * I3
    //       -1/(2*(r1-r0)) * I2 * (I5+I6)
    // with
    //   I1 = int_0^1 f1'(phi) g1'(phi) dphi
    //   I2 = int_0^1 f1(phi) g1(phi) dphi
    //   I3 = int_0^1 f2(s) g2(s) r(s)^{-2} ds
    //   I4 = int_0^1 f2'(s) g2'(s) ds
    //   I5 = int_0^1 f2(s) g2'(s) r(s)^{-1} ds
    //   I6 = int_0^1 f2'(s) g2(s) r(s)^{-1} ds
    // Note that {I1,I2} and {I3,I4,I5,I6} can be computed simultaneously,
    // using the same point evaluations of the integrand.
    // We do this in two helper subroutines.

    double I1, I2, I3, I4, I5, I6;

    basis_.angular_integrals
      (typename WaveletBasis::Index0(lambda.j(), lambda.e()[0], lambda.k()[0]),
       typename WaveletBasis::Index0(mu.j(), mu.e()[0], mu.k()[0]),
       I1, I2);
    
    basis_.radial_integrals
      (typename WaveletBasis::Index1(lambda.j(), lambda.e()[1], lambda.k()[1]),
       typename WaveletBasis::Index1(mu.j(), mu.e()[1], mu.k()[1]),
       I3, I4, I5, I6);
    
    return
      I1*I3/(4*M_PI*M_PI)
      + I2*I4/((basis_.r1()-basis_.r0())*(basis_.r1()-basis_.r0()))
      + I2*I3/4
      - I2*(I5+I6)/(2*(basis_.r1()-basis_.r0()));
  }
   
}

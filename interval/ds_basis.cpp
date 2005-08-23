// implementation for ds_basis

#include <cassert>
#include <cmath>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/refinable.h>
#include <Rd/r_index.h>

namespace WaveletTL
{
  template <int d, int dT>
  DSBasis<d,dT>::DSBasis(const int s0, const int s1, const int sT0, const int sT1,
			 DSBiorthogonalizationMethod bio)
  {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
    
    assert(((s0 == 0 && sT0 > 0) || (s0 > 0 && sT0 == 0))
	   && ((s1 == 0 && sT1 > 0) || (s1 > 0 && sT1 == 0)));
    
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;

    setup_GammaLR();
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::alpha(const int m, const unsigned int r) const
  {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1<d>(); k <= ell2<d>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alpha(0, s);
	  result += cdf.a().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= (double)((1<<(r+1))-2);
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alpha(0, r-i);
      }
    }
    return result;
  }
  
  template <int d, int dT>
  const double
  DSBasis<d,dT>::alphaT(const int m, const unsigned int r) const
  {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1T<d,dT>(); k <= ell2T<d,dT>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alphaT(0, s);
	  result += cdf.aT().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= (double)((1<<(r+1))-2);
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alphaT(0, r-i);
      }
    }
    return result;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaL(const int m, const unsigned int r) const
  {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellTl; q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaLT(const int m, const unsigned int r) const
  {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < elll(); q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaR(const int m, const unsigned int r) const
  {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellTr; q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaRT(const int m, const unsigned int r) const
  {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < ellr(); q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_GammaLR() {
    // IGPMlib reference: I_Mask_Bspline::EvalGammaL(), ::EvalGammaR()

    const unsigned int GammaLsize = std::max(d-s0, dT-sT0);
    GammaL.resize(GammaLsize, GammaLsize);

    const unsigned int GammaRsize = std::max(d-s1, dT-sT1);
    GammaR.resize(GammaRsize, GammaRsize);

    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d,dT> >,2> zmask;
    InfiniteVector<double, MultiIndex<int,2> > zvalues(zmask.evaluate());

    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values

    const int I1Low  = -ell2<d>()+1;
    const int I1Up   = ell_l()-d+dT-1;
    const int I1Up_r = ell_r()-d+dT-1;
    const int I2Low  = -ell2T<d,dT>()+1;
    const int I2Up   = ellT_l()-1;
    const int I2Up_r = ellT_r()-1;

    Matrix<double> I(std::max(I1Up,I1Up_r)-I1Low+1,
		     std::max(I2Up,I2Up_r)-I2Low+1);
    for (int nu = I1Low; nu < -ell1<d>(); nu++)
      for (int mu = I2Low; mu < -ell1T<d,dT>(); mu++) {
  	double help(0);
  	const int diff(mu - nu);
   	const int sLow = std::max(-ell2<d>()+1,-ell2T<d,dT>()+1-diff);
  	for (int s = sLow; s <= nu; s++)
  	  help += zvalues.get_coefficient(MultiIndex<int,2>(s, s+diff));
  	I(-I1Low+nu, -I2Low+mu) = help; // [DKU] (5.1.7)
      }
    for (int nu = I1Low; nu <= std::max(I1Up,I1Up_r); nu++)
      for (int mu = I2Low; mu <= std::max(I2Up,I2Up_r); mu++) {
 	if ((nu >= -ell1<d>()) || ((nu <= ell_l()-1) && (mu >= -ell1T<d,dT>())))
 	  I(-I1Low+nu, -I2Low+mu) = (nu == mu ? 1 : 0); // [DKU] (5.1.6)
      }

//     cout << "I=" << endl << I << endl;
    
    // 3. finally, compute the Gramian GammaL
    for (int r = s0; r < d; r++) {
      for (int k = sT0; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_l(); nu++)
     	  for (int mu = I2Low; mu <= I2Up; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.4)
      }
      for (int k = dT; k < d+sT0-s0; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_l(); nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_l()-dT-s0-sT0+k);
	GammaL(r-s0, k-sT0) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    
    for (int r = d; r < dT+s0-sT0; r++)
      for (int k = sT0; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu <= I2Up; mu++)
    	  help += alphaT(mu, k) * I(-I1Low+ell_l()-d-s0-sT0+r, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.5)
      }

    cout << "GammaL:" << endl << GammaL << endl;

//     Matrix<double> Gammafull(dT);
//     for (int r = 0; r < d; r++)
//       for (int k = 0; k < dT; k++) {
// 	double help = 0;
// 	for (int nu = I1Low; nu < ell_l(); nu++)
// 	  for (int mu = I2Low; mu <= I2Up; mu++)
// 	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
// 	Gammafull(r, k) = help;
//       }
//     cout << "For testing, the full Gramian (without b.c.'s):" << endl
// 	 << Gammafull;

    // The same for GammaR:

    for (int r = s1; r < d; r++) {
      for (int k = sT1; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_r(); nu++)
     	  for (int mu = I2Low; mu <= I2Up_r; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.4)
      }
      for (int k = dT; k < d+sT1-s1; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_r(); nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_r()-dT-s1-sT1+k);
	GammaR(r-s1, k-sT1) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    
    for (int r = d; r < dT+s1-sT1; r++)
      for (int k = sT1; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu <= I2Up_r; mu++)
    	  help += alphaT(mu, k) * I(-I1Low+ell_r()-d-s1-sT1+r, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.5)
      }

    cout << "GammaR:" << endl << GammaL << endl;
    
  }
}

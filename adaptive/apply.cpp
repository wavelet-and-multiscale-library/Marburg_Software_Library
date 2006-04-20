// implementation for APPLY

#include <utils/array1d.h>
#include <list>
#include <map>

#include <adaptive/compression.h>

using MathTL::Array1D;

namespace WaveletTL
{
  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy)
  {
    typedef typename PROBLEM::Index Index;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    if (v.size() > 0) {
      // compute the number of bins V_0,...,V_q
      const double norm_v_sqr = l2_norm_sqr(v);
      const double norm_v = sqrt(norm_v_sqr);
      const double norm_A = P.norm_A();

      const unsigned int q = (unsigned int) std::max(ceil(log(sqrt((double)v.size())*norm_v*norm_A*2/eta)/M_LN2), 0.);
      
      // Setup the bins: The i-th bin contains the entries of v with modulus in the interval
      // (2^{-(i+1)}||v||,2^{-i}||v||], 0 <= i <= q-1, the remaining elements (with even smaller modulus)
      // are collected in the q-th bin.
      Array1D<std::list<std::pair<Index, double> > > bins(q+1);
      for (typename InfiniteVector<double,Index>::const_iterator it(v.begin());
 	   it != v.end(); ++it) {
 	const unsigned int i = std::min(q, (unsigned int)floor(-log(fabs(*it)/norm_v)/M_LN2));
	bins[i].push_back(std::make_pair(it.index(), *it));
      }
      // glue all the bins together
      Array1D<std::pair<Index, double> > v_binned(v.size());
      for (unsigned int bin = 0, id = 0; bin <= q; bin++)
	for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[bin].begin());
	     it != bins[bin].end(); ++it, ++id)
	  v_binned[id] = *it;

      const double theta = 0.5;

      // setup the segments v_{[0]},...,v_{[\ell]},
      // \ell being the smallest number such that
      //   ||A||*||v-\sum_{k=0}^\ell v_{[k]}|| <= theta * eta
      // i.e.
      //   ||v-\sum_{k=0}^\ell v_{[k]}||^2 <= eta^2 * theta^2 / ||A||^2
      // see [S, (3.9)]
      const double threshold = eta*eta*theta*theta/(norm_A*norm_A);
      unsigned int id = 0, k = 0;
      double error_sqr = norm_v_sqr;
      typename std::list<std::list<std::pair<Index, double> > > vks;
      typename std::list<double> vks_norm;
      while (true) {
	// setup the k-th segment v_{[k]}
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	}
	vks.push_back(vk);
	vks_norm.push_back(sqrt(vk_norm_sqr));

	if (error_sqr <= threshold || id >= v.size()) break; // in this case, ell=k
	k++;
      }
      const unsigned int ell = k;
      
      // compute the smallest J >= ell, such that
      //   \sum_{k=0}^{\ell} alpha_{J-k}*2^{-s(J-k)}*||v_{[k]}|| <= (1-theta) * eta
      unsigned int J = ell;
      const double s = P.s_star();
      while (true) {
	double check = 0.0;
	unsigned int k = 0;
	for (std::list<double>::const_iterator it(vks_norm.begin()); k <= ell; ++it, ++k)
	  check += P.alphak(J-k) * pow(ldexp(1.0,J-k),-s) * (*it);
	if (check <= (1-theta)*eta) break;
	J++;
      }


      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
	  add_compressed_column(P, itk->second, itk->first, J-k, w, jmax, strategy);
	}
      }
    }
  }  

//   static int its = 0; // DIRTY HACK, REMOVE THIS SOON!!!

  template <class PROBLEM>
  void RES(const PROBLEM& P,
	   const InfiniteVector<double, typename PROBLEM::Index>& w,
	   const double xi,
	   const double delta,
	   const double epsilon,
	   const int jmax,
	   InfiniteVector<double, typename PROBLEM::Index>& tilde_r,
	   double& nu,
	   const CompressionStrategy strategy)
  {
    unsigned int k = 0;
    double zeta = 2.*xi;
    //double zeta = 2.*epsilon;
    double l2n = 0.;

//     its++;

//     if (its > 100) {
//       InfiniteVector<double, typename PROBLEM::Index> help;
//       zeta *= 1.0 / (1 << 9);
//       P.RHS (zeta/2., tilde_r);
//       cout << "zeta half = " << zeta/2. << endl;
//       APPLY_COARSE(P, w, zeta/2., help, 0.00000001, jmax, strategy);
//       tilde_r -= help;
//       l2n = l2_norm(tilde_r);
//       nu = l2n + zeta;
//     }
//    else
//      {
	do {
	  zeta /= 2.;
	  P.RHS (zeta/2., tilde_r);
	  //cout << tilde_r << endl;
	  InfiniteVector<double, typename PROBLEM::Index> help;
	  //cout << w << endl;
	  cout << "zeta halbe = " << zeta/2. << endl;
	  //cout << "before aply in RES " << endl;
	  //APPLY(P, w, zeta/2., help, jmax, strategy);
	  cout << ++k << "calls of APPlY in RES" << endl;
	  APPLY_COARSE(P, w, zeta/2., help, 0.00000001, jmax, strategy);
	  //cout << "after apply in RES " << endl;
	  tilde_r -= help;
	  l2n = l2_norm(tilde_r);
	  nu = l2n + zeta;
	  //cout << "zeta = " << zeta << endl;
	  //cout << delta*l2n << endl;
	  //cout << "delta = " << delta << endl;
	  //if(k == 10)
	  //	break;
	}
 	while ( (nu > epsilon) && (zeta > delta*l2n) );
	//}
  }

}

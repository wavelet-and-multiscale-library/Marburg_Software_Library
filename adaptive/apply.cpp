// implementation for APPLY

#include <utils/array1d.h>
#include <list>
#include <map>

using MathTL::Array1D;

namespace WaveletTL
{
  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& w,
	     const int jmax)
  {
//     cout << "APPLY called, with v=" << endl << v << endl;

    typedef typename PROBLEM::WaveletBasis::Index Index;

    w.clear();

    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    if (v.size() > 0) {
      // compute the number of bins V_0,...,V_q
      const double norm_v_sqr = l2_norm_sqr(v);
//       cout << "APPLY: norm_v_sqr=" << norm_v_sqr << endl;
      const double norm_v = sqrt(norm_v_sqr);
      const double norm_A = P.norm_A();

      const unsigned int q = (unsigned int) std::max(ceil(log(sqrt((double)v.size())*norm_v*norm_A*2/eta)/M_LN2), 0.);
      
//       cout << "APPLY(): number of bins is q+1=" << q+1 << endl;

      // Setup the bins: The i-th bin contains the entries of v with modulus in the interval
      // (2^{-(i+1)}||v||,2^{-i}||v||], 0 <= i <= q-1, the remaining elements (with even smaller modulus)
      // are collected in the q-th bin.
      Array1D<std::list<std::pair<Index, double> > > bins(q+1);
      for (typename InfiniteVector<double,Index>::const_iterator it(v.begin());
 	   it != v.end(); ++it) {
 	const unsigned int i = std::min(q, (unsigned int)floor(-log(fabs(*it)/norm_v)/M_LN2));
	bins[i].push_back(std::make_pair(it.index(), *it));
      }

//       cout << "APPLY(): v in binned form looks as follows:" << endl;
//       for (unsigned int bin = 0; bin <= q; bin++) {
// 	cout << "bin=" << bin << endl;
// 	for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[bin].begin());
// 	     it != bins[bin].end(); ++it)
// 	  cout << it->first << ", " << it->second << endl;
//       }

      // glue all the bins together
      Array1D<std::pair<Index, double> > v_binned(v.size());
      for (unsigned int bin = 0, id = 0; bin <= q; bin++)
	for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[bin].begin());
	     it != bins[bin].end(); ++it, ++id)
	  v_binned[id] = *it;

//       cout << "APPLY(): the glued bins:" << endl;
//       for (unsigned int i = 0; i < v_binned.size(); i++)
//  	cout << v_binned[i].first << ", " << v_binned[i].second << endl;

      // setup the segments v_{[0]},...,v_{[\ell]},
      // \ell being the smallest number such that
      //   ||A||*||v-\sum_{k=0}^\ell v_{[k]}|| <= eta/2
      // i.e.
      //   ||v-\sum_{k=0}^\ell v_{[k]}||^2 <= eta^2/(4*||A||^2)
      // see [S, (3.9)]
      const double threshold = eta*eta/(4*norm_A*norm_A);
//       cout << "APPLY(): threshold=" << threshold << endl;
      unsigned int id = 0, k = 0;
      double error_sqr = norm_v_sqr;
      typename std::list<std::list<std::pair<Index, double> > > vks;
      typename std::list<double> vks_norm_sqr;
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
	vks_norm_sqr.push_back(vk_norm_sqr);

	if (error_sqr <= threshold || id >= v.size()) break; // in this case, ell=k
	k++;
      }
      const unsigned int ell = k;
//       cout << "APPLY(): all v_{[k]} set up, 0<=k<=ell=" << ell << endl;

      // compute the smallest J >= ell, such that
      //   \sum_{k=0}^{\ell} alpha_{J-k}*2^{-s(J-k)}*||v_{[k]}|| <= eta/2
      unsigned int J = ell;
      const double s = P.s_star();
//       cout << "APPLY(): s=" << s << endl;
      while (true) {
	double check = 0;
	unsigned int k = 0;
	for (std::list<double>::const_iterator it(vks_norm_sqr.begin()); k <= ell; ++it, ++k)
	  check += P.alphak(J-k) * pow(2.0, -s*(J-k)) * sqrt(*it);
	if (check <= eta/2.0) break;
	J++;
      }

//       cout << "APPLY(): J=" << J << endl;
//       unsigned int ncols = 0; k = 0;
//       for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
//  	   k <= ell; ++it, ++k)
//  	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
//  	     itk != it->end(); ++itk)
//  	  ncols++;
//       cout << "APPLY(): number of active columns: " << ncols << endl;

      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
// 	cout << "APPLY(): collecting partial sums, k=" << k << endl;
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
// 	  cout << "APPLY(): requesting column " << itk->first << " with J-k=" << J-k << endl;
	  P.add_column(itk->second, itk->first, J-k, w, jmax);
	}
      }
    }
  }  
}

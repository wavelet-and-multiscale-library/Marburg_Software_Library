// implementation for APPLY

#include <utils/array1d.h>
#include <list>
#include <map>



using MathTL::Array1D;

namespace WaveletTL
{
  template <class PROBLEM>
  void APPLY_TEST(const PROBLEM& P,
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
    cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
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

      //k = 4;

      const unsigned int ell = k;
      // compute the smallest J >= ell, such that
      //   \sum_{k=0}^{\ell} alpha_{J-k}*2^{-s(J-k)}*||v_{[k]}|| <= (1-theta) * eta
      unsigned int J = ell;
      const double s = P.s_star();
#if 1 // old determination of J
      while (true) {
	double check = 0.0;
	unsigned int k = 0;
	for (std::list<double>::const_iterator it(vks_norm.begin()); k <= ell; ++it, ++k)
	  check += 1.0/*P.alphak(J-k)*/ * pow(ldexp(1.0,J-k),-s) * (*it);
	if (check <= (1-theta)*eta) break;
	J++;
      }

      cout << "old J = " << J << endl;
#endif    



      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      cout << "done binning in apply..." << endl;
      cout << "ell = " << ell << endl;
      
      double tol = 10.;
      
      InfiniteVector<double, typename PROBLEM::Index> w_old;
      InfiniteVector<double, typename PROBLEM::Index> w_new;

      Vector<double> ww(P.basis().degrees_of_freedom());
      

      J = ell;

       //cout << *(P.basis().get_wavelet(4000)) << endl;
       // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      unsigned int z = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
	  add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
	  z++;
	}
      }
      
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w_old.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
	}
      }
      

      cout << "size = " << w_old.size() << endl;
      ww.resize(ww.size());
      while ( tol >  eta ) {
	w_new.clear();
	ww.resize(ww.size());
	J++;
	cout << "J = " << J << endl;
	k = 0;
	z = 0;
	for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	     k <= ell; ++it, ++k) {
	  for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	       itk != it->end(); ++itk) {
	    add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
	    z++;
	  }
	}
	
	for (unsigned int i = 0; i < ww.size(); i++) {
	  if (ww[i] != 0.) {
	    w_new.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
	  }
	}

	cout << "size = " << w_new.size() << endl;
	
	cout << "norm diff = " << l2_norm(w_old-w_new) << endl;
	tol = 1./(1.-pow(2.,-s)) * l2_norm(w_old-w_new);
	cout << "tol = " << tol << endl;

	w_old.clear();
	w_old = w_new;

	
      } // end while
      w = w_new;
    }
  }

  template <class PROBLEM>
  void APPLY_OPTIMIZED(const PROBLEM& P,
		       const InfiniteVector<double, typename PROBLEM::Index>& v,
		       const double eta,
		       InfiniteVector<double, typename PROBLEM::Index>& w,
		       double& time,
		       const int jmax,
		       const CompressionStrategy strategy)
  {
    clock_t tstart, tend;
    tstart = clock();


    typedef typename PROBLEM::Index Index;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
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
      Vector<int> maxlevels_in_bucket(q+1);
      int ii = 0;
      while (true) {
	// setup the k-th segment v_{[k]}
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	  // check whether the inserted index has the so far the largest level in the bucket 
	  if (v_binned[id].first.j() > maxlevels_in_bucket[ii])
	    maxlevels_in_bucket[ii] = v_binned[id].first.j();
	}
	ii++;
	vks.push_back(vk);
	vks_norm.push_back(sqrt(vk_norm_sqr));
	if (error_sqr <= threshold || id >= v.size()) break; // in this case, ell=k
	k++;
      }
      cout << "ell = " << k << endl;
      cout << "maximal levels in buckets = " << maxlevels_in_bucket << endl;
      const unsigned int ell = k;
      // compute the smallest J >= ell, such that
      //   \sum_{k=0}^{\ell} alpha_{J-k}*2^{-s(J-k)}*||v_{[k]}|| <= (1-theta) * eta
      unsigned int J = ell;
      const double s = P.s_star();
      while (true) {
	double check = 0.0;
	unsigned int k = 0;
	for (std::list<double>::const_iterator it(vks_norm.begin()); k <= ell; ++it, ++k)
	  check += /*P.alphak(J-k)*/1.0 * pow(ldexp(1.0,J-k),-s) * (*it);
	if (check <= (1-theta)*eta) break;
	J++;
      }

      cout << "J according to commonly used APPLY = " << J << endl;
      cout << "now computing optimized vector of J's... " << endl;
      
      int sum_n_k = 0;
      Vector<double> vks_size(ell+1);
      int i = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it = vks.begin();
	   it != vks.end(); ++it, i++) {
	sum_n_k += (*it).size();
	vks_size[i] = (*it).size();
	//vks_size[i] = ldexp(1.0,i);
      }

      //sum_n_k = ldexp(1.0,ell+1)-1;

      cout << "number of coefficients in any bucket = " << sum_n_k << endl;
      
      Vector<int> J_k(ell+1);
      const double C = 1.;//P.alphak(J-k);
      

#if 0 // 1D case
      //const double lambda = ( 2*C*(ldexp(1.0,ell+1)-1) ) / (eta * s * M_LN2);
      const double lambda = (2*C) / (eta * s * M_LN2) * sum_n_k;
      cout << "lambda = " << lambda << endl;
      i = 0;
      for (std::list<double>::const_iterator it(vks_norm.begin()); i <= ell; ++it, ++i) {
	//J_k[i] = (int)ceil(log((lambda * s * M_LN2 * (*it)) / ldexp(1.0,i)) / (M_LN2 * s));
	J_k[i] = (int)ceil(log((lambda * s * M_LN2 * (*it)) / vks_size[i]) / (M_LN2 * s));
	cout << "J_k" << "[ " << i << "] = " <<  J_k[i] << endl;
      }
#else
      //2D case
      double lambda = pow((2.0*C)/(eta),(s+1)/(s)) / s;
      i = 0;
      double tmp = 0.;
      for (std::list<double>::const_iterator it(vks_norm.begin()); i <= ell; ++it, ++i) {
	tmp += pow(vks_size[i],s/(s+1)) * pow((*it),1./(s+1));
      }
      lambda *= pow(tmp,(s+1)/s);
      i = 0;
      for (std::list<double>::const_iterator it(vks_norm.begin()); i <= ell; ++it, ++i) {
	J_k[i] = (int)ceil(log(lambda * s * (*it) / vks_size[i]) / (M_LN2*(s+1)));
	cout << "J_k" << "[ " << i << "] = " <<  J_k[i] << endl;
      }
#endif

      cout << "done binning in apply..." << endl;

      Vector<double> ww(P.basis().degrees_of_freedom());
      //cout << *(P.basis().get_wavelet(4000)) << endl;
      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      unsigned int z = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
	  //add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
	  add_compressed_column(P, itk->second, itk->first, J_k[k], ww, jmax, strategy);
	  z++;
	}
      }
      //cout << "copying vector" << endl;
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
	}
      }
    }
    time += (double)(clock()-tstart)/CLOCKS_PER_SEC;
    cout << "time = " << time << endl;
  }  


  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy)
  {
      //cout << "AUSGEFÜHRT!: " << P.basis().degrees_of_freedom() << endl; @PHK
    typedef typename PROBLEM::Index Index;
    
    //cout << "bin drin" << endl;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    //cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
    if (v.size() > 0) {
//        cout << "v größer 0" << endl;
      
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
        //cout << i << ", " << it.index() << ", " << *it << endl;
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

//      cout << ell << endl;
//      cout << "J = " << J << endl;
 //     cout << 1 << endl;

      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      //cout << "done binning in apply..." << endl;

      Vector<double> ww(P.basis().degrees_of_freedom());
 //     cout << "AUSGEFÜHRT PART2: " << P.basis().degrees_of_freedom() << endl;//HIER WEITERMACHEN @PHK
      //cout << *(P.basis().get_wavelet(4000)) << endl;
      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      unsigned int z = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
	  //add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
	  //cout << "J-k = " << J-k << endl;
 //           cout << "addcompressed wird ausgeführt" << endl;
//            cout << itk->second << ", " << itk->first << endl;
            //cout << "Beginn cc" << endl;
	  
            add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy, true);
 //         cout << "Ende cc" << endl;
//          //cout << ww << endl;
	  z++;
	}
      }
//      cout << "copying vector" << endl;
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
	}
      }
    }
//    cout << "bin raus" << endl;
  }

  template <class PROBLEM>
  void APPLY_QUARKLET(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy,
             const int pmax,
             const double a,
             const double b)
  {
      //cout << "AUSGEFÜHRT!: " << P.basis().degrees_of_freedom() << endl; @PHK
    typedef typename PROBLEM::Index Index;
    
    //cout << "bin drin" << endl;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    //cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
    if (v.size() > 0) {
//        cout << "v größer 0" << endl;
      
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
        //cout << i << ", " << it.index() << ", " << *it << endl;
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
      std::vector<int> vksize;          //for parallelization
      while (true) {
	// setup the k-th segment v_{[k]}
          vksize.resize(vksize.size()+1);
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
          //vksize[k]++;
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	}
//        cout<<"k="<<k<<endl;
//        cout<<"vksize: "<<vk.size()<<endl;
        vksize[k]=vk.size();
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

//      cout << "ell: " << ell << endl;
//      cout << "J: " << J << endl;
//      J-=5;
//      cout << "J = " << J << endl; //maybe J is too big, and that is the reason for the enormous amount of Degrees of freedom @PHK
 //     cout << 1 << endl;

      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      //cout << "done binning in apply..." << endl;

      Vector<double> ww(P.frame().degrees_of_freedom());
      //Vector<double> wwhelp(P.frame().degrees_of_freedom());
      //cout << ww<<endl;
 //     cout << "AUSGEFÜHRT PART2: " << P.basis().degrees_of_freedom() << endl;//HIER WEITERMACHEN @PHK
      //cout << *(P.basis().get_wavelet(4000)) << endl;
      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
//      for(int i=0;i<vksize.size();i++) cout<<vksize[i]<<endl;
#if PARALLEL==1
      //outer loop parallelization
//#pragma omp parallel
//{       
//          if(omp_get_thread_num()==0)
//            cout<<"number of threads: "<<omp_get_num_threads()<<endl;
//          Vector<double> wwprivate(P.frame().degrees_of_freedom());        
//          #pragma omp for
//      for(k=0;k<ell;k++){
//          typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
//          advance(it, k); 
//          unsigned int z = 0;
//            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
//	     itk != it->end(); ++itk,++z) {
          
          //inner loop parallelization
          k = 0;
          for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	      k <= ell; ++it, ++k) {
//              cout <<vksize[k]<<endl;
//              int my_num_threads = (vksize[k]>=64) ? 8 : 1;
#pragma omp parallel //num_threads(my_num_threads)
{
//            if(omp_get_thread_num()==0)
//                cout<<"number of threads: "<<omp_get_num_threads()<<endl;
            Vector<double> wwprivate(P.frame().degrees_of_freedom()); 
            unsigned int z = 0;
#pragma omp for
            for(int z=0;z<vksize[k];z++){
              typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
              advance(itk,z);          
              
            //collapse parallelization
//#pragma omp parallel
//{
//            Vector<double> wwprivate(P.frame().degrees_of_freedom()); 
//#pragma omp for collapse(2)
//            for(int k=0;k<ell;k++){
//                typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
//                advance(it, k);
//                for(int z=0;z<vksize[k];z++){
//                    typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
//                    advance(itk,z);
#else
      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
           unsigned int z = 0;
            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk,++z) {
#endif
          //  Vector<double> wwhelp(P.frame().degrees_of_freedom());
	  //add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
//            cout << "J-k = " << J-k << endl;
//            cout << "addcompressed wird ausgeführt" << endl;
//            cout << itk->second << ", " << itk->first << endl;
            //cout << "Beginn cc" << endl;
           // add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true, pmax, a, b);
           //Vector<double> wwhelp(P.frame().degrees_of_freedom());
#if PARALLEL==1
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, wwprivate, jmax, strategy, true);
#else
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true);
#endif
//          cout << "Ende cc" << endl;
          //cout << ww << endl;
	  //z++;
	}
//        cout<<"z="<<z<<":"<<vksize[k]<<endl;
           //outer loop or collapse parallelization
//      }
#if PARALLEL==1
#pragma omp critical
{
        for(int i=0;i<ww.size();i++){
            ww(i)+=wwprivate(i);
        }
}
} 
#endif
      //inner loop parallelization
       }
//      cout<<"k= "<<k<<endl;
//          cout<<"vksize="<<vksize<<endl;
//      cout << "copying vector" << endl;
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.frame().get_quarklet(i)), ww[i]);
	}
      }
    }
//    cout << "bin raus" << endl;
  }
      
      template <class PROBLEM>
  void APPLY_QUARKLET_SEQUENTIAL(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy,
             const int pmax,
             const double a,
             const double b)
  {
    typedef typename PROBLEM::Index Index;
    w.clear();
    if (v.size() > 0) {
      const double norm_v_sqr = l2_norm_sqr(v);
      const double norm_v = sqrt(norm_v_sqr);
      const double norm_A = P.norm_A();
      
      const unsigned int q = (unsigned int) std::max(ceil(log(sqrt((double)v.size())*norm_v*norm_A*2/eta)/M_LN2), 0.);
      Array1D<std::list<std::pair<Index, double> > > bins(q+1);
      for (typename InfiniteVector<double,Index>::const_iterator it(v.begin());
 	   it != v.end(); ++it) {
 	const unsigned int i = std::min(q, (unsigned int)floor(-log(fabs(*it)/norm_v)/M_LN2));
	bins[i].push_back(std::make_pair(it.index(), *it));
      }
      
      Array1D<std::pair<Index, double> > v_binned(v.size());
      for (unsigned int bin = 0, id = 0; bin <= q; bin++)
	for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[bin].begin());
	     it != bins[bin].end(); ++it, ++id)
	  v_binned[id] = *it;

      const double theta = 0.5;
      const double threshold = eta*eta*theta*theta/(norm_A*norm_A);
      unsigned int id = 0, k = 0;
      double error_sqr = norm_v_sqr;
      typename std::list<std::list<std::pair<Index, double> > > vks;
      typename std::list<double> vks_norm;
      std::vector<int> vksize;          //for parallelization
      while (true) {
          vksize.resize(vksize.size()+1);
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
          //vksize[k]++;
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	}
//        cout<<"k="<<k<<endl;
//        cout<<"vksize: "<<vk.size()<<endl;
        vksize[k]=vk.size();
	vks.push_back(vk);
	vks_norm.push_back(sqrt(vk_norm_sqr));
	if (error_sqr <= threshold || id >= v.size()) break; // in this case, ell=k
	k++;   
      }
      const unsigned int ell = k;
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
      Vector<double> ww(P.frame().degrees_of_freedom());
//      for(int i=0;i<vksize.size();i++) cout<<vksize[i]<<endl;

      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
           unsigned int z = 0;
            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk,++z) {


            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true);

	}
//        cout<<"z="<<z<<":"<<vksize[k]<<endl;
           //outer loop parallelization
      }

      
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.frame().get_quarklet(i)), ww[i]);
	}
      }
    }
//    cout << "bin raus" << endl;
  }
      
      template <class PROBLEM>
  void APPLY_QUARKLET_PARALLEL_OUTER(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy,
             const int pmax,
             const double a,
             const double b)
  {
      //cout << "AUSGEFÜHRT!: " << P.basis().degrees_of_freedom() << endl; @PHK
    typedef typename PROBLEM::Index Index;
    
    //cout << "bin drin" << endl;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    //cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
    if (v.size() > 0) {
//        cout << "v größer 0" << endl;
      
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
        //cout << i << ", " << it.index() << ", " << *it << endl;
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
      std::vector<int> vksize;          //for parallelization
      while (true) {
	// setup the k-th segment v_{[k]}
          vksize.resize(vksize.size()+1);
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
          //vksize[k]++;
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	}
//        cout<<"k="<<k<<endl;
//        cout<<"vksize: "<<vk.size()<<endl;
        vksize[k]=vk.size();
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

//      cout << "ell: " << ell << endl;
//      cout << "J: " << J << endl;
//      J-=5;
//      cout << "J = " << J << endl; //maybe J is too big, and that is the reason for the enormous amount of Degrees of freedom @PHK
 //     cout << 1 << endl;

      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      //cout << "done binning in apply..." << endl;

      Vector<double> ww(P.frame().degrees_of_freedom());
      //Vector<double> wwhelp(P.frame().degrees_of_freedom());
      //cout << ww<<endl;
 //     cout << "AUSGEFÜHRT PART2: " << P.basis().degrees_of_freedom() << endl;//HIER WEITERMACHEN @PHK
      //cout << *(P.basis().get_wavelet(4000)) << endl;
      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
//      for(int i=0;i<vksize.size();i++) cout<<vksize[i]<<endl;
#if 1
      //outer loop parallelization
#pragma omp parallel
{       
//          if(omp_get_thread_num()==0)
//            cout<<"number of threads: "<<omp_get_num_threads()<<endl;
          Vector<double> wwprivate(P.frame().degrees_of_freedom());        
          #pragma omp for
      for(k=0;k<ell;k++){
          typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
          advance(it, k); 
          unsigned int z = 0;
            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk,++z) {
                       
              
#else
      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
           unsigned int z = 0;
            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk,++z) {
#endif
         
#if 1
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, wwprivate, jmax, strategy, true);
#else
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true);
#endif
//          cout << "Ende cc" << endl;
          //cout << ww << endl;
	  //z++;
	}
//        cout<<"z="<<z<<":"<<vksize[k]<<endl;
           //outer loop parallelization
      }
#if 1
#pragma omp critical
{
        for(int i=0;i<ww.size();i++){
            ww(i)+=wwprivate(i);
        }
}
} 
#endif
      //inner loop parallelization
//       }
//      cout<<"k= "<<k<<endl;
//          cout<<"vksize="<<vksize<<endl;
//      cout << "copying vector" << endl;
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.frame().get_quarklet(i)), ww[i]);
	}
      }
    }
//    cout << "bin raus" << endl;
  }
   
      template <class PROBLEM>
  void APPLY_QUARKLET_PARALLEL_INNER(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::Index>& w,
	     const int jmax,
	     const CompressionStrategy strategy,
             const int pmax,
             const double a,
             const double b)
  {
      //cout << "AUSGEFÜHRT!: " << P.basis().degrees_of_freedom() << endl; @PHK
    typedef typename PROBLEM::Index Index;
    
    //cout << "bin drin" << endl;

    w.clear();
    // Binary Binning variant of APPLY from [S],[B]
    // Remark: it is possible to perform binary binning without actually assembling
    // the bins, however, in this first version we do setup the bins to avoid
    // unnecessary difficulties
    //cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
    if (v.size() > 0) {
//        cout << "v größer 0" << endl;
      
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
        //cout << i << ", " << it.index() << ", " << *it << endl;
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
      std::vector<int> vksize;          //for parallelization
      while (true) {
	// setup the k-th segment v_{[k]}
          vksize.resize(vksize.size()+1);
	std::list<std::pair<Index, double> > vk;
	double vk_norm_sqr = 0;
	for (unsigned int n = 1; error_sqr > threshold && id < v.size() && n <= ldexp(1.0, k)-floor(ldexp(1.0, k-1)); n++, id++) {
	  vk.push_back(v_binned[id]);
          //vksize[k]++;
	  const double help = v_binned[id].second * v_binned[id].second;
	  error_sqr -= help;
	  vk_norm_sqr += help;
	}
//        cout<<"k="<<k<<endl;
//        cout<<"vksize: "<<vk.size()<<endl;
        vksize[k]=vk.size();
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

//      cout << "ell: " << ell << endl;
//      cout << "J: " << J << endl;
//      J-=5;
//      cout << "J = " << J << endl; //maybe J is too big, and that is the reason for the enormous amount of Degrees of freedom @PHK
 //     cout << 1 << endl;

      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      //cout << "done binning in apply..." << endl;

      Vector<double> ww(P.frame().degrees_of_freedom());
      //Vector<double> wwhelp(P.frame().degrees_of_freedom());
      //cout << ww<<endl;
 //     cout << "AUSGEFÜHRT PART2: " << P.basis().degrees_of_freedom() << endl;//HIER WEITERMACHEN @PHK
      //cout << *(P.basis().get_wavelet(4000)) << endl;
      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
//      for(int i=0;i<vksize.size();i++) cout<<vksize[i]<<endl;
#if 1
      //outer loop parallelization
//#pragma omp parallel
//{       
//          if(omp_get_thread_num()==0)
//            cout<<"number of threads: "<<omp_get_num_threads()<<endl;
//          Vector<double> wwprivate(P.frame().degrees_of_freedom());        
//          #pragma omp for
//      for(k=0;k<ell;k++){
//          typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
//          advance(it, k); 
//          unsigned int z = 0;
//            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
//	     itk != it->end(); ++itk,++z) {
          
          //inner loop parallelization
          k = 0;
          for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	      k <= ell; ++it, ++k) {
//              cout <<vksize[k]<<endl;
#pragma omp parallel
{
//            if(omp_get_thread_num()==0)
//                cout<<"number of threads: "<<omp_get_num_threads()<<endl;
            Vector<double> wwprivate(P.frame().degrees_of_freedom());    
#pragma omp for
            for(int z=0;z<vksize[k];z++){
              typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
              advance(itk,z);
              
              
#else
      k = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
           unsigned int z = 0;
            for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk,++z) {
#endif
          //  Vector<double> wwhelp(P.frame().degrees_of_freedom());
	  //add_compressed_column(P, itk->second, itk->first, J-k, ww, jmax, strategy);
//            cout << "J-k = " << J-k << endl;
//            cout << "addcompressed wird ausgeführt" << endl;
//            cout << itk->second << ", " << itk->first << endl;
            //cout << "Beginn cc" << endl;
           // add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true, pmax, a, b);
           //Vector<double> wwhelp(P.frame().degrees_of_freedom());
#if 1
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, wwprivate, jmax, strategy, true);
#else
            add_compressed_column_quarklet(P, itk->second, itk->first, J-k, ww, jmax, strategy, true);
#endif
//          cout << "Ende cc" << endl;
          //cout << ww << endl;
	  //z++;
	}
//        cout<<"z="<<z<<":"<<vksize[k]<<endl;
           //outer loop parallelization
//      }
#if 1
#pragma omp critical
{
        for(int i=0;i<ww.size();i++){
            ww(i)+=wwprivate(i);
        }
}
} 
#endif
      //inner loop parallelization
       }
//      cout<<"k= "<<k<<endl;
//          cout<<"vksize="<<vksize<<endl;
//      cout << "copying vector" << endl;
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (ww[i] != 0.) {
	  w.set_coefficient(*(P.frame().get_quarklet(i)), ww[i]);
	}
      }
    }
//    cout << "bin raus" << endl;
  }
      


  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const int p,
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
    //cout << "entering apply..." << endl;

    //cout << "size = " << v.size() << endl;
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

      // hack: We let 'add_compressed_column' and 'add_level'
      // in cached_problem.cpp/.h work on full vectors. We do this because the call of 
      // 'w.add_coefficient();' in 'add_level' is inefficient. 
      // Below we will then copy ww into the sparse vector w.
      // Probably this will be handled in a more elegant way in the near future.

      //cout << "done binning in apply..." << endl;

      Vector<double> ww(P.basis().degrees_of_freedom());

      // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
      k = 0;
      unsigned int z = 0;
      for (typename std::list<std::list<std::pair<Index, double> > >::const_iterator it(vks.begin());
	   k <= ell; ++it, ++k) {
	for (typename std::list<std::pair<Index, double> >::const_iterator itk(it->begin());
	     itk != it->end(); ++itk) {
	  //cout << "J-k = " << J-k << endl;
	  add_compressed_column(P, p, itk->second, itk->first, J-k, ww, jmax, strategy);
	  z++;
	}
      }
      // copy ww into w
      for (unsigned int i = 0; i < ww.size(); i++) {
	if (fabs(ww[i]) > 1.0e-16) {
	  //if (ww[i] != 0.) {
	  w.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
	}
      }
    }
  
    //cout << "done apply" << endl;
  }  



  //  static int its = 0; // DIRTY HACK, REMOVE THIS SOON!!!

  template <class PROBLEM>
  void RES(const PROBLEM& P,
	   const InfiniteVector<double, typename PROBLEM::Index>& w,
	   const double xi,
	   const double delta,
	   const double epsilon,
	   const int jmax,
	   InfiniteVector<double, typename PROBLEM::Index>& tilde_r,
	   double& nu,
	   unsigned int& niter,
	   const CompressionStrategy strategy
	   )
  {
    //    unsigned int k = 0;
    double zeta = 2.*xi;
    double l2n = 0.;
    do {
      zeta /= 2.;
      P.RHS (zeta/2., tilde_r);
      InfiniteVector<double, typename PROBLEM::Index> help;
      //APPLY(P, w, .0/*zeta/2.*/, help, jmax, strategy);
      APPLY_COARSE(P, w, zeta/2., help, 1.0e-6, jmax, strategy);
      //APPLY_COARSE(P, w, zeta/2., help, 0.5, jmax, strategy);
      tilde_r -= help;
      l2n = l2_norm(tilde_r);
      nu = l2n + zeta;
      ++niter;
    }
    while ( (nu > epsilon) && (zeta > delta*l2n) );
    
  }
  
  template <class PROBLEM>
  void RES_QUARKLET(const PROBLEM& P,
	   const InfiniteVector<double, typename PROBLEM::Index>& w,
	   const double xi,
	   const double delta,
	   const double epsilon,
	   const int jmax,
	   InfiniteVector<double, typename PROBLEM::Index>& tilde_r,
	   double& nu,
	   unsigned int& niter,
	   const CompressionStrategy strategy,
           const int pmax,
           const double a,
           const double b
	   )
  {
    //    unsigned int k = 0;
    double zeta = 2.*xi;
    double l2n = 0.;
    do {
      zeta /= 2.;
      P.RHS (zeta/2., tilde_r);
      InfiniteVector<double, typename PROBLEM::Index> help;
      //APPLY(P, w, .0/*zeta/2.*/, help, jmax, strategy);
      APPLY_QUARKLET_COARSE(P, w, zeta/2., help, jmax, strategy, pmax, a, b, 1.0e-6);
      //APPLY_COARSE(P, w, zeta/2., help, 0.5, jmax, strategy);
      tilde_r -= help;
      l2n = l2_norm(tilde_r);
      nu = l2n + zeta;
      ++niter;
    }
    while ( (nu > epsilon) && (zeta > delta*l2n) );
    
  }

  
           
  
}

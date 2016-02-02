// implementation for apply_tensor.h

namespace WaveletTL
{
    template <class PROBLEM>
    void APPLY_TENSOR(PROBLEM& P,
            const InfiniteVector<double, typename PROBLEM::Index>& v,
            const double eta,
            InfiniteVector<double, typename PROBLEM::Index>& w,
            const int jmax,
            const CompressionStrategy strategy,
            const bool preconditioning)
    {
        // Remark: Remark from APPLY applies here as well, since binary binning part is the similar.
        // linfty norm is used, resulting in the Factor 2 for p.
        
        typedef typename PROBLEM::Index Index;
        w.clear();
        if (v.size() > 0) 
        {
            // compute the number of bins V_0,...,V_q
            const double norm_v = linfty_norm(v);         
#if _APPLY_TENSOR_DEBUGMODE == 1
            if (norm_v > 1e+10)
            {
                cout << "APPLY_TENSOR (Index):: l_infty norm of v is too high!"<<endl;
                cout << "v = " << v << endl;
            }
#endif
            // Explicit sorting into bins is not necessary!
            // We may compute the bin number, norm of all entries in the bins, and number of entries in each bin non the less.
            // this enables us to compute the j_p in [DSS]
            
            // Theory: The i-th bin contains the entries of v with modulus in the interval
            // (2^{-(i+1)/2}||v||_\infty,2^{-i/2}||v||_\infty], 0 <= i <= q-1, the remaining elements (with even smaller modulus)
            // are collected in the q-th bin.
            // q is chosen s.t. elements not in any bin have norm <= threshold (= eta/(2*norm_A))
            // The entries not in any bin are discarded.
            const double norm_A = P.norm_A();
            double threshold = eta/2.0/norm_A;
            //const unsigned int q = (unsigned int) std::max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.);
            const unsigned int q = std::max ((unsigned int)1, 2*log2( (unsigned int) std::max(1, (int) ceil (sqrt(v.size())*norm_v/threshold) ) ) );
#if _APPLY_TENSOR_DEBUGMODE == 1
            assert (q == (unsigned int) std::max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.) );
#endif
            Array1D<double> bin_norm_sqr(q); //square norm of the coeffs in the bins
            Array1D<int> bin_size(q); // number of entries in each bin
            // the bin number corresponding to an entry of v is computed twice in this code. storing it into an InfiniteVector is probably not faster
            //InfiniteVector<int, typename PROBLEM::Index> bin_number; // for each entry of v: store its bin - so the binnumber has to be computed only once - may be unnecessary
#if _APPLY_TENSOR_DEBUGMODE == 1
            for (unsigned int i=0; i< q;++i)
            {
                assert (bin_norm_sqr[i] == 0);
                assert (bin_size[i] == 0);
            }
#endif
            {
                unsigned int temp_i;
                for (typename InfiniteVector<double,Index>::const_iterator it(v.begin());it != v.end(); ++it)
                {
                    //const unsigned int i = std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2));
                    //i = std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2));
                    //temp_i = std::min(q, 2*floor(log2(floor(norm_v/fabs(*it)))) );
#if _APPLY_TENSOR_DEBUGMODE == 1
                    assert (std::min(q, 2*floor(log2(floor(norm_v/fabs(*it)))) ) = std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2)));
#endif
                    temp_i = 2*log2(floor(norm_v/fabs(*it))) -1; // bins: 1,...,q but observe temp_i: 0,...,q-1
                    if (temp_i < q)
                    {
                        //bins[i].push_back(std::make_pair(it.index(), *it));
                        bin_norm_sqr[temp_i] += (*it)*(*it);
                        bin_size[temp_i] += 1;
                    }
                }
            } // end of sorting into bins
            
            // In difference to the isotropic APPLY we use the bins directly as the segments v_{[0]},...,v_{[\ell]}.
            // [DSS]: compute smallest \ell such that
            //   ||v-\sum_{k=0}^\ell v_{[k]}|| <= threshold
            // i.e.
            //   ||v-\sum_{k=0}^\ell v_{[k]}||^2 <= threshold^2,
            // where threshold = eta * theta / ||A||, theta = 1/2
            threshold = threshold*threshold;
            double delta(l2_norm_sqr(v));
            //double bound = temp_d - threshold;
            unsigned int ell=0, num_of_relevant_entries(0); // number of bins needed to approximate v; number_of_coeffs in the first ell buckets
            if (delta <= threshold)
            {
                // v is well approximated by 0
                return;
            }
            else
            {
                while (delta > threshold)
                {
                    num_of_relevant_entries += bin_size[ell];
                    delta -= bin_norm_sqr[ell];
                    //bound -= bin_norm_sqr[ell];
                    ell++;
                    if (ell > q)
                    {
                        break;
                    }
                }
            }
            // we need buckets 0,...,ell-1
            
            
#if 0 // compute ell in a different way to check validity of code. needs explicit computation of bins
            double bound = l2_norm_sqr(v) - threshold;
            //double error_sqr = l2_norm_sqr(v);
            double new_norm_sqr(0);
            //Array1D<double> bin_norm_sqr(q+1); //square norm of the coeffs in the bins
            unsigned int ell2(0);
            unsigned int num_relevant_entries=0;
            if (l2_norm_sqr(v) > threshold)
            {
                unsigned int vsize = v.size();
                while (true)
                {
                    bin_norm_sqr[ell2]=0;
                    for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[ell].begin()), itend(bins[ell].end()); it != itend; ++it)
                    {
                        bin_norm_sqr[ell] += (it->second) * (it->second);
                        //num_relevant_entries++;
                    }
                    if (bin_norm_sqr[ell2] < 1e-15 && bins[ell2].size() > 0) // this bin contains only very small coefs. Thus it and all later bins are not of any importance
                    {
                        //error_sqr = threshold; // we are near the machine limit
                        if (ell2 > 0)
                        {
                            --ell2;
                        }
                        break;
                    }
                    num_relevant_entries += bins[ell2].size();
                    new_norm_sqr += bin_norm_sqr[ell2];
                    //error_sqr -= bin_norm_sqr[ell2];
                    if (new_norm_sqr >= bound)
                        //if (error_sqr <= threshold)
                    {
                        break;
                    }
                    else if (num_relevant_entries == vsize)
                    {
                        //error_sqr = 0;
                        break;
                    }
                    else
                    {
                        ell2++;
                    }
                }
            }
            assert (ell == ell2);
#endif
            // Compute z = \sum_{p=0}^ell A^{jp}v_{[p]},
            // where A^{(jp)} is the optimal approximation matrix for each bin.
            // Optimality is meant in the sense that the computational effort is minimal and the error bound is met.
            // The computation is tailored for biorthogonal anisotropical wavelets:
            // A^{jp} has C*jp^dim nontrivial entries per row/column which need C*jp^dim many operations to be computed, further
            // \| A -A^{jp} \| \leq D 2^{-\rho jp} (A is compressible with s* = \infty)
            // jp is computed such that effort for apply is minimal, i.e.,
            // \sum _{p=0}^ell C*jp^dim * bin_size[p] -> minimal,
            // under the error bound condition
            // \sum_{p=0}^ell D*s^{-\rho jp} \|v_p\|_2 \leq eta-delta,
            // where rho < 1/2 is the compressibility of A (valid for biorthogonl anisotropical wavelets)
            // eta is the target accuracy of APPLY,
            // delta = \|A\|\|v - \sum_{p=1}^\ell v_p\|_2 (<= eta/2)
            delta = sqrt(delta) * norm_A;
            assert (delta <= eta/2);
            // we use the formula for jp tilde from [DSS] for orthogonal wavelets, i.e., effort C*jp^1.
            // effect: error \eta-\delta is met, but effort may be suboptimal
            // => with D = P.alphak(p)
            //jp_tilde = ceil (log(sqrt(bin_norm_sqr[p])*num_relevant_entries*P.alphak(p)/(bin_size[p]*(eta-delta)))/M_LN2 );
            
            Array1D<int> jp_tilde(ell);
            for (unsigned int i = 0; i < ell; ++i)
            {
                //jp_tilde[i] = ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ); 
                jp_tilde[i] = ceil (log2(ceil (sqrt(bin_norm_sqr[i]) * num_of_relevant_entries * P.alphak(i) / (bin_size[i]*(eta-delta))))  );
#if _APPLY_TENSOR_DEBUGMODE == 1
                assert (jp_tilde[i] ==  ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ) );
#endif
            }
            // hack: We work with full vectors (of size degrees_of_freedom).
            // We do this because adding sparse vectors seems to be inefficient.
            // Below we will then copy ww into the sparse vector w.
            // Probably this will be handled in a more elegant way someday.
            // ww should be of a similar structure to the cache in P
            Vector<double> ww(P.basis().degrees_of_freedom());
            //cout << *(P.basis().get_wavelet(4000)) << endl;
            // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}
            for (typename InfiniteVector<double,Index>::const_iterator it(v.begin()), itend(v.end());it != itend; ++it)
            {
#if _APPLY_TENSOR_DEBUGMODE == 1
                assert ( (2*log2(floor(norm_v/fabs(*it))) -1) >= 0 );
                assert ( (2*log2(floor(norm_v/fabs(*it))) -1) < ell );
#endif
                //temp_i = 2*log2(floor(norm_v/fabs(*it))) -1
                //add_compressed_column_tensor(P, *it, it, jp_tilde[2*log2(floor(norm_v/fabs(*it))) -1], ww, jmax, strategy, preconditioning);
                P.add_ball(it.index(),ww,jp_tilde[2*log2(floor(norm_v/fabs(*it))) -1],*it,jmax,strategy,preconditioning);
            }
            // copy ww into w
            for (unsigned int i = 0; i < ww.size(); i++) 
            {
                if (ww[i] != 0.)
                {
                    w.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
                }
            }
        }
    }
    
    template <class PROBLEM>
    void APPLY_TENSOR(PROBLEM& P,
            const InfiniteVector<double, int>& v,
            const double eta,
            InfiniteVector<double, int>& w,
            const int jmax,
            const CompressionStrategy strategy,
            const bool preconditioning)
    {
        // Remark: Remark from APPLY applies here as well, since binary binning part is the similar.
        // linfty norm is used, resulting in the Factor 2 for p.
        w.clear();
        if (v.size() > 0) 
        {
            // compute the number of bins V_0,...,V_q
            const double norm_v = linfty_norm(v);  
#if _APPLY_TENSOR_DEBUGMODE == 1
            if (norm_v > 1e+10)
            {
                cout << "APPLY_TENSOR (int):: l_infty norm of v is too high!"<<endl;
                cout << "v = " << v << endl;
            }
#endif
            // Explicit sorting into bins is not necessary!
            // We may compute the bin number, norm of all entries in the bins, and number of entries in each bin non the less.
            // this enables us to compute the j_p in [DSS]
            
            // Theory: The i-th bin contains the entries of v with modulus in the interval
            // (2^{-(i+1)/2}||v||_\infty,2^{-i/2}||v||_\infty], 0 <= i <= q-1, the remaining elements (with even smaller modulus)
            // are collected in the q-th bin (discarded).
            // q is chosen s.t. elements not in any bin have norm <= threshold (= eta/(2*norm_A))
            // The entries not in any bin are discarded.
            const double norm_A = P.norm_A();
            double threshold = eta/2.0/norm_A;
            const unsigned int q = (unsigned int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.);
            //const unsigned int q = (unsigned int)(log2( (unsigned int) ceil (v.size()*norm_v*norm_v/threshold/threshold) ) +1);
            //const unsigned int q = (unsigned int)(2* log2( (unsigned int) ceil (sqrt(v.size())*norm_v/threshold) ) +1);
#if _APPLY_TENSOR_DEBUGMODE == 1
            assert (ceil (v.size()*norm_v*norm_v/threshold/threshold) > 0);
// CLEANUP            
                    
            if (abs (q - (unsigned int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) ) > 1 )
            {
                cout << "q = " << q << endl;
                cout << "v.size()*norm_v*norm_v/threshold/threshold = " << v.size()*norm_v*norm_v/threshold/threshold << endl;
                cout << "log (v.size()*norm_v*norm_v/threshold/threshold) / M_LN2 = " << log (v.size()*norm_v*norm_v/threshold/threshold) / M_LN2 << endl<< endl;
                
                cout << "sqrt(v.size())*norm_v/threshold = " << sqrt(v.size())*norm_v/threshold << endl;
                cout << "log (sqrt(v.size())*norm_v/threshold) / M_LN2 = " << log (sqrt(v.size())*norm_v/threshold) / M_LN2 << endl << endl;
                
                cout << "(unsigned int) ceil (sqrt(v.size())*norm_v/threshold) = " << (unsigned int) ceil (sqrt(v.size())*norm_v/threshold) << endl;
                cout << "log2( (unsigned int) ceil (sqrt(v.size())*norm_v/threshold) ) = " << log2( (unsigned int) ceil (sqrt(v.size())*norm_v/threshold) ) << endl << endl;
                
                cout << "(unsigned int) ceil (v.size()*norm_v*norm_v/threshold/threshold) = " << (unsigned int) ceil (v.size()*norm_v*norm_v/threshold/threshold) << endl;
                cout << "(unsigned int) ceil (v.size()*norm_v*norm_v*norm_A*2.0/eta*norm_A*2.0/eta) = " << (unsigned int) ceil (v.size()*norm_v*norm_v*norm_A*2.0/eta*norm_A*2.0/eta)<< endl;
                cout << "(sqrt((double)v.size())*norm_v*norm_A*2.0/eta) = " << (sqrt((double)v.size())*norm_v*norm_A*2.0/eta) << endl;
                cout << "ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2)  = " << ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) << endl;
                cout << "v.size() = " << v.size() << "; norm_v = " << norm_v << "; norm_A = " << norm_A << "; eta = " << eta << endl;
                
                cout << "q = " << q << "; (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) = " << (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) << endl;
                cout << "abs ((int)q - (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) ) = " << abs ((int)q - (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) ) << endl;
            }
                
            assert (abs ((int)q - (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) ) <= 1 );
            assert ((int)q<= (int) std::max(ceil(2*log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) , 0.) );
#endif
            Array1D<double> bin_norm_sqr(q); //square norm of the coeffs in the bins
            Array1D<int> bin_size(q); // number of entries in each bin
            // the bin number corresponding to an entry of v is computed twice in this code. storing it into an InfiniteVector is probably not faster! E.g.,
            //InfiniteVector<int, typename PROBLEM::Index> bin_number; // for each entry of v: store its bin - so the binnumber has to be computed only once - may be unnecessary
            for (unsigned int i=0; i<q;++i)
            {
                bin_norm_sqr[i] = 0;
                bin_size[i] = 0;
            }
            unsigned int temp_i;
            for (typename InfiniteVector<double,int>::const_iterator it(v.begin());it != v.end(); ++it)
            {
                temp_i = (unsigned int)floor(2*log(norm_v/fabs(*it))/M_LN2); // bins: 1,...,q but observe temp_i: 0,...,q-1
                // temp_i = 2*log2( (unsigned int) floor(norm_v/fabs(*it))); // bins: 1,...,q but observe temp_i: 0,...,q-1
#if _APPLY_TENSOR_DEBUGMODE == 1
                if (std::min(q, temp_i) != std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2)))
                {
                    cout << endl << "Line 307" << endl;
                    cout << "q = " << q << "; norm_v = " << norm_v << "; fabs(*it) = " << fabs(*it) << endl;
                    cout << "floor(norm_v/fabs(*it)) = " << floor(norm_v/fabs(*it)) << endl;
                    cout << "(unsigned int) floor(norm_v/fabs(*it)) = " << (unsigned int) floor(norm_v/fabs(*it)) << endl;
                    cout << "2*log2( (unsigned int) floor(norm_v/fabs(*it))) = " << 2*log2( (unsigned int) floor(norm_v/fabs(*it))) << endl << endl;;
                    cout << "fabs(*it)/norm_v = " << fabs(*it)/norm_v << endl;
                    cout << "floor(-2*log(fabs(*it)/norm_v)/M_LN2) = " << floor(-2*log(fabs(*it)/norm_v)/M_LN2) << endl;
                    cout << "(unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2) = " << (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2) << endl;
                }
                assert (std::min(q, temp_i) == std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2)));
#endif
                if (temp_i < q)
                {
                    //bins[i].push_back(std::make_pair(it.index(), *it));
                    bin_norm_sqr[temp_i] += (*it)*(*it);
                    bin_size[temp_i] += 1;
                }
            } // end of sorting into bins
            
#if _APPLY_TENSOR_DEBUGMODE == 1
            for (unsigned int i=0; i<q; ++i)
            {
                assert ( (bin_size[i] == 0) == (bin_norm_sqr[i] == 0));
            }
#endif
            // q is pessimistic!
            
            int number_of_buckets(q);
            while (number_of_buckets > 0)
            {
                if (bin_size[number_of_buckets-1] != 0)
                    break;
                --number_of_buckets;
            }
            
            
            //cout << "bin_size = " << bin_size << endl;
            
            // In difference to the isotropic APPLY we use the bins directly as the segments v_{[0]},...,v_{[\ell]}.
            // [DSS]: compute smallest \ell such that
            //   ||v-\sum_{k=0}^\ell v_{[k]}|| <= threshold
            // i.e.
            //   ||v-\sum_{k=0}^\ell v_{[k]}||^2 <= threshold^2,
            // where threshold = eta * theta / ||A||, theta = 1/2
           
            // if threshold^2 > 1e-14 the computation of the squared l2 norm fails.
            // we will always need all entries of v in this setting
            
            threshold = threshold*threshold;
            double delta(l2_norm_sqr(v));
            if (delta <= threshold)
            {
                // v is well approximated by 0
                return;
            }
            unsigned int ell=0, num_of_relevant_entries(0); // number of bins needed to approximate v; number_of_coeffs in the first ell buckets
            if (threshold > 1e-12)
            {
                //double bound = temp_d - threshold;
                double temp_d(0);
                cout << "delta = " << delta << endl;
                while (delta > threshold)
                {
                    if (bin_size[ell] != 0)
                    {
                        num_of_relevant_entries += bin_size[ell];
                        delta -= bin_norm_sqr[ell];
                        temp_d += bin_norm_sqr[ell];
                        cout << "bin_norm_sqr[" << ell << "]= " << bin_norm_sqr[ell] << "; temp_d = " << temp_d << "; delta = " << delta << "; threshold = " << threshold << endl;
    //                    bound -= bin_norm_sqr[ell];
                    }
                    ell++;
                    if (ell >= number_of_buckets)
                    {
                        cout << "we need all buckets!" << endl;
                        break;
                    }
                }
                if (delta < 0)
                {
                    delta = 0;
                    assert (num_of_relevant_entries == v.size());
                }
                // we need buckets 0,...,ell-1
            }
            else
            {
                ell = number_of_buckets;
                num_of_relevant_entries = v.size();
                delta = 0;
            }

            // Compute z = \sum_{p=0}^ell A^{jp}v_{[p]},
            // where A^{(jp)} is the optimal approximation matrix for each bin.
            // Optimality is meant in the sense that the computational effort is minimal and the error bound is met.
            // The computation is tailored for biorthogonal anisotropical wavelets:
            // A^{jp} has C*jp^dim nontrivial entries per row/column which need C*jp^dim many operations to be computed, further
            // \| A -A^{jp} \| \leq D 2^{-\rho jp} (A is compressible with s* = \infty)
            // jp is computed such that effort for apply is minimal, i.e.,
            // \sum _{p=0}^ell C*jp^dim * bin_size[p] -> minimal,
            // under the error bound condition
            // \sum_{p=0}^ell D*s^{-\rho jp} \|v_p\|_2 \leq eta-delta,
            // where rho < 1/2 is the compressibility of A (valid for biorthogonal anisotropical wavelets)
            // eta is the target accuracy of APPLY,
            // delta = \|A\|\|v - \sum_{p=1}^\ell v_p\|_2 (<= eta/2)
//            assert (delta <= threshold);
            delta = sqrt(delta) * norm_A;
            //if (!( (! (delta <= eta/2) ) == (delta > eta/2)))
            if (!(delta <= eta/2))
            {
                cout << "was soll das?" << endl;
                cout << "delta = " << delta << "; eta = " << eta << "; eta/2 = " << eta/2 << endl;
                cout << "(delta <= eta/2) = " << (delta <= eta/2) << "; (delta > eta/2) = " << (delta > eta/2) << endl;
                cout << "l2_norm_sqr(v) = " << l2_norm_sqr(v) << endl;
                cout << "bin_size = " << bin_size << endl;
                cout << "bin_norm_sqr = " << bin_norm_sqr << endl;
                cout << "ell = " << ell << endl;
            }
            assert (delta <= eta/2);
            // we use the formula for jp tilde from [DSS] for orthogonal wavelets, i.e., effort C*jp^1.
            // effect: error \eta-\delta is met, but effort may be suboptimal
            // => with D = P.alphak(p)
            //jp_tilde = ceil (log(sqrt(bin_norm_sqr[p])*num_relevant_entries*P.alphak(p)/(bin_size[p]*(eta-delta)))/M_LN2 );
            
            Array1D<int> jp_tilde(ell); // jp_tilde[i] == approximation that is needed for the ith bin
            for (unsigned int i = 0; i < ell; ++i)
            {
                //jp_tilde[i] = ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ); 
                if (bin_size[i] > 0)
                {
                    jp_tilde[i] = (int)ceil (log(sqrt(bin_norm_sqr[i]) * num_of_relevant_entries * P.alphak(i) / (bin_size[i]*(eta-delta))) / M_LN2 );
                    // jp_tilde[i] = (int)ceil (log2(ceil (sqrt(bin_norm_sqr[i]) * num_of_relevant_entries * P.alphak(i) / (bin_size[i]*(eta-delta))))  );
#if _APPLY_TENSOR_DEBUGMODE == 1
                    if (jp_tilde[i] < 0)
                    {
                        cout << "i = " << i << "; jp_tilde["<< i << "] = "  << jp_tilde[i] << endl;
                        cout << "bin_norm_sqr[" << i<< "] = " << bin_norm_sqr[i] << "; num_of_relevant_entries = " << num_of_relevant_entries << "; P.alphak(" << i<< ") = " << P.alphak(i) << "; bin_size[" << i<< "] = " << bin_size[i] << ";  eta = " << eta << "; delta = " << delta << endl;
                    }
                    assert(jp_tilde[i] >= 0);
                    if (jp_tilde[i] !=  ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ) )
                    {
                        cout << "i = " << i << "; jp_tilde["<< i << "] = "  << jp_tilde[i] << endl;
                        cout << "bin_norm_sqr[" << i<< "] = " << bin_norm_sqr[i] << "; num_of_relevant_entries = " << num_of_relevant_entries << "; P.alphak(" << i<< ") = " << P.alphak(i) << "; bin_size[" << i<< "] = " << bin_size[i] << ";  eta = " << eta << "; delta = " << delta << endl;
                        cout << "sqrt(bin_norm_sqr[" << i<< "])*num_of_relevant_entries*P.alphak(" << i<< ")/(bin_size[" << i <<" ]*(eta-delta)) = " << sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)) << endl;
                        cout << "ceil (log(" << sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)) << ")/M_LN2 ) = " << ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ) << endl;
                    }
                    assert (jp_tilde[i] ==  ceil (log(sqrt(bin_norm_sqr[i])*num_of_relevant_entries*P.alphak(i)/(bin_size[i]*(eta-delta)))/M_LN2 ) );
#endif
                }
                else
                {
                    jp_tilde[i] = 0;
                }
            }
            // hack: We work with full vectors (of size degrees_of_freedom).
            // We do this because adding sparse vectors seems to be inefficient.
            // Below we will then copy ww into the sparse vector w.
            // Probably this will be handled in a more elegant way someday.
            // ww should be of a similar structure to the cache in P
            Vector<double> ww(P.basis()->degrees_of_freedom());
            // compute w = \sum_{k=0}^(\ell-1) A_{J-k}v_{[k]}
            for (typename InfiniteVector<double,int>::const_iterator it(v.begin()), itend(v.end());it != itend; ++it)
            {
#if _APPLY_TENSOR_DEBUGMODE == 1
                if (log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) < 0)
                {
                    cout << "norm_v = " << norm_v << "; fabs(*it) = " << fabs(*it) << "; log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) = " << log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) << endl;
                }
                assert ( log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) >= 0 );
//                if (log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) >= q)
//                {
//                    cout << "ell = " << ell << endl;
//                    cout << "norm_v = " << norm_v << "; fabs(*it) = " << fabs(*it) << "; log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) = " << log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) << endl;
//                }
//                assert ( log2((unsigned int)(norm_v*norm_v/fabs(*it)/fabs(*it))) < q );
#endif
                //temp_i = 2*log2(floor(norm_v/fabs(*it))) -1
                //add_compressed_column_tensor(P, *it, it, jp_tilde[2*log2(floor(norm_v/fabs(*it))) -1], ww, jmax, strategy, preconditioning);
                temp_i = (unsigned int)floor(2*log(norm_v/fabs(*it))/M_LN2); // bins: 1,...,q but observe temp_i: 0,...,q-1
                //temp_i = log2((unsigned int)(floor(norm_v*norm_v/fabs(*it)/fabs(*it))));
                if ( temp_i < ell )
                {
                    P.add_ball(it.index() ,ww,jp_tilde[temp_i],*it   ,jmax,strategy,preconditioning);
                }
            }
            // copy ww into w
            for (unsigned int i = 0; i < ww.size(); i++) 
            {
                if (ww[i] != 0.)
                {
                    w.set_coefficient(i, ww[i]);
                }
            }
        }
    }    

}

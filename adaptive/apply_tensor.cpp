// implementation for apply_tensor.h

#include <utils/array1d.h>
#include <utils/tiny_tools.h>
#include <list>
#include <map>

#include <limits> // to use -inf as a literal

using MathTL::Array1D;

namespace WaveletTL
{
    /* For sorting */
    template<class T> struct index_cmp {
        index_cmp(const T arr) : arr(arr) {}
        bool operator()(const size_t a, const size_t b) const
        { return arr[a] > arr[b]; } // decreasing order
        const T arr;
    };

    template <class PROBLEM>
    void APPLY(const PROBLEM& P,
            const InfiniteVector<double, typename PROBLEM::Index>& v,
            const double eta,
            InfiniteVector<double, typename PROBLEM::Index>& w,
            const int jmax,
            const CompressionStrategy strategy)
    {
        Array1D<int> jp_guess(0);
        APPLY(P,v,eta,jp_guess,w,jmax,strategy);
    }

    template <class PROBLEM>
    void APPLY(const PROBLEM& P,
               const InfiniteVector<double, typename PROBLEM::Index>& v,
               const double eta,
               //const set<unsigned int> levelwindow,
               Array1D<int>& jp_guess,
               InfiniteVector<double, typename PROBLEM::Index>& w,
               const int jmax,
               const CompressionStrategy strategy)
    {
      typedef typename PROBLEM::Index Index;
      w.clear();
      // Binary Binning variant of APPLY from [DSS]
      // Remark: Remark from APPLY applies here as well, since binary binning part is the similar.
      // linfty norm is used, resulting in the Factor 2 for p.
      //cout << "entering apply_tensor..." << endl;
      //cout << "size = " << v.size() << endl;
      if (v.size() > 0) {
          // compute the number of bins V_0,...,V_q
          const double norm_v = linfty_norm(v);
          const double norm_A = P.norm_A();
          const unsigned int q = (unsigned int) std::max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.);
          // Setup the bins: The i-th bin contains the entries of v with modulus in the interval
          // (2^{-(i+1)/2}||v||,2^{-i/2}||v||], 0 <= i <= q-1, the remaining elements (with even smaller modulus)
          // are collected in the q-th bin.
          // q is chosen s.t. elements in the last bin have norm <= eta/(2*norm_A)
          Array1D<std::list<std::pair<Index, double> > > bins(q+1);
          {
              unsigned int i;
              for (typename InfiniteVector<double,Index>::const_iterator it(v.begin());it != v.end(); ++it)
              {

                  //const unsigned int i = std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2));
                  i = std::min(q, (unsigned int)floor(-2*log(fabs(*it)/norm_v)/M_LN2));
                  bins[i].push_back(std::make_pair(it.index(), *it));
              }
          }
          //cout << "done binning in apply..." << endl;

          // In difference to the isotropic APPLY we use the bins directly as the segments v_{[0]},...,v_{[\ell]}.
          // We are interested in \ell being the smallest number such that
          //   ||A||*||v-\sum_{k=0}^\ell v_{[k]}|| <= theta * eta
          // i.e.
          //   ||v-\sum_{k=0}^\ell v_{[k]}||^2 <= eta^2 * theta^2 / ||A||^2

          // const double theta = 0.5; // this explains the use of twotothejhalf from tinytools
          // const double threshold = eta*eta*theta*theta/(norm_A*norm_A);
          double threshold = eta/2.0/norm_A;
          threshold = threshold*threshold;
          double error_sqr = l2_norm_sqr(v);
          Array1D<double> bin_norm_sqr(q+1); //square norm of the coeffs in the bins
          unsigned int num_relevant_entries=0;
          unsigned int ell=0; // number of bins needed to approximate v
          if (l2_norm_sqr(v) > threshold)
          {
              unsigned int vsize = v.size();
              while (true)
              {
                  bin_norm_sqr[ell]=0;
                  for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[ell].begin()), itend(bins[ell].end()); it != itend; ++it)
                  {
                      bin_norm_sqr[ell] += (it->second) * (it->second);
                      num_relevant_entries++;
                  }
                  error_sqr -= bin_norm_sqr[ell];
                  if (error_sqr <= threshold)
                  {
                      break;
                  }
                  else if (num_relevant_entries == vsize)
                  {
                      error_sqr = 0;
                      break;
                  }
                  else
                  {
                      ell++;
                  }
              }
          }

          // Compute z = \sum_{p=0}^ell A^{(jp_tilde)}v_{[p]}
          // where jp_tilde is the solution of a minimization problem
          // (No projection to supp_s as in [DSS])
          // We assume the special case where
          // s^*=infty, alphak=D_const*2^{-rho*j}.
          // This is tailored for the tbasis basis construction.

          // the optimal effort for the pth bin is approximated by jptilde.
          // We have for p=1,...,q
          // 2^(-rho*jp)*D*||w_p|| / ((2jp-1)*size(w_p)) \leq (epsilon-delta) / (sum_{l=1}^q size(w_l)*(2j_l-1))
          //
          // We compute (jp_tilde)_{p=1}^q in a greedy way. The vector is initialized with jp_tilde from the L2orthogonal case as an initial guess.
          // If the initial error is to high the (jp) vector is succesisfly increased, choosing the most efficient entry in each step. (1/ (error_decrease*cost_increase) maximal)
          // If the initial error is already below epsilon-delta we decrease jp_tilde as often as possible.

          // In [DSS] the matrix has linear many nontrivial entries per radius J.
          // For a matrix stemming from a biorthogonal basis we have polynomial many.
          // This leads to the above minimization Problem for jp_tilde.
          // For a L2 orthogonal basis we would have
          // D_const = 1; // Value unknown, but -> infty for rho -> 1
          // rho = 1/2; // this is actually the limit value
          // D is given by P.alphak(k). Here a little hack (implicit assumtion) is used since alphak is konstant and also hardcoded

          const double delta = (error_sqr < 0)?(0):(norm_A*sqrt(error_sqr));

// CLEANUP
          assert (delta < eta);

          // =>
          // jp_tilde = ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);
          Array1D<int> jp_tilde(ell+1);
          Array1D<double> log_2_efficiency_sqr(ell+1);
          int jpmax(0);
          unsigned int jpmaxpos(0);

          //initialize

          if (P.space_dimension == 1 || jp_guess.size() == 0)
          {
              int int_temp;
              for (unsigned int k = 0; k<= ell; k++)
              {
                  if (bins[k].size() > 0)
                  {
                      int_temp = (int) ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);
                      jp_tilde[k] = (int_temp > 0)?int_temp:0;
                  }
                  else
                  {
                      jp_tilde[k] = std::numeric_limits<int>::quiet_NaN();
                  }
                  if (jp_tilde[k] > jpmax)
                  {
                      jpmax=jp_tilde[k];
                      jpmaxpos=k;
                  }
              }
          }
          else
          {
              unsigned int k;
              for (k = 0; k<= ((ell < jp_guess.size())?ell:jp_guess.size()-1); k++)
              {
                  if (bins[k].size() > 0)
                  {
                      jp_tilde[k] = jp_guess[k];
                      if (jp_tilde[k] > jpmax)
                      {
                          jpmax=jp_tilde[k];
                          jpmaxpos=k;
                      }
                  }
                  else
                  {
                      jp_tilde[k] = std::numeric_limits<int>::quiet_NaN();
                  }
              }
              for (; k<= ell; k++)
              {
                  if (bins[k].size() > 0)
                  {
                      jp_tilde[k] = ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);
                      if (jp_tilde[k] > jpmax)
                      {
                          jpmax=jp_tilde[k];
                          jpmaxpos=k;
                      }
                  }
                  else
                  {
                      jp_tilde[k] = std::numeric_limits<int>::quiet_NaN();
                  }
              }
          }
// CLEANUP
          if (jpmax > 70)
          {
              cout << "  jpmax is really high!" << endl;
              cout << jp_guess << " <-- jp_guess (initialization)"<< endl;
              cout << jp_tilde << " <-- jp_tilde after initialization"<< endl;
          }
          if (P.space_dimension > 1) // for dim=1 there is no minimization problem to be solved
          {
              std::vector<size_t> sort_index;
              // error = \sum_{p=1}^l 2^{-\rho*j_p}||w_p|| <= (eta-delta)  <=>  \sum 2^{\rho*(jpmax-jp)}*||w_p|| <= (eta-delta)*2^{-rho*jpmax}..
              double error(0);
              for (unsigned int k = 0; k<= ell; k++)
              {
                  //error +=  pow(2,(jpmax-jp_tilde[k])/2.0)*sqrt(bin_norm_sqr[k]);
                  if (bin_norm_sqr[k] > 0)
                  {
                      error +=  sqrt(bin_norm_sqr[k])* twotothejhalf(jpmax-jp_tilde[k]);
                      assert (jpmax - jp_tilde[k] < 60);
                  }
              }

              // to increase the convergence to the solution jp_tilde we introduce a step size N.
              // increase mode: find the most efficient bin and increase its jp value by N
              // if as a result we have error < eta-delta we switch to "decrease mode" and N/2, otherwise we continue with the now most efficient bin
              // decrease mode: search for the least efficient mode and decrease its jp value by N. switch to increase mode with N/2 if needed, otherwise continue decreasing.
              // break if N would drop below 1
              // Alternativly we could start with a stepsize of 1, jp_tilde[k] == 0 for all k and increase till the desired value is reached. Hopefully the code here is faster.

              unsigned int stepsize = 4, N; // the actual decrease in the current iteration may be less than stepsize since jp_tilde[k] cannot be reduced below 0
              double update;

              // compute efficiencies of buckets, i.e., error reduction / cost_increase if jp_tilde would be increased by 1
              for (unsigned int k=0; k<= ell; k++)
              {
                  sort_index.push_back(k);
                  if (bins[k].size() == 0)
                  {
                      log_2_efficiency_sqr[k] = -std::numeric_limits<double>::infinity();
                  }
                  else
                  {
                      if (P.space_dimension == 2)
                      {
                          update = bins[k].size() * (2*jp_tilde[k]+1);
                          update = update*update;
                          log_2_efficiency_sqr[k] = -jp_tilde[k]+log(bin_norm_sqr[k]/update)/M_LN2;
                      }
                      else
                      {
                          //update = (pow(jp_tilde[k]+1,P.space_dimension)-pow(jp_tilde[k],P.space_dimension))*bins[k].size();
                          update = (intpower(jp_tilde[k]+1,P.space_dimension)-intpower(jp_tilde[k],P.space_dimension))*bins[k].size();
                          update = update * update;
                          log_2_efficiency_sqr[k] = -jp_tilde[k]+log(bin_norm_sqr[k]/update)/M_LN2;
                      }
                  }
              }

              sort(sort_index.begin(), sort_index.end(), index_cmp<Array1D<double>&>(log_2_efficiency_sqr)); // entries are sorted in decreasing order

              bool increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 ); // we expect high values for jp_tilde, thats why I dont like to divide by (eps-delta)
              unsigned int pos;
              unsigned int lauf(0);
              int temp_int;
              while (true)
              {
                  if (increase_mode)
                  {
                      pos = sort_index[0]; // Position of most efficient entry
                      lauf = 0; // Position in sort_index
                      N = stepsize; // there is no problem if we increase jp_tilde
                  }
                  else
                  {
                      lauf = ell; // Position in sort_index
                      while (lauf > 0)
                      {
                          if ( (bins[sort_index[lauf]].size() == 0) || (jp_tilde[sort_index[lauf]] == 0) )
                          {
                              --lauf;
                          }
                          else
                          {
                              break;
                          }
                      }
                      // in case jp_tilde[i] = 0 for all i and bins[0] is empty: (would crash if all bins were empty what cannot happen)
                      if (lauf == 0)
                      {
                          while (bins[sort_index[lauf]].size() == 0)
                          {
                              lauf++;
                          }
                      }
                      pos = sort_index[lauf]; // Position of the least efficient entry
                      N = ( stepsize <jp_tilde[pos] ) ? stepsize:jp_tilde[pos]; // jp_tilde[k] may never become less than 0
                  }

                  // compute update for log_2_efficiency_sqr
                  if (P.space_dimension == 2)
                  {
                      update = (2*jp_tilde[pos]+1) / (double)(2*(jp_tilde[pos]+N)+1);

                  }
                  else if (P.space_dimension == 3)
                  {
                      update = (3*jp_tilde[pos]*jp_tilde[pos]+3*jp_tilde[pos]+1) / (double)(3*(jp_tilde[pos]+N)*(jp_tilde[pos]+N)+3*(jp_tilde[pos]+N)+1);

                  }
                  else
                  {
                      //update = (pow(jp_tilde[pos]+1,P.space_dimension)-pow(jp_tilde[pos],P.space_dimension)) / (pow(jp_tilde[pos]+N+1,P.space_dimension)-pow(jp_tilde[pos]+N,P.space_dimension));
                      update = (double)(intpower(jp_tilde[pos]+1,P.space_dimension)-intpower(jp_tilde[pos],P.space_dimension)) / (double)(intpower(jp_tilde[pos]+N+1,P.space_dimension)-intpower(jp_tilde[pos]+N,P.space_dimension));
                  }
                  //update = update*update*pow(2,-jp_tilde[pos]); // 2*rho = 1
                  update = update*update;
                  update = -(double)N+log(update)/M_LN2;

                  if (increase_mode)
                  {
                      jp_tilde[pos] = jp_tilde[pos] + N; // update value of jp_tilde
                      log_2_efficiency_sqr[pos] = log_2_efficiency_sqr[pos] + update; // efficiency after the increase/decrease of jp_tilde
                      if (jp_tilde[pos] > jpmax)
                      {
                          error = error *twotothejhalf(jp_tilde[pos]-jpmax);
                          jpmax=jp_tilde[pos];
                          jpmaxpos = pos;
                      }
                      error = error + sqrt(bin_norm_sqr[pos]) * (twotothejhalf(jpmax-jp_tilde[pos])-twotothejhalf(jpmax-(jp_tilde[pos]-(int)N)));
                      // if diffrence between lowest and highest jp_tilde value gets (or was in this case) to big, twotothejhalf will fail:
                      assert (jpmax-(jp_tilde[pos]-(int)N) < 60);
                  }
                  else
                  {
                      jp_tilde[pos] = jp_tilde[pos] - (int)N; // update value of jp_tilde
                      log_2_efficiency_sqr[pos] = log_2_efficiency_sqr[pos] - update; // efficiency after the increase/decrease of jp_tilde
                      error = error + sqrt(bin_norm_sqr[pos]) * (twotothejhalf(jpmax-jp_tilde[pos]) - twotothejhalf(jpmax-(jp_tilde[pos]+(int)N)));
                      // if diffrence between lowest and highest jp_tilde value gets to big, twotothejhalf will fail:
                      assert (jpmax-jp_tilde[pos] < 60);
                      if (pos == jpmaxpos)
                      {
                          jpmax = jp_tilde[pos];
                          // find the biggest jp_tilde
                          for (unsigned int k=0;k<=ell;k++)
                          {
                              if (jp_tilde[k] > jpmax)
                              {
                                  jpmax=jp_tilde[k];
                                  jpmaxpos=k;
                              }
                          }
                          error = error / twotothejhalf((jp_tilde[pos]+(int)N)-jpmax);
                      }
                  }
// CLEANUP
                  if (jpmax > 70)
                  {
                      cout << "jpmax is really high!" << endl;
                      cout << jp_guess << " <-- jp_guess (initialization)"<< endl;
                      cout << jp_tilde << " <-- jp_tilde (current)"<< endl;
                  }
                  // sort efficiency vector
                  if (increase_mode)
                  {
                      N = sort_index[lauf]; // observe the different meaning of N
                      while ( lauf < ell)
                      {
                          if (log_2_efficiency_sqr[N] < log_2_efficiency_sqr[sort_index[lauf+1]])
                          {
                              //sort_index.swap(lauf,lauf+1);
                              ++lauf;
                              sort_index[lauf-1] = sort_index[lauf];
                              sort_index[lauf] = N;
                          }
                          else
                          {
                              break;
                          }
                      }
                  }
                  else
                  {
                      N = sort_index[lauf]; // observe the different meaning of N
                      while ( lauf > 0)
                      {
                          if (log_2_efficiency_sqr[N] > log_2_efficiency_sqr[sort_index[lauf-1]])
                          {
                              //sort_index.swap(lauf-1,lauf);
                              --lauf;
                              sort_index[lauf+1] = sort_index[lauf];
                              sort_index[lauf] = N;
                          }
                          else
                          {
                              break;
                          }
                      }
                  }
                  // switch mode if necessary
                  if (increase_mode)
                  {
                      if (log(error)/M_LN2 < jpmax/2.0+log(eta-delta)/M_LN2 )
                      {
                          increase_mode = false;
                          stepsize = stepsize / 2;;
                          if (stepsize < 1)
                          {
                              break;
                          }
                      }
                  }
                  else
                  {
                      if (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 )
                      {
                          increase_mode = true;
                          if ((stepsize / 2) < 1)
                          {
                              jp_tilde[pos] = jp_tilde[pos] + stepsize; // undo the last step
                              break;
                          }
                          stepsize = stepsize / 2;
                      }
                  }
              } // end of while(true)
          } // end of if (P.space_dimension > 1)

          // storing jp_tilde for the next application of APPLY
          jp_guess = jp_tilde;
          
// CLEANUP
          assert (delta < eta);
          assert (!isnan(delta));

          //cout << "J = " << J << endl;
          // hack: We work with full vectors (of size degrees_of_freedom).
          // We do this because adding sparse vectors seems to be inefficient.
          // Below we will then copy ww into the sparse vector w.
          // Probably this will be handled in a more elegant way in the near future.
          Vector<double> ww(P.basis().degrees_of_freedom());
          //cout << *(P.basis().get_wavelet(4000)) << endl;
          // compute w = \sum_{k=0}^\ell A_{J-k}v_{[k]}

          for (unsigned int k=0; k<=ell; k++)
          {
              for (typename std::list<std::pair<Index, double> >::const_iterator bin_it(bins[k].begin()), bin_end_it(bins[k].end()); bin_it != bin_end_it; ++bin_it)
              {
                      add_compressed_column(P, bin_it->second, bin_it->first, jp_tilde[k], ww, jmax, strategy);
              }
          }
          // copy ww into w
          for (unsigned int i = 0; i < ww.size(); i++) {
              if (ww[i] != 0.)
              {
                  w.set_coefficient(*(P.basis().get_wavelet(i)), ww[i]);
              }
          }
      }
    }

}

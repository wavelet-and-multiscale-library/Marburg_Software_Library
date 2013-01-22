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

    template<class T> struct reverse_cmp {
        reverse_cmp(const T arr) : arr(arr) {}
        bool operator()(const size_t a, const size_t b) const
        { return arr[a] < arr[b]; } // increasing order
        const T arr;
    };

    template <class PROBLEM>
    void APPLY(const PROBLEM& P,
            const InfiniteVector<double, typename PROBLEM::Index>& v,
            const double eta,
            InfiniteVector<double, typename PROBLEM::Index>& w,
            const int jmax,
            const CompressionStrategy strategy,
            const bool preconditioning)
    {
        Array1D<int> jp_guess(0);
        APPLY(P,v,eta,jp_guess,w,jmax,strategy, preconditioning);
    }

    template <class PROBLEM>
    void APPLY(const PROBLEM& P,
               const InfiniteVector<double, typename PROBLEM::Index>& v,
               const double eta,
               //const set<unsigned int> levelwindow,
               Array1D<int>& jp_guess,
               InfiniteVector<double, typename PROBLEM::Index>& w,
               const int jmax,
               const CompressionStrategy strategy,
               const bool preconditioning)
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
// CLEANUP
          if (norm_v > 1e+10)
          {
              cout << "apply_tensor:: l_infty norm of v is too high!"<<endl;
              cout << "v = " << v << endl;
          }


          const double norm_A = P.norm_A();
          const unsigned int q = (unsigned int) std::max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.);
/*
          cout << "APPLY called with v = " << endl << v << endl;
          cout << "compute q." << endl << " v.size() = " << v.size() << " norm_A = " << norm_A << " eta = " << eta << endl;
          cout << "(sqrt((double)v.size())*norm_v*norm_A*2.0/eta) = " << (sqrt((double)v.size())*norm_v*norm_A*2.0/eta) << endl;
          cout << "log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2 = " << log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2 << endl;
          cout << "2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) = " << 2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2) << endl;
          cout << "max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.) = " << std::max(2*ceil(log(sqrt((double)v.size())*norm_v*norm_A*2.0/eta)/M_LN2), 0.) << endl;
*/
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
          //double error_sqr = l2_norm_sqr(v);
          double bound = l2_norm_sqr(v) - threshold;
          double new_norm_sqr(0);

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
                      //num_relevant_entries++;
                  }
                  if (bin_norm_sqr[ell] < 1e-15 && bins[ell].size() > 0) // this bin contains only very small coefs. Thus it and all later bins are not of any importance
                  {
                      //error_sqr = threshold; // we are near the machine limit
                      if (ell > 0)
                      {
                          --ell;
                      }
                      break;
                  }
                  num_relevant_entries += bins[ell].size();
                  new_norm_sqr += bin_norm_sqr[ell];
                  //error_sqr -= bin_norm_sqr[ell];
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
                      ell++;
                  }
              }
          }
          else
          {
              // v is well approximated by 0
              return;
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
          // D is given by P.alphak(k). Here a little hack (implicit assumtion) is used since alphak is constant and also hardcoded

          const double delta = ((l2_norm_sqr(v)-bound) < 0)?(0):(norm_A*sqrt((l2_norm_sqr(v)-bound)));

// CLEANUP
          assert (delta < eta);

          // =>
          // jp_tilde = ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);
          Array1D<int> jp_tilde(ell+1);


#if _SKIP_MINIMIZATION_PROBLEM == 1
          for (unsigned int i = 0; i <= ell; ++i)
          {
              jp_tilde[i] = jmax;
          }
#else

          //Array1D<double> log_2_efficiency_sqr(ell+1);
          //int pos1, pos2;

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
              int int_temp;
              for (k = 0; k<= ((ell < jp_guess.size())?ell:jp_guess.size()-1); k++)
              {
                  if (bins[k].size() > 0)
                  {
                      int_temp = jp_guess[k];
                      if (int_temp == 0) // a guess of 0 is considered bad. This is because jmax-jmin may become very big in this case
                      {
                          int_temp = (int) ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);
                          int_temp = (int_temp > 0)?int_temp:0;
                      }
                      jp_tilde[k] = int_temp;
                      if (int_temp > jpmax)
                      {
                          jpmax=int_temp;
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
#if _NO_JPMAX_grt_70_WARNING
          //cout << "JPMAX_WARNING_A ";
#else
          if (jpmax > 70)
          {
              //cout << "  jpmax is really high!" << endl;
              cout << jp_guess << " <-- jp_guess (initialization)"<< endl;
              cout << jp_tilde << " <-- jp_tilde after initialization"<< endl;
              cout << "i jp_tilde bins[i].size() bin_norm_sqr[i]" << endl;
              for (unsigned int i=0; i<=ell ; ++i)
              {
                  cout << i << " " << jp_tilde[i] << " " << bins[i].size() << " " << bin_norm_sqr[i] << endl;
              }
              for (unsigned int i=ell+1; i< bins.size(); ++i)
              {
                  cout << i << " * " << bins[i].size() << " " << bin_norm_sqr[i] << endl;
              }
          }
#endif
          if (P.space_dimension > 1) // for dim=1 there is no minimization problem to be solved
          {
              // estimate current error of the approximation.

              std::vector<size_t> sort_index_inc_mode, sort_index_dec_mode;
              // error = \sum_{p=1}^l D*2^{-\rho*j_p}||w_p|| <= (eta-delta)  <=>  \sum 2^{\rho*(jpmax-jp)}*||w_p|| <= 1/D * (eta-delta)*2^{rho*jpmax}
              double error(0);
              for (unsigned int k = 0; k<= ell; k++)
              {
                  //error +=  pow(2,(jpmax-jp_tilde[k])/2.0)*sqrt(bin_norm_sqr[k]); // choose rho = 1/2
                  if (bin_norm_sqr[k] > 0)
                  {
                      if ((jpmax-jp_tilde[k]) >= 60)
                      {
                          cout << "craisy high level differences detected Line 292" << endl;
                          error += pow(2,(jpmax-jp_tilde[k])/2.0)*sqrt(bin_norm_sqr[k]);
                      }
                      else
                      {
                          error +=  sqrt(bin_norm_sqr[k])* twotothejhalf(jpmax-jp_tilde[k]);
                      }
// CLEANUP
#if 0
                      if (jpmax - jp_tilde[k] >= 60)
                      {
                          cout << "apply_tensor:: jpmax - jp_tilde[k] > 60. Some debugging output:" << endl;
                          cout << "apply_tensor:: jp_tilde = " << endl << jp_tilde << endl
                               << "apply_tensor:: bin_norm.size() = " << endl << "[";
                          for (unsigned int i = 0; i<bins.size()-1; ++i)
                          {
                              cout << bins[i].size() << " ";
                          }
                          cout << bins[bins.size()-1].size() << "]" << endl;
                          //cout << "apply_tensor:: bin_norm_sqr = " << endl << bin_norm_sqr << endl;

                          cout << "i jp_tilde log2_efficiency_sqr bins[i].size() bin_norm_sqr[i]" << endl;
                          for (unsigned int i=0; i<=ell ; ++i)
                          {
                              cout << i << " " << jp_tilde[i] << " " << log_2_efficiency_sqr[i] << " " << bins[i].size() << " " << bin_norm_sqr[i] << endl;
                          }
                          for (unsigned int i=ell+1; i< bins.size(); ++i)
                          {
                              cout << i << " * * " << bins[i].size() << " " << bin_norm_sqr[i] << endl;
                          }
                          cout << "stored coeficients"<< endl;
                          for (unsigned int i=0; i < bins[i].size(); i++)
                          {
                              cout << "bins["<<i<<"] = { ";
                              for (typename std::list<std::pair<Index, double> >::const_iterator it(bins[i].begin()), itend(bins[i].end()); it != itend; ++it)
                              {
                                  cout << (it->second) << " ";
                              }
                              cout << "}" << endl;
                          }
                          cout << "jpmax = " << jpmax << " k = " << k << " ell = " << ell << endl;
                          cout << "computation of some jp_tilde values" << endl;
                          cout << "jp_tilde[k] = ceil (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2);" << endl;
                          unsigned int int_temp = 0;
                          cout << " k = "<< int_temp << " bin_norm_sqr[k] = " << bin_norm_sqr[int_temp] << " num_relevant_entries = " << num_relevant_entries << " P.alphak(k) = " << P.alphak(int_temp) << " bins[k].size() = " << bins[int_temp].size() << " (eta-delta) = " << (eta-delta) << endl;
                          cout << " (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2) = " << (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          cout << " jp_tilde[k] = " << ceil (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          int_temp = 39;
                          cout << " k = "<< int_temp << " bin_norm_sqr[k] = " << bin_norm_sqr[int_temp] << " num_relevant_entries = " << num_relevant_entries << " P.alphak(k) = " << P.alphak(int_temp) << " bins[k].size() = " << bins[int_temp].size() << " (eta-delta) = " << (eta-delta) << endl;
                          cout << " (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2) = " << (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          cout << " jp_tilde[k] = " << ceil (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          int_temp = 46;
                          cout << " k = "<< int_temp << " bin_norm_sqr[k] = " << bin_norm_sqr[int_temp] << " num_relevant_entries = " << num_relevant_entries << " P.alphak(k) = " << P.alphak(int_temp) << " bins[k].size() = " << bins[int_temp].size() << " (eta-delta) = " << (eta-delta) << endl;
                          cout << " (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2) = " << (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          cout << " jp_tilde[k] = " << ceil (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          int_temp = 47;
                          cout << " k = "<< int_temp << " bin_norm_sqr[k] = " << bin_norm_sqr[int_temp] << " num_relevant_entries = " << num_relevant_entries << " P.alphak(k) = " << P.alphak(int_temp) << " bins[k].size() = " << bins[int_temp].size() << " (eta-delta) = " << (eta-delta) << endl;
                          cout << " (log(sqrt(bin_norm_sqr[k])*num_relevant_entries*P.alphak(k)/(bins[k].size()*(eta-delta)))/M_LN2*2) = " << (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                          cout << " jp_tilde[k] = " << ceil (log(sqrt(bin_norm_sqr[int_temp])*num_relevant_entries*P.alphak(int_temp)/(bins[int_temp].size()*(eta-delta)))/M_LN2*2) << endl;
                      }
#endif
                      assert (jpmax - jp_tilde[k] < 60);
                  }
              }

              assert (error > 0);
              // assert error <= 2^{rho*jpmax}*(eps-delta)
              // the minimization problem is solved in a greedy way:
              // if initial guess has error <= eps-delta start in decrease mode, else in increase mode
              // in each step: find most efficient bucket and increase or decrease jp.
              // continue until this bucket is no longer the most efficient one or until the error passes over the point (eps-delte)
              // in that case we have to switch the mode and try to find a possible bucket to continue (because the initial guess may be stupid)
              // if the best bucket to decrease would drop the error above (eps-delta) stop.
              // if the best bucket to decrease is the same that we just increased (in the step before) stop. (loop detection)

              // cost to increase jp->jp+1 :: C*#w_P*(2jp+1)
              // error reduction  jp->jp+1 :: D||w_p||_2*2^{-rho*jp}(2^{-rho}-1)
              // jp-> jp+1::
              // log_2 (-(error(jp+1)-error(jp)) / (cost(jp+1)-cost(jp))) = const -rho*jp-log_2(2jp+1)
              // jp+1 ->jp
              // log_2 (-(error(jp)-error(jp-1)) / (cost(jp)-cost(jp-1))) = const -rho*(jp-1)-log_2(2(jp-1)+1) //with the same constant
              // replace 2jp+1 with (jp+1)^DIM-jp^DIM if needed

              Array1D<double> log_2_inc_efficiency(ell+1);
              Array1D<double> log_2_dec_efficiency(ell+1);
              Array1D<double> relative_importance(ell+1);
              bool increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 ); // we expect high values for jp_tilde, thats why I dont like to divide by (eps-delta)
//CLEANUP
              //cout << "rel_importance log2_inc_eff log2_dec_eff" << endl;
              for (unsigned int k=0; k<=ell; k++)
              {
                  //initialize log_2 efficiency vector
                  sort_index_inc_mode.push_back(k);
                  sort_index_dec_mode.push_back(k);
                  if (bins[k].size() == 0)
                  {
                      relative_importance[k] = 0;
                      log_2_inc_efficiency[k] = -std::numeric_limits<double>::infinity();
                      log_2_dec_efficiency[k] = std::numeric_limits<double>::infinity();
//CLEANUP
                      //cout << relative_importance[k] << "  " << log_2_inc_efficiency[k] << "  " << log_2_dec_efficiency[k] << endl;
                  }
                  else
                  {
                      relative_importance[k] = 0.5*log(bin_norm_sqr[k]/(bins[k].size()*bins[k].size()))/M_LN2;
                      if (P.space_dimension == 2)
                      {
                          log_2_inc_efficiency[k] = relative_importance[k] -0.5*jp_tilde[k]-log(2*jp_tilde[k]+1)/M_LN2;
                          // since it is impossible to decrease jp = 0.
                          // in decrease_mode we search for the bucket with the highest value of - delta (error/cost)
                          // thats why this operation gets the efficiency std::numeric_limits<double>::infinity()
                          if (jp_tilde[k] == 0)
                          {
                              log_2_dec_efficiency[k] = std::numeric_limits<double>::infinity();
                          }
                          else
                          {
                              log_2_dec_efficiency[k] = relative_importance[k] -0.5*(jp_tilde[k]-1)-log(2*jp_tilde[k]-1)/M_LN2;
                          }
//CLEANUP
                          //cout << relative_importance[k] << "  " << log_2_inc_efficiency[k] << "  " << log_2_dec_efficiency[k] << endl;
                      }
                      else
                      {
                          log_2_inc_efficiency[k] = relative_importance[k] -0.5*jp_tilde[k]-log(intpower(jp_tilde[k]+1,P.space_dimension)-intpower(jp_tilde[k],P.space_dimension))/M_LN2;
                          // since it is impossible to decrease jp = 0.
                          // in decrease_mode we search for the bucket with the highest value of - delta (error/cost)
                          // thats why this operation gets the efficiency std::numeric_limits<double>::infinity()
                          if (jp_tilde[k] == 0)
                          {
                              log_2_dec_efficiency[k] = std::numeric_limits<double>::infinity();
                          }
                          else
                          {
                              log_2_dec_efficiency[k] = relative_importance[k] -0.5*(jp_tilde[k]-1)-log(intpower(jp_tilde[k],P.space_dimension)-intpower(jp_tilde[k]-1,P.space_dimension))/M_LN2;
                          }
// CLEANUP
                          //cout << relative_importance[k] << "  " << log_2_inc_efficiency[k] << "  " << log_2_dec_efficiency[k] << endl;
                      }
                  }
              }

              //unsigned int stepsize = 4, N; // the actual decrease in the current iteration may be less than stepsize since jp_tilde[k] cannot be reduced below 0
              //double update;

              // compute efficiencies of buckets, i.e., error reduction / cost_increase if jp_tilde would be increased by 1
              /*
              for (unsigned int k=0; k<= ell; k++)
              {
                  sort_index.push_back(k);
                  if (bins[k].size() == 0)
                  {
                      //log_2_efficiency_sqr[k] = -std::numeric_limits<double>::infinity();
                      log_2_efficiency[k] = -std::numeric_limits<double>::infinity();
                  }
                  else
                  {
                      if (P.space_dimension == 2)
                      {
                          update = bins[k].size() * (2*jp_tilde[k]+1);
                          update = update*update;
                          //log_2_efficiency_sqr[k] = -jp_tilde[k]+log(bin_norm_sqr[k]/update)/M_LN2;
                          relative_importance[k] = 0.5*log_2(bin_norm_sqr[k]/(bins[k].size()*bins[k].size()));
                          log_2_efficiency[k] = relative_importance[k] -0.5*jp_tilde[k]+log(bin_norm_sqr[k]/update)/M_LN2;
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
              */

              // sort efficiency vectors to determine where to start increasesing/decreasing jp_tilde (the first entry will hold the interesting bucket)
              sort(sort_index_inc_mode.begin(), sort_index_inc_mode.end(), index_cmp<Array1D<double>&>(log_2_inc_efficiency)); // entries are sorted in decreasing order
              sort(sort_index_dec_mode.begin(), sort_index_dec_mode.end(), reverse_cmp<Array1D<double>&>(log_2_dec_efficiency)); // entries are sorted in increasing order

              //bool increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 ); // we expect high values for jp_tilde, thats why I dont like to divide by (eps-delta)
              unsigned int last_altered_pos = ell+1; // initialized with an unreachable position
              unsigned int pos;
              bool done = false;
              bool mode_switch = false;
              // loops are still possible! (inc at 0, dec at 1, dec at 0, inc at 1 , ...)
              // I assume that jp_tilde values will only have two values in a circle, this makes for a maximum loop length of 2^(ell+1). I use the number 2ell+1 instead
              int change_counter(0), max_number_of_switches(2*ell+1);
              while (!done)
              {
// CLEANUP
                  /*
                  cout << "error = " << error << endl;
                  cout << "increase_mode == " << increase_mode << " log2(error) = " << log(error)/M_LN2 << " jpmax/2.0+log2(eta-delta) = " << jpmax/2.0+log(eta-delta)/M_LN2 << endl;
                  //increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 );
                  cout << "jp_tilde = " << endl << jp_tilde << endl;
                  cout << "log_2_inc_efficiency = " << endl << log_2_inc_efficiency << endl;
                  //cout << "sort_index_inc_mode = " << sort_index_inc_mode << endl;
                  cout << "[ ";
                  for (unsigned int k=0; k <= ell; ++k)
                  {
                      cout << sort_index_inc_mode[k] << " ";
                  }
                  cout << "]" << endl;
                  cout << "log_2_dec_efficiency = " << endl << log_2_dec_efficiency << endl;
                  //cout << "sort_index_dec_mode = " << sort_index_dec_mode << endl;
                  cout << "[ ";
                  for (unsigned int k=0; k <= ell; ++k)
                  {
                      cout << sort_index_dec_mode[k] << " ";
                  }
                  cout << "]" << endl;
                   */

                  pos = (increase_mode? sort_index_inc_mode[0]:sort_index_dec_mode[0]); // position to alter jp_tilde
                  if (increase_mode)
                  {
                      //update jp_tilde
                      jp_tilde[pos] = jp_tilde[pos] + 1; // update value of jp_tilde
                      if ( ((last_altered_pos == pos) && mode_switch) || (change_counter >= max_number_of_switches ) )
                      {
                          // we have switched the mode, but not the bucket, thus we are finished after we have increased jp_tilde a final time
                          done = true;
                      }
                      last_altered_pos = pos;
                      // update efficiency vectors
                      log_2_dec_efficiency[pos] = log_2_inc_efficiency[pos]; 
                      log_2_inc_efficiency[pos] = relative_importance[pos] -0.5*jp_tilde[pos] + ((P.space_dimension == 2)? (-log(2*jp_tilde[pos]+1)/M_LN2) : -log(intpower(jp_tilde[pos]+1,P.space_dimension)-intpower(jp_tilde[pos],P.space_dimension))/M_LN2);
                      // update error
                      if (jp_tilde[pos] > jpmax)
                      {
                          error = error * twotothejhalf(jp_tilde[pos]-jpmax);
                          jpmax=jp_tilde[pos];
                          jpmaxpos = pos;
                      }                      
                      error = error + sqrt(bin_norm_sqr[pos]) * (twotothejhalf(jpmax-jp_tilde[pos])-twotothejhalf(jpmax-(jp_tilde[pos]-1)));
                      // if diffrence between lowest and highest jp_tilde value gets (or was in this case) to big, twotothejhalf will fail:
                      assert (jpmax-(jp_tilde[pos]-1) < 60);
                      // decide what to do next:
                      increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 );
                      if (increase_mode)
                      {
                          mode_switch = false;
                          if (ell >= 1)
                          {
                              if (log_2_inc_efficiency[sort_index_inc_mode[0]] < log_2_inc_efficiency[sort_index_inc_mode[1]])
                              {
                                  // update sort vectors
                                  // this may done in a more clever way, since only the first element is (possibly) not in order
                                  sort(sort_index_inc_mode.begin(), sort_index_inc_mode.end(), index_cmp<Array1D<double>&>(log_2_inc_efficiency)); // entries are sorted in decreasing order
                                  sort(sort_index_dec_mode.begin(), sort_index_dec_mode.end(), reverse_cmp<Array1D<double>&>(log_2_dec_efficiency)); // entries are sorted in increasing order
                              }
                          }
                      }
                      else
                      {
                          // we are switching the mode!
                          mode_switch = true; change_counter++;
                          if (ell >= 1)
                          {
                                  // update sort vectors
                                  // this may done in a more clever way, since only the first element is (possibly) not in order
                              if ( ((pos == sort_index_dec_mode[0]) && (log_2_dec_efficiency[sort_index_inc_mode[0]] > log_2_dec_efficiency[sort_index_inc_mode[1]]) )
                                      || (jp_tilde[sort_index_dec_mode[0]] == 0) ) // in the second case we have altered an element somewhere in the vector AND since the first one wont be decreased, we end up with at the very same position
                              {
                                  sort(sort_index_inc_mode.begin(), sort_index_inc_mode.end(), index_cmp<Array1D<double>&>(log_2_inc_efficiency)); // entries are sorted in decreasing order
                                  sort(sort_index_dec_mode.begin(), sort_index_dec_mode.end(), reverse_cmp<Array1D<double>&>(log_2_dec_efficiency)); // entries are sorted in increasing order
                              }
                          }
                      }
                  }
                  else
                  {
                      if ((last_altered_pos == pos) && mode_switch)
                      {
                          // we have switched the mode, but not the bucket, thus we are finished!
                          done = true;
                      }
                      else
                      {
                          int position_in_sort = 0; // if this is a legal decreasing position, ie jp_tilde[0] > 0, sort_index_dec_mode[0] == pos will be valid. if jp_tilde[pos] == 0 we search and store here.
                          //update jp_tilde
                          if (jp_tilde[pos] > 0)
                          {
                              jp_tilde[pos] = jp_tilde[pos] - 1; // update value of jp_tilde
                          }
                          else
                          {
                              // search for bucket to decrease jp_tilde
                              bool temp_bool = false;
                              for (unsigned int k=1;k<=ell;++k)
                              {
                                  if (jp_tilde[sort_index_dec_mode[k]] > 0)
                                  {
                                      position_in_sort = k;
                                      pos = sort_index_dec_mode[k];
                                      jp_tilde[pos] = jp_tilde[pos] - 1;
                                      temp_bool = true;
                                  }
                              }
                              if (temp_bool == false)
                              {
                                  // all jp_tilde are minimal and error is still low enough? then we are done!
                                  done = true;
                              }
                          }
                          if (!done)
                          {
                              last_altered_pos = pos;
                              // update efficiency vectors
                              log_2_inc_efficiency[pos] = log_2_dec_efficiency[pos];
                              if (jp_tilde[pos] == 0)
                              {
                                  log_2_dec_efficiency[pos] = std::numeric_limits<double>::infinity();
                              }
                              else
                              {
                                  log_2_dec_efficiency[pos] = relative_importance[pos] -0.5*(jp_tilde[pos]-1) 
                                                              + ((P.space_dimension == 2)? (-log(2*jp_tilde[pos]-1)/M_LN2)
                                                                                         : (-log(intpower(jp_tilde[pos],P.space_dimension)-intpower(jp_tilde[pos]-1,P.space_dimension))/M_LN2)
                                                                );
                              }
                              // update error
                              error = error + sqrt(bin_norm_sqr[pos]) * (twotothejhalf(jpmax-jp_tilde[pos]) - twotothejhalf(jpmax-(jp_tilde[pos]+1)));
                              // if difference between lowest and highest jp_tilde value gets to big, twotothejhalf will fail:
                              assert (jpmax-jp_tilde[pos] < 60);
                              //adjust jpmax
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
                                  if (jpmaxpos == pos)
                                  {
                                      // jpmax has decreased by 1
                                      error = error / sqrt(2);
                                  }
                              }
                              
                              // decide what to do next:
                              increase_mode = (log(error)/M_LN2 > jpmax/2.0+log(eta-delta)/M_LN2 );
                              if (increase_mode)
                              {
                                  // we are switching the mode!
                                  mode_switch = true; change_counter++;
                                  if (ell >= 1)
                                  {
                                      if  ( (pos != sort_index_inc_mode[0])
                                              && ( (pos != sort_index_inc_mode[1])
                                                  || ( (pos == sort_index_inc_mode[1]) && (log_2_inc_efficiency[sort_index_inc_mode[0]] < log_2_inc_efficiency [sort_index_inc_mode[1]]) )
                                                 )
                                          )
                                      {
                                          // update sort vectors
                                          // this may done in a more clever way, since only the first element is (possibly) not in order
                                          sort(sort_index_inc_mode.begin(), sort_index_inc_mode.end(), index_cmp<Array1D<double>&>(log_2_inc_efficiency)); // entries are sorted in decreasing order
                                          sort(sort_index_dec_mode.begin(), sort_index_dec_mode.end(), reverse_cmp<Array1D<double>&>(log_2_dec_efficiency)); // entries are sorted in increasing order
                                      }
                                  }
                              }
                              else
                              {
                                  mode_switch = false;
                                  if (ell >= 1)
                                  {
                                      if ( (position_in_sort == ell)
                                           || (jp_tilde[sort_index_dec_mode[position_in_sort+1]] == 0)
                                           || (log_2_dec_efficiency[pos] > log_2_dec_efficiency[sort_index_inc_mode[position_in_sort+1]]) )
                                      {
                                          // update sort vectors
                                          // this may done in a more clever way, since only the first element is (possibly) not in order
                                          sort(sort_index_inc_mode.begin(), sort_index_inc_mode.end(), index_cmp<Array1D<double>&>(log_2_inc_efficiency)); // entries are sorted in decreasing order
                                          sort(sort_index_dec_mode.begin(), sort_index_dec_mode.end(), reverse_cmp<Array1D<double>&>(log_2_dec_efficiency)); // entries are sorted in increasing order
                                      }
                                  }
                              }
                          }
                      }
                  }
              }
#if 0 // old code with unterminated loop (sometimes)
              unsigned int lauf(0);
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
#endif
#if _NO_JPMAX_grt_70_WARNING
                  //cout << "JPMAX_WARNING_B ";
#else
                  if (jpmax > 70)
                  {
                      cout << "jpmax is really high!" << endl;
                      cout << jp_guess << " <-- jp_guess (initialization)"<< endl;
                      cout << jp_tilde << " <-- jp_tilde (current)"<< endl;
                  }
#endif
          } // end of if (P.space_dimension > 1)
#endif // end of _SKIP_MINIMIZATION_PROBLEM
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
                      add_compressed_column(P, bin_it->second, bin_it->first, jp_tilde[k], ww, jmax, strategy, preconditioning);
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

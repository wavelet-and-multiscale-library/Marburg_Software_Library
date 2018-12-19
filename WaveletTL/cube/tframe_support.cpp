// implementation for tframe_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>
#include <iostream>
#include <cube/tframe_index.h>
#include <interval/pq_support.h>

using MathTL::multi_degree;
using MathTL::FixedArray1D;

namespace WaveletTL
{
    template <class IFRAME, unsigned int DIM>
    inline
    void
    support(const TensorFrame<IFRAME,DIM>& frame,
            const typename TensorFrame<IFRAME,DIM>::Index& lambda,
            typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        frame.support(lambda, supp);
    }
    
    template <class IFRAME, unsigned int DIM>
    inline
    void
    support(const TensorFrame<IFRAME,DIM>& frame,
            const int& lambda_num,
            typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        frame.support(lambda_num, supp);
    }
    
    template <class IFRAME, unsigned int DIM>
    bool
    intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
                       const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                       const typename TensorFrame<IFRAME,DIM>::Index& mu)
    {
        typename TensorFrame<IFRAME,DIM>::Support supp_lambda;
        frame.support(lambda, supp_lambda);
        typename TensorFrame<IFRAME,DIM>::Support supp_mu;
        frame.support(mu, supp_mu);
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        double a=0;
        double b=0;
        for (unsigned int i=0;i<DIM;i++)
        {
            if (supp_lambda.j[i] > supp_mu.j[i]) {
//                supp.j[i] = supp_lambda.j[i];
                const int adjust = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
                supp_mu.a[i] *= adjust;
                supp_mu.b[i] *= adjust;
            } else {
//                supp.j[i] = supp_mu.j[i];
                const int adjust = 1<<(supp_mu.j[i]-supp_lambda.j[i]);
                supp_lambda.a[i] *= adjust;
                supp_lambda.b[i] *= adjust;
            }
            a = std::max(supp_lambda.a[i],supp_mu.a[i]);
            b = std::min(supp_lambda.b[i],supp_mu.b[i]);
            if (a >= b)
                return false;
        }
        return true;
    }

    template <class IFRAME, unsigned int DIM>
    bool
    intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
                       const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                       const typename TensorFrame<IFRAME,DIM>::Support& supp_mu,
                       typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        typename TensorFrame<IFRAME,DIM>::Support supp_lambda;
        frame.support(lambda, supp_lambda);
//        typename TensorFrame<IFRAME,DIM>::Support supp_mu;
//        frame.support(mu, supp_mu);
//        typename TensorFrame<IFRAME,DIM>::Support supp_mu(supp2);       //todo: delete copy constructor
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        for (unsigned int i=0;i<DIM;i++)
        {
            int adjustlambda=1;
            int adjustmu=1;
            if (supp_lambda.j[i] > supp_mu.j[i]) {
                supp.j[i] = supp_lambda.j[i];
                adjustmu = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
//                supp_mu.a[i] *= adjust;
//                supp_mu.b[i] *= adjust;
            } else {
                supp.j[i] = supp_mu.j[i];
                adjustlambda = 1<<(supp_mu.j[i]-supp_lambda.j[i]);
//                supp_lambda.a[i] *= adjust;
//                supp_lambda.b[i] *= adjust;
            }
            supp.a[i] = std::max(supp_lambda.a[i]*adjustlambda,supp_mu.a[i]*adjustmu);
            supp.b[i] = std::min(supp_lambda.b[i]*adjustlambda,supp_mu.b[i]*adjustmu);
            if (supp.a[i] >= supp.b[i])
                return false;
        }
        return true;
    }
    
    template <class IFRAME, unsigned int DIM>
    bool
    intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
                       const int& lambda_num,
                       const typename TensorFrame<IFRAME,DIM>::Support& supp_mu,
                       typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        typename TensorFrame<IFRAME,DIM>::Support supp_lambda;
        frame.support(lambda_num, supp_lambda);
//        typename TensorFrame<IFRAME,DIM>::Support supp_mu;
//        frame.support(mu, supp_mu);
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        for (unsigned int i=0;i<DIM;i++)
        {
            int adjustlambda=1;
            int adjustmu=1;
            if (supp_lambda.j[i] > supp_mu.j[i]) {
                supp.j[i] = supp_lambda.j[i];
                adjustmu = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
//                supp_mu.a[i] *= adjust;
//                supp_mu.b[i] *= adjust;
            } else {
                supp.j[i] = supp_mu.j[i];
                adjustlambda = 1<<(supp_mu.j[i]-supp_lambda.j[i]);
//                supp_lambda.a[i] *= adjust;
//                supp_lambda.b[i] *= adjust;
            }
            supp.a[i] = std::max(supp_lambda.a[i]*adjustlambda,supp_mu.a[i]*adjustmu);
            supp.b[i] = std::min(supp_lambda.b[i]*adjustlambda,supp_mu.b[i]*adjustmu);
            if (supp.a[i] >= supp.b[i])
                return false;
        }
        return true;
    }
    
    template <class IFRAME, unsigned int DIM>
    inline
    bool intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
			  const typename TensorFrame<IFRAME,DIM>::Index& lambda1,
			  const typename TensorFrame<IFRAME,DIM>::Index& lambda2,
			  typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        typedef typename TensorFrame<IFRAME,DIM>::Support Support;

        Support supp2(frame.get_support(lambda2.number()));
//    frame.support(lambda2, supp2);

        return intersect_supports(frame, lambda1, supp2, supp);
    }
    
    template <class IFRAME, unsigned int DIM>
    inline
    bool intersect_supports(const TensorFrame<IFRAME,DIM>& frame,
			  const int& lambda_num,
			  const int& mu_num,
			  typename TensorFrame<IFRAME,DIM>::Support& supp)
    {
        typedef typename TensorFrame<IFRAME,DIM>::Support Support;

        Support supp2(frame.get_support(mu_num));
    

        return intersect_supports(frame, lambda_num, supp2, supp);
    }
  


    template <class IFRAME, unsigned int DIM>
    void intersecting_quarklets(const TensorFrame<IFRAME,DIM>& frame,
                               const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                               const typename TensorFrameIndex<IFRAME,DIM>::level_type& j,
                               std::list<int>& intersecting,
                               const MultiIndex<int,DIM> p)
    {
        typedef typename TensorFrame<IFRAME,DIM>::Index Index;
        typedef typename TensorFrame<IFRAME,DIM>::Support Support;
        intersecting.clear();
        
        Support supp, supp_lambda(frame.get_support(lambda.number()));
//        cout<<"bin hier"<<endl;

    // a brute force solution
    if (j==frame.j0()) {
        
      Index last_qua(last_quarklet<IFRAME,DIM>(&frame, j, p));
      int lastquarkletnumber(last_qua.number());

      Index mu (first_generator<IFRAME,DIM>(&frame, j, p));
//      cout << "First generator: " << first_generator<IFRAME>(&frame, j, p) << endl;
//      cout << "Last quarklet: " << last_qua << endl;
      for (int it = mu.number();; ++it) {
          frame.get_support(it);
        if (intersect_supports(frame, lambda, *(frame.get_quarklet(it)))){

          intersecting.push_back(it);
        }
	if (it == lastquarkletnumber) break;
      }
    
    } 

    else {
      Index last_qua(last_quarklet<IFRAME,DIM >(&frame, j, p));
      int lastquarkletnumber(last_qua.number());
      Index mu = first_quarklet<IFRAME,DIM>(&frame, j, p);
//      cout << "First quarklet: " << first_quarklet<IFRAME>(&frame, j, p) << endl;
//      cout << "Last quarklet: " << last_qua << endl;
      for (int it = mu.number() /*frame.get_first_wavelet_numbers()[0]*/;; ++it) {
//          cout << "Mu: " << mu << ", " << mu.number() << endl;
          frame.get_support(it);
        if (intersect_supports(frame, lambda, *(frame.get_quarklet(it)))){
          intersecting.push_back(it);
        }
	if (it == lastquarkletnumber) break;
      }
    }

    }
    
    template <class IFRAME, unsigned int DIM>
  void intersecting_quarklets(const TensorFrame<IFRAME,DIM>& frame,
			     const int& lambda_num,
			     const typename TensorFrameIndex<IFRAME,DIM>::level_type& j,
			     std::list<int>& intersecting,
                             const typename TensorFrameIndex<IFRAME,DIM>::polynomial_type& p)
  {
//    typedef typename TensorFrame<IFRAME,DIM>::Index Index;
    typedef typename TensorFrame<IFRAME,DIM>::Support Support;

    intersecting.clear();

    Support supp, supp_lambda;
    frame.support(lambda_num, supp_lambda);
    int lastquarkletnumber(frame.get_last_wavelet_numbers()[0]);
//    intersect_supports(frame, 0, supp_lambda, supp);
    
    // a brute force solution
    if (j==frame.j0()) {
        
//      Index last_qua(last_quarklet<IFRAME>(&frame, j, p));
//      Index mu (first_generator<IFRAME>(&frame, j, p));
      for (int it = /*mu.number()*/ 0;; ++it) {
	if (intersect_supports(frame, it, supp_lambda, supp))
            intersecting.push_back(it);
	if (/* *frame.get_quarklet(it)*/ it == /*last_qua*/ lastquarkletnumber) break;
      }
    } 
    else {
//      Index last_qua(last_quarklet<IFRAME>(&frame, j, p));
//      Index mu = first_quarklet<IFRAME>(&frame, j, p);
//      cout << "First quarklet: " << first_quarklet<IFRAME>(&frame, j, p) << endl;
//      cout << "Last quarklet: " << last_qua << endl;
      for (int it = /*mu.number()*/ frame.get_first_wavelet_numbers()[0];; ++it) {
//          cout << "Mu: " << mu << ", " << mu.number() << endl;
	if (intersect_supports(frame, it, supp_lambda, supp))
	  intersecting.push_back(it);
	if (/* *frame.get_quarklet(it)*/ it == /*last_qua*/lastquarkletnumber) break;
      }
    }
  }
    


    template <class IFRAME, unsigned int DIM>
    bool intersect_singular_support(const TensorFrame<IFRAME,DIM>& frame,
                                    const typename TensorFrame<IFRAME,DIM>::Index& lambda,
                                    const typename TensorFrame<IFRAME,DIM>::Index& mu)
    {
//        cout<<"bin hier"<<endl;
        int j, k1, k2;
        for (unsigned int i = 0; i < DIM; i++) {
//            cout<<i<<endl;
            if (!(intersect_singular_support(*frame.frames()[i],                  //doch nicht logisch falsch
                                             lambda.j()[i], lambda.e()[i], lambda.k()[i],
                                             mu.j()[i], mu.e()[i], mu.k()[i],
                                             j, k1, k2) ))
                return false;
        }
        return true;
//            if ((intersect_singular_support(*frame.frames()[i],        
//                                             lambda.j()[i], lambda.e()[i], lambda.k()[i],
//                                             mu.j()[i], mu.e()[i], mu.k()[i],
//                                             j, k1, k2) ))
//                return true;
//        }
//        return false;
    }
}

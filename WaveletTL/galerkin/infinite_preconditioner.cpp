// implementation for infinite_preconditioner.h

#include <cmath>

namespace WaveletTL
{
  template <class INDEX>
  InfinitePreconditioner<INDEX>::~InfinitePreconditioner() {}

  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::apply_left_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    apply_preconditioner(y, x);
  }
    
  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::reverse_left_preconditioner
  (const InfiniteVector<double,INDEX>& x,
   InfiniteVector<double,INDEX>& y) const
  {
    reverse_preconditioner(x, y);
  }
    
  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::apply_right_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    apply_preconditioner(y, x);
  }
  
  template <class INDEX>
  inline
  void
  InfiniteSymmetricPreconditioner<INDEX>::reverse_right_preconditioner
  (const InfiniteVector<double,INDEX>& x,
   InfiniteVector<double,INDEX>& y) const
  {
    reverse_preconditioner(x, y);
  }
  
  template <class INDEX>
  inline
  void
  FullyDiagonalPreconditioner<INDEX>::apply_preconditioner
  (const InfiniteVector<double,INDEX>& y,
   InfiniteVector<double,INDEX>& x) const
  {
    x = y;
    x.scale(this, -1);
  };
  
  template <class INDEX>
  inline
  void
  FullyDiagonalPreconditioner<INDEX>::reverse_preconditioner
  (const InfiniteVector<double,INDEX>& x,
   InfiniteVector<double,INDEX>& y) const
  {
    y = x;
    y.scale(this, 1);
  };
  
  template <class INDEX>
  inline
  double
  FullyDiagonalDyadicPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return pow(ldexp(1.0, lambda.j()), operator_order());
  }

 
  template <class INDEX>
  inline
  double
  FullyDiagonalQuarkletPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
#ifdef DYADIC      
#if _WAVELETTL_USE_TFRAME==1
    double hspreconditioner(0), l2preconditioner(1);
//    H^s weights, cf. Diss Keding Formula (6.1.20):
//    (\sum_{i=1}^d (p_i+1)^{4s+\delta_2}4^{s j_i})^{1/2}
//    \prod_{i=1}^d \left(p_i+1\right)^{\delta_1/2}, \delta_1>0,\,\delta_2>1
    
    for (int i=0; i<_DIM; i++){
        hspreconditioner+=pow(1+lambda.p()[i],DELTA2+4*operator_order())*(1<<(2*lambda.j()[i])*(int) operator_order());
        l2preconditioner*=pow(1+lambda.p()[i],DELTA1*0.5);
    }
    double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
    return preconditioner;
    
#else
//    H^s weights, cf. Diss Keding Formula (5.3.9):
//    2^{js}*(p+1)^(2s+\delta_1/2+\delta_2/2), falls operator_order()>0 
    
    if (operator_order()==0) 
        return pow(1+lambda.p(),DELTA1*0.5);
    else
        return (1<<lambda.j()* (int) operator_order())*pow(1+lambda.p(),DELTA1*0.5+DELTA2*0.5+2*operator_order()); 
#endif
#else
    return 1;
#endif    
  }

  
  template <class INDEX>
  inline
  double
  FullyDiagonalQuarkletEnergyNormPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
//      double polynomialpreconditioner(1);
//      for (int i=0; i<_DIM; i++){
//        polynomialpreconditioner*=pow(1+lambda.p()[i],2);        
//      }
//      return sqrt(a(lambda, lambda)) *polynomialpreconditioner;
//      
      return sqrt(a(lambda, lambda));
      
  };

  template <class INDEX>
  inline
  double
  FullyDiagonalEnergyNormPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return sqrt(a(lambda, lambda));    
  };
  
  
  template <class INDEX>
  inline
  double
  FullyDiagonalDyPlusEnNormPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
   double hspreconditioner(0), l2preconditioner(1);

   //     weights used in the experiments for Diss Philipp Keding:
//    (\sum_{i=1}^d (p_i+1)^{4s+\delta_2})^{1/2}
//    \prod_{i=1}^d \left(p_i+1\right)^{\delta_1/2}*sqrt(a(\psi_\lambda,\psi_\lambda))/sqrt(dimension), \delta_1>0,\,\delta_2>1
    for (int i=0; i<_DIM; i++){
        hspreconditioner+=pow(1+lambda.p()[i],DELTA2+4*operator_order());
        l2preconditioner*=pow(1+lambda.p()[i],DELTA1*0.5);
    }
    double preconditioner = sqrt(hspreconditioner)*l2preconditioner*sqrt(a(lambda, lambda))/sqrt(_DIM);
        
    return preconditioner;      
   
    
  };

  template <class INDEX>
  inline
  double
  TrivialPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return 1;
  };
}

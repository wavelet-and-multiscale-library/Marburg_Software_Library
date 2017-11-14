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
      
#if _WAVELETTL_USE_TFRAME==1
    double hspreconditioner(0), l2preconditioner(1);
//    int space_dimension = (*(lambda.frame())).space_dimension;
//    int space_dimension = 2;
    for (int i=0; i<_DIM; i++){
        hspreconditioner+=pow(1+lambda.p()[i],8)*ldexp(1.0,2*lambda.j()[i]);
        l2preconditioner*=pow(1+lambda.p()[i],2);
    }
    double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
    return preconditioner;
    
#else
    return pow((1<<lambda.j())*pow(1+lambda.p(),4),operator_order())*pow(1+lambda.p(),2); //2^j*(p+1)^(2+\delta), falls operator_order()=1 (\delta=4)
#endif
    
  }


  template <class INDEX>
  inline
  double
  FullyDiagonalEnergyNormPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return sqrt(a(lambda, lambda));
    //return ldexp(1.0, lambda.j()); //ATTENTION!!! HAS TO BE CHANGED BACK; ONLY FOR EXPERIMENTING
  };

  template <class INDEX>
  inline
  double
  TrivialPreconditioner<INDEX>::diag(const INDEX& lambda) const
  {
    return 1;
  };
}

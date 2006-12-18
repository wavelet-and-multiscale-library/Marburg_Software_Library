// implementation for ldomain_jl_support.h

namespace WaveletTL
{
  inline
  void support(const LDomainJLBasis& basis,
 	       const Index& lambda,
 	       LDomainJLBasis::Support& supp)
  {
    if (lambda.e()[0]+lambda.e()[1] == 0) {
      // generator
      supp.j = lambda.j();
      supp.xmin = lambda.k()[0]-1;
      supp.xmax = lambda.k()[0]+1;
      supp.ymin = lambda.k()[1]-1;
      supp.ymax = lambda.k()[1]+1;
    } else {
      // wavelet
      supp.j = lambda.j()+1;
      supp.xmin = 2*(lambda.k()[0]-1);
      supp.xmax = 2*(lambda.k()[0]+1);
      supp.ymin = 2*(lambda.k()[1]-1);
      supp.ymax = 2*(lambda.k()[1]+1);
    }
  }
  
}

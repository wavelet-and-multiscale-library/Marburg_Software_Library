// implementation for ldomain_support.h

namespace WaveletTL
{
  template <class IBASIS>
  void support(const LDomainBasis<IBASIS>& basis,
	       const typename LDomainBasis<IBASIS>::Index& lambda,
	       typename LDomainBasis<IBASIS>::Support& supp)
  {
    typedef typename LDomainBasis<IBASIS>::Index Index;

    // by convention, we initialize the support with an "empty" set
    for (int p = 0; p <= 2; p++) {
      supp.xmin[p] = supp.xmax[p] = supp.ymin[p] = supp.ymax[p] = -1;
    }
    
    const int ecode = lambda.e()[0]+2*lambda.e()[1];

    if (ecode == 0) {
      // psi_lambda is a generator. Here we know by construction of the
      // composite basis that per patch, psi_lambda looks like a single
      // tensor product of 1D generators (possibly weighted by a factor).

      supp.j = lambda.j();

      switch (lambda.p()) {
      case 0:
 	// psi_lambda completely lives on patch 0

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[0],
		supp.xmax[0]);
	
	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[0],
		supp.ymax[0]);

 	break;
      case 1:
 	// psi_lambda completely lives on patch 1

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	break;
      case 2:
 	// psi_lambda completely lives on patch 2

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[2],
		supp.xmax[2]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[2],
		supp.ymax[2]);

 	break;
      case 3:
  	// psi_lambda lives on patches 0 and 1
	
	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[0],
		supp.xmax[0]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       basis.basis1d().DeltaLmin(),
				       &basis.basis1d()),
		supp.ymin[0],
		supp.ymax[0]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       basis.basis1d().DeltaRmax(lambda.j()),
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	break;
      case 4:
 	// psi_lambda lives on patches 1 and 2

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       basis.basis1d().DeltaRmax(lambda.j()),
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       basis.basis1d().DeltaLmin(),
				       &basis.basis1d()),
		supp.xmin[2],
		supp.xmax[2]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambda.j(),
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[2],
		supp.ymax[2]);

	break;
      }

    } else {
      // wavelet
    }
  }
}

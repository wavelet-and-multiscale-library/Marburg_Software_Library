// implementation for sturm_equation.h

#include <cmath>
#include <algorithm>
#include <list>
#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/gauss_data.h>



namespace WaveletTL
{
  template <class WBASIS>
  SturmEquation<WBASIS>::SturmEquation(const SimpleSturmBVP& bvp,
				       const bool precompute_f)
    : bvp_(bvp), basis_(bvp.bc_left(), bvp.bc_right()), normA(0.0), normAinv(0.0)
  {
#ifdef ENERGY
//      compute_diagonal();
#endif
    if (precompute_f) precompute_rhs();
    
    
    //const int jmax = 12;
    //basis_.set_jmax(jmax);
  }

  template <class WBASIS>
  SturmEquation<WBASIS>::SturmEquation(const SimpleSturmBVP& bvp,
				       const WBASIS& basis,
				       const bool precompute_f)
    : bvp_(bvp), basis_(basis), normA(0.0), normAinv(0.0)
  {
#ifdef ENERGY
      compute_diagonal();
#endif
    if (precompute_f) precompute_rhs();
    
    
    //const int jmax = 12;
    //basis_.set_jmax(jmax);
  }

  template <class WBASIS>
  void
  SturmEquation<WBASIS>::precompute_rhs() const
  {
    typedef typename WaveletBasis::Index Index;
    cout << "precompute rhs.." << endl;
    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    InfiniteVector<double,int> fhelp_int;
    
#ifdef FRAME
//    cout << basis_.degrees_of_freedom() << endl;

    for (int i=0; i<basis_.degrees_of_freedom();i++) {
//        cout << "hallo" << endl;
//        cout << *(basis_.get_quarklet(i)) << endl;
        const double coeff = f(*(basis_.get_quarklet(i)))/D(*(basis_.get_quarklet(i)));
        fhelp.set_coefficient(*(basis_.get_quarklet(i)), coeff);
        fhelp_int.set_coefficient(i, coeff);
//        cout << *(basis_.get_quarklet(i)) << endl;
    }
//    cout << "bin hier1" << endl;
#else
    for (int i=0; i<basis_.degrees_of_freedom();i++) {
//        cout << "bin hier: " << i << endl;
//        cout << D(*(basis_.get_wavelet(i))) << endl;
//        cout << *(basis_.get_wavelet(i)) << endl;
        const double coeff = f(*(basis_.get_wavelet(i)))/D(*(basis_.get_wavelet(i)));
//        cout << f(*(basis_.get_wavelet(i))) << endl;
//        cout << coeff << endl;
        fhelp.set_coefficient(*(basis_.get_wavelet(i)), coeff);
        fhelp_int.set_coefficient(i, coeff);
//        cout << *(basis_.get_wavelet(i)) << endl;
    }
//    const int j0   = basis().j0();
//    for (Index lambda(basis_.first_generator(j0));;++lambda)
//      {
//	const double coeff = f(lambda)/D(lambda);
//	if (fabs(coeff)>1e-15)
//	  fhelp.set_coefficient(lambda, coeff);
//          fhelp_int.set_coefficient(i, coeff);
//	if (lambda == basis_.last_wavelet(jmax))
//	  break;
//        
//            
//      }
#endif
    fnorm_sqr = l2_norm_sqr(fhelp);
    
    // sort the coefficients into fcoeffs
    fcoeffs.resize(fhelp.size());
    fcoeffs_int.resize(fhelp_int.size());
    unsigned int id(0), id2(0);
    for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
	 it != itend; ++it, ++id)
      fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
    sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
    
    for (typename InfiniteVector<double,int>::const_iterator it(fhelp_int.begin()), itend(fhelp_int.end());
	 it != itend; ++it, ++id2)
      fcoeffs_int[id2] = std::pair<int,double>(it.index(), *it);
    sort(fcoeffs_int.begin(), fcoeffs_int.end(), typename InfiniteVector<double,int>::decreasing_order());
    
    rhs_precomputed = true;
    cout << "end precompute rhs.." << endl;
//    cout << fhelp << endl;
//    cout << fcoeffs << endl;
  }

  template <class WBASIS>
  inline
  double
  SturmEquation<WBASIS>::D(const typename WBASIS::Index& lambda) const
  {
#ifdef FRAME
  #ifdef DYADIC
      return mypow((1<<lambda.j())*mypow(1+lambda.p(),6),operator_order())*mypow(1+lambda.p(),2); //2^j*(p+1)^6, falls operator_order()=1 (\delta=4)
//      return 1<<(lambda.j()*(int) operator_order());
#endif
    #ifdef TRIVIAL
      return 1;
#endif

#ifdef ENERGY 
      return stiff_diagonal[lambda.number()]*(lambda.p()+1);
#endif
#endif


#ifdef BASIS
  #ifdef DYADIC
    return 1<<(lambda.j()*(int) operator_order());
//    return pow(ldexp(1.0, lambda.j()),operator_order());
  #else
    #ifdef TRIVIAL
      return 1;
    #else
      #ifdef ENERGY
//        return sqrt(a(lambda, lambda));
        return stiff_diagonal[lambda.number()];
      #else
        return sqrt(a(lambda, lambda));
      #endif  
    #endif
  #endif
#else
  return 1;
#endif

//    return 1;
//     return lambda.e() == 0 ? 1.0 : ldexp(1.0, lambda.j()); // do not scale the generators
//     return lambda.e() == 0 ? 1.0 : sqrt(a(lambda, lambda)); // do not scale the generators

  }

  template <class WBASIS>
  inline
  double
  SturmEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
			   const typename WBASIS::Index& nu) const
  {
    return a(lambda, nu, 2*WBASIS::primal_polynomial_degree());
  }
  
  template <class WBASIS>
  double
  SturmEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
			   const typename WBASIS::Index& nu,
			   const unsigned int p) const
  {
    // a(u,v) = \int_0^1 [p(t)u'(t)v'(t)+q(t)u(t)v(t)] dt

    double r = 0;

    // Remark: There are of course many possibilities to evaluate
    // a(u,v) numerically.
    // In this implementation, we rely on the fact that the primal functions in
    // WBASIS are splines with respect to a dyadic subgrid.
    // We can then apply an appropriate composite quadrature rule.
    // In the scope of WBASIS, the routines intersect_supports() and evaluate()
    // must exist, which is the case for DSBasis<d,dT>.

    // First we compute the support intersection of \psi_\lambda and \psi_\nu:
    typedef typename WBASIS::Support Support;

    Support supp;

    if (intersect_supports(basis_, lambda, nu, supp))
      {
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
        

        #ifdef FRAME
	const unsigned int N_Gauss = std::min((unsigned int)10,(p+1)/2+ (lambda.p()+nu.p()+1)/2);
//        const unsigned int N_Gauss = 10;
//        const unsigned int N_Gauss = (p+1)/2;
         #else 
        const unsigned int N_Gauss = (p+1)/2;
        #endif
        const double h = ldexp(1.0, -supp.j);
	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values, der1values, der2values;
	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

	// - compute point values of the integrands
	evaluate(basis_, lambda, gauss_points, func1values, der1values);        
	evaluate(basis_, nu, gauss_points, func2values, der2values);
//        if((lambda.number()==19 && nu.number()==19) || (lambda.number()==26 && nu.number()==26)){
//            cout << lambda << endl;
//            cout << gauss_points << endl;
//            cout << func1values << endl;
//            cout << func2values << endl;
//        }

	// - add all integral shares
	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    const double t = gauss_points[id];
	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	    const double pt = bvp_.p(t);
  	    if (pt != 0)
	      r += pt * der1values[id] * der2values[id] * gauss_weight;
	    
	    const double qt = bvp_.q(t);
  	    if (qt != 0)
	      r += qt * func1values[id] * func2values[id] * gauss_weight;
	  }
      }

    return r;
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;
#ifdef FRAME
     
      const int pmax = std::min(basis().get_pmax_(),2);
      //const int pmax = 0;
      int p = 0;
      
      for (Index lambda = basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == basis().last_wavelet(jmax,p)){
            ++p;
            lambda = basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
#else
      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == last_wavelet(&basis(), jmax)) break;
      }
#endif      
      
      
      
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
    }

    return normA;
  }
   
  template <class WBASIS>
  double
  SturmEquation<WBASIS>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;
      
#ifdef FRAME
      
      const int pmax = std::min(basis().get_pmax_(),2);
      //const int pmax = 0;
      int p = 0;
      
      for (Index lambda = basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == basis().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == basis().last_wavelet(jmax,p)){
            ++p;
            lambda = basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
#else
      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == last_wavelet(&basis(), jmax)) break;
      }
#endif
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);

#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
    }

    return normAinv;
  }

  template <class WBASIS>
  double
  SturmEquation<WBASIS>::f(const typename WBASIS::Index& lambda) const
  {
    // f(v) = \int_0^1 g(t)v(t) dt
//      cout << "bin in f" << endl;
    double r = 0;

    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis_, lambda, k1, k2);

    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 7; //perhaps we need +lambda.p()/2 @PHK
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1)), vvalues;
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;

    // - compute point values of the integrand
    evaluate(basis_, 0, lambda, gauss_points, vvalues);
//    cout << "bin immer noch in f" << endl;
    // - add all integral shares
    for (int patch = k1, id = 0; patch < k2; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	const double t = gauss_points[id];
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    
	const double gt = bvp_.g(t);
	if (gt != 0)
	  r += gt
	    * vvalues[id]
	    * gauss_weight;
      }
    
#ifdef DELTADIS
//    double tmp = 1;
//    Point<1> p1;
//    p1[0] = 0.5;
//    Point<1> p2;
//    chart->map_point_inv(p1,p2);
//    tmp =  evaluate(basis_, 0,
//			       typename WBASIS::Index(lambda.j(),
//						      lambda.e()[0],
//						      lambda.k()[0],
//						      basis_),
//			       p2[0]);
//    tmp /= chart->Gram_factor(p2);
//  
//  
//    return 4.0*tmp + r;    
#ifdef NONZERONEUMANN
    return r + 4*basis_.evaluate(0, lambda, 0.5)+3*M_PI*(basis_.evaluate(0, lambda, 1)+basis_.evaluate(0, lambda, 0));
#else
    return r+ 4*basis_.evaluate(0, lambda, 0.5);
#endif    
#else
    return r;
#endif
  }
  
  template <class WBASIS>
  inline
  void
  SturmEquation<WBASIS>::RHS(const double eta,
			     InfiniteVector<double, typename WBASIS::Index>& coeffs) const
  {
    if (!rhs_precomputed) precompute_rhs();

    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    typedef typename WBASIS::Index Index;
    typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs.end() && coarsenorm < bound);
  }
  
  template <class WBASIS>
  inline
  void
  SturmEquation<WBASIS>::RHS(const double eta,
			     InfiniteVector<double,int>& coeffs) const
  {
    if (!rhs_precomputed) precompute_rhs();

    coeffs.clear();
    double coarsenorm(0);
    double bound(fnorm_sqr - eta*eta);
    typename Array1D<std::pair<int, double> >::const_iterator it(fcoeffs_int.begin());
    do {
      coarsenorm += it->second * it->second;
      coeffs.set_coefficient(it->first, it->second);
      ++it;
    } while (it != fcoeffs_int.end() && coarsenorm < bound);
  }



  template <class WBASIS>
  void
  SturmEquation<WBASIS>::compute_diagonal()
  {
    cout << "SturmEquation(): precompute diagonal of stiffness matrix..." << endl;

    SparseMatrix<double> diag(1,basis_.degrees_of_freedom());
    char filename[50];
    char matrixname[50];
#ifdef ONE_D
    int d = WBASIS::primal_polynomial_degree();
    int dT = WBASIS::primal_vanishing_moments();
#else
#ifdef TWO_D
    int d = WBASIS::primal_polynomial_degree();
    int dT = WBASIS::primal_vanishing_moments();
#endif
#endif
    
    // prepare filenames for 1D and 2D case
#ifdef ONE_D
    sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_interval_lap07_d", d, "_dT", dT);
    sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_1D_lap07_d", d, "_dT", dT);
#endif
#ifdef TWO_D
    sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_lshaped_lap1_d", d, "_dT", dT);
    sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_2D_lap1_d", d, "_dT", dT);
#endif
#ifndef PRECOMP_DIAG
    std::list<Vector<double>::size_type> indices;
    std::list<double> entries;
#endif
#ifdef PRECOMP_DIAG
    cout << "reading in diagonal of unpreconditioned stiffness matrix from file "
	 << filename << "..." << endl;
    diag.matlab_input(filename);
    cout << "...ready" << endl;
#endif
    stiff_diagonal.resize(basis_.degrees_of_freedom());
    for (int i = 0; i < basis_.degrees_of_freedom(); i++) {
#ifdef PRECOMP_DIAG
      stiff_diagonal[i] = diag.get_entry(0,i);
#endif 
#ifndef PRECOMP_DIAG
#ifdef FRAME
      stiff_diagonal[i] = sqrt(a(*(basis_.get_quarklet(i)),*(basis_.get_quarklet(i))));
#endif
#ifdef BASIS
      stiff_diagonal[i] = sqrt(a(*(basis_.get_wavelet(i)),*(basis_.get_wavelet(i))));
#endif
      
      indices.push_back(i);
      entries.push_back(stiff_diagonal[i]);
#endif
      //cout << stiff_diagonal[i] << " " << *(basis_->get_wavelet(i)) << endl;
    }
#ifndef PRECOMP_DIAG
    diag.set_row(0,indices, entries);
    diag.matlab_output(filename, matrixname, 1);
#endif

    cout << "... done, diagonal of stiffness matrix computed" << endl;
  }
}

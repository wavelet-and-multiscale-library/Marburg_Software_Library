// implementation for pq_expansion.h

#include <set>
#include <list>

#include <utils/array1d.h>
#include <utils/function.h>
#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <galerkin/gramian.h>
#include <galerkin/cached_problem.h>
#include <interval/interval_bspline.h>
#include <adaptive/cdd1.h>

#include "i_q_index.h"
#include "pq_frame.h"

namespace WaveletTL
{
  template <int d, int dT>
  double integrate(const Function<1>* f,
		   const PQFrame<d,dT>& basis,
		   const typename PQFrame<d,dT>::Index& lambda)
  {
    double r = 0;
    
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis, lambda, k1, k2);
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5/*+(lambda.p())*/;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1));
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
    
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	for (int patch = k1; patch < k2; patch++)
	  {
	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
	    const double ft = f->value(Point<1>(t));
	    if (ft != 0)
	      r += ft
		* evaluate(basis, 0, lambda, t)
		* gauss_weight;
            //cout << "QStelle: " << t << ", Wert: " << evaluate(basis, 0, lambda, t) << endl;
	  }
      }
    
    return r;
  }

  template <int d, int dT>
  double integrate(const PQFrame<d,dT>& basis,
		   const typename PQFrame<d,dT>::Index& lambda,
		   const typename PQFrame<d,dT>::Index& mu, const int derivative)
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef typename PQFrame<d,dT>::Support Support;
    Support supp;

    if (intersect_supports(basis, lambda, mu, supp))
      {
 	// Set up Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = d+(lambda.p()+mu.p())/2;
 	const double h = ldexp(1.0, -supp.j);
 	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
 	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
 	// - compute point values of the integrands
  	evaluate(basis, derivative, lambda.p(), lambda.j(), lambda.e(), lambda.k(), gauss_points, func1values);
 	evaluate(basis, derivative, mu.p(), mu.j(), mu.e(), mu.k(), gauss_points, func2values);
	
 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += func1values[id] * func2values[id] * gauss_weight;
 	  }
      }
    
    return r;
  }
  
  
  
  
  
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PQFrame<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename PQFrame<d,dT>::Index>& coeffs,
         const int pmax)
  {
    typedef typename PQFrame<d,dT>::Index Index;
    const int j0 = basis.j0();
    assert(jmax >= j0);
    
    coeffs.clear();
    int p = 0;
    for (Index lambda = basis.first_generator(j0);;)
      {
 	coeffs.set_coefficient(lambda, integrate(f, basis, lambda));
 	if (lambda == basis.last_wavelet(jmax,pmax))
	  break;
        if (lambda == basis.last_wavelet(jmax,p)){
            ++p;
            lambda = basis.first_generator(j0, p);
        }
        else
            ++lambda;
        
      }
    if (!primal) {
#if 0
      IntervalGramian<PQFrame<d,dT> > G(basis, coeffs);
      CachedProblem<IntervalGramian<PQFrame<d,dT> > > GC(&G);
      InfiniteVector<double, typename PQFrame<d,dT>::Index> x;
      CDD1_SOLVE(GC, 1e-6, x, jmax);
      coeffs.swap(x);
#else
      
      //cout << "Beginn Indexset" << endl;
      // setup active index set
      p = 0;
      std::set<Index> Lambda;
      for (Index lambda = basis.first_generator(j0);;) {
 	Lambda.insert(lambda);
 	if (lambda == basis.last_wavelet(jmax,pmax))
	  break;
        if (lambda == basis.last_wavelet(jmax,p)){
            ++p;
            lambda = basis.first_generator(j0, p);
        }
        else
            ++lambda;
        
      }
      //cout << "Indexset fertig" << endl;
      // setup Gramian A_Lambda
      //cout << "Größe Lambda: " << Lambda.size() << endl;
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	     it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(basis, *it2, *it1);

	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 
      
      //cout << "Gramsche fertig" << endl;
      // solve A_Lambda*x = b
      Vector<double> b(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	b[row] = coeffs.get_coefficient(*it);

      Vector<double> x(b);
      unsigned int iterations;
      
      //cout << "Beginn CG" << endl;
      
      CG(A_Lambda, b, x, 1e-12, 200, iterations);
      //Richardson(A_Lambda, b, x, 0.8, 1e-12, 200, iterations);
      //SSOR(A_Lambda, b, x, 1, 1e-12, 200, iterations);
      
      //cout << "Ende CG" << endl;
     // cout << "Iterationen: " << iterations << endl;
      
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	coeffs.set_coefficient(*it, x[row]);
#endif
    }
  }

  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PQFrame<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs)
  {
    typedef typename PQFrame<d,dT>::Index Index;
    const int j0 = basis.j0();
    assert(jmax >= j0);
    
    coeffs.resize(basis.Deltasize(jmax+1));
    int id = 0;
    for (Index lambda = basis.first_generator(j0);;++lambda, ++id)
      {
 	coeffs[id] = integrate(f, basis, lambda);
	if (lambda == basis.last_wavelet(jmax))
	  break;
      }
    
    if (!primal) {
      // setup active index set
      std::set<Index> Lambda;
      for (Index lambda = basis.first_generator(j0);; ++lambda) {
 	Lambda.insert(lambda);
	if (lambda == basis.last_wavelet(jmax)) break;
      }
      
      // setup Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	       it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(basis, *it2, *it1);
	      
	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 
      
      // solve A_Lambda*x = b
      Vector<double> b(coeffs);
      
      unsigned int iterations;
      CG(A_Lambda, b, coeffs, 1e-15, 500, iterations);
    }
  }
  
  template <int d, int dT>
  double factor(const PQFrame<d,dT>& basis, const int l, const int k, const int r, const int q, const int j){
      double v=0;
      typedef typename PQFrame<d,dT>::Index Index;
      Index lambda(0,j+1,0,l,&basis);
      Index mu(0,j,1,k,&basis);
      
      unsigned int degree=q+r;
      InfiniteVector<double, Index> c;
      basis.reconstruct_1(mu, j+1, c);
//      cout << "mu: " << mu << endl;
//      cout << c << endl;
      
      MonomeFunction f(degree);
      v=c[lambda]*integrate(&f, basis, lambda);
      return v;
  }
  
  template <int d, int dT>
  double factor2(const PQFrame<d,dT>& basis, const int r, const int p, const int j, const int l){
      double v=0;
      typedef typename PQFrame<d,dT>::Index Index;
      Index lambda(p,j+1,0,l,&basis);
      MonomeFunction f(r);
      v=integrate(&f, basis, lambda);
      return v;
  }
  
  
  template <int d, int dT>
  void setup_factor_matrix(const PQFrame<d,dT>& basis, SparseMatrix<double>& A_Lambda, const int j, const int p)
  {
      //cout << "Hallo" << endl;
      A_Lambda.resize(dT*basis.Nablasize(j),(p+1)*basis.Deltasize(j+1));
    //cout << A_Lambda.row_dimension() << endl;
    
    typedef typename SparseMatrix<double>::size_type size_type;
    typedef typename PQFrame<d,dT>::Index Index;
    
    std::set<Index> Lambda;
      for (Index lambda = basis.first_wavelet(j);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == basis.last_wavelet(j)) break;
        //if (i==7) break;
      }
    
    std::set<Index> Mu;
      for (Index mu = basis.first_generator(j+1);; ++mu) {
	Mu.insert(mu);
	if (mu == basis.last_generator(j+1)) break;
        //if (i==7) break;
      }
    
    size_type row = 0;
    for(int r=0; r<dT; ++r){
        for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
            it1 != itend; ++it1, ++row)
        {
            //cout << (*it1).j() << ", " << row << endl;
            Index nu(*it1);
            //cout << nu.k()<< endl;
            std::list<size_type> indices;
            std::list<double> entries;

            size_type column = 0;
            for (typename std::set<Index>::const_iterator it2(Mu.begin()), itend(Mu.end());
                it2 != itend; ++it2){
                    for(int q=0; q<=p; ++q, ++column){
                        double entry = factor(basis, (*it2).k(), (*it1).k(), r, q, j);
                            if (fabs(entry) > 1e-15) {
                                indices.push_back(column);
                                entries.push_back(entry);
                            }
                    }
            }
	A_Lambda.set_row(row, indices, entries);
      }
    }
  }
  
  template <int d, int dT>
  double rightside(const PQFrame<d,dT>& basis, typename PQFrame<d,dT>::Index lambda, const int r, const bool leftborder){
      
      double wert = 0;
      typedef typename PQFrame<d,dT>::Index Index;
      InfiniteVector<double, Index> c;
      basis.reconstruct_1(lambda, lambda.j()+1, c);
      MonomeFunction f(r);
      
      
      if(leftborder==1){
        typename InfiniteVector<double, Index>::const_iterator it(c.begin()), it2(c.end());
        for(int k=0; k<dT; k++, it++){
            
        }
        
      
      
        for(it; it<it2; it++){
          wert-=*it*integrate(&f, basis, it.index());           
        }
      }
      else{
        typename InfiniteVector<double, Index>::const_reverse_iterator it3(c.rbegin()), it4(c.rend());
        for(int k=0; k<dT; k++, it3++){
            
        }            
        
      
      
        for(it3; it3<it4; it3++){
          wert-=*it3*integrate(&f, basis, it3.index());           
        }
      }
      
      
      
      return wert;
      
      
  }
  
  template <int d, int dT>
  void rightsidevector(const PQFrame<d,dT>& basis, typename PQFrame<d,dT>::Index lambda, Vector<double>& coeffs, const bool leftborder){
      coeffs.resize(dT);
      for(int i=0; i<dT; i++){
          coeffs[i]=rightside(basis, lambda, i, leftborder);
      }
      
  }
  
  template <int d, int dT>
  void system_matrix(const PQFrame<d,dT>& basis, Matrix<double>& A, typename PQFrame<d,dT>::Index lambda, const bool leftborder){
      typedef typename PQFrame<d,dT>::Index Index;
      
      A.resize(dT,dT);
      InfiniteVector<double, Index> c;
      std::set<Index> supp;
      basis.reconstruct_1(lambda, lambda.j()+1, c);
      
      
      if(leftborder == 1) {
        typename InfiniteVector<double, Index>::const_iterator it(c.begin());
        for(int l=0; l<dT; l++, it++){
            for(int k=0; k<dT; k++){
                MonomeFunction f(k);
                A(k,l)=integrate(&f, basis, it.index());
            }
          
        }
      }
      else{
          
        typename InfiniteVector<double, Index>::const_reverse_iterator it2(c.rbegin());
        for(int l=dT-1; l>=0; l--, it2++){
            for(int k=0; k<dT; k++){
                MonomeFunction f(k);
                A(k,l)=integrate(&f, basis, it2.index());
            }
        }
      }
  }
}

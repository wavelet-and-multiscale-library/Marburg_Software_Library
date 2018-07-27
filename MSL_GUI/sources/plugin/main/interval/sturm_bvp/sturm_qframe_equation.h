/*  -*- c++ -*-

   +-----------------------------------------------------------------------+
   | MSL GUI - A Graphical User Interface for the Marburg Software Library |
   |                                                                       |
   | Copyright (C) 2018 Henning Zickermann                                 |
   | Contact: <zickermann@mathematik.uni-marburg.de>                       |
   +-----------------------------------------------------------------------+

     This file extensively uses code from the Marburg Software Library,
     WaveletTL, which is Copyright (C) 2002-2009 Thorsten Raasch, Manuel Werner.


     This file is part of MSL GUI.

     MSL GUI is free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     MSL GUI is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with MSL GUI.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef WAVELETTL_STURM_QFRAME_EQUATION_H
#define WAVELETTL_STURM_QFRAME_EQUATION_H


#include <set>
#include <MathTL/utils/array1d.h>
#include <MathTL/numerics/sturm_bvp.h>
#include <WaveletTL/galerkin/galerkin_utils.h>
#include <WaveletTL/galerkin/infinite_preconditioner.h>
#include <cmath>
#include <algorithm>
#include <list>
#include <MathTL/algebra/vector.h>
#include <MathTL/algebra/sparse_matrix.h>
#include <MathTL/numerics/eigenvalues.h>
#include <MathTL/numerics/gauss_data.h>

#include "WaveletTL/interval/i_q_index.h"
#include "WaveletTL/galerkin/cached_quarklet_problem.h"

using namespace MathTL;

namespace WaveletTL
{

/*!
 * Overload of add_compressed_column_quarklet(...) for CachedQuarkletProblem.
 * We overload this method for being able to use the non-TFRAME-branch of the
 * method independently of the _WAVELETTL_USE_TFRAME macro.
 */
template <class PROBLEM>
void
add_compressed_column_quarklet(const CachedQuarkletProblem<PROBLEM>& P,
                               const double factor,
                               const typename PROBLEM::Index& lambda,
                               const int J,
                               //InfiniteVector<double, typename PROBLEM::Index>& w,
                               Vector<double>& w,
                               const int jmax,
                               const CompressionStrategy strategy,
                               const bool preconditioning,
                               const int pmax,
                               const double a,
                               const double b) //a and b prefactors in strategy DKOR
{

  //typedef typename PROBLEM::QuarkletFrame QuarkletFrame;
  //    typedef typename PROBLEM::Index Index;
  //typedef typename WaveletBasis::Support Support;


  //     if (P.local_operator())

  // differential operators

  if (strategy == DKR) {
    ////
    // Quarklet strategy:
    // active row indices nu have to fulfill ||nu|-|lambda|| <= J/(d*b) and
    // the supports of psi_lambda and psi_nu have to intersect


    //            cout << "bin in DKOR drin" << endl;
    //ATTENTION: does not seem to work correctly with b and a
    const int maxlevel = std::min(lambda.j()+ (int) (J/(P.space_dimension)), jmax);
    //            cout << maxlevel << endl;
    //            cout << std::max(P.basis().j0()-1, lambda.j()- (int) (J/(P.space_dimension * 2))) <<endl;
    //              cout << lambda << endl;


    for (int level = std::max(P.frame().j0()-1, lambda.j()- (int) (J/(P.space_dimension * b)));
         level <= maxlevel; level++)
    {
      const int minplevel = std::max(0, lambda.p() + 1  - (int) pow(2,(J-b*abs(level-lambda.j())/a)));
      const int maxplevel = std::min(lambda.p() - 1  + (int) pow(2,(J-b*abs(level-lambda.j())/a)), pmax);

      for (int polynomial = minplevel;
           polynomial <= maxplevel; polynomial++)
      {
        //                        cout << "adding level: " << level << endl;
        P.add_level(lambda,w,polynomial,level,factor,J,strategy,jmax,pmax,a,b);
        //                        P.add_level(lambda,w,0,level,factor,J);
        //                cout << w << endl;
        //                cout << "Stop" << endl;
      }
    }

  }

  if (strategy == CDD1) {
    // [CDD1] strategy:
    // active row indices nu have to fulfill ||nu|-|lambda|| <= J/d and
    // the supports of psi_lambda and psi_nu have to intersect

    const int maxlevel = std::min(lambda.j()+(J/P.space_dimension), jmax);
    for (int level = std::max(P.frame().j0()-1, lambda.j()-(J/P.space_dimension));
         level <= maxlevel; level++)
    {
      //cout << "adding level: " << level << endl;
      P.add_level(lambda,w,0,level,factor,J);
    }
  }



  //     else
  //       {
  // 	// integral operators: branch is not implemented so far
  //       }
}



/*!
 * Partial specialization of FullyDiagonalQuarkletPreconditioner<INDEX, DIM>
 * for INDEX == IntervalQIndex<IFRAME>, DIM == 1 (same as original with _WAVELETTL_USE_TFRAME==0)
 */
template <class IFRAME>
class FullyDiagonalQuarkletPreconditioner<IntervalQIndex<IFRAME>, 1>
    : public FullyDiagonalPreconditioner< IntervalQIndex<IFRAME> >
{
public:
    FullyDiagonalQuarkletPreconditioner(double delta1 = 6, double delta2 = 2)
        : delta1_(delta1), delta2_(delta2)
    {

    }

    void set_delta1(double delta1) { delta1_ = delta1; }
    void set_delta2(double delta2) { delta2_ = delta2; }

  /*!
    (half) operator order t
   */
  virtual double operator_order() const = 0;

  virtual double a(const IntervalQIndex<IFRAME>& lambda,
                   const IntervalQIndex<IFRAME>& nu) const = 0;

  /*!
    evaluate the diagonal preconditioner D
  */
  double diag(const IntervalQIndex<IFRAME>& lambda) const
  {
    //return pow((1<<lambda.j())*pow(1+lambda.p(),4),operator_order())*pow(1+lambda.p(),2); //2^j*(p+1)^(2+\delta), falls operator_order()=1 (\delta=4)

    //    H^s weights, cf. Diss Keding Formula (5.3.9):
    //    2^{js}*(p+1)^(2s+\delta_1/2+\delta_2/2), falls operator_order()>0

        if (operator_order()==0)
            return pow(1+lambda.p(),delta1_*0.5);
        else
            return (1<<lambda.j()* (int) operator_order())*pow(1+lambda.p(),delta1_*0.5+delta2_*0.5+2*operator_order());
  }

private:
  double delta1_, delta2_;
};


/*!
 * The class SturmQFrameEquation<WBASIS> is identical to SturmEquation<WBASIS> with
 * macros FRAME and DYADIC defined, but is independent of these macros.
 */
template <class WBASIS>
class SturmQFrameEquation
    : public FullyDiagonalQuarkletPreconditioner<typename WBASIS::Index, 1>

{
public:
  /*!
      make template argument accessible
    */
  typedef WBASIS WaveletBasis;

  SturmQFrameEquation(const SimpleSturmBVP& bvp,
                      const bool precompute_rhs = true);

  SturmQFrameEquation(const SimpleSturmBVP& bvp,
                      const WaveletBasis& basis,
                      const bool precompute_rhs = true);

  /*!
      wavelet index class
    */
  typedef typename WaveletBasis::Index Index;

  /*!
      read access to the basis
    */
  const WBASIS& basis() const { return basis_; }

  /*!
      space dimension of the problem
    */
  static const int space_dimension = 1;

  /*!
      differential operators are local
    */
  static bool local_operator() { return true; }

  /*!
      (half) order t of the operator
      (inherited from FullyDiagonalDyadicPreconditioner)
    */
  double operator_order() const { return (bvp_.p(0.0)==1 ? 1. : 0.); }

  /*!
      evaluate the diagonal preconditioner D
    */
  double D(const Index& lambda) const;

  /*!
      evaluate the (unpreconditioned) bilinear form a
      (inherited from FullyDiagonalEnergyNormPreconditioner)
    */
  double a(const Index& lambda,
           const Index& nu) const;

  /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
  double a(const Index& lambda,
           const Index& nu,
           const unsigned int p) const;

  /*!
      estimate the spectral norm ||A||
    */
  double norm_A() const;

  /*!
      estimate the spectral norm ||A^{-1}||
    */
  double norm_Ainv() const;

  /*!
      estimate compressibility exponent s^*
    */
  double s_star() const {
    return 1.0 + WBASIS::primal_vanishing_moments(); // [St04a], Th. 2.3 for n=1
    //return 1.5;
  }

  /*!
      estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
  double alphak(const unsigned int k) const {
    return 2*norm_A(); // suboptimal
  }

  /*!
      evaluate the (unpreconditioned) right-hand side f
    */
  double f(const Index& lambda) const;

  /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
  void RHS(const double eta, InfiniteVector<double,Index>& coeffs) const;

  //int version
  void RHS(const double eta, InfiniteVector<double,int>& coeffs) const;

  /*!
      compute (or estimate) ||F||_2
    */
  double F_norm() const { return sqrt(fnorm_sqr); }

protected:
  const SimpleSturmBVP& bvp_;
  WBASIS basis_;

  // flag whether right-hand side has already been precomputed
  mutable bool rhs_precomputed;

  /*!
      precomputation of the right-hand side
      (constness is not nice but necessary to have RHS a const function)
    */
  void precompute_rhs() const;

  // right-hand side coefficients on a fine level, sorted by modulus
  mutable Array1D<std::pair<Index,double> > fcoeffs;
  mutable Array1D<std::pair<int,double> > fcoeffs_int;

  //! Square root of coefficients on diagonal of stiffness matrix.
  Array1D<double> stiff_diagonal;

  // (squared) \ell_2 norm of the precomputed right-hand side
  mutable double fnorm_sqr;

  // estimates for ||A|| and ||A^{-1}||
  mutable double normA, normAinv;
};
}


/*##################################################################################################
    Implementation
##################################################################################################*/

namespace WaveletTL
{
  template <class WBASIS>
  SturmQFrameEquation<WBASIS>::SturmQFrameEquation(const SimpleSturmBVP& bvp,
                                       const bool precompute_f)
    : bvp_(bvp), basis_(bvp.bc_left(), bvp.bc_right()), normA(0.0), normAinv(0.0)
  {
    if (precompute_f) precompute_rhs();


    //const int jmax = 12;
    //basis_.set_jmax(jmax);
  }

  template <class WBASIS>
  SturmQFrameEquation<WBASIS>::SturmQFrameEquation(const SimpleSturmBVP& bvp,
                                       const WBASIS& basis,
                                       const bool precompute_f)
    : bvp_(bvp), basis_(basis), normA(0.0), normAinv(0.0)
  {
    if (precompute_f) precompute_rhs();


    //const int jmax = 12;
    //basis_.set_jmax(jmax);
  }

  template <class WBASIS>
  void
  SturmQFrameEquation<WBASIS>::precompute_rhs() const
  {
    typedef typename WaveletBasis::Index Index;
    cout << "precompute rhs.." << endl;
    // precompute the right-hand side on a fine level
    InfiniteVector<double,Index> fhelp;
    InfiniteVector<double,int> fhelp_int;


    cout << basis_.degrees_of_freedom() << endl;
    for (int i=0; i<basis_.degrees_of_freedom();i++) {
//        cout << "hallo" << endl;
//        cout << *(basis_.get_quarklet(i)) << endl;
        const double coeff = f(*(basis_.get_quarklet(i)))/D(*(basis_.get_quarklet(i)));
        fhelp.set_coefficient(*(basis_.get_quarklet(i)), coeff);
        fhelp_int.set_coefficient(i, coeff);
        cout << *(basis_.get_quarklet(i)) << endl;
    }
//    cout << "bin hier1" << endl;

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
  SturmQFrameEquation<WBASIS>::D(const typename WBASIS::Index& lambda) const
  {
      return mypow((1<<lambda.j())*mypow(1+lambda.p(),4),operator_order())*mypow(1+lambda.p(),2); //2^j*(p+1)^6, falls operator_order()=1 (\delta=4)
//      return 1<<(lambda.j()*(int) operator_order());
  }

  template <class WBASIS>
  inline
  double
  SturmQFrameEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
                           const typename WBASIS::Index& nu) const
  {
    return a(lambda, nu, 2*WBASIS::primal_polynomial_degree());
  }

  template <class WBASIS>
  double
  SturmQFrameEquation<WBASIS>::a(const typename WBASIS::Index& lambda,
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


        const unsigned int N_Gauss = std::min((unsigned int)10,(p+1)/2+ (lambda.p()+nu.p()+1)/2);
//        const unsigned int N_Gauss = 10;
//        const unsigned int N_Gauss = (p+1)/2;

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
  SturmQFrameEquation<WBASIS>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;

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
  SturmQFrameEquation<WBASIS>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = j0+3;

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
  SturmQFrameEquation<WBASIS>::f(const typename WBASIS::Index& lambda) const
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
  SturmQFrameEquation<WBASIS>::RHS(const double eta,
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
  SturmQFrameEquation<WBASIS>::RHS(const double eta,
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

}



#endif // WAVELETTL_STURM_QFRAME_EQUATION_H

// implementation for row_stage_equation.h

#include <adaptive/apply.h>

namespace WaveletTL
{
  template <class ELLIPTIC_EQ, class GRAMIAN>
  void
  ROWStageEquation<ELLIPTIC_EQ,GRAMIAN>::setup_rhs
  (unsigned int i,
   const double tolerance,
   const double t_n,
   const double h,
   const InfiniteVector<double,Index>& D_un,
   const std::list<InfiniteVector<double,Index> >& Dalpha_uj,
   const int jmax)
  {
    typedef InfiniteVector<double,Index> V;
    V help, w;
    
    const unsigned int stages = row_method_->A.row_dimension(); // for readability

    // start with zero rhs
    y.clear();

    // first summand DD^{-1}<A Psi,Psi>^T D^{-1}w, where
    // w = Du^{(n)} + DD_alpha^{-1}sum_{j=1}^{i-1} a_{i,j}D_alpha u_j
    if (i > 0) {
      typename std::list<V>::const_iterator it = Dalpha_uj.begin();
      for (unsigned int j = 0; j < i; j++, ++it) {
	w.add(row_method_->A(i,j), *it);
      }
      w.scale(this, -1);    // w *= D_alpha^{-1}
      w.scale(elliptic_, 1); // w *= D
    }
    w.add(D_un);
    APPLY(*elliptic_, w, tolerance, help, jmax, St04a); // yields -D^{-1}<A Psi,Psi>D^{-1}w
    help.scale(elliptic_, 1); // help *= D
    help.scale(-1.0); // result = -D(-D^{-1}<A Psi,Psi>^T D^{-1}w)
    y = help;
//     cout << "ROWStageEquation::setup_rhs() done, y1=" << endl << y << endl;
    
    // second summand <f(t_n+alpha_i*h),Psi>^T
    if (f_ != 0) {
      f_->set_time(t_n+h*row_method_->alpha_vector[i]);
      w.clear();
      expand(f_, elliptic_->basis(), true, jmax, w); // expand in the dual (!) basis
//       cout << "ROWStageEquation::setup_rhs() done, y2=" << endl << w << endl;
      y.add(w);
    }

    // third summand <Psi,Psi>^T D_alpha^{-1}sum_{j=1}^{i-1}(c_{i,j}/h)D_alpha u_j
    if (i > 0) {
      w.clear();
      typename std::list<V>::const_iterator it = Dalpha_uj.begin();
      for (unsigned int j = 0; j < i; j++, ++it) {
	w.add(row_method_->C(i,j)/h, *it);
      }
      w.scale(this, -1); // w *= D_alpha^{-1}
      APPLY(*gramian_, w, tolerance/(4*stages), help, jmax, St04a); // yields <Psi,Psi>^T D_alpha^{-1}sum(...)
//       cout << "ROWStageEquation::setup_rhs() done, y3=" << endl << help << endl;
      y.add(help);
    }
    
    // fourth summand h*gamma_i*<f'(t_n),Psi>^T
    if (ft_ != 0) {
      ft_->set_time(t_n);
      w.clear();
      expand(ft_, elliptic_->basis(), true, jmax, w); // expand in the dual (!) basis
//       cout << "ROWStageEquation::setup_rhs() done, y4=" << endl << w << endl;
      y.add(h*row_method_->gamma_vector[i], w);
    }
    
//     cout << "ROWStageEquation::setup_rhs() done, y=" << endl << y << endl;
  }

  template <class ELLIPTIC_EQ, class GRAMIAN>
  void
  ROWStageEquation<ELLIPTIC_EQ,GRAMIAN>
  ::add_level (const Index& lambda,
	       InfiniteVector<double, Index>& w, const int j,
	       const double factor,
	       const int J,
	       const CompressionStrategy strategy) const
  {
    // We have to compute a (level-)part of the lambda-th column of 
    //   D_alpha^{-1}<(alpha*I-A)Psi,Psi>^T D_alpha^{-1}
    //   = alpha*D_alpha^{-1}<Psi,Psi>^T D_alpha^{-1}
    //     - D_alpha^{-1}D ( D^{-1}<APsi,Psi>^T D^{-1} ) DD_alpha^{-1}

    // Gramian part
    InfiniteVector<double,Index> help;
    gramian_->add_level(lambda, help, j, factor * alpha_/D(lambda), J, strategy);
    help.scale(this, -1); // help *= D_alpha^{-1}
    w.add(help);
   
    // elliptic part
    help.clear();
    elliptic_->add_level(lambda, help, j, factor*elliptic_->D(lambda)/D(lambda), J, strategy);
    help.scale(elliptic_, 1); // help *= D
    help.scale(this, -1);     // help *= D_alpha^{-1}
    w.add(help);
  }

//   template <int d, int dT>
//   void
//   ROWStageEquationFull1D<d,dT>::setup_rhs(unsigned int i,
// 					  const double tolerance,
// 					  const double t_n,
// 					  const double h,
// 					  const double alpha,
// 					  const Vector<double>& D_un,
// 					  const std::list<Vector<double> >& Dalpha_uj,
// 					  const int jmax)
//   {
//   }


}

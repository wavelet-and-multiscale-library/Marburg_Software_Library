// implementation for spline_basis.h

#include <Rd/cdf_utils.h>

namespace WaveletTL
{
  template <int d, int dT>
  SplineBasis<d,dT>::SplineBasis(const char* flavor,
				 const char* options,
				 const int s0, const int s1, const int sT0, const int sT1)
    : SplineBasisData<d,dT>(flavor,options,s0,s1,sT0,sT1)
  {
    if (SplineBasisData<d,dT>::flavor_ == "P") {
      // cf. PBasis<d,dT> ...
      DeltaLmin_       = 1-d-ell1<d>()+s0;
      DeltaRmax_offset = -1-ell1<d>()-s1;
    }
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Mj(const int j, const Vector<double>& x, Vector<double>& y) const
  {
    SplineBasisData<d,dT>::Mj0_->set_level(j);
    SplineBasisData<d,dT>::Mj1_->set_level(j);

    // decompose x=(x1 x2) appropriately
    Vector<double> x1(SplineBasisData<d,dT>::Mj0_->column_dimension());
    for (unsigned int i = 0; i < x1.size(); i++) x1[i] = x[i];
    Vector<double> x2(SplineBasisData<d,dT>::Mj1_->column_dimension());
    for (unsigned int i = 0; i < x2.size(); i++) x2[i] = x[i+x1.size()];
    SplineBasisData<d,dT>::Mj0_->apply(x1, y);
    Vector<double> yhelp(y.size());
    SplineBasisData<d,dT>::Mj1_->apply(x2,yhelp);
    y.add(yhelp);
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Mj_transposed(const int j, const Vector<double>& x, Vector<double>& y) const
  {
    SplineBasisData<d,dT>::Mj0_->set_level(j);
    SplineBasisData<d,dT>::Mj1_->set_level(j);

    // y is a block vector
    Vector<double> y1(SplineBasisData<d,dT>::Mj0_->column_dimension());
    SplineBasisData<d,dT>::Mj0_->apply_transposed(x, y1);
    Vector<double> y2(SplineBasisData<d,dT>::Mj1_->column_dimension());
    SplineBasisData<d,dT>::Mj1_->apply_transposed(x, y2);
    for (unsigned int i = 0; i < y1.size(); i++) y[i] = y1[i];
    for (unsigned int i = 0; i < y2.size(); i++) y[i+y1.size()] = y2[i];
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Gj(const int j, const Vector<double>& x, Vector<double>& y) const
  {
    SplineBasisData<d,dT>::Mj0T_->set_level(j);
    SplineBasisData<d,dT>::Mj1T_->set_level(j);

    // y is a block vector
    Vector<double> y1(SplineBasisData<d,dT>::Mj0T_->column_dimension());
    SplineBasisData<d,dT>::Mj0T_->apply_transposed(x, y1);
    Vector<double> y2(SplineBasisData<d,dT>::Mj1T_->column_dimension());
    SplineBasisData<d,dT>::Mj1T_->apply_transposed(x, y2);
    for (unsigned int i = 0; i < y1.size(); i++) y[i] = y1[i];
    for (unsigned int i = 0; i < y2.size(); i++) y[i+y1.size()] = y2[i];
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Gj_transposed(const int j, const Vector<double>& x, Vector<double>& y) const
  {
    SplineBasisData<d,dT>::Mj0T_->set_level(j);
    SplineBasisData<d,dT>::Mj1T_->set_level(j);

    // decompose x=(x1 x2) appropriately
    Vector<double> x1(SplineBasisData<d,dT>::Mj0T_->column_dimension());
    for (unsigned int i = 0; i < x1.size(); i++) x1[i] = x[i];
    Vector<double> x2(SplineBasisData<d,dT>::Mj1T_->column_dimension());
    for (unsigned int i = 0; i < x2.size(); i++) x2[i] = x[i+x1.size()];
    SplineBasisData<d,dT>::Mj0T_->apply(x1, y);
    Vector<double> yhelp(y.size());
    SplineBasisData<d,dT>::Mj1T_->apply(x2,yhelp);
    y.add(yhelp);
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Tj(const int j, const Vector<double>& x, Vector<double>& y) const
  { 
    y = x;
    for (int k = j0(); k <= j; k++) {
      // decompose current vector appropriately
      Vector<double> y1(Deltasize(k+1));
      for (unsigned int i = 0; i < y1.size(); i++) y1[i] = y[i];
      Vector<double> Mky1(Deltasize(k+1));
      apply_Mj(k, y1, Mky1);
      for (unsigned int i = 0; i < Mky1.size(); i++) y[i] = Mky1[i];
    }
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Tjinv(const int j, const Vector<double>& x, Vector<double>& y) const
  { 
    y = x;

    // T_j^{-1}=diag(G_{j0},I)*...*diag(G_{j-1},I)*G_j
    for (int k = j; k >= j0(); k--) {
      // decompose current vector appropriately
      Vector<double> y1(Deltasize(k+1));
      for (unsigned int i = 0; i < y1.size(); i++) y1[i] = y[i];
      Vector<double> Gky1(Deltasize(k+1));
      apply_Gj(k, y1, Gky1);
      for (unsigned int i = 0; i < Gky1.size(); i++) y[i] = Gky1[i];
    }
  }

  template <int d, int dT>
  void
  SplineBasis<d,dT>::apply_Tj_transposed(const int j, const Vector<double>& x, Vector<double>& y) const
  { 
    y = x;

    for (int k = j; k >= j0(); k--) {
      // decompose current vector appropriately
      Vector<double> y1(Deltasize(k+1));
      for (unsigned int i = 0; i < y1.size(); i++) y1[i] = y[i];
      Vector<double> MkTy1(Deltasize(k+1));
      apply_Mj_transposed(k, y1, MkTy1);
      for (unsigned int i = 0; i < MkTy1.size(); i++) y[i] = MkTy1[i];
    }
  }
}

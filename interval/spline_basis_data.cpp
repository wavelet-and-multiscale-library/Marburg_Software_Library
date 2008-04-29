// implementation for spline_basis_data.h

#include <sstream>
#include <cstring>

namespace WaveletTL
{
  // the different values of j0 for the available spline wavelet bases

  template <>
  const int
  SplineBasisData_j0<2,2,DS_construction_bio5,0,0,0,0>::j0 = 3; // 2 does not work with composite bases

  template <>
  const int
  SplineBasisData_j0<2,2,DS_construction_bio5e,0,0,0,0>::j0 = 3; // 2 does not work with composite bases

  template <>
  const int
  SplineBasisData_j0<3,3,DS_construction_bio5,0,0,0,0>::j0 = 4;

  template <>
  const int
  SplineBasisData_j0<3,3,DS_construction_bio5e,0,0,0,0>::j0 = 4;

  template <>
  const int
  SplineBasisData_j0<1,1,P_construction,0,0,0,0>::j0 = 1; // for technical reasons, QuasiStationaryMatrix should have at least 2 columns

  template <>
  const int
  SplineBasisData_j0<1,3,P_construction,0,0,0,0>::j0 = 3;

  template <>
  const int
  SplineBasisData_j0<2,2,P_construction,0,0,0,0>::j0 = 2;

  template <>
  const int
  SplineBasisData_j0<2,2,P_construction,1,0,0,0>::j0 = 2;

  template <>
  const int
  SplineBasisData_j0<2,2,P_construction,0,1,0,0>::j0 = 2;

  template <>
  const int
  SplineBasisData_j0<2,2,P_construction,1,1,0,0>::j0 = 3;

  template <>
  const int
  SplineBasisData_j0<3,3,P_construction,0,0,0,0>::j0 = 3;

  template <>
  const int
  SplineBasisData_j0<3,3,P_construction,1,0,0,0>::j0 = 3;

  template <>
  const int
  SplineBasisData_j0<3,3,P_construction,0,1,0,0>::j0 = 3;

  template <>
  const int
  SplineBasisData_j0<3,3,P_construction,1,1,0,0>::j0 = 3;


  //
  // first the "generic" stuff for SplineBasisData: destructor and check routine for all flavors

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::~SplineBasisData()
  {
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::~SplineBasisData()
  {
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  void
  SplineBasisData<d,dT,flavor,s0,s1,sT0,sT1>::check() const
  {
    cout << endl
 	 << "-------------------------------" << endl
 	 << "SplineBasisData::check() called" << endl
 	 << "-------------------------------" << endl << endl;
    
    cout << "* some basic data:" << endl;
    cout << "  d=" << d << endl;
    cout << "  dT=" << dT << endl;
    cout << "  flavor: " << flavor << endl;
    cout << "  s0=" << s0 << endl;
    cout << "  s1=" << s1 << endl;
    cout << "  sT0=" << sT0 << endl;
    cout << "  sT1=" << sT1 << endl;

    cout << "Mj0=" << endl << Mj0_;
    cout << "Mj0T=" << endl << Mj0T_;
    cout << "Mj1=" << endl << Mj1_;
    cout << "Mj1T=" << endl << Mj1T_;
    cout << "Mj1c=" << endl << Mj1c_;
    
    cout << "CLA=" << endl << CLA_;
    cout << "CRA=" << endl << CRA_;
    cout << "CLAT=" << endl << CLAT_;
    cout << "CRAT=" << endl << CRAT_;

    cout << "check biorthogonality of the generators on different levels:" << endl;
    for (int j = SplineBasisData_j0<d,dT,flavor,s0,s1,sT0,sT1>::j0; j <= 8; j++) {
      cout << "* j=" << j;
      Mj0_.set_level(j);
      Mj0T_.set_level(j);
      Vector<double> x(Mj0T_.column_dimension()),
	y(Mj0T_.row_dimension());
      double maxerr = 0;
      for (unsigned int i = 0; i < Mj0T_.column_dimension(); i++) {
	x.scale(0);
	x[i] = 1.0;
	Mj0T_.apply(x,y);
	Mj0_.apply_transposed(y,x);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj0^T*Mj0T-I||_1: " << maxerr;
      
      maxerr = 0;
      for (unsigned int i = 0; i < Mj0_.column_dimension(); i++) {
	x.scale(0);
	x[i] = 1.0;
	Mj0_.apply(x,y);
	Mj0T_.apply_transposed(y,x);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj0T^T*Mj0-I||_1: " << maxerr << endl;
    }

    cout << endl << "check biorthogonality of the full systems on different levels:" << endl;
    for (int j = SplineBasisData_j0<d,dT,flavor,s0,s1,sT0,sT1>::j0; j <= 8; j++) {
      cout << "* j=" << j;
      Mj0_.set_level(j);
      Mj1_.set_level(j);
      Mj0T_.set_level(j);
      Mj1T_.set_level(j);
      Vector<double> x(Mj0T_.row_dimension()),
	y1(Mj0T_.column_dimension()),
	y2(Mj1T_.column_dimension()),
	xhelp(Mj0T_.row_dimension());
      
      double maxerr = 0;
      for (unsigned int i = 0; i < Mj0T_.row_dimension(); i++) {
	x.scale(0);
	x[i] = 1;
	Mj0T_.apply_transposed(x,y1);
	Mj1T_.apply_transposed(x,y2);
	Mj0_.apply(y1,x);
	Mj1_.apply(y2,xhelp);
	x.add(xhelp);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj*Gj-I||_1: " << maxerr;
      
      maxerr = 0;
      for (unsigned int i = 0; i < Mj0_.row_dimension(); i++) {
	x.scale(0);
	x[i] = 1;
	Mj0_.apply_transposed(x,y1);
	Mj1_.apply_transposed(x,y2);
	Mj0T_.apply(y1,x);
	Mj1T_.apply(y2,xhelp);
	x.add(xhelp);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Gj*Mj-I||_1: " << maxerr << endl;
    }
  }

  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  void
  SplineBasisData<d,dT,P_construction,s0,s1,sT0,sT1>::check() const
  {
    cout << endl
 	 << "-------------------------------" << endl
 	 << "SplineBasisData::check() called" << endl
 	 << "-------------------------------" << endl << endl;
    
    cout << "* some basic data:" << endl;
    cout << "  d=" << d << endl;
    cout << "  dT=" << dT << endl;
    cout << "  flavor: " << P_construction << endl;
    cout << "  s0=" << s0 << endl;
    cout << "  s1=" << s1 << endl;
    cout << "  sT0=" << sT0 << endl;
    cout << "  sT1=" << sT1 << endl;

    cout << "Mj0=" << endl << Mj0_;
    cout << "Mj0T=" << endl << Mj0T_;
    cout << "Mj1=" << endl << Mj1_;
    cout << "Mj1T=" << endl << Mj1T_;

    cout << "check biorthogonality of the generators on different levels:" << endl;
    for (int j = SplineBasisData_j0<d,dT,P_construction,s0,s1,sT0,sT1>::j0; j <= 8; j++) {
      cout << "* j=" << j;
      Mj0_.set_level(j);
#if 0
      {
	SparseMatrix<double> S;
	Mj0_.to_sparse(S);
	cout << endl << "Mj0(j)=" << endl << S << endl;
      }
#endif 
      Mj0T_.set_level(j);
#if 0
      {
	SparseMatrix<double> S;
	Mj0T_.to_sparse(S);
	cout << endl << "Mj0T(j)=" << endl << S << endl;
      }
#endif 
      Vector<double> x(Mj0T_.column_dimension()),
	y(Mj0T_.row_dimension());
      double maxerr = 0;
      for (unsigned int i = 0; i < Mj0T_.column_dimension(); i++) {
	x.scale(0);
	x[i] = 1.0;
	Mj0T_.apply(x,y);
	Mj0_.apply_transposed(y,x);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj0^T*Mj0T-I||_1: " << maxerr;
      
      maxerr = 0;
      for (unsigned int i = 0; i < Mj0_.column_dimension(); i++) {
	x.scale(0);
	x[i] = 1.0;
	Mj0_.apply(x,y);
	Mj0T_.apply_transposed(y,x);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj0T^T*Mj0-I||_1: " << maxerr << endl;
    }

    cout << endl << "check biorthogonality of the full systems on different levels:" << endl;
    for (int j = SplineBasisData_j0<d,dT,P_construction,s0,s1,sT0,sT1>::j0; j <= 8; j++) {
      cout << "* j=" << j;
      Mj0_.set_level(j);
      Mj1_.set_level(j);
      Mj0T_.set_level(j);
      Mj1T_.set_level(j);
      Vector<double> x(Mj0T_.row_dimension()),
	y1(Mj0T_.column_dimension()),
	y2(Mj1T_.column_dimension()),
	xhelp(Mj0T_.row_dimension());
      
      double maxerr = 0;
      for (unsigned int i = 0; i < Mj0T_.row_dimension(); i++) {
	x.scale(0);
	x[i] = 1;
	Mj0T_.apply_transposed(x,y1);
	Mj1T_.apply_transposed(x,y2);
	Mj0_.apply(y1,x);
	Mj1_.apply(y2,xhelp);
	x.add(xhelp);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Mj*Gj-I||_1: " << maxerr;
      
      maxerr = 0;
      for (unsigned int i = 0; i < Mj0_.row_dimension(); i++) {
	x.scale(0);
	x[i] = 1;
	Mj0_.apply_transposed(x,y1);
	Mj1_.apply_transposed(x,y2);
	Mj0T_.apply(y1,x);
	Mj1T_.apply(y2,xhelp);
	x.add(xhelp);
	x[i] = x[i] - 1.0;
	maxerr = std::max(maxerr, linfty_norm(x));
      }
      cout << ", ||Gj*Mj-I||_1: " << maxerr << endl;
    }
  }

  //
  //
  // some precomputed data for DS bases

  template <>
  SplineBasisData<2,2,DS_construction_bio5,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,DS_construction_bio5,0,0,0,0>::j0;
    std::ostringstream entries;
    entries << "0.707106781186547 0 "
	    << "0.288735268984507 0.353553390593274 "
	    << "-0.117851130197758 1.01015254455221 "
	    << "-0.0589255650988791 0.505076272276106";
    Matrix<double> Mj0_l(4, 2, entries.str().c_str());
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(3, "0.35355339059327 0.70710678118655 0.35355339059327");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
//     Mj0_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
    entries.str("");
    entries << "0.707106781186547 0 "
	    << "1.51522881682832 0.353553390593274 "
	    << "-0.53033008588991 0.742462120245874 "
	    << "0 0.247487373415291 "
	    << "0 -0.123743686707646";
    Matrix<double> Mj0T_l(5, 2, entries.str().c_str());
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(5, "-0.17677669529664 0.35355339059327 1.06066017177982 0.35355339059327 -0.17677669529664");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
//     Mj0T_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
    entries.str("");
    entries << "0.375 0 "
	    << "-0.109375 -0.0875 "
	    << "0.1875 -0.25 "
	    << "-0.34375 0.75 "
	    << "0.125 -0.25 "
	    << "0.0625 -0.125";
    Matrix<double> Mj1_l(6, 2, entries.str().c_str());
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
//     Mj1_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
    entries.str("");
    entries << "1.33333333333333 0.666666666666667 "
	    << "-2.85714285714286 -1.42857142857143 "
	    << "1 0 "
	    << "0 1 "
	    << "0 -0.5";
    Matrix<double> Mj1T_l(5, 2, entries.str().c_str());
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_lr(3, "-0.5 1 -0.5");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
//     Mj1T_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
    entries.str("");
    entries << "-0.2625 -0.0875 "
	    << "0.25 -0.25 "
	    << "-0.3125 0.75 "
	    << "0.125 -0.25 "
	    << "0.0625 -0.125";
    Matrix<double> Mj1c_l(5, 2, entries.str().c_str());
    Matrix<double> Mj1c_r; Mj1c_l.reflect(Mj1c_r);
    Mj1c_ = QuasiStationaryMatrix<double>(j0, 15, 8, Mj1c_l, Mj1c_r, Mj1_band_lr, Mj1_band_lr, 2, 2, 1.0);
//     Mj1c_ = QuasiStationaryMatrix<double>(j0, 7, 4, Mj1c_l, Mj1c_r, Mj1_band_lr, Mj1_band_lr, 2, 2, 1.0);
    CLA_ = Matrix<double>(2, 2, "1 0 -0.166666666666667 1.42857142857143");
    CRA_ = CLA_;
    CLAT_ = Matrix<double>(3, 2, "6 -0.7 3 0 0 0.7");
    CRAT_ = CLAT_;
  }

  template <>
  SplineBasisData<2,2,DS_construction_bio5e,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,DS_construction_bio5e,0,0,0,0>::j0;
    std::ostringstream entries;
    entries << "0.707106781186547 0 "
	    << "0.288735268984507 0.353553390593274 "
	    << "-0.117851130197758 1.01015254455221 "
	    << "-0.0589255650988791 0.505076272276106";
    Matrix<double> Mj0_l(4, 2, entries.str().c_str());
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(3, "0.35355339059327 0.70710678118655 0.35355339059327");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
//     Mj0_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
    entries.str("");
    entries << "0.707106781186547 0 "
	    << "1.51522881682832 0.353553390593274 "
	    << "-0.53033008588991 0.742462120245874 "
	    << "0 0.247487373415291 "
	    << "0 -0.123743686707646";
    Matrix<double> Mj0T_l(5, 2, entries.str().c_str());
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(5, "-0.17677669529664 0.35355339059327 1.06066017177982 0.35355339059327 -0.17677669529664");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
//     Mj0T_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
    entries.str("");
    entries << "0.375 0 "
	    << "-0.109375 -0.00380794872403089 "
	    << "0.1875 -0.0108798534972311 "
	    << "-0.34375 0.0326395604916934 "
	    << "0.125 -0.0108798534972311 "
	    << "0.0625 -0.00543992674861556";
    Matrix<double> Mj1_l(6, 2, entries.str().c_str());
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
//     Mj1_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
    entries.str("");
    entries << "1.33333333333333 15.3188337241014 "
	    << "-2.85714285714286 -32.8260722659316 "
	    << "1 0 "
	    << "0 22.9782505861521 "
	    << "0 -11.4891252930761";
    Matrix<double> Mj1T_l(5, 2, entries.str().c_str());
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_lr(3, "-0.5 1 -0.5");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
//     Mj1T_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
    entries.str("");
    entries << "-0.2625 -0.00380794872403089 "
	    << "0.25 -0.0108798534972311 "
	    << "-0.3125 0.0326395604916934 "
	    << "0.125 -0.0108798534972311 "
	    << "0.0625 -0.00543992674861556 ";
    Matrix<double> Mj1c_l(5, 2, entries.str().c_str());
    Matrix<double> Mj1c_r; Mj1c_l.reflect(Mj1c_r);
    Mj1c_ = QuasiStationaryMatrix<double>(j0, 15, 8, Mj1c_l, Mj1c_r, Mj1_band_lr, Mj1_band_lr, 2, 2, 1.0);
//     Mj1c_ = QuasiStationaryMatrix<double>(j0, 7, 4, Mj1c_l, Mj1c_r, Mj1_band_lr, Mj1_band_lr, 2, 2, 1.0);
    CLA_ = Matrix<double>(2, 2, "1 0 -0.166666666666667 1.42857142857143");
    CRA_ = CLA_;
    CLAT_ = Matrix<double>(3, 2, "6 -0.7 3 0 0 0.7");
    CRAT_ = CLAT_;
  }

  template <>
  SplineBasisData<3,3,DS_construction_bio5,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,DS_construction_bio5,0,0,0,0>::j0;
    std::ostringstream entries;
    entries << "0.707106781186398 0 0 "
	    << "0.690720459146663 0.353553390593505 0 "
	    << "-0.32389524472094 0.353553390593443 0.176776695296556 "
	    << "-0.142570613712759 0.343481380196646 0.258256676839446 "
	    << "-0.100706298224861 0.276334644218305 0.387385015259226 "
	    << "0.0242033624305832 0.15753657287204 0.54233902136297 "
	    << "0.064993644568729 0.073603152899181 0.464862018311131 "
	    << "0.0216645481895764 0.024534384299727 0.15495400610371";
    Matrix<double> Mj0_l(8, 3, entries.str().c_str());
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(4, "0.176776695296637 0.530330085889911 0.530330085889911 0.17677669529663");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 30, 14, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 6, 6, 1.0);
    entries.str("");
    entries << "0.707106781185991 0 0 "
	    << "0.538529735380284 0.874897243093867 -0.307592872975334 "
	    << "-0.0387385015263574 1.18325516575686 -0.3445671572034 "
	    << "-0.344714555828851 0.71808950353182 -0.0606618632922345 "
	    << "-0.428683486094793 0.609367134416539 0.0160888811452143 "
	    << "0.260745625562884 -0.597886052237867 1.08129556773563 "
	    << "0.278423295092582 -0.70452095186605 1.09631288121043 "
	    << "-0.0574524259714819 0.13603342823686 -0.165970639943909 "
	    << "-0.0397747564418011 0.111157750183301 -0.224400360517013 "
	    << "0.013258252147267 -0.0370525833944338 0.0748001201723378";
    Matrix<double> Mj0T_l(10, 3, entries.str().c_str());
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    entries.str("");
    entries << "0.0662912607362388 -0.198873782208717 "
	    << "-0.154679608384557 0.994368911043582 "
	    << "0.994368911043582 -0.154679608384557 "
	    << "-0.198873782208717 0.0662912607362388";
    Vector<double> Mj0T_band_lr(8, entries.str().c_str());
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 30, 14, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 4, 4, 1.0);
    entries.str("");
    entries << "0.24375 0.303125 -0.184375 -0.20625 -0.025 "
	    << "-0.0157815824466431 0.0806566378547729 0.0312823027480479 0.0607154255317074 0.0105128546099001 "
	    << "-0.35481078041896 -0.357136498365546 0.104691052055713 0.149634617686268 0.0111238220301554 "
	    << "0.719869791666666 -0.274578993055555 -0.036714409722222 -0.00614583333333299 -0.012534722222222 "
	    << "-0.2096484375 0.782207031249999 -0.227402343750001 -0.219375000000001 -0.04640625 "
	    << "-0.0718828125 -0.09434765625 0.50144921875 -0.525375 -0.09346875 "
	    << "-0.00225 -0.0361875 -0.4404375 0.5810625 0.0789166666666665 "
	    << "-0.00075 -0.0433125 -0.0530625 0.0999375 0.47075 "
	    << "0 -0.03515625 0.10546875 -0.10546875 -0.46875 "
	    << "0 -0.01171875 0.03515625 -0.03515625 -0.0729166666666667 "
	    << "0 0 0 0 0.09375 "
	    << "0 0 0 0 0.03125";
    Matrix<double> Mj1_l(12, 5, entries.str().c_str());
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    entries.str("");
    entries << "-0.03125 -0.09375 "
	    << "0.0729166666666666 0.46875 "
	    << "-0.46875 -0.0729166666666667 "
	    << "0.09375 0.03125";
    Vector<double> Mj1_band_l(8, entries.str().c_str());
    entries.str("");
    entries << "0.03125 0.09375 "
	    << "-0.0729166666666666 -0.46875 "
	    << "0.46875 0.0729166666666667 "
	    << "-0.09375 -0.03125";
    Vector<double> Mj1_band_r(8, entries.str().c_str());
    Mj1_ = QuasiStationaryMatrix<double>(j0, 30, 16, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 6, 6, 1.0);
    entries.str("");
    entries << "-0.945626477541817 -2.23847517730569 -4.00109929078116 -3.25494089834593 "
	    << "0.489408327246703 1.40978816654599 2.62235208181327 2.1524226929646 "
	    << "-1.46092037983956 -2.19138056975949 -3.06793279766344 -2.33747260774361 "
	    << "1 0 0 0 "
	    << "0 1 0 0 "
	    << "0 0 1 0 "
	    << "0 0 0 1 "
	    << "0 0 0 -0.333333333333333";
    Matrix<double> Mj1T_l(8, 4, entries.str().c_str());
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_l(4, "-0.375 1.125 -1.125 0.375");
    Vector<double> Mj1T_band_r(4, "0.375 -1.125 1.125 -0.375");
    Mj1T_ = QuasiStationaryMatrix<double>(j0,30, 16, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 6, 6, 1.0);
    entries.str("");
    entries << "-0.253882978723144 -0.215443816489063 0.211384640957205 0.262185837765684 0.0349335106382616 "
	    << "-0.24315937499993 -0.218287955729068 0.0202367838541345 0.055160351562466 -0.000327604166669463 "
	    << "0.769015957446905 -0.213461325354493 -0.0738890735816355 -0.0477310505319997 -0.0175753546099395 "
	    << "-0.174933510638256 0.825378158244729 -0.253661070478757 -0.248749168883016 -0.0499667553191536 "
	    << "-0.0802260638298284 -0.104723238031966 0.507760139627693 -0.518315325797839 -0.0926130319148893 "
	    << "-0.0246542553192109 -0.064049202127735 -0.423490691489316 0.600019946808566 0.0812145390070988 "
	    << "-0.00821808510640363 -0.0525997340425783 -0.0474135638297715 0.106256648936188 0.471515957446811 "
	    << "0 -0.03515625 0.10546875 -0.10546875 -0.46875 "
	    << "0 -0.01171875 0.03515625 -0.03515625 -0.0729166666666666 "
	    << "0 0 0 0 0.09375 "
	    << "0 0 0 0 0.03125";
    Matrix<double> Mj1c_l(11, 5, entries.str().c_str());
    Matrix<double> Mj1c_r; Mj1c_l.reflect(Mj1c_r);
    Mj1c_ = QuasiStationaryMatrix<double>(j0, 28, 16, Mj1c_l, Mj1c_r, Mj1_band_l, Mj1_band_r, 5, 5, 1.0);
    entries.str("");
    entries << "1.52801418439628 -0.270270270270127 0 "
	    << "0.471985815602211 0.270270270270539 0 "
	    << "-0.114267139480271 0.518626734843214 0.146092037983819 "
	    << "-0.230744680851163 0.474799123447897 0.438276113951726 "
	    << "0.122553191489537 0.13878743608459 0.876552227903642";
    CLA_ = Matrix<double>(5, 3, entries.str().c_str());
    CRA_ = CLA_;
    entries.str("");
    entries << "11.6000000000102 -6.88829787232354 3.49367907800638 "
	    << "7.20000000000617 -3.47170212764972 1.7060542553146 "
	    << "3.80000000000308 -0.999787234037989 0.475791134749664 "
	    << "1.40000000000094 0.527446808511657 -0.197110283688414 "
	    << "0 1.10999999999922 -0.312649999999638 "
	    << "-0.400000000000452 0.747872340424688 0.129171985815993 "
	    << "0.200000000000291 -0.558936170211928 1.12835567375848";
    CLAT_ = Matrix<double>(7, 3, entries.str().c_str());
    CRAT_ = CLAT_;
  }

  template <>
  SplineBasisData<3,3,DS_construction_bio5e,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,DS_construction_bio5e,0,0,0,0>::j0;
    std::ostringstream entries;
    entries << "0.707106781186398 0 0 "
	    << "0.690720459146663 0.353553390593505 0 "
	    << "-0.32389524472094 0.353553390593443 0.176776695296556 "
	    << "-0.142570613712759 0.343481380196646 0.258256676839446 "
	    << "-0.100706298224861 0.276334644218305 0.387385015259226 "
	    << "0.0242033624305832 0.15753657287204 0.54233902136297 "
	    << "0.064993644568729 0.073603152899181 0.464862018311131 "
	    << "0.0216645481895764 0.024534384299727 0.15495400610371";
    Matrix<double> Mj0_l(8, 3, entries.str().c_str());
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(4, "0.176776695296637 0.530330085889911 0.530330085889911 0.17677669529663");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 30, 14, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 6, 6, 1.0);
    entries.str("");
    entries << "0.707106781185991 0 0 "
	    << "0.538529735380284 0.874897243093867 -0.307592872975334 "
	    << "-0.0387385015263574 1.18325516575686 -0.3445671572034 "
	    << "-0.344714555828851 0.71808950353182 -0.0606618632922345 "
	    << "-0.428683486094793 0.609367134416539 0.0160888811452143 "
	    << "0.260745625562884 -0.597886052237867 1.08129556773563 "
	    << "0.278423295092582 -0.70452095186605 1.09631288121043 "
	    << "-0.0574524259714819 0.13603342823686 -0.165970639943909 "
	    << "-0.0397747564418011 0.111157750183301 -0.224400360517013 "
	    << "0.013258252147267 -0.0370525833944338 0.0748001201723378";
    Matrix<double> Mj0T_l(10, 3, entries.str().c_str());
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    entries.str("");
    entries << "0.0662912607362388 -0.198873782208717 "
	    << "-0.154679608384557 0.994368911043582 "
	    << "0.994368911043582 -0.154679608384557 "
	    << "-0.198873782208717 0.0662912607362388";
    Vector<double> Mj0T_band_lr(8, entries.str().c_str());
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 30, 14, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 4, 4, 1.0);
    entries.str("");
    entries << "0.24375 0.000513179450639216 0.0110376597715015 0.000834764817451311 0.0175720391268678 "
	    << "-0.0157815824466431 0.000915259105101272 0.00123569370492313 0.000259761184582998 -0.0118511099734895 "
	    << "-0.35481078041896 -0.000503380263708669 -0.0113304911719982 -0.00160180091504651 -0.00290194552658878 "
	    << "0.719869791666666 -0.000698057837431133 -0.00672020092585003 -0.00167785717416517 0.0181771060099316 "
	    << "-0.2096484375 0.00271451003664066 0.0238002244077162 0.00135890841065896 -0.00232059603222293 "
	    << "-0.0718828124999998 -0.0179463859220804 0.00021907041986547 0.00624722305865691 -0.00110767705768198 "
	    << "-0.00225 0.0172690055281071 -0.00443602034104831 -0.007412749888786 -0.00218678975094893 "
	    << "-0.00075 0.0060159614050828 -0.00355744700822455 0.0175012339473863 9.61814633350482e-05 "
	    << "0 -0.00724505933141968 0.00128023008453561 -0.0171734824977549 -0.000485103251330258 "
	    << "0 -0.0017666637479702 0.000115444407456376 -0.00231441658339277 -4.01222834164913e-05 "
	    << "0 0.000729400532815903 -0.00035021132331243 0.00383633728034121 0.000136776150405294 "
	    << "0 0.000243133510938634 -0.00011673710777081 0.00127877909344707 4.55920501350979e-05";
    Matrix<double> Mj1_l(12, 5, entries.str().c_str());
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    entries.str("");
    entries << "-0.03125 -0.09375 "
	    << "0.0729166666666666 0.46875 "
	    << "-0.46875 -0.0729166666666667 "
	    << "0.09375 0.03125";
    Vector<double> Mj1_band_l(8, entries.str().c_str());
    entries.str("");
    entries << "0.03125 0.09375 "
	    << "-0.0729166666666666 -0.46875 "
	    << "0.46875 0.0729166666666667 "
	    << "-0.09375 -0.03125";
    Vector<double> Mj1_band_r(8, entries.str().c_str());
    Mj1_ = QuasiStationaryMatrix<double>(j0, 30, 16, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 6, 6, 1.0);
    entries.str("");
    entries << "-0.945626477541817 1.44619498325877 -1.25892531388003 -1.9019389395207 42.4104198710965 "
	    << "0.489408327246703 -0.796571535511957 -1.23153895840456 1.02701738908654 -27.7114414458721 "
	    << " -1.46092037983956 2.36231004319434 -18.0058010150959 -3.27727100838957 33.2207868750493 "
	    << "1 0 0 0 0 "
	    << "0 6.44604880527398 31.7402843965424 1.77743764143808 -2.95954169217808 "
	    << "0 -26.1398430757433 -7.3489082723521 4.48898186486304 -5.32624950014288 "
	    << "0 22.5045317214334 -10.6754664638928 -14.2576710910886 -4.47813106503936 "
	    << "0 5.16573154667264 -1.06151465505845 26.356802651387 1.57573964945565 "
	    << "0 -14.2506473855442 5.19750391090056 -24.3047763236522 -0.093407956247842 "
	    << "0 4.75021579518141 -1.73250130363352 8.10159210788406 0.0311359854159473";
    Matrix<double> Mj1T_l(10, 5, entries.str().c_str());
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_l(4, "-0.375 1.125 -1.125 0.375");
    Vector<double> Mj1T_band_r(4, "0.375 -1.125 1.125 -0.375");
    Mj1T_ = QuasiStationaryMatrix<double>(j0,30, 16, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 8, 8, 1.0);
    entries.str("");
    entries << "-0.253882978723144 0.000413971951306335 -0.00954618200058477 -0.000555658994278584 -0.029015938902841 "
	    << "-0.24315937499993 -0.000268314799519983 -0.00627461332164192 -0.00121943100709685 0.00514705084096629 "
	    << "0.769015957446905 -0.0005945878790471 -0.00449472951269848 -0.00150954747116643 0.0217200735915737 "
	    << "-0.174933510638256 0.00278759716318785 0.0253722102608124 0.00147779579356996 0.000182017518958457 "
	    << "-0.0802260638298284 -0.0179639513995531 -0.000158734583485131 0.00621865012461265 -0.00170914552460688 "
	    << "-0.0246542553192109 0.0172218366934951 -0.00545054566472882 -0.00748947720817742 -0.00380192185793247 "
	    << "-0.00821808510640363 0.00600023846021212 -0.00389562211611805 0.0174756581742558 -0.000442195905659466 "
	    << "0 -0.00724505933141968 0.00128023008453561 -0.0171734824977548 -0.00048510325133026 "
	    << "0 -0.0017666637479702 0.000115444407456376 -0.00231441658339276 -4.0122283416492e-05 "
	    << "0 0.000729400532815903 -0.00035021132331243 0.00383633728034121 0.000136776150405294 "
	    << "0 0.000243133510938634 -0.00011673710777081 0.00127877909344707 4.5592050135098e-05";
    Matrix<double> Mj1c_l(11, 5, entries.str().c_str());
    Matrix<double> Mj1c_r; Mj1c_l.reflect(Mj1c_r);
    Mj1c_ = QuasiStationaryMatrix<double>(j0, 28, 16, Mj1c_l, Mj1c_r, Mj1_band_l, Mj1_band_r, 5, 5, 1.0);
    entries.str("");
    entries << "1.52801418439628 -0.270270270270127 0 "
	    << "0.471985815602211 0.270270270270539 0 "
	    << "-0.114267139480271 0.518626734843214 0.146092037983819 "
	    << "-0.230744680851163 0.474799123447897 0.438276113951726 "
	    << "0.122553191489537 0.13878743608459 0.876552227903642";
    CLA_ = Matrix<double>(5, 3, entries.str().c_str());
    CRA_ = CLA_;
    entries.str("");
    entries << "11.6000000000102 -6.88829787232354 3.49367907800638 "
	    << "7.20000000000617 -3.47170212764972 1.7060542553146 "
	    << "3.80000000000308 -0.999787234037989 0.475791134749664 "
	    << "1.40000000000094 0.527446808511657 -0.197110283688414 "
	    << "0 1.10999999999922 -0.312649999999638 "
	    << "-0.400000000000452 0.747872340424688 0.129171985815993 "
	    << "0.200000000000291 -0.558936170211928 1.12835567375848";
    CLAT_ = Matrix<double>(7, 3, entries.str().c_str());
    CRAT_ = CLAT_;
  }

  //
  //
  // some precomputed data for P bases

  template <>
  SplineBasisData<1,1,P_construction,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<1,1,P_construction,0,0,0,0>::j0;
    Matrix<double> Mj0_lr(0); // empty corner blocks
    Vector<double> Mj0_band_lr(2, "1.0 1.0");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 4, 2, Mj0_lr, Mj0_lr, Mj0_band_lr, Mj0_band_lr, 0, 0, M_SQRT1_2);
    Matrix<double> Mj0T_lr(0); // empty corner blocks
    Vector<double> Mj0T_band_lr(2, "1.0 1.0");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 4, 2, Mj0T_lr, Mj0T_lr, Mj0T_band_lr, Mj0T_band_lr, 0, 0, M_SQRT1_2);
    Matrix<double> Mj1_lr(0); // empty corner blocks
    Vector<double> Mj1_band_l(2, "1.0 -1.0");
    Vector<double> Mj1_band_r(2, "-1.0 1.0"); // note: Haar wavelets are here antisymmetric w.r.t 0.5
    Mj1_ = QuasiStationaryMatrix<double>(j0, 4, 2, Mj1_lr, Mj1_lr, Mj1_band_l, Mj1_band_r, 0, 0, M_SQRT1_2);
    Matrix<double> Mj1T_lr(0); // empty corner blocks
    Vector<double> Mj1T_band_l(2, "1.0 -1.0");
    Vector<double> Mj1T_band_r(2, "-1.0 1.0");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 4, 2, Mj1T_lr, Mj1T_lr, Mj1T_band_l, Mj1T_band_r, 0, 0, M_SQRT1_2);
    CDF_factor = 1.0;
  }

  template <>
  SplineBasisData<1,3,P_construction,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<1,3,P_construction,0,0,0,0>::j0;
    Matrix<double> Mj0_lr(0); // empty corner blocks
    Vector<double> Mj0_band_lr(2, "1.0 1.0");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj0_lr, Mj0_lr, Mj0_band_lr, Mj0_band_lr, 0, 0, M_SQRT1_2);
    std::ostringstream entries;
    entries << "1.375 -0.5 0.125 "
	    << "0.625 0.5 -0.125 "
	    << "0.125 1 -0.125 "
	    << "-0.125 1 0.125 "
	    << "0 0.125 1 "
	    << "0 -0.125 1 "
	    << "0 0 0.125 "
	    << "0 0 -0.125";
    Matrix<double> Mj0T_l(8, 3, entries.str().c_str());
    Matrix<double> Mj0T_r(8, 3); Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(6, "-0.125 0.125 1.0 1.0 0.125 -0.125");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 4, 4, M_SQRT1_2);
    entries.str("");
    entries << "0.3125 -0.0625 "
	    << "-0.6875 -0.0625 "
	    << "0.25 0.5 "
	    << "0.25 -0.5 "
	    << "-0.0625 0.0625 "
	    << "-0.0625 0.0625";
    Matrix<double> Mj1_l(6, 2, entries.str().c_str());
    Matrix<double> Mj1_r(6, 2); Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_l(6, "-0.125 -0.125 1.0 -1.0 0.125 0.125");
    Vector<double> Mj1_band_r(6, "0.125 0.125 -1.0 1.0 -0.125 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 2, 2, M_SQRT1_2);
    entries.str("");
    entries << "2 0 "
	    << "-2 0 "
	    << "0 2 "
	    << "0 -2";
    Matrix<double> Mj1T_l(4, 2, entries.str().c_str());
    Matrix<double> Mj1T_r(4, 2); Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_l(2, "1.0 -1.0");
    Vector<double> Mj1T_band_r(2, "-1.0 1.0");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 4, 4, M_SQRT1_2);
    CDF_factor = 1.0;
  }

  template <>
  SplineBasisData<2,2,P_construction,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,P_construction,0,0,0,0>::j0;
    Matrix<double> Mj0_l(2, 1, "1 0.5");
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 1, M_SQRT1_2);
    Matrix<double> Mj0T_l(5, 2, "1.25 -0.125 1.5 0.25 -0.75 1.625 0 0.5 0 -0.25");
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, M_SQRT1_2);
    Matrix<double> Mj1_l(4, 1, "-0.75 0.5625 -0.125 -0.0625");
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 1, 1, M_SQRT1_2);
    Matrix<double> Mj1T_lr(0); // empty corner blocks
    Vector<double> Mj1T_band_lr(3, "-1 2 -1");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 9, 4, Mj1T_lr, Mj1T_lr, Mj1T_band_lr, Mj1T_band_lr, 0, 0, M_SQRT1_2);
    CDF_factor = -0.5;
  }

  template <>
  SplineBasisData<2,2,P_construction,1,1,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,P_construction,1,1,0,0>::j0;
    Matrix<double> Mj0_lr(0); // empty corner blocks
    Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 15, 7, Mj0_lr, Mj0_lr, Mj0_band_lr, Mj0_band_lr, 0, 0, M_SQRT1_2);
    Matrix<double> Mj0T_l(6, 2, "1.5 -0.5 1 0 0.5 0.5 -0.25 1.5 0 0.5 0 -0.25");
    Matrix<double> Mj0T_r; Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 15, 7, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 3, M_SQRT1_2);
    Matrix<double> Mj1_l(5, 1, "0.625 -0.75 -0.25 0.25 0.125");
    Matrix<double> Mj1_r; Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 15, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 0, 0, M_SQRT1_2);
    Matrix<double> Mj1T_l(2, 1, "2 -1");
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_lr(3, "-1 2 -1");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 15, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 1, 1, M_SQRT1_2);
    CDF_factor = -0.5;
  }

  template <>
  SplineBasisData<2,2,P_construction,1,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,P_construction,1,0,0,0>::j0;
    Matrix<double> Mj0_l(0); // empty corner block
    Matrix<double> Mj0_r(2, 1, "0.5 1");
    Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 0, 1, M_SQRT1_2);
    Matrix<double> Mj0T_l(6, 2, "1.5 -0.5 1 0 0.5 0.5 -0.25 1.5 0 0.5 0 -0.25");
    Matrix<double> Mj0T_r(5, 2, "-0.25 0 0.5 0 1.625 -0.75 0.25 1.5 -0.125 1.25");
    Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 2, M_SQRT1_2);
    Matrix<double> Mj1_l(5, 1, "0.625 -0.75 -0.25 0.25 0.125");
    Matrix<double> Mj1_r(4, 1, "-0.0625 -0.125 0.5625 -0.75");
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 0, 1, M_SQRT1_2);
    Matrix<double> Mj1T_l(2, 1, "2 -1");
    Matrix<double> Mj1T_r(0); // empty corner block
    Vector<double> Mj1T_band_lr(3, "-1 2 -1");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 1, 0, M_SQRT1_2);
    CDF_factor = -0.5;
  }

  template <>
  SplineBasisData<2,2,P_construction,0,1,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<2,2,P_construction,0,1,0,0>::j0;
    Matrix<double> Mj0_l(2, 1, "1 0.5");
    Matrix<double> Mj0_r(0); // empty corner block
    Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 0, M_SQRT1_2);
    Matrix<double> Mj0T_l(5, 2, "1.25 -0.125 1.5 0.25 -0.75 1.625 0 0.5 0 -0.25");
    Matrix<double> Mj0T_r(6, 2, "-0.25 0 0.5 0 1.5 -0.25 0.5 0.5 0 1 -0.5 1.5");
    Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 3, M_SQRT1_2);
    Matrix<double> Mj1_l(4, 1, "-0.75 0.5625 -0.125 -0.0625");
    Matrix<double> Mj1_r(5, 1, "0.125 0.25 -0.25 -0.75 0.625");
    Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
    Mj1_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 1, 0, M_SQRT1_2);
    Matrix<double> Mj1T_l(0); // empty corner block
    Matrix<double> Mj1T_r(2, 1, "-1 2");
    Vector<double> Mj1T_band_lr(3, "-1 2 -1");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 8, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 0, 1, M_SQRT1_2);
    CDF_factor = -0.5;
  }

  template <>
  SplineBasisData<3,3,P_construction,0,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,P_construction,0,0,0,0>::j0;
    Matrix<double> Mj0_l(4, 2, "1 0 0.5 0.5 0 0.75 0 0.25");
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(4, "0.25 0.75 0.75 0.25");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 18, 10, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 2, 2, M_SQRT1_2);
    std::ostringstream entries;
    entries << "13.5 -3.28125 0.828125 -0.140625 "
	    << "9 6.5625 -1.65625 0.28125 "
	    << "-6.75 19.96875 -2.484375 0.421875 "
	    << "2.25 -1.03125 10.765625 -1.828125 "
	    << "0 -6.328125 14.9765625 -2.6015625 "
	    << "0 2.109375 -2.7421875 12.8671875 "
	    << "0 0 -2.53125 12.65625 "
	    << "0 0 0.84375 -1.96875 "
	    << "0 0 0 -2.53125 "
	    << "0 0 0 0.84375";
    Matrix<double> Mj0T_l(10, 4, entries.str().c_str()); Mj0T_l.scale(1./9.);
    Matrix<double> Mj0T_r(10, 4); Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(8, "0.09375 -0.28125 -0.21875 1.40625 1.40625 -0.21875 -0.28125 0.09375");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 18, 10, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 4, 4, M_SQRT1_2);
    entries.str("");
    entries << "1.125 0 "
	    << "-1.27734375 -0.46875 "
	    << "0.5634765625 -0.1171875 "
	    << "0.0498046875 1.5234375 "
	    << "-0.146484375 -1.171875 "
	    << "-0.025390625 -0.203125 "
	    << "0.0263671875 0.2109375 "
	    << "0.0087890625 0.0703125";
    Matrix<double> Mj1_l(8, 2, entries.str().c_str()); Mj1_l.scale(1./3.);
    entries.str("");
    entries << "0.0703125 0.0087890625 "
	    << "0.2109375 0.0263671875 "
	    << "-0.203125 -0.025390625 "
	    << "-1.171875 -0.146484375 "
	    << "1.5234375 0.0498046875 "
	    << "-0.1171875 0.5634765625 "
	    << "-0.46875 -1.27734375 "
	    << "0 1.125";
    Matrix<double> Mj1_r(8, 2, entries.str().c_str()); Mj1_r.scale(1./3.);
    Vector<double> Mj1_band_l(8, "-0.09375 -0.28125 0.21875 1.40625 -1.40625 -0.21875 0.28125 0.09375"); Mj1_band_l.scale(1./3.);
    Vector<double> Mj1_band_r(8, "0.09375 0.28125 -0.21875 -1.40625 1.40625 0.21875 -0.28125 -0.09375"); Mj1_band_r.scale(1./3.);
    Mj1_ = QuasiStationaryMatrix<double>(j0, 18, 8, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 2, 2, M_SQRT1_2);
    Matrix<double> Mj1T_l(4, 1, "4 -8 6 -2"); Mj1T_l.scale(1./3.);
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_l(4, "-0.75 2.25 -2.25 0.75");
    Vector<double> Mj1T_band_r(4, "0.75 -2.25 2.25 -0.75");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 18, 8, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 2, 2, M_SQRT1_2);
    CDF_factor = 1./3.;
  }

  template <>
  SplineBasisData<3,3,P_construction,1,0,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,P_construction,1,0,0,0>::j0;
    Matrix<double> Mj0_l(3, 1, "0.5 0.75 0.25");
    Matrix<double> Mj0_r(4, 2, "0.25 0 0.75 0 0.5 0.5 0 1");
    Vector<double> Mj0_band_lr(4, "0.25 0.75 0.75 0.25");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 2, M_SQRT1_2);
    std::ostringstream entries;
    entries << "5.4375 -2.09375 0.46875 "
	    << "4.21875 0.328125 -0.140625 "
	    << "0.46875 3.203125 -0.515625 "
	    << "-2.109375 4.9921875 -0.8671875 "
	    << "0.703125 -0.9140625 4.2890625 "
	    << "0 -0.84375 4.21875 "
	    << "0 0.28125 -0.65625 "
	    << "0 0 -0.84375 "
	    << "0 0 0.28125";
    Matrix<double> Mj0T_l(9, 3, entries.str().c_str()); Mj0T_l.scale(1./3.);
    entries.str("");
    entries << "0.84375 0 0 0 "
	    << "-2.53125 0 0 0 "
	    << "-1.96875 0.84375 0 0 "
	    << "12.65625 -2.53125 0 0 "
	    << "12.8671875 -2.7421875 2.109375 0 "
	    << "-2.6015625 14.9765625 -6.328125 0 "
	    << "-1.828125 10.765625 -1.03125 2.25 "
	    << "0.421875 -2.484375 19.96875 -6.75 "
	    << "0.28125 -1.65625 6.5625 9 "
	    << "-0.140625 0.828125 -3.28125 13.5";
    Matrix<double> Mj0T_r(10, 4, entries.str().c_str()); Mj0T_r.scale(1./9.);
    Vector<double> Mj0T_band_lr(8, "0.09375 -0.28125 -0.21875 1.40625 1.40625 -0.21875 -0.28125 0.09375");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 4, M_SQRT1_2);
    entries.str("");
    entries << "-1.23046875 -0.46875 "
	    << "1.3330078125 -0.1171875 "
	    << "-0.0791015625 1.5234375 "
	    << "-0.544921875 -1.171875 "
	    << "-0.064453125 -0.203125 "
	    << "0.1318359375 0.2109375 "
	    << "0.0439453125 0.0703125";
    Matrix<double> Mj1_l(7, 2, entries.str().c_str()); Mj1_l.scale(1./3.);
    entries.str("");
    entries << "-0.6328125 0.0703125 "
	    << "-1.8984375 0.2109375 "
	    << "1.828125 -0.203125 "
	    << "10.546875 -1.171875 "
	    << "-13.7109375 0.3984375 "
	    << "1.0546875 4.5078125 "
	    << "4.21875 -10.21875 "
	    << "0 9";
    Matrix<double> Mj1_r(8, 2, entries.str().c_str()); Mj1_r.scale(1./27.);
    Vector<double> Mj1_band_lr(8, "-0.09375 -0.28125 0.21875 1.40625 -1.40625 -0.21875 0.28125 0.09375"); Mj1_band_lr.scale(1./3.);
    Mj1_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 1, 2, M_SQRT1_2);
    Matrix<double> Mj1T_l(3, 1, "-8 6 -2"); Mj1T_l.scale(1./3.);
    Matrix<double> Mj1T_r(4, 1, "-2.25 6.75 -9 4.5"); Mj1T_r.scale(1./3.);
    Vector<double> Mj1T_band_lr(4, "-0.75 2.25 -2.25 0.75");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 1, 2, M_SQRT1_2);
    CDF_factor = 1./3.;
  }

  template <>
  SplineBasisData<3,3,P_construction,0,1,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,P_construction,0,1,0,0>::j0;
    Matrix<double> Mj0_l(4, 2, "1 0 0.5 0.5 0 0.75 0 0.25");
    Matrix<double> Mj0_r(3, 1, "0.25 0.75 0.5");
    Vector<double> Mj0_band_lr(4, "0.25 0.75 0.75 0.25");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 2, 1, M_SQRT1_2);
    std::ostringstream entries;
    entries << "13.5 -3.28125 0.828125 -0.140625 "
	    << "9 6.5625 -1.65625 0.28125 "
	    << "-6.75 19.96875 -2.484375 0.421875 "
	    << "2.25 -1.03125 10.765625 -1.828125 "
	    << "0 -6.328125 14.9765625 -2.6015625 "
	    << "0 2.109375 -2.7421875 12.8671875 "
	    << "0 0 -2.53125 12.65625 "
	    << "0 0 0.84375 -1.96875 "
	    << "0 0 0 -2.53125 "
	    << "0 0 0 0.84375";
    Matrix<double> Mj0T_l(10, 4, entries.str().c_str()); Mj0T_l.scale(1./9.);
    entries.str("");
    entries << "0.28125 0 0 "
	    << "-0.84375 0 0 "
	    << "-0.65625 0.28125 0 "
	    << "4.21875 -0.84375 0 "
	    << "4.2890625 -0.9140625 0.703125 "
	    << "-0.8671875 4.9921875 -2.109375 "
	    << "-0.515625 3.203125 0.46875 "
	    << "-0.140625 0.328125 4.21875 "
	    << "0.46875 -2.09375 5.4375";
    Matrix<double> Mj0T_r(9, 3, entries.str().c_str()); Mj0T_r.scale(1./3.);
    Vector<double> Mj0T_band_lr(8, "0.09375 -0.28125 -0.21875 1.40625 1.40625 -0.21875 -0.28125 0.09375");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 4, 3, M_SQRT1_2);
    entries.str("");
    entries << "10.125 0 "
	    << "-11.49609375 -4.21875 "
	    << "5.0712890625 -1.0546875 "
	    << "0.4482421875 13.7109375 "
	    << "-1.318359375 -10.546875 "
	    << "-0.228515625 -1.828125 "
	    << "0.2373046875 1.8984375 "
	    << "0.0791015625 0.6328125";
    Matrix<double> Mj1_l(8, 2, entries.str().c_str()); Mj1_l.scale(1./27.);
    entries.str("");
    entries << "-0.2109375 0.1171875 "
	    << "-0.6328125 0.3515625 "
	    << "0.609375 -0.171875 "
	    << "3.515625 -1.453125 "
	    << "-4.5703125 -0.2109375 "
	    << "0.3515625 3.5546875 "
	    << "1.40625 -3.28125";
    Matrix<double> Mj1_r(7, 2, entries.str().c_str()); Mj1_r.scale(1./9.);
    Vector<double> Mj1_band_lr(8, "-0.09375 -0.28125 0.21875 1.40625 -1.40625 -0.21875 0.28125 0.09375"); Mj1_band_lr.scale(1./3.);
    Mj1_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 2, 1, M_SQRT1_2);
    Matrix<double> Mj1T_l(4, 1, "4 -8 6 -2"); Mj1T_l.scale(1./3.);
    Matrix<double> Mj1T_r(3, 1, "-0.75 2.25 -3");
    Vector<double> Mj1T_band_lr(4, "-0.75 2.25 -2.25 0.75");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 2, 1, M_SQRT1_2);
    CDF_factor = 1./3.;
  }

  template <>
  SplineBasisData<3,3,P_construction,1,1,0,0>::SplineBasisData()
  {
    const int j0 = SplineBasisData_j0<3,3,P_construction,1,1,0,0>::j0;
    Matrix<double> Mj0_l(3, 1, "0.5 0.75 0.25");
    Matrix<double> Mj0_r; Mj0_l.reflect(Mj0_r);
    Vector<double> Mj0_band_lr(4, "0.25 0.75 0.75 0.25");
    Mj0_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 1, M_SQRT1_2);
    std::ostringstream entries;
    entries << "5.4375 -2.09375 0.46875 "
	    << "4.21875 0.328125 -0.140625 "
	    << "0.46875 3.203125 -0.515625 "
	    << "-2.109375 4.9921875 -0.8671875 "
	    << "0.703125 -0.9140625 4.2890625 "
	    << "0 -0.84375 4.21875 "
	    << "0 0.28125 -0.65625 "
	    << "0 0 -0.84375 "
	    << "0 0 0.28125";
    Matrix<double> Mj0T_l(9, 3, entries.str().c_str()); Mj0T_l.scale(1./3.);
    Matrix<double> Mj0T_r(9, 3); Mj0T_l.reflect(Mj0T_r);
    Vector<double> Mj0T_band_lr(8, "0.09375 -0.28125 -0.21875 1.40625 1.40625 -0.21875 -0.28125 0.09375");
    Mj0T_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 3, M_SQRT1_2);
    entries.str("");
    entries << "-1.23046875 -0.46875 "
	    << "1.3330078125 -0.1171875 "
	    << "-0.0791015625 1.5234375 "
	    << "-0.544921875 -1.171875 "
	    << "-0.064453125 -0.203125 "
	    << "0.1318359375 0.2109375 "
	    << "0.0439453125 0.0703125";
    Matrix<double> Mj1_l(7, 2, entries.str().c_str()); Mj1_l.scale(1./3.);
    Matrix<double> Mj1_r(7, 2); Mj1_l.reflect(Mj1_r);
    Vector<double> Mj1_band_l(8, "-0.09375 -0.28125 0.21875 1.40625 -1.40625 -0.21875 0.28125 0.09375"); Mj1_band_l.scale(1./3.);
    Vector<double> Mj1_band_r(8, "0.09375 0.28125 -0.21875 -1.40625 1.40625 0.21875 -0.28125 -0.09375"); Mj1_band_r.scale(1./3.);
    Mj1_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 1, 1, M_SQRT1_2);
    Matrix<double> Mj1T_l(3, 1, "-8 6 -2"); Mj1T_l.scale(1./3.);
    Matrix<double> Mj1T_r; Mj1T_l.reflect(Mj1T_r);
    Vector<double> Mj1T_band_l(4, "-0.75 2.25 -2.25 0.75");
    Vector<double> Mj1T_band_r(4, "0.75 -2.25 2.25 -0.75");
    Mj1T_ = QuasiStationaryMatrix<double>(j0, 16, 8, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 1, 1, M_SQRT1_2);
    CDF_factor = 1./3.;
  }

}

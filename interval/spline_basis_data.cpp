// implementation for spline_basis_data.h

#include <sstream>

namespace WaveletTL
{
  template <>
  SplineBasisData<2,2,P_construction>::SplineBasisData
  (const char* options,
   const int s0, const int s1, const int sT0, const int sT1)
    : s0_(s0), s1_(s1), sT0_(sT0), sT1_(sT1),
      Mj0_(0), Mj1_(0), Mj0T_(0), Mj1T_(0)
  {
    assert(sT0==0 && sT1==0);

    if (s0==0 && s1==0) { // no Dirichlet b.c.'s
      j0_ = 2;
      Matrix<double> Mj0_l(2, 1, "1 0.5");
      Matrix<double> Mj0_r; Mj0_l.mirror(Mj0_r);
      Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
      Mj0_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 1, M_SQRT1_2);
      Matrix<double> Mj0T_l(5, 2, "1.25 -0.125 1.5 0.25 -0.75 1.625 0 0.5 0 -0.25");
      Matrix<double> Mj0T_r; Mj0T_l.mirror(Mj0T_r);
      Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
      Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, M_SQRT1_2);
      Matrix<double> Mj1_l(4, 1, "-0.75 0.5625 -0.125 -0.0625");
      Matrix<double> Mj1_r; Mj1_l.mirror(Mj1_r);
      Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
      Mj1_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 1, 1, M_SQRT1_2);
      Matrix<double> Mj1T_lr(0); // empty corner blocks
      Vector<double> Mj1T_band_lr(3, "-1 2 -1");
      Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1T_lr, Mj1T_lr, Mj1T_band_lr, Mj1T_band_lr, 0, 0, M_SQRT1_2);
    }
      
    if (s0==1 && s1==0) { // homogeneous Dirichlet b.c.'s at x=0
      j0_ = 2;
      Matrix<double> Mj0_l(0); // empty corner block
      Matrix<double> Mj0_r(2, 1, "0.5 1");
      Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
      Mj0_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 0, 1, M_SQRT1_2);
      Matrix<double> Mj0T_l(6, 2, "1.5 -0.5 1 0 0.5 0.5 -0.25 1.5 0 0.5 0 -0.25");
      Matrix<double> Mj0T_r(5, 2, "-0.25 0 0.5 0 1.625 -0.75 0.25 1.5 -0.125 1.25");
      Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
      Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 2, M_SQRT1_2);
      Matrix<double> Mj1_l(5, 1, "0.625 -0.75 -0.25 0.25 0.125");
      Matrix<double> Mj1_r(4, 1, "-0.0625 -0.125 0.5625 -0.75");
      Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
      Mj1_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 0, 1, M_SQRT1_2);
      Matrix<double> Mj1T_l(2, 1, "2 -1");
      Matrix<double> Mj1T_r(0); // empty corner block
      Vector<double> Mj1T_band_lr(3, "-1 2 -1");
      Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 1, 0, M_SQRT1_2);
    }

    if (s0==0 && s1==1) { // homogeneous Dirichlet b.c.'s at x=1
      j0_ = 2;
      Matrix<double> Mj0_l(2, 1, "1 0.5");
      Matrix<double> Mj0_r(0); // empty corner block
      Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
      Mj0_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 0, M_SQRT1_2);
      Matrix<double> Mj0T_l(5, 2, "1.25 -0.125 1.5 0.25 -0.75 1.625 0 0.5 0 -0.25");
      Matrix<double> Mj0T_r(6, 2, "-0.25 0 0.5 0 1.5 -0.25 0.5 0.5 0 1 -0.5 1.5");
      Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
      Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 3, M_SQRT1_2);
      Matrix<double> Mj1_l(4, 1, "-0.75 0.5625 -0.125 -0.0625");
      Matrix<double> Mj1_r(5, 1, "0.125 0.25 -0.25 -0.75 0.625");
      Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
      Mj1_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 1, 0, M_SQRT1_2);
      Matrix<double> Mj1T_l(0); // empty corner block
      Matrix<double> Mj1T_r(2, 1, "-1 2");
      Vector<double> Mj1T_band_lr(3, "-1 2 -1");
      Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 8, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 0, 1, M_SQRT1_2);
    }

    if (s0==1 && s1==1) { // homogeneous Dirichlet b.c.'s at x=0 and x=1
      j0_ = 3;
      Matrix<double> Mj0_lr(0); // empty corner blocks
      Vector<double> Mj0_band_lr(3, "0.5 1.0 0.5");
      Mj0_ = new QuasiStationaryMatrix<double>(j0_, 15, 7, Mj0_lr, Mj0_lr, Mj0_band_lr, Mj0_band_lr, 0, 0, M_SQRT1_2);
      Matrix<double> Mj0T_l(6, 2, "1.5 -0.5 1 0 0.5 0.5 -0.25 1.5 0 0.5 0 -0.25");
      Matrix<double> Mj0T_r; Mj0T_l.mirror(Mj0T_r);
      Vector<double> Mj0T_band_lr(5, "-0.25 0.5 1.5 0.5 -0.25");
      Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 15, 7, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 3, M_SQRT1_2);
      Matrix<double> Mj1_l(5, 1, "0.625 -0.75 -0.25 0.25 0.125");
      Matrix<double> Mj1_r; Mj1_l.mirror(Mj1_r);
      Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
      Mj1_ = new QuasiStationaryMatrix<double>(j0_, 15, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 0, 0, M_SQRT1_2);
      Matrix<double> Mj1T_l(2, 1, "2 -1");
      Matrix<double> Mj1T_r; Mj1T_l.mirror(Mj1T_r);
      Vector<double> Mj1T_band_lr(3, "-1 2 -1");
      Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 15, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 1, 1, M_SQRT1_2);
    }
  }

  template<>
  SplineBasisData<3,3,P_construction>::SplineBasisData
  (const char* options,
   const int s0, const int s1, const int sT0, const int sT1)
    : s0_(s0), s1_(s1), sT0_(sT0), sT1_(sT1),
      Mj0_(0), Mj1_(0), Mj0T_(0), Mj1T_(0)
  {
    assert(sT0 == 0 && sT1 == 0);
      
    if (s0==1 && s1==1) { // homogeneous Dirichlet b.c.'s at x=0 and x=1
      j0_ = 3;
      Matrix<double> Mj0_l(3, 1, "0.5 0.75 0.25");
      Matrix<double> Mj0_r; Mj0_l.mirror(Mj0_r);
      Vector<double> Mj0_band_lr(4, "0.25 0.75 0.75 0.25");
      Mj0_ = new QuasiStationaryMatrix<double>(j0_, 16, 8, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 1, 1, M_SQRT1_2);
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
      Matrix<double> Mj0T_r(9, 3); Mj0T_l.mirror(Mj0T_r);
      Vector<double> Mj0T_band_lr(8, "0.09375 -0.28125 -0.21875 1.40625 1.40625 -0.21875 -0.28125 0.09375");
      Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 16, 8, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 3, 3, M_SQRT1_2);
      entries.str("");
      entries << "-1.23046875 -0.46875 "
	      << "1.3330078125 -0.1171875 "
	      << "-0.0791015625 1.5234375 "
	      << "-0.544921875 -1.171875 "
	      << "-0.064453125 -0.203125 "
	      << "0.1318359375 0.2109375 "
	      << "0.0439453125 0.0703125";
      Matrix<double> Mj1_l(7, 2, entries.str().c_str()); Mj1_l.scale(1./3.);
      Matrix<double> Mj1_r(7, 2); Mj1_l.mirror(Mj1_r);
      Vector<double> Mj1_band_l(8, "-0.09375 -0.28125 0.21875 1.40625 -1.40625 -0.21875 0.28125 0.09375"); Mj1_band_l.scale(1./3.);
      Vector<double> Mj1_band_r(8, "0.09375 0.28125 -0.21875 -1.40625 1.40625 0.21875 -0.28125 -0.09375"); Mj1_band_r.scale(1./3.);
      Mj1_ = new QuasiStationaryMatrix<double>(j0_, 16, 8, Mj1_l, Mj1_r, Mj1_band_l, Mj1_band_r, 1, 1, M_SQRT1_2);
      Matrix<double> Mj1T_l(3, 1, "-8 6 -2"); Mj1T_l.scale(1./3.);
      Matrix<double> Mj1T_r; Mj1T_l.mirror(Mj1T_r);
      Vector<double> Mj1T_band_l(4, "-0.75 2.25 -2.25 0.75");
      Vector<double> Mj1T_band_r(4, "0.75 -2.25 2.25 -0.75");
      Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 16, 8, Mj1T_l, Mj1T_r, Mj1T_band_l, Mj1T_band_r, 1, 1, M_SQRT1_2);
    }
  }

  template <>
  SplineBasisData<2,2,DS_construction>::SplineBasisData
  (const char* options,
   const int s0, const int s1, const int sT0, const int sT1)
    : s0_(s0), s1_(s1), sT0_(sT0), sT1_(sT1),
      Mj0_(0), Mj1_(0), Mj0T_(0), Mj1T_(0)
  {
    if (s0==0 && s1==0 && sT0==0 && sT1==0) {
      if (options == "bio5") {
// 	j0_ = 3;
	j0_ = 2;
	std::ostringstream entries;
	entries << "0.707106781186547 0 "
		<< "0.288735268984507 0.353553390593274 "
		<< "-0.117851130197758 1.01015254455221 "
		<< "-0.0589255650988791 0.505076272276106";
	Matrix<double> Mj0_l(4, 2, entries.str().c_str());
	Matrix<double> Mj0_r; Mj0_l.mirror(Mj0_r);
	Vector<double> Mj0_band_lr(3, "0.35355339059327 0.70710678118655 0.35355339059327");
// 	Mj0_ = new QuasiStationaryMatrix<double>(j0_, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
	Mj0_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
	entries.str("");
	entries << "0.707106781186547 0 "
		<< "1.51522881682832 0.353553390593274 "
		<< "-0.53033008588991 0.742462120245874 "
		<< "0 0.247487373415291 "
		<< "0 -0.123743686707646";
	Matrix<double> Mj0T_l(5, 2, entries.str().c_str());
	Matrix<double> Mj0T_r; Mj0T_l.mirror(Mj0T_r);
	Vector<double> Mj0T_band_lr(5, "-0.17677669529664 0.35355339059327 1.06066017177982 0.35355339059327 -0.17677669529664");
// 	Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
	Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
	entries.str("");
	entries << "0.375 0 "
		<< "-0.109375 -0.0875 "
		<< "0.1875 -0.25 "
		<< "-0.34375 0.75 "
		<< "0.125 -0.25 "
		<< "0.0625 -0.125";
 	Matrix<double> Mj1_l(6, 2, entries.str().c_str());
	Matrix<double> Mj1_r; Mj1_l.mirror(Mj1_r);
	Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
// 	Mj1_ = new QuasiStationaryMatrix<double>(j0_, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
	Mj1_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
	entries.str("");
	entries << "1.33333333333333 0.666666666666667 "
		<< "-2.85714285714286 -1.42857142857143 "
		<< "1 0 "
		<< "0 1 "
		<< "0 -0.5";
	Matrix<double> Mj1T_l(5, 2, entries.str().c_str());
	Matrix<double> Mj1T_r; Mj1T_l.mirror(Mj1T_r);
	Vector<double> Mj1T_band_lr(3, "-0.5 1 -0.5");
// 	Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
	Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
      }

      if (options == "bio5-energy") {
// 	j0_ = 3;
	j0_ = 2;
	std::ostringstream entries;
	entries << "0.707106781186547 0 "
		<< "0.288735268984507 0.353553390593274 "
		<< "-0.117851130197758 1.01015254455221 "
		<< "-0.0589255650988791 0.505076272276106";
	Matrix<double> Mj0_l(4, 2, entries.str().c_str());
	Matrix<double> Mj0_r; Mj0_l.mirror(Mj0_r);
	Vector<double> Mj0_band_lr(3, "0.35355339059327 0.70710678118655 0.35355339059327");
// 	Mj0_ = new QuasiStationaryMatrix<double>(j0_, 17, 9, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
	Mj0_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0_l, Mj0_r, Mj0_band_lr, Mj0_band_lr, 3, 3, 1.0);
	entries.str("");
	entries << "0.707106781186547 0 "
		<< "1.51522881682832 0.353553390593274 "
		<< "-0.53033008588991 0.742462120245874 "
		<< "0 0.247487373415291 "
		<< "0 -0.123743686707646";
	Matrix<double> Mj0T_l(5, 2, entries.str().c_str());
	Matrix<double> Mj0T_r; Mj0T_l.mirror(Mj0T_r);
	Vector<double> Mj0T_band_lr(5, "-0.17677669529664 0.35355339059327 1.06066017177982 0.35355339059327 -0.17677669529664");
// 	Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 17, 9, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
	Mj0T_ = new QuasiStationaryMatrix<double>(j0_, 9, 5, Mj0T_l, Mj0T_r, Mj0T_band_lr, Mj0T_band_lr, 2, 2, 1.0);
	entries.str("");
	entries << "0.375 0 "
		<< "-0.109375 -0.00380794872403089 "
		<< "0.1875 -0.0108798534972311 "
		<< "-0.34375 0.0326395604916934 "
		<< "0.125 -0.0108798534972311 "
		<< "0.0625 -0.00543992674861556";
 	Matrix<double> Mj1_l(6, 2, entries.str().c_str());
	Matrix<double> Mj1_r; Mj1_l.mirror(Mj1_r);
	Vector<double> Mj1_band_lr(5, "-0.125 -0.25 0.75 -0.25 -0.125");
// 	Mj1_ = new QuasiStationaryMatrix<double>(j0_, 17, 8, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
	Mj1_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1_l, Mj1_r, Mj1_band_lr, Mj1_band_lr, 3, 3, 1.0);
	entries.str("");
	entries << "1.33333333333333 15.3188337241014 "
		<< "-2.85714285714286 -32.8260722659316 "
		<< "1 0 "
		<< "0 22.9782505861521 "
		<< "0 -11.4891252930761";
	Matrix<double> Mj1T_l(5, 2, entries.str().c_str());
	Matrix<double> Mj1T_r; Mj1T_l.mirror(Mj1T_r);
	Vector<double> Mj1T_band_lr(3, "-0.5 1 -0.5");
// 	Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 17, 8, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
	Mj1T_ = new QuasiStationaryMatrix<double>(j0_, 9, 4, Mj1T_l, Mj1T_r, Mj1T_band_lr, Mj1T_band_lr, 4, 4, 1.0);
      }
    }
  }


  template <int d,int dT,SplineBasisFlavor flavor>
  SplineBasisData<d,dT,flavor>::~SplineBasisData()
  {
    if (Mj0_) delete Mj0_;
    if (Mj1_) delete Mj1_;
    if (Mj0T_) delete Mj0T_;
    if (Mj1T_) delete Mj1T_;
  }

  template <int d,int dT,SplineBasisFlavor flavor>
  void
  SplineBasisData<d,dT,flavor>::check() const
  {
    cout << endl
	 << "-------------------------------" << endl
	 << "SplineBasisData::check() called" << endl
	 << "-------------------------------" << endl << endl;

    cout << "* some basic data:" << endl;
    cout << "  flavor: " << flavor << endl;
    cout << "  d=" << d << endl;
    cout << "  dT=" << dT << endl;
    cout << "  s0=" << s0_ << endl;
    cout << "  s1=" << s1_ << endl;
    cout << "  sT0=" << sT0_ << endl;
    cout << "  sT1=" << sT1_ << endl;
    cout << "  j0=" << j0_ << endl;
    if (Mj0_) cout << "Mj0=" << endl << *Mj0_;
    if (Mj1_) cout << "Mj1=" << endl << *Mj1_;
    if (Mj0T_) cout << "Mj0T=" << endl << *Mj0T_;
    if (Mj1T_) cout << "Mj1T=" << endl << *Mj1T_;

    if (Mj0_ && Mj0T_) {
      cout << "check biorthogonality of the generators on different levels:" << endl;
      for (int j = j0_; j <= 8; j++) {
	cout << "* j=" << j;
	Mj0_->set_level(j);
	Mj0T_->set_level(j);
	Vector<double> x(Mj0T_->column_dimension()),
	  y(Mj0T_->row_dimension());
	double maxerr = 0;
	for (unsigned int i = 0; i < Mj0T_->column_dimension(); i++) {
	  x.scale(0);
	  x[i] = 1.0;
	  Mj0T_->apply(x,y);
	  Mj0_->apply_transposed(y,x);
	  x[i] = x[i] - 1.0;
	  maxerr = std::max(maxerr, linfty_norm(x));
	}
	cout << ", ||Mj0^T*Mj0T-I||_1: " << maxerr;
	
	maxerr = 0;
	for (unsigned int i = 0; i < Mj0_->column_dimension(); i++) {
	  x.scale(0);
	  x[i] = 1.0;
	  Mj0_->apply(x,y);
	  Mj0T_->apply_transposed(y,x);
	  x[i] = x[i] - 1.0;
	  maxerr = std::max(maxerr, linfty_norm(x));
	}
	cout << ", ||Mj0T^T*Mj0-I||_1: " << maxerr << endl;
      }
    }

    if (Mj0_ && Mj0T_ && Mj1_ && Mj1T_) {
      cout << endl << "check biorthogonality of the full systems on different levels:" << endl;
      for (int j = j0_; j <= 8; j++) {
	cout << "* j=" << j;
	Mj0_->set_level(j);
	Mj1_->set_level(j);
	Mj0T_->set_level(j);
	Mj1T_->set_level(j);
	Vector<double> x(Mj0T_->row_dimension()),
	  y1(Mj0T_->column_dimension()),
	  y2(Mj1T_->column_dimension()),
	  xhelp(Mj0T_->row_dimension());
	
	double maxerr = 0;
	for (unsigned int i = 0; i < Mj0T_->row_dimension(); i++) {
	  x.scale(0);
	  x[i] = 1;
	  Mj0T_->apply_transposed(x,y1);
	  Mj1T_->apply_transposed(x,y2);
	  Mj0_->apply(y1,x);
	  Mj1_->apply(y2,xhelp);
	  x.add(xhelp);
	  x[i] = x[i] - 1.0;
	  maxerr = std::max(maxerr, linfty_norm(x));
	}
	cout << ", ||Mj*Gj-I||_1: " << maxerr;
	
	maxerr = 0;
	for (unsigned int i = 0; i < Mj0_->row_dimension(); i++) {
	  x.scale(0);
	  x[i] = 1;
	  Mj0_->apply_transposed(x,y1);
	  Mj1_->apply_transposed(x,y2);
	  Mj0T_->apply(y1,x);
	  Mj1T_->apply(y2,xhelp);
	  x.add(xhelp);
	  x[i] = x[i] - 1.0;
	  maxerr = std::max(maxerr, linfty_norm(x));
	}
	cout << ", ||Gj*Mj-I||_1: " << maxerr << endl;
      }
    }

 
  }
}

// implementation for w_method.h

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

namespace MathTL
{
  template <class VECTOR>
  WMethodStageEquationHelper<VECTOR>::~WMethodStageEquationHelper()
  {
  }

  template <class VECTOR>
  WMethodPreprocessRHSHelper<VECTOR>::~WMethodPreprocessRHSHelper()
  {
  }

  template <class VECTOR>
  WMethod<VECTOR>::WMethod(const Method method,
			   const WMethodStageEquationHelper<VECTOR>* s)
    : stage_equation_helper(s), preprocessor(0)
  {
    LowerTriangularMatrix<double> Alpha, Gamma;
    Vector<double> b, bhat;
    double gamma;

    switch(method)
      {
      case ROS2:
	Alpha.resize(2,2);
	Alpha.set_entry(1, 0, 1.0);

	Gamma.resize(2,2);
	gamma = 1.0 + M_SQRT1_2;
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -2*gamma);
	Gamma.set_entry(1, 1, gamma);
	
	b.resize(2);
	b[0] = b[1] = 0.5;

	bhat.resize(2);
	bhat[0] = 1.0;

	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);
	
 	p = 2;
 	break;
      case RODAS3:
	Alpha.resize(4,4);
	Alpha.set_entry(2, 0,  1.);
	Alpha.set_entry(3, 0,  3./4.);
	Alpha.set_entry(3, 1, -1./4.);
	Alpha.set_entry(3, 2,  1./2.);

	Gamma.resize(4,4);
	gamma = 0.5;
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0,  1.);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0, -1./4.);
	Gamma.set_entry(2, 1, -1./4.);
	Gamma.set_entry(2, 2, gamma);
	Gamma.set_entry(3, 0,  1./12.);
	Gamma.set_entry(3, 1,  1./12.);
	Gamma.set_entry(3, 2, -2./3.);
	Gamma.set_entry(3, 3, gamma);

	b.resize(4);
	b[0] = 5./6.;
	b[1] = b[2] = -1./6.;
	b[3] = 1./2.;
	
	bhat.resize(4);
	bhat[0] =  3./4.;
	bhat[1] = -1./4.;
	bhat[2] = 1./2.;

	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);
	
 	p = 3;
 	break;
      case ROS3P:
	Alpha.resize(3,3);
	Alpha.set_entry(1, 0, 1.);
	Alpha.set_entry(2, 0, 1.);

	Gamma.resize(3,3);
	gamma = 0.5 + sqrt(3.)/6.;
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -1.);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0, -gamma);
	Gamma.set_entry(2, 1, 0.5-2*gamma);
	Gamma.set_entry(2, 2, gamma);

	b.resize(3);
	b[0] = 2./3.;
	b[2] = 1./3.;

	bhat.resize(3);
	bhat[0] = bhat[1] = bhat[2] = 1./3.;

	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);

	p = 3;
	break;
      case ROWDA3:
	Alpha.resize(3,3);
	Alpha.set_entry(1, 0, 0.7);
	Alpha.set_entry(2, 0, 0.7);
	
	Gamma.resize(3,3);
	gamma = 0.43586652150845899941601945119356; // gamma^3-3*gamma^2+1.5*gamma-1/6=0
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, 0.1685887625570998);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0, 4.943922277836421);
	Gamma.set_entry(2, 1, 1.);
	Gamma.set_entry(2, 2, gamma);

	b.resize(3);
	b[0] =  0.3197278911564624;
	b[1] =  0.7714777906171382;
	b[2] = -0.09120568177360061;
	
	bhat.resize(3);
	bhat[0] = 0.926163587124091;
	bhat[1] = 0.073836412875909;

	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);
	
 	p = 3;
 	break;
      case ROS3:
	gamma = 0.43586652150845899941601945119356; // gamma^3-3*gamma^2+1.5*gamma-1/6=0

	Alpha.resize(3,3);
	Alpha.set_entry(1, 0, gamma);
	Alpha.set_entry(2, 0, gamma);
	
	Gamma.resize(3,3);
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -0.19294655696029095575009695436041);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 1, 1.74927148125794685173529749738960);
	Gamma.set_entry(2, 2, gamma);
	
	b.resize(3);
	b[0] = -0.75457412385404315829818998646589;
	b[1] =  1.94100407061964420292840123379419;
	b[2] = -0.18642994676560104463021124732829;
	
	bhat.resize(3);
	bhat[0] = -1.53358745784149585370766523913002;
	bhat[1] =  2.81745131148625772213931745457622;
	bhat[2] = -0.28386385364476186843165221544619;
	
	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);
	
 	p = 3;
 	break;
      case ROS3Pw:
	gamma = 0.5 + sqrt(3.)/6.;

	Alpha.resize(3,3);
	Alpha.set_entry(1, 0, 2*gamma);
	Alpha.set_entry(2, 0, 0.5);

	Gamma.resize(3,3);
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -2*gamma);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0, -.67075317547305480);
	Gamma.set_entry(2, 1, -.17075317547305482);
	Gamma.set_entry(2, 2, gamma);

	b.resize(3);
	b[0] = 0.10566243270259355;
	b[1] = 0.049038105676657971;
	b[2] = 0.84529946162074843;

	bhat.resize(3);
	bhat[0] = -0.17863279495408180;
	bhat[1] = 1./3.;
	bhat[2] = 0.84529946162074843;

	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);

	p = 3;
	break;
      case ROSI2P2:
	gamma = 0.43586652150845899941601945119356; // gamma^3-3*gamma^2+1.5*gamma-1/6=0

	Alpha.resize(4,4);
	Alpha.set_entry(1, 0, 0.5);
	Alpha.set_entry(2, 0, -0.51983699657507165);
	Alpha.set_entry(2, 1,  1.5198369965750715);
	Alpha.set_entry(3, 0, -0.51983699657507165);
	Alpha.set_entry(3, 1,  1.5198369965750715);

	Gamma.resize(4,4);
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -0.5);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0, -0.40164172503011392);
	Gamma.set_entry(2, 1,  1.174271852697665);
	Gamma.set_entry(2, 2, gamma);
	Gamma.set_entry(3, 0,  1.1865036632417383);
	Gamma.set_entry(3, 1, -1.5198369965750715);
	Gamma.set_entry(3, 2, -0.10253318817512568);
	Gamma.set_entry(3, 3, gamma);

	b.resize(4);
	b[0] = 2./3.;
	b[2] = -0.10253318817512568;
	b[3] =  0.435866521508459;

	bhat.resize(4);
	bhat[0] = -0.95742384859111473;
	bhat[1] =  2.9148476971822297;
	bhat[2] =  0.5;
	bhat[3] = -1.4574238485911146;
	
	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);
	
	p = 3;
	break;
      case ROSI2PW:
	gamma = 0.43586652150845899941601945119356; // gamma^3-3*gamma^2+1.5*gamma-1/6=0

	Alpha.resize(4,4);
	Alpha.set_entry(1, 0,  8.7173304301691801e-1);
	Alpha.set_entry(2, 0, -7.9937335839852708e-1);
	Alpha.set_entry(2, 1, -7.9937335839852708e-1);
	Alpha.set_entry(3, 0,  7.0849664917601007e-1);
	Alpha.set_entry(3, 1,  3.1746327955312481e-1);
	Alpha.set_entry(3, 2, -2.5959928729134892e-2);

	Gamma.resize(4,4);
	Gamma.set_entry(0, 0, gamma);
	Gamma.set_entry(1, 0, -8.7173304301691801e-1);
	Gamma.set_entry(1, 1, gamma);
	Gamma.set_entry(2, 0,  3.0647867418622479);
	Gamma.set_entry(2, 1,  3.0647867418622479);
	Gamma.set_entry(2, 2, gamma);
	Gamma.set_entry(3, 0, -1.0424832458800504e-1);
	Gamma.set_entry(3, 1, -3.1746327955312481e-1);
	Gamma.set_entry(3, 2, -1.4154917367329144e-2);
	Gamma.set_entry(3, 3, gamma);

	b.resize(4);
	b[0] =  6.0424832458800504e-1;
	b[2] = -4.0114846096464034e-2;
	b[3] =  4.3586652150845900e-1;

	bhat.resize(4);
	bhat[0] = bhat[1] = 4.4315753191688778e-1;
	bhat[3] = 1.1368493616622447e-1;
	
	transform_coefficients(Alpha, Gamma, b, bhat,
			       A, C, m, e, alpha_vector, gamma_vector);

	p = 3;
	break;
      case GRK4T:
	A.resize(4,4); // corresponds to the [HW] notation
	A(1,0) = 0.2000000000000000e+01;
        A(2,0) = A(3,0) = 0.4524708207373116e+01;
        A(2,1) = A(3,1) = 0.4163528788597648e+01;
	
	C.resize(4,4); // corresponds to the [HW] notation
	C(0,0) = C(1,1) = C(2,2) = C(3,3) = 0.231; // gamma
	C(1,0) = -0.5071675338776316e+01;
	C(2,0) = 0.6020152728650786e+01;
 	C(2,1) = 0.1597506846727117e+00;
 	C(3,0) = -0.1856343618686113e+01;
 	C(3,1) = -0.8505380858179826e+01;
 	C(3,2) = -0.2084075136023187e+01;
	
	gamma_vector.resize(4); // corresponds to Di in [HW]
	gamma_vector[0] = 0.2310000000000000e+00;
        gamma_vector[1] = -0.3962966775244303e-01;
	gamma_vector[2] = 0.5507789395789127e+00;
        gamma_vector[3] = -0.5535098457052764e-01;

	alpha_vector.resize(4); // corresponds to Ci in [HW]
	alpha_vector[1] = 0.4620000000000000e+00;
	alpha_vector[2] = alpha_vector[3]
	  = 0.8802083333333334e+00;

 	m.resize(4); // corresponds to Bi in [HW] (the m_i in the [HW II] book)
 	m[0] = 0.3957503746640777e+01;
	m[1] = 0.4624892388363313e+01;
	m[2] = 0.6174772638750108e+00;
	m[3] = 0.1282612945269037e+01;

 	e.resize(4);
	e[0] = 0.2302155402932996e+01;
	e[1] = 0.3073634485392623e+01;
	e[2] = -0.8732808018045032e+00;
	e[3] = -0.1282612945269037e+01;

	p = 4;
	break;
      case RODAS:
	A.resize(6,6);

	// values from KARDOS:
	A(1,0) =   1.5440000000000000;
	A(2,0) =   0.9466785281022800;
	A(2,1) =   0.2557011699000000;
	A(3,0) =   3.3148251870787937;
	A(3,1) =   2.8961240159798773;
	A(3,2) =   0.9986419140000000;
	A(4,0) =   1.2212245090707980;
	A(4,1) =   6.0191344810926299;
	A(4,2) =  12.5370833291149566;
	A(4,3) =  -0.6878860361200000;
	A(5,0) =   1.2212245092986209;
	A(5,1) =   6.0191344813485754;
	A(5,2) =  12.5370833293196604;
	A(5,3) =  -0.6878860360800001;
	A(5,4) =   1.0000000000000000;
	
	C.resize(6, 6);
 	C(0,0) = C(1,1) = C(2,2) = C(3,3)
	  = C(4,4) = C(5,5) = 0.25; // = gamma

	C(1,0) =  -5.6688000000000000;
	C(2,0) =  -2.4300933568670464;
	C(2,1) =  -0.2063599157120000;
	C(3,0) =  -0.1073529065055983;
	C(3,1) =  -9.5945622510667228;
	C(3,2) = -20.4702861487999996;
	C(4,0) =   7.4964433159050206;
	C(4,1) = -10.2468043146053738;
	C(4,2) = -33.9999035259299589;
	C(4,3) =  11.7089089319999999;
	C(5,0) =   8.0832467990118602;
	C(5,1) =  -7.9811329880455499;
	C(5,2) = -31.5215943254324245;
	C(5,3) =  16.3193054312706352;
	C(5,4) =  -6.0588182388799998;

	// values from KARDOS:
	gamma_vector.resize(6);
	gamma_vector[0] =  0.2500000000000000;
	gamma_vector[1] = -0.1043000000000000;
	gamma_vector[2] =  0.1034999999980000;
	gamma_vector[3] = -0.0362000000000000;
	gamma_vector[4] =  0.0000000000000000;
	gamma_vector[5] =  0.0000000000000000;

	// values from KARDOS:
	alpha_vector.resize(6);
	alpha_vector[0] = 0.0000000000000000;
	alpha_vector[1] = 0.3860000000000000;
	alpha_vector[2] = 0.2100000000000000;
	alpha_vector[3] = 0.6300000000000000;
	alpha_vector[4] = 1.0000000000000000;
	alpha_vector[5] = 1.0000000000000000;

	// values from KARDOS:
	m.resize(6);
	m[0] =  1.2212245092981915;
	m[1] =  6.0191344813101981;
	m[2] = 12.5370833292377792;
	m[3] = -0.6878860360960002;
	m[4] =  1.0000000000000000;
	m[5] =  1.0000000000000000;

	e.resize(6);
  	e[0] = 0.0;
  	e[1] = 0.0;
  	e[2] = 0.0;
  	e[3] = 0.0;
  	e[4] = 0.0;
  	e[5] = m[5];
	
	p = 4;
	break;
      case RODASP:
	A.resize(6, 6);

	// values from KARDOS:
	A(1,0) =   3.0000000000000000;
	A(2,0) =   1.8310367935359999;
	A(2,1) =   0.4955183967600000;
	A(3,0) =   2.3043765826379414;
	A(3,1) =  -0.0524927524844542;
	A(3,2) =  -1.1767987618400000;
	A(4,0) =  -7.1704549640449367;
	A(4,1) =  -4.7416366720041934;
	A(4,2) = -16.3100263134518535;
	A(4,3) =  -1.0620040441200000;
	A(5,0) =  -7.1704549641649340;
	A(5,1) =  -4.7416366720441925;
	A(5,2) = -16.3100263134518570;
	A(5,3) =  -1.0620040441200000;
	A(5,4) =   1.0000000000000000;

	C.resize(6, 6);
 	C(0,0) = C(1,1) = C(2,2) = C(3,3)
	  = C(4,4) = C(5,5) = 0.25; // = gamma

	// values from KARDOS:
	C(1,0) = -12.0000000000000000;
	C(2,0) =  -8.7917951740800000;
	C(2,1) =  -2.2078655870400000;
	C(3,0) =  10.8179305689176530;
	C(3,1) =   6.7802706116824574;
	C(3,2) =  19.5348594463999987;
	C(4,0) =  34.1909500739412096;
	C(4,1) =  15.4967115394459682;
	C(4,2) =  54.7476087604061235;
	C(4,3) =  14.1600539214399994;
	C(5,0) =  34.6260583162319335;
	C(5,1) =  15.3008497633150125;
	C(5,2) =  56.9995557863878588;
	C(5,3) =  18.4080700977581699;
	C(5,4) =  -5.7142857142399999;

	// values from KARDOS:
	gamma_vector.resize(6);
	gamma_vector[0] =  0.2500000000000000;
	gamma_vector[1] = -0.5000000000000000;
	gamma_vector[2] = -0.0235040000000000;
	gamma_vector[3] = -0.0362000000000000;
	gamma_vector[4] =  0.0               ;
	gamma_vector[5] =  0.0               ;

	// values from KARDOS:
	alpha_vector.resize(6);
	alpha_vector[0] = 0.0;
	alpha_vector[1] = 0.75;
	alpha_vector[2] = 0.21;
	alpha_vector[3] = 0.63;
	alpha_vector[4] = 1.0;
	alpha_vector[5] = 1.0;

	// values from KARDOS:
	m.resize(6);
	m[0] =  -7.1704549641649322;
	m[1] =  -4.7416366720441925;
	m[2] = -16.3100263134518570;
	m[3] =  -1.0620040441200000;
	m[4] =   1.0000000000000000;
	m[5] =   1.0000000000000000;
	
	e.resize(6);
  	e[0] = 0.0;
  	e[1] = 0.0;
  	e[2] = 0.0;
  	e[3] = 0.0;
  	e[4] = 0.0;
  	e[5] = m[5];

	p = 4;
	break;
      default:
	break;
      }
  }
  
  template <class VECTOR>
  void
  WMethod<VECTOR>::transform_coefficients(const LowerTriangularMatrix<double>& Alpha,
					  const LowerTriangularMatrix<double>& Gamma,
					  const Vector<double>& b,
					  const Vector<double>& bhat,
					  LowerTriangularMatrix<double>& A,
					  LowerTriangularMatrix<double>& C,
					  Vector<double>& m,
					  Vector<double>& e,
					  Vector<double>& alpha_vector,
					  Vector<double>& gamma_vector)
  {
    const unsigned int s = Alpha.row_dimension();

    alpha_vector.resize(s, true);
    for (unsigned int i = 0; i < s; i++) {
      double alphai = 0;
      for (unsigned int j = 0; j < i; j++) alphai += Alpha.get_entry(i, j);
      alpha_vector[i] = alphai;
    }

    gamma_vector.resize(s, true);
    for (unsigned int i = 0; i < s; i++) {
      double gammai = 0;
      // note that here, we use "j<=i" (compare order condition check below):
      for (unsigned int j = 0; j <= i; j++) gammai += Gamma.get_entry(i, j);
      gamma_vector[i] = gammai;
    }

    Gamma.inverse(C);
    A.resize(s, s);
    A = Alpha * C; // A = Alpha * Gamma^{-1}

    m.resize(s, false);
    C.apply_transposed(b, m); // m^T = b^T * Gamma^{-1}

    e.resize(s, true);
    C.apply_transposed(bhat, e);
    e.sadd(-1.0, m); // e = m-mhat

    C.scale(-1.0);
    for (unsigned int i = 0; i < s; i++)
      C.set_entry(i, i, Gamma.get_entry(0,0)); // gamma

//     cout << "* coefficients after transform_coefficients():" << endl;
//     cout << "alpha_vector=" << alpha_vector << endl;
//     cout << "gamma_vector=" << gamma_vector << endl;
//     cout << "A=" << endl << A;
//     cout << "C=" << endl << C;
//     cout << "m=" << m << endl;
//     cout << "e=" << e << endl;

//     check_order_conditions(Alpha, Gamma, b, bhat, false);
  }

  template <class VECTOR>
  void
  WMethod<VECTOR>::check_order_conditions(const LowerTriangularMatrix<double>& Alpha,
					  const LowerTriangularMatrix<double>& Gamma,
					  const Vector<double>& b,
					  const Vector<double>& bhat,
					  const bool wmethod)
  {
    cout << "* checking algebraic order conditions for a ROW method:" << endl;

    const unsigned int s = Alpha.row_dimension();
    double help, gamma = Gamma.get_entry(0,0);

    LowerTriangularMatrix<double> Beta(s, s);
    for (unsigned int i = 0; i < s; i++)
      for (unsigned int j = 0; j < i; j++)
	Beta.set_entry(i, j, Alpha.get_entry(i, j) + Gamma.get_entry(i, j));

    Vector<double> alpha_vector(s);
    for (unsigned int i = 0; i < s; i++) {
      double alphai = 0;
      for (unsigned int j = 0; j < i; j++) alphai += Alpha.get_entry(i, j);
      alpha_vector[i] = alphai;
    }
    cout << "  alpha_vector=" << alpha_vector << endl;

    Vector<double> gamma_vector(s);
    for (unsigned int i = 0; i < s; i++) {
      double gammai = 0;
      // for the order conditions, we use gamma_i=sum_{j<i}gamma_j (not "<=" as in transform_coefficients())
      for (unsigned int j = 0; j < i; j++) gammai += Gamma.get_entry(i, j);
      gamma_vector[i] = gammai;
    }
    cout << "  gamma_vector=" << gamma_vector << endl;

    Vector<double> beta_vector(alpha_vector+gamma_vector);
    cout << "  beta_vector=" << beta_vector << endl;

    cout << "  (A1)  (sum_i b_i)-1=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += b[i];
    cout << help-1 << endl;

    cout << "  (A2)  (sum_i b_i beta_i)-(1/2-gamma)=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += b[i] * beta_vector[i];
    cout << help-(0.5-gamma) << endl;

    cout << "  (A3a) (sum_i b_i alpha_i^2)-1/3=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += b[i] * alpha_vector[i] * alpha_vector[i];
    cout << help-1./3. << endl;

    cout << "  (A3b) (sum_{i,j} b_i beta_{i,j} beta_j)-(1/6-gamma+gamma^2)=";
    help = 0;
    for (unsigned int i = 0; i < s; i++)
      for (unsigned int j = 0; j < i; j++)
	help += b[i] * Beta.get_entry(i,j) * beta_vector[j];
    cout << help-1./6.+gamma-gamma*gamma << endl;

    cout << "  (A1)  (sum_i bhat_i)-1=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += bhat[i];
    cout << help-1 << endl;

    cout << "  (A2)  (sum_i bhat_i beta_i)-(0.5-gamma)=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += bhat[i] * beta_vector[i];
    cout << help-(0.5-gamma) << endl;

    cout << "  (A3a) (sum_i bhat_i alpha_i^2)-1/3=";
    help = 0;
    for (unsigned int i = 0; i < s; i++) help += bhat[i] * alpha_vector[i] * alpha_vector[i];
    cout << help-1./3. << endl;

    cout << "  (A3b) (sum_{i,j} bhat_i beta_{i,j} beta_j)-(1/6-gamma+gamma^2)=";
    help = 0;
    for (unsigned int i = 0; i < s; i++)
      for (unsigned int j = 0; j < i; j++)
	help += bhat[i] * Beta.get_entry(i,j) * beta_vector[j];
    cout << help-1./6.+gamma-gamma*gamma << endl;
  }

  template <class VECTOR>
  void
  WMethod<VECTOR>::increment(const AbstractIVP<VECTOR>* ivp,
			     const double t_m, const VECTOR& u_m,
			     const double tau,
			     VECTOR& u_mplus1,
			     VECTOR& error_estimate,
			     const double tolerance) const
  {
    const unsigned int stages = A.row_dimension(); // for readability

    Array1D<VECTOR> u(stages);
    u[0] = u_m; u[0].scale(0.0); // ensures correct size
    for (unsigned int i = 1; i < stages; i++)
      u[i] = u[0];

    VECTOR rhs(u[0]), help(u[0]); // ensures correct size
    
    // solve stage equations (TODO: adjust the tolerances appropriately)
    for (unsigned int i(0); i < stages; i++) {
      // setup i-th right-hand side
      //   f(t_m + \tau * \alpha_i, u^{(m)} + \sum_{j=1}^{i-1} a_{i,j} * u_j)
      //   + \sum_{j=1}^{i-1} \frac{c_{i,j}}{\tau} * u_j
      //   + \tau * \gamma_i * g

      help = u_m;
      for (unsigned int j(0); j < i; j++)
  	help.add(A(i,j), u[j]);
      ivp->evaluate_f(t_m+tau*alpha_vector[i], help, tolerance/(4*stages), rhs);

      if (preprocessor == 0) { // no preprocessing necessary
	for (unsigned int j(0); j < i; j++)
	  rhs.add(C(i,j)/tau, u[j]);
      } else {
	help.scale(0.0);
	for (unsigned int j(0); j < i; j++)
	  help.add(C(i,j)/tau, u[j]);
	preprocessor->preprocess_rhs_share(help, tolerance/(4*stages));
	rhs.add(help);
      }
	
      stage_equation_helper->approximate_ft(ivp, t_m, u_m, tolerance/(4*stages), help);
      rhs.add(tau*gamma_vector[i], help);
      
      // solve i-th stage equation
      // (\tau*\gamma_{i,i})^{-1}I - T) u_i = rhs
      stage_equation_helper->solve_W_stage_equation(ivp, t_m, u_m, 1./(tau*C(i,i)), rhs, tolerance/(4*stages), u[i]);
    }
    
    // update u^{(m)} -> u^{(m+1)} by the k_i
    u_mplus1 = u_m;
    for (unsigned int i(0); i < stages; i++)
      u_mplus1.add(m[i], u[i]);
    
    // error estimate
    error_estimate = u_m; error_estimate.scale(0.0); // ensures correct size
    for (unsigned int i = 0; i < stages; i++)
      error_estimate.add(e[i], u[i]);
  }
}

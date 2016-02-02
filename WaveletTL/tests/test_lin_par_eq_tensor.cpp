// more output for the cached problemm (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0
#define _NO_JPMAX_grt_70_WARNING 1 // turn off warning in APPLY_TENSOR
// switch between isotropic and anisotropic Wavelets (in cdd1.h)
#define _WAVELETTL_USE_TBASIS 1
#define _DIMENSION 2
#define _PROBLEM_NO 5
#define _RANDBEDINGUNGEN false
#define _D 3
#define _DT 3
#define _OFFSET 0
#define _TESTTHEHELPER 0
#define _HELPERSALPHA 1.0
#define _VORBERECHNETE_RHS 0
#define _TESTPARABOLIC_FUNCTIONALITY 0 // test evaluate_f, evaluate_ft, solve_ROW_stage_equation, preprocess_rhs_share
#define _TESTPARABOLIC_RUN 1 // solve the parabolic problem (with ROS2)
#define _EXPANSIONTYPE_F_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_f; 1 liefert Unterschied falls f konstant ist und man TIME_CONSTANT =0,1 wählt
#define _EXPANSIONTYPE_FT_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_ft
#define _TIME_CONSTANT_DRIVING_TERM 0 // may choose 1 for 1D problems 4,6,7, 2D problem 1. should change nothing!
// CDD1 or CDD2 ?
#define _CDD_TYPE 1
// use delta-distribution as RHS
// #define _DELTA 0

#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
//#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <adaptive/cdd1.h>
//#include <adaptive/cdd2.h>
#include <geometry/sampled_mapping.h>

#if _WAVELETTL_USE_TBASIS == 1

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tbasis_indexplot.h>
//#include <cube/tbasis_expansion.h>

#include <parabolic/lin_par_eq_tensor.h>
/*
 * for comparisn with the isotropic case
 */
#else
#include <cube/cube_basis.h>
#include <galerkin/cube_equation.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>
#include <cube/cube_indexplot.h>

#include <cube/cube_expansion.h>

#include <parabolic/lin_par_equation.h>
#endif


// ---- from test_lin_par_equation


#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <list>

#include <time.h>

// for comparison in the 1d case
#include <galerkin/gramian.h>
#include <galerkin/cached_problem.h>

//#define _MATHTL_ONESTEPSCHEME_VERBOSITY 1
//#define _WAVELETTL_CDD1_VERBOSITY 0
// --------------
//- --- - - - - - - -- - -- - ---
#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>
//#include <geometry/grid.h>
//#include <geometry/sampled_mapping.h>

//#include <interval/i_index.h>
//#include <interval/ds_basis.h>
//#include <interval/ds_expansion.h>
//#include <interval/p_basis.h>

//#include <galerkin/sturm_equation.h>




//#include <algebra/infinite_vector.h>
//#include <numerics/ivp.h>
//#include <numerics/w_method.h>
//#include <utils/function.h>
//#include <galerkin/cached_tproblem.h>
//#include <adaptive/apply_tensor.h>
//#include <adaptive/cdd1.h>


using namespace std;
using namespace MathTL;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;


/* -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega */
template <unsigned int DIM>
class PertubedBVP
: public EllipticBVP<DIM>
{
public:
    /*!
      constructor with given right-hand side
    */
    PertubedBVP(const Function<DIM>* f, const double eta)
    : eta_(eta), EllipticBVP<DIM> (f, f, f)
    {
    }

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return eta_; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 1.0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
  protected:
      const double eta_;
};

/* 0 Operator*/
template <unsigned int DIM>
class ZeroBVP
: public EllipticBVP<DIM>
{
public:
    /*!
      constructor with given right-hand side
    */
    ZeroBVP(const Function<DIM>* f)
    : EllipticBVP<DIM> (f, f, f)
    {
    }

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return 0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
};

  /* Some 1D righthand sides for different problems
   * 1: -delta u = f; hom Dirichlet BC; u=x*(1-x); f = 2;
   * 2: (I-delta) u = f; hom Dirichlet BC; u=x*(1-x); f=x*(1-x)+2;
   * 3: (I-delta) u = f; hom Neumann BC; u = x^3-1.5x^2+d; f=x^3-1.5x^2-6x+3+d;
   * 4: u'=delta u + f; hom Dirichlet BC, f time independent; u=(1-exp(-pi^2t))sin(pix); f = pi^2sin(pix)
   * 5: u'=delta u + f; u = tcos(pix); f = (1+pi^2)cos(pix)
   * 6: u'=delta u + f; u = ??, u_0=hat; f=0;
   * 7: u_0 = \Xi_[0,0.5]; f=0;
   * 8: u'=delta u + f; hom Neumann BC; u= exp(alpha*(1-t)) * sin^2(2pix); f= u_t-u_xx=exp(alpha(1-t))*(-alpha*sin^2+8pi^2sin^2-8pi^2cos^2)
   * 9: u'=delta u + f; hom Neumann BC; u_x = 0; u=p(t); f= p'(t); p(t) = exp(2t)
   *10: u'=delta u + f; hom Neumann BC; u_x = 0; u= t; f= 1
   *11: u'=delta u + f; hom Neumann BC; u_x = 0; f= 0, u= const
   */
  template <unsigned int N>
  class Driving_f1D
  : public Function<1,double>

  {
  public:
      virtual ~Driving_f1D() {};
      double value(const Point<1>& p, const unsigned int component = 0) const
      {
          switch(N) {
              case 1:
                  return 2.0;
                  break;
              case 2:
                  return p[0]*(1-p[0])+2.0;
                  break;
              case 3:
                  return 3.0 + p[0]*(-6.0+p[0]*(-1.5+p[0]));
                  break;
              case 4:
                  return M_PI*M_PI*sin(M_PI*p[0]);
                  break;
              case 5:
                  return (1+M_PI*M_PI*get_time())*cos(M_PI*p[0]);
                  break;
              case 6:
              case 7:
                  return 0;
                  break;
              case 8:
                  return exp(1*(1-get_time()))*( (8*M_PI*M_PI-1)*sin(2*M_PI*p[0])*sin(2*M_PI*p[0]) -8*M_PI*M_PI*cos(2*M_PI*p[0])*cos(2*M_PI*p[0]) ) ;
                  break;
              case 9:
                  return 2*exp(2*get_time());
                  break;
              case 10:
                  return 1;
                  break;
              case 11:
                  return 0;
                  break;
              default:
                  abort();
                  break;
          }
      }
      void vector_value(const Point<1>& p, Vector<double>& values) const
      {
          values[0] = value(p);
      }
  };

  /* Some 1D solutions for different problems
   * 1: -delta u = f; hom Dirichlet BC; u=x*(1-x); f = 2;
   * 2: (I-delta) u = f; hom Dirichlet BC; u=x*(1-x); f=x*(1-x)+2;
   * 3: (I-delta) u = f; hom Neumann BC; u = x^3-1.5x^2+d; f=x^3-1.5x^2-6x+3+d;
   * 4: u'=delta u + f; hom Dirichlet BC, f time independent; u=(1-exp(-pi^2t))sin(pix); f = pi^2sin(pix)
   * 5: u'=delta u + f; u = tcos(pix); f = (1+pi^2)cos(pix)
   * 6: u'=delta u + f; u = ??, u_0=hat; f=0;
   * 7: u_0 = \Xi_[0,0.5]; f=0;
   * 8: u'=delta u + f; u= exp(alpha*(1-t)) * sin^2(2pix); f= u_t-u_xx=exp(alpha(1-t))*(-alpha*sin^2+8pi^2sin^2-8pi^2cos^2);
   * 9: u'=delta u + f; u=p(t); f= p'(t);
   */
  template <unsigned int N>
  class Exact_Sol1D
  : public Function<1,double>
  {
  public:
      virtual ~Exact_Sol1D() {};
      double value(const Point<1>& p, const unsigned int component = 0) const {
          switch(N)
          {
              case 1:
                  return p[0]*(1-p[0]);
                  break;
              case 2:
                  return p[0]*(1-p[0]);
                  break;
              case 3:
                  return p[0]*p[0]*(-1.5+p[0]);
                  break;
              case 4:
                  return (1-exp(-M_PI*M_PI*get_time()))*sin(M_PI*p[0]);
                  break;
              case 5:
                  return get_time()*cos(M_PI*p[0]);
                  break;
              case 6:
                  return 0.5-abs(p[0]-0.5);
                  break;
              case 7:
                  return (p[0] <= 0.5)?1.0:0;
                  break;
              case 8:
                  return exp(1*(1-get_time()))*sin(2*M_PI*p[0])*sin(2*M_PI*p[0]);
                  break;
              case 9:
                  return exp(2*get_time());
                  break;
              case 10:
                  return get_time();
                  break;
              case 11:
                  return 1;
                  break;
              default:
                  abort();
                  break;
          }
      }
      void vector_value(const Point<1>& p, Vector<double>& values) const
      {
          values[0] = value(p);
      }
  };

/*
  Some driving terms on the cube
  1: u(x,y) = (1-exp(-pi^2 t))*sin(pi*x)*sin(pi*y), f = 2*pi*sin(pi*x)*sin(pi*y); u_t = u'' + f; homogeneous Dirichlet BC
  2: u(x,y) = x*(1-x)*y*(1-y), f = -2(x*(1-x)+y*(1-y)); -u'' = f; homogeneous Dirichlet BC
  3: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y); inhomogeneous BC ... not supported
  //4: -delta u + u = f : u = (2x^3-3x^2+1)*(2y^3-3y^2+2), -delta u + u = f =-(12x-6)(2y^3-3y^2+2) -(12y-6)(2x^3-3x^2+2) +  (2x^3-3x^2+2)(2y^3-3y^2+2)
  4: -delta u + u = f : u = (x^2*(x-1)^2*y^2*(y-1)^2, -delta u + u = f = - (12x^2-12*x+2)*y^2*(y-1)^2 - (12y^2-12*y+2)*x^2*(x-1)^2 + x^2*(x-1)^2*y^2*(y-1)^2 homogeneous Neumann BC
  5: u = t*x^2(x-1)^2y^2(y-1)^2+10; u_t = delta u+f; f = - (12x^2-12*x+2)*y^2*(y-1)^2*t - (12y^2-12*y+2)*x^2*(x-1)^2*t + x^2*(x-1)^2*y^2*(y-1)^2
*/
template <unsigned int N>
class Driving_f2D
  : public Function<2,double>
{
public:
  virtual ~Driving_f2D() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 2:
      return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
      break;
    case 3:
      return (200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
              * exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
      break;
   case 4:
      //return -(12*p[0]-6)*(2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)
      //       -(12*p[1]-6)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2)
      //       +(2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2);
      return -(12*p[0]*p[0]-12*p[0]+2)*p[1]*p[1]*(p[1]-1)*(p[1]-1)
             -(12*p[1]*p[1]-12*p[1]+2)*p[0]*p[0]*(p[0]-1)*(p[0]-1)
             + p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1); //
      break;
    case 5:
      return -(12*p[0]*p[0]-12*p[0]+2)*p[1]*p[1]*(p[1]-1)*(p[1]-1)*get_time()
             -(12*p[1]*p[1]-12*p[1]+2)*p[0]*p[0]*(p[0]-1)*(p[0]-1)*get_time()
             + p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1);
      break;
    case 6:
      return 0;
      break;
    default:
      abort();
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

/*
  Some driving terms for the parabolic problem on the cube
  1: u(x,y) = (1-exp(-pi^2 t))*sin(pi*x)*sin(pi*y), f = 2*pi*sin(pi*x)*sin(pi*y); u_t = u'' + f
  6: u= sum of 2 gauss kernels
*/
template <unsigned int N>
class Exact_Sol2D
  : public Function<2,double>
{
public:
  virtual ~Exact_Sol2D() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return (1-exp(-2*M_PI*M_PI*get_time()))*sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 2:
      return p[0]*(1-p[0])*p[1]*(1-p[1]);
      break;
    case 3:
      return exp(-50*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
      break;
    case 4:
      //return (2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2);
      return p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1);
      break;
    case 5:
      return get_time()*p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1)+10.0;
      break;
    case 6:
      return 2*exp(-50*((p[0]-0.2)*(p[0]-0.2)+(p[1]-0.3)*(p[1]-0.3))) + 3*exp(-40*((p[0]-0.8)*(p[0]-0.8)+(p[1]-0.6)*(p[1]-0.6)));
      break;
    default:
      abort();
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
    cout << "Testing linear parabolic equations..." << endl;
    const int d  = _D;
    const int dT = _DT;
    const int dim = _DIMENSION;
    const int offset = _OFFSET; // Radius for the index set Lambda
    typedef PBasis<d,dT> Basis1d; string basis_str = "P";
    //typedef DSBasis<d,dT> Basis1d; string basis_str = "DS";

    typedef TensorBasis<Basis1d,dim> Basis;
    typedef Basis::Index Index;
    typedef Index::level_type index_lt;

    FixedArray1D<bool,(2*dim)> bc;
    if (dim==1) {bc[0] = bc[1] = _RANDBEDINGUNGEN;}
    else { bc[0] = bc[1] = bc[2] = bc[3] = _RANDBEDINGUNGEN;}

    clock_t tstart, tend;
    double time;
    unsigned int jmax;

    cout << "main:: setting up equation" << endl;

    // uncoment for constant RHS ...
    //ConstantFunction<dim> constant_rhs(Vector<double>(1, "2.0"));
    //PoissonBVP<dim> elliptic(&constant_rhs);

    const unsigned int N=_PROBLEM_NO;
#if _DIMENSION == 1
    Driving_f1D<N> cf;
    //Driving_f1D<N> uexact;
    Exact_Sol1D<N> uexact;
#else
    Driving_f2D<N> cf;
    Exact_Sol2D<N> uexact;
#endif
    cf.set_time(0);
    uexact.set_time(0);

    IdentityBVP<dim> identity(&cf); // rechte Seite muss gesetzt werden, um dualen coeff vektor berechnen zu koennen
    double eta = 1;
    //PertubedBVP<dim> elliptic(&cf,eta); // -eta u'' + u
    //ZeroBVP<dim> elliptic(&cf); // 0 Operator. Just for testing Problems where u'' should be 0. That are for instance 9,10
    PoissonBVP<dim> elliptic(&cf); // -u''
    //PoissonBVP<dim> poisson(&cf);
#if _VORBERECHNETE_RHS // initialize the RHS vectors of the two cached problems for testing
    TensorEquation<Basis1d,dim,Basis> eq(&elliptic, bc, true);
    TensorEquation<Basis1d,dim,Basis> gram(&identity, bc, true);
    //TensorEquation<Basis1d,dim,Basis> poi(&poisson, bc, true);
    jmax = multi_degree(eq.basis_.j0())+offset;
    eq.set_jmax(jmax, true);
    gram.set_jmax(jmax, true);
    //poi.set_jmax(jmax, true);
#else // Do not initialize the RHS vectors of the two cached problems, since they arent used in the parabolic problem
    TensorEquation<Basis1d,dim,Basis> eq(&elliptic, bc, false);
    TensorEquation<Basis1d,dim,Basis> gram(&identity, bc, false);
    //TensorEquation<Basis1d,dim,Basis> poi(&poisson, bc, false);
    jmax = multi_degree(eq.basis_.j0())+offset;
    eq.set_jmax(jmax, false);
    gram.set_jmax(jmax, false);
    //poi.set_jmax(jmax, false);
#endif

    cout << "main:: setting up problems" << endl;
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem(&eq,5.5,35.0);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctgramian(&gram,5.5,35.0);
    //CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctpoisson(&poi,5.5,35.0);

    cout << "main:: compute time-constant driving term f and initial value u0" << endl;
    
    double tol = 1e-3;
    unsigned int maxlevel = ctproblem.basis().get_jmax();
    InfiniteVector<double,Index> f, u0;
    // Entwicklung bezüglich der primalen Basis
    ctproblem.basis().expand(&cf, false, jmax, f);
    ctproblem.basis().expand(&uexact, false, jmax, u0);


    // Entwicklung bezüglich der dualen Basis (weil das so sein muss für die Anfangswerte)
    // gleichbedeutend mit :: compute ctproblem.basis().expand(&cf, TRUE, jmax, f)
    InfiniteVector<double,Index> fT, u0T;
    ctgramian.set_f(&cf);
    CDD1_SOLVE(ctgramian, tol, fT, maxlevel, tensor_simple);
    fT.scale(&ctgramian,-1);
    fT.compress(1e-14); // Entwicklungskoeffizienten von cf bezüglich der primalen basis, dh <cf,\tilde psi_lambda>=fT_lambda

    // gleichbedeutend mit :: compute ctproblem.basis().expand(&uexact, TRUE, jmax, u0);
    ctgramian.set_f(&uexact);
    CDD1_SOLVE(ctgramian, tol, u0T, maxlevel, tensor_simple);
    u0T.scale(&ctgramian,-1);
    u0T.compress(1e-14);


#if 1 // Kontrollplot und Kontrollauswertung der rechten Seite
    SampledMapping<dim> s_u0T(evaluate(ctproblem.basis(), u0T, true, d+dim+offset+1));
    std::ofstream ustream_u0T;
    ustream_u0T.open("u0_init.m");
    s_u0T.matlab_output(ustream_u0T);
    ustream_u0T.close();

    SampledMapping<dim> s_fT(evaluate(ctproblem.basis(), fT, true, d+dim+offset+1));
    std::ofstream ustream_fT;
    ustream_fT.open("f_init.m");
    s_fT.matlab_output(ustream_fT);
    ustream_fT.close();


    const double h = 1e-6;
    InfiniteVector<double,Index> fhelp;
    cf.set_time(h);
    ctgramian.set_f(&cf);
    CDD1_SOLVE(ctgramian, tol, fhelp, maxlevel, tensor_simple); // => fhelp = (<f,\tilde\psi_\lambda>_\lambda)
    fhelp.scale(&ctgramian,-1);
    fhelp.compress(1e-14);
    SampledMapping<dim> s_2(evaluate(ctproblem.basis(), fhelp, true, d+dim+offset+1));
    std::ofstream ustream_2;
    ustream_2.open("fhelp.m");
    s_2.matlab_output(ustream_2);
    ustream_2.close();

    //zum Vergleich
    ctproblem.basis().expand(&cf, false, maxlevel, fhelp); // expand in the primal basis
    fhelp.compress(1e-14);
    SampledMapping<dim> s_3(evaluate(ctproblem.basis(), fhelp, true, d+dim+offset+1));
    std::ofstream ustream_3;
    ustream_3.open("fhelp2.m");
    s_3.matlab_output(ustream_3);
    ustream_3.close();
    cf.set_time(0);
    ctgramian.set_f(&cf);
#endif

    tstart = clock();

    cout << "setting up the parabolic problem"<<endl;
#if _TIME_CONSTANT_DRIVING_TERM
    // with time independent RHS given by an InfiniteVector
    LinearParabolicEquationTensor<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > parabolic(&ctproblem, &ctgramian, u0T, f, jmax);
#else
    // with time dependent RHS given by a function:
    LinearParabolicEquationTensor<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > parabolic(&ctproblem, &ctgramian, u0T, &cf, jmax);
#endif

    //expand(&func_u0, ctproblem.basis(), false, jmax, u0);
    //ctproblem.basis().expand(&func_u0, false, jmax, u0);

  // handle different test cases:
  // 1: u0 = hat function, f(t)=0
  // 2: u0 = haar function, f(t)=0
  // 3: u0 = 0, f(t,x) = pi^2*sin(pi*x), u(t,x)=(1-exp(-pi^2*t))*sin(pi*x)
  // 4: u(t,x)=t*x*(1-x)^3+(1-t)*x^3*(1-x), f(t,x)=u_t(t,x)-u_{xx}(t,x)
  // 5: u(t,x)=t*x*(1-x), f(t,x) = x*(1-x)+2*t
  // 6: u(t,x)=x*(1-x), f(t,x) = 2
  // 7: u(t,x)=exp(-100*(x-0.6+0.2*t)^2), f(t,x)=(224-40x-8t-(120-200x-40t)^2)*u(t,x)
  // 8: u(t,x)=sin(pi*t)*x*(1-x)^3+(1-sin(pi*t))*x^3*(1-x), f=u_t-u_{xx}
  // 9: u0 = haar function, f(t)=chi_{[0,1/2)}
  // 10: u0 = 0, f(t)=chi_{[0,1/2)}(t)*chi_{[1/4,3/4]}(x)



#if _TESTTHEHELPER
    /* test the helper.
     * output:
     * u_helper = (aI+A) u = f
     * u_ctproblem = A u = f
     * u_helper2 = aI+I = f
     * u_gramian = I u = f
     */
    cout << "test LinParEqTenROWStageEquationHelper" << endl;
    double alpha = _HELPERSALPHA; // transform a laplacian into a pertubed problem and vice versa
    
    LinParEqTenROWStageEquationHelper<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > helper(alpha, &ctproblem, &ctgramian, f);
    
    //CDD1_SOLVE(helper, tol, result, jmax_); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
    InfiniteVector<double,Index> u_epsilon;
    CDD1_SOLVE(helper, tol, u_epsilon,maxlevel,tensor_simple);
    u_epsilon.scale(&helper, -1); // Dx -> x
    cout << "  ... done. u_epsilon.size() = " << u_epsilon.size() << endl;
    cout << "saving computed solution to u_helper.m" << endl;
    SampledMapping<dim> s(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
    std::ofstream u_stream("u_helper.m");
    s.matlab_output(u_stream);
    u_stream.close();
    cout << "  ... done."<<endl;


    cout << "CDD1_SOLVE with identity problem" << endl;
    LinParEqTenROWStageEquationHelper<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > helper2(alpha, &ctgramian, &ctgramian, f);
    CDD1_SOLVE(helper2, tol, u_epsilon,maxlevel,tensor_simple);
    u_epsilon.scale(&helper2, -1); // Dx -> x
    SampledMapping<dim> s_h2(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
    std::ofstream u_stream_h2("u_helper2.m");
    s_h2.matlab_output(u_stream_h2);
    u_stream_h2.close();
    CDD1_SOLVE(ctgramian, tol, u_epsilon,maxlevel,tensor_simple);
    u_epsilon.scale(&ctgramian, -1); // Dx -> x
    SampledMapping<dim> s_gramian(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
    std::ofstream u_stream_gramian("u_ctgramian.m");
    s_gramian.matlab_output(u_stream_gramian);
    u_stream_gramian.close();
    cout << "  ... done."<<endl;

    cout << "CDD1_SOLVE with cached problem" << endl;
    CDD1_SOLVE(ctproblem, tol, u_epsilon,maxlevel,tensor_simple);
    u_epsilon.scale(&ctproblem, -1); // Dx -> x
    cout << "  ... done. u_epsilon.size() = " << u_epsilon.size() << endl;
    cout << "saving computed solution to u_ctproblem.m" << endl;
    SampledMapping<dim> s2(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
    std::ofstream u_stream2("u_ctproblem.m");
    s2.matlab_output(u_stream2);
    u_stream2.close();
    cout << "  ... done."<<endl;
/*
    cout << "CDD1_SOLVE with poisson problem" << endl;
    CDD1_SOLVE(ctpoisson, tol, u_epsilon,maxlevel,tensor_simple);
    u_epsilon.scale(&ctpoisson, -1); // Dx -> x
    SampledMapping<dim> s_poi(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
    std::ofstream u_stream_poi("u_poisson.m");
    s_poi.matlab_output(u_stream_poi);
    u_stream_poi.close();
*/

    //for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda) 
    //{
//        Lambda.insert(lambda);
//        if (lambda == itend) break;
//    }
#endif

#if _TESTPARABOLIC_FUNCTIONALITY
    cout << "Testing functionalitry of LinearParabolicEquationTensor"<<endl;

#if _DIMENSION == 1
    // why isn't the diagonal on I == 1 ??
    double value;
    cout << "Testing diagonal entries of ctgramian."<<endl;
    cout << "ctgramian.basis().degrees_of_freedom() = " << ctgramian.basis().degrees_of_freedom() << endl;
    cout << "(ctproblem.basis().degrees_of_freedom() = " << ctproblem.basis().degrees_of_freedom() << ")" << endl;
    Basis1d einebasis(false,false);
    //! Gramian
    IntervalGramian<Basis1d> G(einebasis, InfiniteVector<double,Basis1d::Index>());
    //! cached Gramian
    CachedProblem<IntervalGramian<Basis1d> > GC(&G);

    einebasis.set_jmax(jmax);
    cout << "einebasis.degrees_of_freedom() = " << einebasis.degrees_of_freedom() << endl;
    cout << " the following diagonal entry of I is not equal to 1: " << endl;

    for (unsigned int i=0; i< ctgramian.basis().degrees_of_freedom();++i)
    {
        // diagonal entries
        value = ctgramian.D(ctgramian.basis().get_wavelet(i));
        if (value != 1.0)
        {
            cout << " i = " << i << " " << *(ctgramian.basis().get_wavelet(i));
            //cout << " value = " << ctgramian.D(ctgramian.basis().get_wavelet(i)) << endl;
            //cout << " compare = " << GC.D(einebasis.get_wavelet(i)) << endl;
            cout << " value = " << gram.a(gram.basis_.get_wavelet(i),gram.basis_.get_wavelet(i)) << endl;
            cout << " i = " << i << " " << *(einebasis.get_wavelet(i));
            cout << " compare = " << G.a(einebasis.get_wavelet(i),einebasis.get_wavelet(i)) << endl;
        }
        for (unsigned int j=0;j<ctgramian.basis().degrees_of_freedom();++j)
        {
            if (abs( gram.a(gram.basis_.get_wavelet(i),gram.basis_.get_wavelet(j))- G.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j)) ) > 1e-14)
            {
                cout << "PROBLEM!" << endl;
                cout << "i = " << i << " j = " << j << " gram.a = " << gram.a(gram.basis_.get_wavelet(i),gram.basis_.get_wavelet(j)) << " G.a = " << G.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j)) << " difference = " << abs( gram.a(gram.basis_.get_wavelet(i),gram.basis_.get_wavelet(j))- G.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j)) ) << endl;
            }
        }
    }

    //

        for (unsigned int i=0; i< ctgramian.basis().degrees_of_freedom();++i)
    {
        for (unsigned int j=0;j<ctgramian.basis().degrees_of_freedom();++j)
        {
            if (abs( ctgramian.a(ctgramian.basis().get_wavelet(i),ctgramian.basis().get_wavelet(j))- GC.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j)) ) > 1e-14)
            {
                cout << "PROBLEM mit ctgramian und GC !" << endl;
                cout << "i = " << i << " j = " << j 
                     << " ctgramian.a = " << ctgramian.a(ctgramian.basis().get_wavelet(i),ctgramian.basis().get_wavelet(j))
                     << " CG.a = " << GC.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j))
                     << " difference = " << abs( ctgramian.a(ctgramian.basis().get_wavelet(i),ctgramian.basis().get_wavelet(j))- GC.a(einebasis.get_wavelet(i),einebasis.get_wavelet(j)) ) << endl;
            }
        }
    }
#endif
    

    // some transforming from primal to dual side and some plots
    int steps = 1<<4;
    double h2 = 1.0/steps;
    InfiniteVector<double,Index> vec,res;
    cout << "  h = " << h2 << ":" << endl;
    cout << "Testing evaluate_f and evaluate_ft for different t and v=0"<<endl;
    std::ofstream resultstream2;
    for (unsigned int i=0; i<=steps; ++i)
    {
        ostringstream output_f, output_ft;
        output_f << "f" << i << ".m";
        output_ft << "ft" << i << ".m";

        parabolic.evaluate_f(i*h2,vec,1e-6,res);
        resultstream2.open(output_f.str().c_str());
        SampledMapping<dim> eval_f_plot(evaluate(ctproblem.basis(), res, false, d+dim+offset+1));
        eval_f_plot.matlab_output(resultstream2);
        resultstream2.close();

        parabolic.evaluate_ft(i*h2,vec,1e-6,res);
        resultstream2.open(output_ft.str().c_str());
        SampledMapping<dim> eval_ft_plot(evaluate(ctproblem.basis(), res, false, d+dim+offset+1));
        eval_ft_plot.matlab_output(resultstream2);
        resultstream2.close();
    }
    cout << "test preprocess_rhs_share" << endl;
    InfiniteVector<double,Index> u0Tp(u0T);
    parabolic.preprocess_rhs_share(u0Tp,1e-6);

    SampledMapping<dim> u0true(evaluate(ctproblem.basis(),    u0, true, d+dim+offset+1));
    SampledMapping<dim> u0false(evaluate(ctproblem.basis(),   u0, false, d+dim+offset+1));
    SampledMapping<dim> u0Ttrue(evaluate(ctproblem.basis(),   u0T, true, d+dim+offset+1));
    SampledMapping<dim> u0Tfalse(evaluate(ctproblem.basis(),  u0T, false, d+dim+offset+1));
    SampledMapping<dim> u0Tptrue(evaluate(ctproblem.basis(),  u0Tp, true, d+dim+offset+1));
    SampledMapping<dim> u0Tpfalse(evaluate(ctproblem.basis(), u0Tp, false, d+dim+offset+1));
    resultstream2.open("u0true.m");    u0true.matlab_output(resultstream2);    resultstream2.close();
    resultstream2.open("u0false.m");   u0false.matlab_output(resultstream2);   resultstream2.close();
    resultstream2.open("u0Ttrue.m");   u0Ttrue.matlab_output(resultstream2);   resultstream2.close();
    resultstream2.open("u0Tfalse.m");  u0Tfalse.matlab_output(resultstream2);  resultstream2.close();
    resultstream2.open("u0Tptrue.m");  u0Tptrue.matlab_output(resultstream2);  resultstream2.close();
    resultstream2.open("u0Tpfalse.m"); u0Tpfalse.matlab_output(resultstream2); resultstream2.close();

#if _DIMENSION == 1
    // use the preprocessing from the original 1d linear parabolic problem:
    //u0Tp = u0T;
    InfiniteVector<double, CachedProblem<IntervalGramian<Basis1d> >::Index> help1, help2;
    Vector<double> help3, help4;
    // umspeichern:

    help3.resize(einebasis.degrees_of_freedom());
    for (InfiniteVector<double,Index>::const_iterator it(u0T.begin()), itend(u0T.end()); it!=itend; ++it )
    {
        help3[it.index().number()]=*it;
        //help1.set_coefficient(it
        //        it->second())
    }
    //APPLY(GC, help1, 1e-6, help2, jmax,St04a); // geht nicht weil das ein nichttensoraufruf wäre

    set<CachedProblem<IntervalGramian<Basis1d> >::Index> Lambda;
    for (unsigned int i=0; i< einebasis.degrees_of_freedom(); ++i)
    {
        Lambda.insert(einebasis.get_wavelet(i));
    }

    SparseMatrix<double> A_Lambda;
    setup_stiffness_matrix(GC, Lambda, A_Lambda);
    help4.resize(help3.size());
    A_Lambda.apply(help3,help4);

    //for (InfiniteVector<double,CachedProblem<IntervalGramian<Basis1d> >::Index>::const_iterator it(einebasis.get_wavelet(0), itend(einebasis.get_wavelet(einebasis.degrees_of_freedom()-1)); it!=itend; ++it )
    for (unsigned int i=0; i< help4.size(); ++i)
    {
        help1.set_coefficient(einebasis.get_wavelet(i),help4[i]);
    }


    SampledMapping<dim> u0Tp2true(evaluate(einebasis,  help1, true, d+dim+offset+1));
    SampledMapping<dim> u0Tp2false(evaluate(einebasis, help1, false, d+dim+offset+1));
    resultstream2.open("u0Tp2true.m");  u0Tp2true.matlab_output(resultstream2);  resultstream2.close();
    resultstream2.open("u0Tp2false.m"); u0Tp2false.matlab_output(resultstream2); resultstream2.close();
#endif

    cout << "should be 0 = " << l2_norm(u0-u0Tp) << endl;

    // are ctgramian and GC equal?
#endif
#if _TESTPARABOLIC_RUN
  // einzelner Testlauf mit konstanter Schrittweite, Ausgabe der ell_2-Fehler
    /* output:
     * u0.m,...,u2^expo.m
     */
  cout << "* testing linearly-implicit scheme (with constant stepsizes):" << endl;

  ROWMethod<InfiniteVector<double,Index> > method(WMethod<InfiniteVector<double,Index> >::ROS2);
  method.set_preprocessor(&parabolic);
  OneStepScheme<InfiniteVector<double,Index> >* scheme = &method;
  IVPSolution<InfiniteVector<double,Index> > results;

  InfiniteVector<double,Index> temp, result, error_estimate, temp_f, temp_ft;

  for (int expo = 3; expo <= 3; expo++) {
    temp = parabolic.u0;

    const int resolution = d+dim+offset+1;
    SampledMapping<dim> u0_plot(evaluate(ctproblem.basis(), temp, true, resolution));

    std::ofstream resultstream;
    resultstream.open("u0.m");
    u0_plot.matlab_output(resultstream);
    resultstream.close();

    // Werte evaluate_f oder Teile davon oder von evaluate_ft aus
#if 0
    // das ganze evaluate_f
    parabolic.evaluate_f(0,temp,1e-6,temp_f);
#else
#if 0
    // Der zweite Teil von evaluate_f, konstante rechte seite
    if (_TIME_CONSTANT_DRIVING_TERM == 1)
    {
        temp_f = f;
    }
#else
#if 1

#if _EXPANSIONTYPE_F_
    // hier stand verwirrrenderweise: elliptic->basis().expand(f_, true, jmax_, w); // expand in the dual (!) basis
    // Zum berechnen der Entwicklungskoeffizienten bezüglich der primalen Basis:
    /*
    ctgramian.set_f(&cf);
    cout << "cf.get_time = "<< cf.get_time() << endl;
    CDD1_SOLVE(ctgramian, 1e-6, temp_f, maxlevel, tensor_simple);
    temp_f.scale(&ctgramian,-1);
     */
#else
    /*
    //identity->set_f(f_);
    cout << "cf.get_time = "<< cf.get_time() << endl;
    ctproblem.basis().expand(&cf, false, maxlevel, temp_f); // expand in the primal basis
     * */
#endif

#endif
#endif
#endif
    

    // Werte evaluate_ft oder Teile von evaluate_f aus
#if 0
    // das ganze evaluate_ft:
    parabolic.evaluate_ft(0,temp,1e-6,temp_ft);
#else
#if 1

    // der erste Teil aus evaluate_f, der Aufruf von APPLY:
    temp_f=temp;
    temp_f.scale(&ctproblem, 1); // w = Dv
    APPLY_TENSOR(ctproblem, temp_f, 1e-6, temp_ft, maxlevel, tensor_simple); // yields -D^{-1}AD^{-1}w
    temp_ft.scale(&ctproblem, 1);
    temp_ft.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
#endif
#endif

    SampledMapping<dim> f0(evaluate(ctproblem.basis(), temp_f, false, resolution)); // output of evaluate_f should be dual coefficients
    resultstream.open("f0.m");
    f0.matlab_output(resultstream);
    resultstream.close();
/*
    SampledMapping<dim> ft0(evaluate(ctproblem.basis(), temp_ft, false, resolution)); // output of evaluate_ft should be dual coefficients
    resultstream.open("ft0.m");
    ft0.matlab_output(resultstream);
    resultstream.close();
*/
    int N = 1<<expo;
    double h = 1.0/N;
    cout << "  h = " << h << ":" << endl;

    results.t.push_back(0);
    results.u.push_back(temp);

    for (int i = 1; i <= N; i++) {
//       cout << "---------------- before increment() -----------------------" << endl;

        cout << "main:: increment with i = " << i << endl;
      method.increment(&parabolic, (i-1)*h, temp, h, result, error_estimate, 1e-5);
      //scheme->increment(&parabolic, (i-1)*h, temp, h, result, error_estimate, 1e-5);
      results.t.push_back(i*h);
      results.u.push_back(result);

      ostringstream output_filename, output_filename2, output_filename3;
      output_filename << "u" << i << ".m";
      output_filename2 << "f" << i << ".m";
      output_filename3 << "ft" << i << ".m";

      resultstream.open(output_filename.str().c_str());
      SampledMapping<dim> ui_plot(evaluate(ctproblem.basis(), result, true, resolution));
      ui_plot.matlab_output(resultstream);
      resultstream.close();
      //eval_f or parts thereof
            /*
      if (_TIME_CONSTANT_DRIVING_TERM == 1)
      {
          temp_f = f;
      }
       * */
#if _EXPANSIONTYPE_F_
/*
      // hier stand verwirrrenderweise: elliptic->basis().expand(f_, true, jmax_, w); // expand in the dual (!) basis
      // Zum berechnen der Entwicklungskoeffizienten bezüglich der primalen Basis:
      ctgramian.set_f(&cf);
      cout << "cf.get_time = "<< cf.get_time() << endl;
      CDD1_SOLVE(ctgramian, 1e-6, temp_f, maxlevel, tensor_simple);
      temp_f.scale(&ctgramian,-1);
      */
#else
      //identity->set_f(f_);
      /*
      cout << "cf.get_time = "<< cf.get_time() << endl;
      ctproblem.basis().expand(&cf, false, maxlevel, temp_f); // expand in the primal basis
       * */
#endif

      //eval_ft or parts thereof

      //parabolic.evaluate_ft(i*h,result,1e-6,temp_ft);
      /*
      temp_f=result;
      temp_f.scale(&ctproblem, 1); // w = Dv
      APPLY(ctproblem, temp_f, 1e-6, temp_ft, maxlevel, tensor_simple); // yields -D^{-1}AD^{-1}w
      temp_ft.scale(&ctproblem, 1);
      temp_ft.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
      */

      parabolic.evaluate_f(i*h,result,1e-6,temp_f);
      
      resultstream.open(output_filename2.str().c_str());
      SampledMapping<dim> fi_plot(evaluate(ctproblem.basis(), temp_f, true, resolution));
      fi_plot.matlab_output(resultstream);
      resultstream.close();
      /*
      resultstream.open(output_filename3.str().c_str());
      SampledMapping<dim> fti_plot(evaluate(ctproblem.basis(), temp_ft, false, resolution));
      fti_plot.matlab_output(resultstream);
      resultstream.close();
      */

//       cout << "---------------- after increment() -----------------------" << endl;

      InfiniteVector<double,Index> uexact_coeffs;
      uexact.set_time(i*h);
      ctproblem.basis().expand(&uexact, false, jmax, uexact_coeffs);
      
      ctgramian.set_f(&uexact);
      CDD1_SOLVE(ctgramian, tol, uexact_coeffs, maxlevel, tensor_simple);
      uexact_coeffs.scale(&ctgramian,-1);
      uexact_coeffs.compress(1e-14);

      cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result - uexact_coeffs) << endl;
      temp = result;
    }
  }
#endif

  return 0;
}

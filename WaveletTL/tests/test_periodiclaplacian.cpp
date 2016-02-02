#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/periodiclp.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <utils/function.h>

#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/ds_basis.h>
#include <interval/periodic.h>
#include <interval/periodic_support.h>
#include <Rd/cdf_basis.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>
#include <galerkin/periodic_gramian.h>
#include <galerkin/periodic_laplacian.h>

#include <adaptive/cdd2.h>

#define DSBASIS
#undef PERIODIC_CDFBASIS
#undef PERIODIC_QUARKLETFRAME




using namespace std;
using namespace WaveletTL;
using namespace MathTL;

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=t*(1-t), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};


/*1:y(t)=t^3-1.5*t^2+0.5*t, -y''(t)=-6*t+3*/
template <unsigned int N>
class PeriodicTestProblem
  : public PeriodicLaplacianProblem
{
public:
  double g(const double t) const {
    switch(N) {
    case 1:
      return (0.3<t<0.5 ? t*t*t-1.5*t*t+0.5*t : 0);
      break;
    case 2:
      //return t*t*t-1.5*t*t+0.5*t; 
        return 3.14;
    default:
      return 0;
      break;
    }
  }
};

class MyFunction : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            //return 6*p[0]-3;
            return p[0]*p[0]*p[0]-1.5*p[0]*p[0]+0.5*p[0];
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};


int main()
{
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD2_SOLVE ..." << endl;

  TestProblem<1> T;
  PeriodicTestProblem<2> TPer;

  const int d  = 3;
  const int dT = 3;
#ifdef DSBASIS
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;
#endif
  
#ifdef PERIODIC_CDFBASIS
  typedef PeriodicBasis<CDFBasis<d,dT> >Basis;
  typedef Basis::Index Index;
#endif
  
#ifdef PERIODIC_QUARKLETFRAME
  typedef PeriodicFrame<Quarklet_Frame<d,dT> >Basis;
  typedef Basis::Index Index;
#endif
  
//#ifdef PERIODIC_CDFBASIS
#if 0
  Basis basis;
  Basis::Support mysupp;
  int k1, k2, k3, k4, k5, k6, k7, k8, j;
  const int j1 = 3;
  const int e1 = 0;
  const int l1 = 2;
  const int j2 = 3;
  const int e2 = 1;
  const int l2 = 4;
  const int m = 3;
  const int a = 0;
  const int b = 3;
  Index lambda(j1,e1,l1);
  Index nu(j2,e2,l2);
  basis.support(lambda, k1, k2);
  cout << ldexp(1.0,-j1-e1)*k1 << ", " << ldexp(1.0,-j1-e1)*k2 << endl;
  basis.support(nu, k7, k8);
  cout << ldexp(1.0,-j2-e2)*k7 << ", " << ldexp(1.0,-j2-e2)*k8 << endl;
  //support(basis, lambda, k3, k4);
  //cout << ldexp(1.0,-4)*k3 << ", " << ldexp(1.0,-4)*k4 << endl;
  InfiniteVector<double,Index> y;
  PeriodicIntervalGramian<CDFBasis<d,dT> > problem(basis, y);
  cout << ldexp(1.0,-m)*a << ", " << ldexp(1.0,-m)*b << endl;
  cout << intersect_supports(basis, lambda, m, a, b, j, k5, k6) << endl;
  cout << ldexp(1.0,-j)*k5 << ", " << ldexp(1.0,-j)*k6 << endl;
  cout << intersect_supports(basis, lambda, nu, mysupp)<< endl;
  cout << ldexp(1.0,-mysupp.j)*mysupp.k1 << ", " << ldexp(1.0,-mysupp.j)*mysupp.k2 << endl;
  
  
#endif
  
  
  Basis basis;
  //basis.set_jmax(12);
  //cout << basis.get_wavelet(0) << endl;
  //basis.full_collection[0];
  /*
  Index lambda(4,0,0, &basis);
  const int j = 5;
  const bool generators = 0;
  std::list<Index> intersecting;
  
  
  intersecting_wavelets(basis, lambda, j, generators, intersecting);
  //cout << lambda << endl;
  cout << intersecting.size() << endl;*/
  
#ifdef PERIODIC_CDFBASIS 
  
  /*InfiniteVector<double,Index> y;
  y[Index(4,0,2)]=0.072262;
  y[Index(4,0,3)]=0.0570544;
  y[Index(4,0,4)]=0.03125;
  y[Index(4,0,5)]=0.03125;
  y[Index(4,0,6)]=0.03125;
  y[Index(4,0,7)]=0.03125;
  y[Index(4,0,8)]=0.03125;
  y[Index(4,0,9)]=0.03125;
  y[Index(4,0,10)]=0.03125;
  y[Index(4,0,11)]=0.03125;
  y[Index(4,0,12)]=0.0570544;
  y[Index(4,0,13)]=0.072262;
  
  cout << y;*/
  
  
  PeriodicIntervalLaplacian<CDFBasis<d,dT> > problem(TPer, basis);
  typedef CachedProblem<PeriodicIntervalLaplacian<CDFBasis<d,dT> > > PROBLEM;
  PROBLEM cproblem(&problem);
  InfiniteVector<double, Index> F_eta, F_nu;
  Point<1, double> mypoint(0.1);
  Vector<double> value(1);
  value[0] = 2;
  ConstantFunction<1> f(value);
  MyFunction myfunction;
  //cout << "Funktionswert: " << myfunction.value(mypoint) << endl;
  const bool primal = 0;
  const int jmax = 5;
  
  basis.expand(&myfunction, primal, jmax, F_eta);
  
  
  InfiniteVector<double, Index> F1_eta;
  cproblem.RHS(1e-6, F1_eta);
  cout << F1_eta << endl;
  const double nu1 = cproblem.norm_Ainv() * l2_norm(F1_eta);
  cout << "nu = " << nu1 << endl;
  InfiniteVector<double, Index> solution1_epsilon;
  double epsilon1 = 1e-6;
  CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon);
  //F_eta[Index(3,1,2)] = 1;
  //F_eta[Index(3,0,3)] = 1;
  //cproblem.RHS(1e-6, F_nu);

  //cout << F_eta << endl;
  
  //cout << F_nu << endl;
  //F_eta..scale(&cproblem, -1); /* scaling because ... */
  SampledMapping<1> sm1(basis.evaluate(solution1_epsilon, 8)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
  std::ofstream u_stream1("../../Arbeitsfläche/plotthis.m");
  sm1.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "mygraph" << endl;
  u_stream1.close();
  
  
  //cout << problem.a(lambda,nu) << endl;
#endif 
#ifdef DSBASIS
  SturmEquation<Basis> problem(T);
  
  typedef CachedProblem<SturmEquation<Basis> > PROBLEM;
  
  PROBLEM cproblem(&problem);
  InfiniteVector<double, Index> F_eta;
  cproblem.RHS(1e-6, F_eta);
  cout << F_eta << endl;
  
  
  InfiniteVector<double, Index> F1_eta;
  cproblem.RHS(1e-6, F1_eta);
  cout << F1_eta << endl;
  const double nu1 = cproblem.norm_Ainv() * l2_norm(F1_eta);
  cout << "nu = " << nu1 << endl;
  InfiniteVector<double, Index> solution1_epsilon;
  double epsilon1 = 0.8;
  CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon);
  cout << "Fertig" << endl;
  
  solution1_epsilon.scale(&cproblem, -1);
  //cout << F_eta << endl;
  SampledMapping<1> sm1(evaluate(cproblem.basis(), solution1_epsilon, true, 16)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
  std::ofstream u_stream1("../../Arbeitsfläche/plotthis.m");
  sm1.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "mygraph" << endl;
  u_stream1.close();
#endif

  //InfiniteVector<double, Index> F_eta;
  //cproblem.RHS(1e-6, F_eta);
  //cout << F_eta;
  //const double nu = cproblem.norm_Ainv(); //* l2_norm(F_eta);
  //cout << "nu = " << nu << endl;
#if 0
  InfiniteVector<double, Index> u_epsilon;
  
  Matrix<double> mymat(3,3,"1 0 0 0 3 0 0 0 -4");
  
  cout << mymat;
  
  const double tol = 0.1;
  double lambdamin, lambdamax;
  unsigned int iterations;
  Vector<double> xk(3, false);
  xk = 1;
  
  
  double ew = PowerIteration(mymat, xk, tol, 100, iterations);
  cout << ew;
#endif  
  //LanczosIteration(mymat, tol, lambdamin, lambdamax, 100, iterations);
  
  //cout << lambdamin << ", " << lambdamax << endl;
  
  
  //CDD2_SOLVE(cproblem, nu, 1e-1, u_epsilon);
  
  /*const double eta = 0.1;
  InfiniteVector<double, Index> v,w;
  v[Index(0,0,0)] = 1;
  
  APPLY(cproblem, v, eta, w, 8, CDD1);*/

  //cout << w << endl;
  
  
  
  return 0;
}





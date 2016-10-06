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
#include <interval/i_q_index.h>
#include <interval/ds_basis.h>
#include <interval/periodic.h>
#include <interval/periodic_support.h>
#include <interval/periodic_frame_support.h>
#include <interval/periodic_frame.h>
#include <Rd/cdf_basis.h>
#include <Rd/quarklet_frame.h>
#include <galerkin/sturm_equation.h>
#undef DSBASIS
#define PERIODIC_CDFBASIS
#undef PERIODIC_QUARKLETFRAME

#ifdef PERIODIC_CDFBASIS
#include <galerkin/cached_problem.h>
#endif
#ifdef PERIODIC_QUARKLETFRAME
#include <galerkin/cached_quarklet_problem.h>
#endif
#include <galerkin/periodic_gramian.h>
#include <galerkin/periodic_laplacian.h>
#include <galerkin/periodic_frame_laplacian.h>
#include <galerkin/periodic_frame_gramian.h>

#include <adaptive/cdd2.h>







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
    case 2:
      return 0;
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
    case 2:
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
    case 2:
      return 1;
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
      //return 6*t-3;
      break;
    case 2:
      return t*(1-t);
        //return 2.15;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};


/*1:y(t)=t^3-1.5*t^2+0.5*t, -y''(t)=-6*t+3 */
template <unsigned int N>
class PeriodicTestProblem
  : public PeriodicLaplacianProblem
{
public:
  double g(const double t) const {
    switch(N) {
    case 1:
      //u(t) = 100*(t*t*t-1.5*t*t+0.5*t)
        //return 100*(t*t*t-1.5*t*t+0.5*t);
        //return -6*t + 3;
        //return -12*t*t+12*t-2;
      if(t>0.5-0.0001 && t<0.5+0.0001){
          return 1/0.0001;
      }
      else return 0;
      break;
    case 2:
       return -2+12*t-12*t*t;
    case 3:
      //return 2*t;
      return -9*M_PI*M_PI*cos(3*M_PI*t);
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
            double t = 0;
            if(0<=p[0] && p[0]<0.5){
                if(p[0]<0.25)
                    t=8*p[0];
                else
                    t=4-8*p[0];
            }
            return t;
            //return p[0]*p[0]*p[0]-1.5*p[0]*p[0]+0.5*p[0];
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};

class MySinus : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            double t = 0;
            t = sin(2*M_PI*p[0]);
            
            return t;
            //return p[0]*p[0]*p[0]-1.5*p[0]*p[0]+0.5*p[0];
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
  const int jmax = 8;
  const int pmax = 3;

#ifdef DSBASIS
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  Basis basis(true,true);
  basis.set_jmax(jmax);
  SturmEquation<Basis> problem(T,basis);
  typedef CachedProblem<SturmEquation<Basis> > PROBLEM;
  PROBLEM cproblem(&problem);
  const char* basis_type = "DS";
#endif
  
#ifdef PERIODIC_CDFBASIS
  typedef PeriodicBasis<CDFBasis<d,dT> >Basis;
  typedef Basis::Index Index;
  Basis basis;
  basis.set_jmax(jmax);
  
  
  //PeriodicIntervalGramian<CDFBasis<d,dT> > problem(TPer, basis);
  PeriodicIntervalLaplacian<CDFBasis<d,dT> > problem(TPer, basis);
  //typedef CachedProblem<PeriodicIntervalGramian<CDFBasis<d,dT> > > PROBLEM;
  typedef CachedProblem<PeriodicIntervalLaplacian<CDFBasis<d,dT> > > PROBLEM;

  PROBLEM cproblem(&problem);
  const char* basis_type = "Periodic";
#endif
  
#ifdef PERIODIC_QUARKLETFRAME
  typedef PeriodicFrame<QuarkletFrame<d,dT> >Basis;
  typedef Basis::Index Index;
  Basis basis;
  
  const int J = 3;
  const int j = 5;
  const double a = 2;
  const double b = 2;
  
  
  basis.set_jpmax(jmax, pmax);
  PeriodicFrameIntervalLaplacian<QuarkletFrame<d,dT> > problem(TPer, basis);
  //PeriodicFrameIntervalGramian<QuarkletFrame<d,dT> > problem(TPer, basis);
  //typedef CachedQuarkletProblem<PeriodicFrameIntervalGramian<QuarkletFrame<d,dT> > > PROBLEM;
  typedef CachedQuarkletProblem<PeriodicFrameIntervalLaplacian<QuarkletFrame<d,dT> > > PROBLEM;
  PROBLEM cproblem(&problem);
  //cout << "Index: " << *basis.get_wavelet(1200) << endl;
  const char* basis_type = "Periodic Frame";
  
#endif
  
#if 0
 
  cout << cproblem.norm_A() << endl;
  cout << cproblem.norm_Ainv() << endl;
  const double cond_A = cproblem.norm_A() * cproblem.norm_Ainv();
  cout << cond_A << endl;
  const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "CDD2_SOLVE: rho=" << rho << endl;
  //const double theta = 0.333;
  const double theta = 2.0/7.0;
  cout << theta << endl;
  const double a = (log10(theta/3.0) / log10(rho));
  const int K = (int) ceil(log10(theta/3.0) / log10(rho));
  cout << a << endl;
  cout << K << endl;
  cout << 3*pow(rho,K) << endl;
  
#endif
  
  //cproblem.norm_A();
//  InfiniteVector<double,Index> f;
//  double eta = 0.01;
//  cproblem.RHS(eta, f);
//  cout << f << endl;
  
#if 1
  InfiniteVector<double, Index> F1_eta, v, w, u;
  cproblem.RHS(1e-6, F1_eta);
//  MySinus testfkt;
//  //Index lambda(2,3,0,0);
//  double eta = 0.01;
  
//  basis.expand(&testfkt, 1, jmax, v);
//  cout << v << endl;
//  basis.expand(&testfkt, 0, jmax, w);
//  cout << w << endl;
//  APPLY(cproblem, w, eta, u, jmax, CDD1);
//  cout << u << endl;
  
//  basis.expand(&testfkt, 1, jmax, pmax, v);
//  //cout << v << endl;
//  basis.expand(&testfkt, 0, jmax, pmax, w);
//  APPLY(cproblem, w,eta,u, jmax, DKOR, pmax, a,b);
  
  //w.compress(1e-5);
  //cout << w << endl;
  //cout << "Integral: " << basis.integrate(&testfkt, lambda) << endl;
  
  

  
  cproblem.norm_A();
  const double nu1 = cproblem.norm_Ainv() * l2_norm(F1_eta);
  cout << "nu = " << nu1 << endl;
  InfiniteVector<double, Index> solution1_epsilon;
  
  double epsilon1 = 1e-6;
#ifdef PERIODIC_QUARKLETFRAME
  CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon, jmax, DKR, pmax, a, b);
#endif
#ifdef PERIODIC_CDFBASIS
  //CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon, jmax, CDD1);
#endif
  cout << solution1_epsilon << endl;
  //cout << solution1_epsilon.size() << endl;
  cout << "Fertig" << endl;
  
#endif  
  
#if 0
  
  InfiniteVector<double,Index> f, Av;
  double teta = 0.0001;
  cproblem.RHS(eta, f);
  cout << f << endl;
  
  APPLY(cproblem, solution1_epsilon, teta, Av, jmax, CDD1);
  Av.compress(1e-3);
  cout << Av << endl;
  cout << l2_norm(f-Av) << endl;
  
#endif
  
#if 0  
  InfiniteVector<double, Index> u, v, Av, w, Atemp;
  double eta = 1;
  Point<1, double> mypoint(0.1);
  Vector<double> value(1);
  value[0] = 1;
  ConstantFunction<1> f(value);
  MyFunction myfunction;
  //cout << myfunction.value(mypoint)<< endl;
  //cout << "Funktionswert: " << myfunction.value(mypoint) << endl;
#ifdef PERIODIC_CDFBASIS
  Index lambda1(3,0,0), lambda2(3,0,0), lambda3(2,0,2), lambda4(2,0,3);
  basis.expand(&myfunction, 1, jmax, v);
  basis.expand(&myfunction, 0, jmax, w);
  //cout << v << endl;
  //cout << w << endl;
  /*v.COARSE(0.01,w);
  cout << w << endl;*/
#endif
#ifdef DSBASIS
  Index lambda1(4,0,12,&basis), lambda2(4,1,15, &basis);/* lambda3(2,0,2), lambda4(2,0,3);*/
#endif
  //cout << "Integral: " << basis.integrate(&myfunction, lambda1) << endl;
  
  
  
  
  w[lambda1]=1;
  //APPLY(cproblem, w, eta, Av, 8, CDD1);
  //cout << Av << endl;
  //cproblem.RHS(1e-6, u);
  //cout << cproblem.f(lambda1) << endl;
  //cout << u << endl;
  //cout << "Integral: " << cproblem.a(lambda1, lambda2) << endl;
  //cout << "Integral: " << basis.integrate(0, lambda1, lambda2) << endl;
  
  /*;
  v[lambda2]=0.1;
  v[lambda3]=0.01;
  v[lambda4]=0.001;*/
  
  //v.COARSE(0.356,w);
  
  /*cproblem.RHS(1e-6, w);
  
  cout << w << endl;
  APPLY(cproblem, w, eta, Av, 8, CDD1);
  cout << Av << endl;
  Av.COARSE(eta, Atemp);
  cout << Atemp << endl;*/
  
  //evaluate(cproblem.basis(), v, true, 16);
  
#endif
  
//  InfiniteVector<double, Index> solution1_epsilon;
//  const Index lambda(2,3,0,1);
//  const Index mu(1,4,1,10);
//  solution1_epsilon[lambda]= 1;
//  solution1_epsilon[mu]= 10;
//  cout << "Größe: " << solution1_epsilon.size()<< endl;
//  solution1_epsilon.compress(5);
//  cout << "Größe: " << solution1_epsilon.size()<< endl;

 
#if 0
      InfiniteVector<double, Index> testvkt;
//      cproblem.RHS(0.011, testvkt);
//      testvkt[Index(2,0,1)] = 0.891203;
//      testvkt[Index(2,0,3)] = -0.891203;
//      testvkt[Index(2,1,0)] = 0.0309778;
//      testvkt[Index(2,1,1)] = 0.0309778;
//      testvkt[Index(2,1,2)] = -0.0309778;
//      testvkt[Index(2,1,3)] = -0.0309778;
      
//      testvkt[Index(2,0,1)] = 1;
//      testvkt[Index(2,0,3)] = -0.580683;
//      testvkt[Index(2,1,0)] = -0.0544891;
//      testvkt[Index(2,1,1)] = -0.0544891;
//      testvkt[Index(2,1,2)] = 0.0544891;
//      testvkt[Index(2,1,3)] = 0.0544891;

      
      cout << cproblem.f(Index(3,0,0)) << endl;
      cout << cproblem.f(Index(3,0,1)) << endl;
      cout << cproblem.f(Index(3,0,2)) << endl;
      cout << cproblem.f(Index(3,0,3)) << endl << endl;
      cout << cproblem.f(Index(3,0,4)) << endl;
      cout << cproblem.f(Index(3,0,5)) << endl;
      cout << cproblem.f(Index(3,0,6)) << endl;
      cout << cproblem.f(Index(3,0,7)) << endl << endl;
      
      
      
      
      cout << testvkt << endl;
      testvkt.clear();
      testvkt[Index(3,0,1)] = 1;
      SampledMapping<1> sm1(basis.evaluate(testvkt, 8));

#endif
  
  
#if 1
  solution1_epsilon.scale(&cproblem, -1);
  //solution1_epsilon.compress(1e-4);
  //solution1_epsilon.COARSE(0.01, w);
  //SampledMapping<1> sm1(basis.evaluate(solution1_epsilon, 8)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
  //SampledMapping<1> sm1(basis.evaluate(F1_eta, 8));
  //SampledMapping<1> sm1(basis.evaluate(mu, 8));
  //basis.evaluate(solution1_epsilon, 8);
  //SampledMapping<1> sm1(basis.evaluate(w, 8));
 SampledMapping<1> sm1(basis.evaluate(solution1_epsilon, 8));
  //SampledMapping<1> sm1(basis.evaluate(v, 8));
#endif 

#ifdef DSBASIS 
  solution1_epsilon.scale(&cproblem, -1);
  w.scale(&cproblem, -1);
  //cout << solution1_epsilon << endl;
  //SampledMapping<1> sm1(evaluate(cproblem.basis(), solution1_epsilon, true, 16)); // Increase last parameter, if "Assertion `resolution >= 0' failed."
  SampledMapping<1> sm1(evaluate(cproblem.basis(), solution1_epsilon, true, 16));
#endif
  
#if 1  
  
  std::ofstream u_stream1("../../Desktop/plotthis.m");
  sm1.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('mygraph " << basis_type << "');" << endl;
  u_stream1.close();
  
#endif 
  
 //cout << solution1_epsilon << endl; 
  
  
  

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
  
//  const double eta = 0.1;
//  InfiniteVector<double, Index> v,w;
//  v[Index(0,3,0,0)] = 1;
//  
//  //APPLY(cproblem, v, eta, w, 8, CDD1); //hier weitermachen, APPLY testen mit CDD1 und DKOR @PHK
//  APPLY(cproblem, v, eta, w, jmax, DKOR, pmax, a, 1);
//  
//  //cproblem.norm_A();
//  
//  cout << w << endl;
//  
//  
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
  
  /*InfiniteVector<double, Index> F_eta, F_nu;
  Point<1, double> mypoint(0.1);
  Vector<double> value(1);
  value[0] = 2;
  ConstantFunction<1> f(value);
  MyFunction myfunction;
  //cout << "Funktionswert: " << myfunction.value(mypoint) << endl;
  const bool primal = 0;
  const int jmax = 5;
  
  basis.expand(&myfunction, primal, jmax, F_eta);*/
  
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
  
  
  return 0;
}






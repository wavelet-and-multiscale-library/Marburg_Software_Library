/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE

#undef BASIS
#define FRAME

#ifdef FRAME
#define _WAVELETTL_USE_TFRAME 1
#else
#define _WAVELETTL_USE_TBASIS 1
#endif
#define _DIM 2
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 2

#include <iostream>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <cube/tbasis_index.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tframe_evaluate.h>
#include <interval/pq_frame.h>
#include <cube/tframe.h>
#include <cube/tframe_evaluate.h>
#include <cube/tframe_indexplot.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <cube/tbasis_indexplot.h>
#include <galerkin/tframe_equation.h>
#include <galerkin/cached_quarklet_tproblem.h>
//#include "TestFunctions2d.h"




using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

//test functions and derivatives
// 1: f(x,y)=x(1-x)y(1-y), Laplace(x,y)=2x(1-x)+2y(1-y)
template <unsigned int N>
class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return p[0]*(1-p[0])*p[1]*(1-p[1]);
      break;
    case 2:
      return sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 3:
      return
	(p[0]<=0.5 ? 4*p[0]*p[0]*(3-4*p[0]) : (2-2*p[0])*(2-2*p[0])*(4*p[0]-1))
	* (p[1]<=0.5 ? 4*p[1]*p[1]*(3-4*p[1]) : (2-2*p[1])*(2-2*p[1])*(4*p[1]-1));
      break;
    default:
      return 0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

template <unsigned int N>
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
      return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
      break;
    case 2:
      return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
      break;
    case 3:
      return
	-((p[0]<=0.5 ? 24-96*p[0] : 96*p[0]-72)
	  * (p[1]<=0.5 ? 4*p[1]*p[1]*(3-4*p[1]) : (2-2*p[1])*(2-2*p[1])*(4*p[1]-1))
	  +(p[0]<=0.5 ? 4*p[0]*p[0]*(3-4*p[0]) : (2-2*p[0])*(2-2*p[0])*(4*p[0]-1))
	  * (p[1]<=0.5 ? 24-96*p[1] : 96*p[1]-72));
      break;
    default:
      return 1;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
  string toString() const
        {
            switch(N)
            {
                case 1:
                    return "L(x,y)=2x(1-x)+2y(1-y)";
                    break;
                case 2:
                    return "";
                    break;
                case 3:
                    return "1";
                    break;
                case 4:
                    return "";
                    break;
                   
                default:
                    return "TestProblem: N not defined.";
                    break;
            }
        }
};

int main()
{
    cout << "Testing tbasis adaptive" << endl;
    const int d  = 3;
    const int dT = 3;
    const int a=2;
    const int b=2;
    const unsigned int dim = 2; 
    const int jmax=8;
    const int pmax = 0;

        FixedArray1D<int,2*dim> s;          //set order of boundary conditions
    s[0]=0, s[1]=0; s[2]=0, s[3]=0; 
    //cout << s << endl;
    FixedArray1D<bool, 2*dim> bc;
    bc[0]=true, bc[1]=true, bc[2]=true, bc[3]=true;
    
    //Basis1d basis1d(true,true);
    //basis1d.set_jmax(jmax1);
    
    //Basis basis(bc); 
    //basis.set_jmax(jmax);
    
    mySolution<2> solution1;
    myRHS<2> rhs1;
    //Function2 function2;
    
    PoissonBVP<dim> poisson1(&rhs1);
    
#ifdef BASIS
    typedef PBasis<d,dT> Basis1d;
//    typedef Basis1d::Index Index1d;
    typedef TensorBasis<Basis1d,dim> Basis;
    typedef Basis::Index Index;
    TensorEquation<Basis1d,dim,Basis> eq(&poisson1, bc);
    eq.set_jmax(jmax);
    
    
#endif
    

    
#ifdef FRAME
    
    typedef PQFrame<d,dT> Frame1d;
//    typedef Basis1d::Index Index1d;
    typedef TensorFrame<Frame1d,dim> Frame;
    typedef Frame::Index Index;
    TensorFrameEquation<Frame1d,dim,Frame> eq(&poisson1, bc);

    eq.set_jpmax(jmax, pmax);
    Frame frame = eq.frame();

#endif  
    
    
    
    
#if 0
    abort();
#endif    
    
#if 1    
    
    
    
#if 0
    //plot stiffness matrix index version 
    int n=eq.basis().degrees_of_freedom();
    Matrix<double> L(n,n);
    int i=0;
    for (Index ind1 = eq.basis().first_generator(), itend = eq.basis().last_wavelet(jmax); ind1 <= itend; ++ind1, i++){
        int j=0;
        for (Index ind2 = eq.basis().first_generator(), itend = eq.basis().last_wavelet(jmax); ind2 <= itend; ++ind2, j++){
            //Index ind1=*it1;
            //Index ind2=*it2;
            L(i,j)=eq.a(ind2,ind1)/(eq.D(ind1)*eq.D(ind2));
        }     
    }
    L.matlab_output("L", "L",0);
    cout << "plot stiffness matrix done" << endl;
#endif
   
#ifdef NONADAPTIVE
    //setup index set
    set<Index> Lambda;  
    MultiIndex<int,dim> p;p[0]=0;p[1]=0;
    Index lambda = eq.frame().first_generator(eq.frame().j0(), p);
    int zaehler=0;
    for (int l = 0; l < eq.frame().degrees_of_freedom(); l++) {
        Lambda.insert(lambda);
        cout << lambda << ", Number: " << lambda.number() << endl;
        if(lambda==eq.frame().last_quarklet(jmax, p)){
            ++p;
            lambda=eq.frame().first_generator(eq.frame().j0(), p);
        }
        else
        ++lambda;
        ++zaehler;
        
    }
    //for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it){
        //cout << *it << endl; 
          //Index ind = *it;
          //cout << ind.number() << endl;
    //}
    
    cout << eq.norm_A() << endl;
    cout << eq.norm_Ainv() << endl;
    
    cout << "setup stiffness matrix" << endl;
    SparseMatrix<double> A;
    setup_stiffness_matrix(eq,Lambda,A); 
    cout << "setup stiffness matrix done" << endl;
    Vector<double> F;
    setup_righthand_side(eq, Lambda, F);
    
    Vector<double> x(Lambda.size());
    x =0;
    unsigned int iterations;
    const int maxiterations = 99;
    //const double omega = 2.0 / (eq.norm_A() + 1.0/eq.norm_Ainv());
    //cout << omega << endl;
    
    CG(A,F,x,1e-8,maxiterations,iterations); 
    cout << "iterative solution computed" << endl;

    //plot solution
    cout << "plotting solution" << endl;
    InfiniteVector<double,Index> u,v;
    unsigned int i2 = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i2){
        u.set_coefficient(*it, x[i2]);
    }
    u.scale(&eq, -1);
    SampledMapping<2> sm1(evaluate(eq.frame(), u , true, 6));
    std::ofstream stream1("solution_na.m");
    sm1.matlab_output(stream1);
    stream1.close();
    cout << "solution plotted" << endl;
    

    //plot coefficients
    u.scale(&eq,1);
    std::ofstream coeff_stream2;
    coeff_stream2.open("coefficients_na.m");
    MultiIndex<int,dim> pstart;//pstart[0]=0;pstart[1]=0;
    MultiIndex<int,dim> jstart;//jstart[0]=3;jstart[1]=3;
    MultiIndex<int,dim> estart;//estart[0]=0;estart[1]=0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it){
        Index ind = *it;
        if(!(ind.j()==jstart && ind.e()==estart && ind.p()==pstart)){
            cout << ind.p()[0]<<ind.p()[1]<<ind.j()[0]-1+ind.e()[0]<<ind.j()[1]-1+ind.e()[1] << endl;
            jstart=ind.j();
            estart=ind.e();
            pstart=ind.p();
            plot_indices_tframe2(&eq.frame(), u,  coeff_stream2, ind.p(), ind.j(), ind.e(), "(flipud(gray))", false, true, -5);
            coeff_stream2<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<ind.p()[0]<<","<<ind.p()[1]<<","<<ind.j()[0]-1+ind.e()[0]<<","<<ind.j()[1]-1+ind.e()[1]<<"))"<<endl;
        }    
    }    
    coeff_stream2.close();
    cout << "coefficients plotted" << endl;
    
#if 0    


    //alternative plot solutionkeding
    cout << "plotting solution 2" << endl;
    int N=100;
    double h = 1e-2;
    Matrix<double> u_values(N+1,N+1);
    for (unsigned int i = 0; i <= N; i++) {
        //cout << i << endl;
        for(unsigned int j=0; j<=N; j++) {
            Point<2> point(i*h,j*h);
            int id = 0;
            for (set<Index>::const_iterator it(Lambda.begin()); it != Lambda.end(); ++it, ++id) 
                u_values(i,j) += x[id] * eq.basis().evaluate(0, *it, point)*1./eq.D(*it);
        }
    }
    u_values.matlab_output("solution_na_alt","solution_na_alt",0);      //richtiger output für quarklets
    cout << "solution 2 plotted" << endl;
#endif
#endif
#ifdef ADAPTIVE
#ifdef BASIS
    cout << "hier1" << endl;
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > cproblem1(&eq);
#endif
#ifdef FRAME
    CachedQuarkletTProblem<TensorFrameEquation<Frame1d,dim,Frame> > cproblem1(&eq);
    
    Frame1d frame1d;
//    cout << "integral: " << eq.integrate(frame1d, 0,5,1,31,0,5,1,31,0) << endl;
    set<Index> Lambda;
  for (int i=0; i<frame.degrees_of_freedom();i++) {
    Lambda.insert(*frame.get_quarklet(i));
        cout << *frame.get_quarklet(i) << endl;
  }
    
    cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(cproblem1, Lambda, rh);
//  cout << rh << endl;
  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;
  
  clock_t tstart, tend;
  double time;
  tstart = clock();

//    WaveletTL::setup_stiffness_matrix(cproblem1, Lambda, stiff, false);
//    WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  WaveletTL::setup_stiffness_matrix(eq, Lambda, stiff, false);
//  WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff);
  

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  stiff.matlab_output("stiff_frame_2D_out", "stiff",1); 
  
  abort();
#endif
//    Index mu=eq.frame().first_generator();
//    Index lambda=eq.frame().first_generator();
    
    
#if 1
// cout << "NORMTEST" << endl;   
//    cout.precision(20);
//    cout << eq.a(lambda,mu) << endl;;
//    cout << cproblem1.a(lambda,mu) << endl;;
    
//    cproblem1.normtest(1,0);
    
    double time2=0.0, tstart2; 
    
    // time measurment variables
 
//     tstart1 = clock();              // start 
//
//    
//
//    
//    double time1=0.0, tstart1; 
//    cout << "Uncached: " << endl;
//    cout << "norm_A: " << eq.norm_A() <<endl;
//    cout << "norm_Ainv: " << eq.norm_Ainv() <<endl;
//    time1 += (clock() - tstart1)/CLOCKS_PER_SEC;
//    cout << "time = " << time1 << " sec." << endl;
    
    
#endif
    
    InfiniteVector<double, Index> F_eta; 
    cproblem1.RHS(1e-6, F_eta);
////    cout << F_eta << endl;
    const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
//    cout << "TEST: " << l2_norm(F_eta) << ", " << cproblem1.norm_Ainv() << endl;
    double epsilon = 1e-3;
    InfiniteVector<double, Index> u_epsilon, v, test, test2;
    
    //cout << eq.s_star() << endl;
    //cout << cproblem1.s_star() << endl;
    //cout << eq.space_dimension << endl;
    //cout << cproblem1.space_dimension << endl;
//    //CDD1_SOLVE(cproblem1, epsilon, u_epsilon, 2*jmax1, tensor_simple);
    
    tstart2 = clock();
    
#ifdef FRAME

    CDD2_QUARKLET_SOLVE(cproblem1, nu, epsilon, u_epsilon, jmax, tensor_simple, pmax, a, b); 
    
    
    
#if 0
    
    typedef typename Index::level_type level_type;
    typedef typename Index::polynomial_type polynomial_type;
    level_type level;
    polynomial_type p;
    level[0]=4, level[1]=3;
    p[0]=1,p[1]=0;
    cout << "Nummer: " << p.number() << endl;
    Index lambda = cproblem1.frame().first_quarklet(level, p);
    Index mu = cproblem1.frame().first_quarklet(level);
    cout << mu << endl;
    ++lambda;
    test[lambda]=1;
    cout << test << endl;
    cout << "Integral: " << cproblem1.a(lambda,lambda) << endl;
    cout << "Preconditioner: " << cproblem1.D(lambda) << ", " << cproblem1.D(mu) <<endl;
    APPLY_QUARKLET(cproblem1, test, 1e-10, test2, jmax, tensor_simple, pmax, a, b);
    cout << test2 << endl;
    
//    cproblem1.a(lambda,lambda);
//    
//    Vector<double> w(cproblem1.frame().degrees_of_freedom());
//    
//    cproblem1.add_ball(lambda,w,4,1,14,tensor_simple,true,7,a,b);
#endif

#else
    CDD2_SOLVE(cproblem1, nu, epsilon, u_epsilon, jmax, tensor_simple);
#endif
    time2 += (clock() - tstart2)/CLOCKS_PER_SEC;
    cout << "time = " << time2 << " sec." << endl;
    
    cout << u_epsilon << endl;
    
//    
//    //APPLY(cproblem1, u_epsilon, 1e-3, v, 2*jmax1, tensor_simple);
//    //APPLY_TEST(cproblem1, v, 10^-3, u_epsilon, 8, tensor_simple);
//    //CompressionStrategy strategy1=tensor_simple;
//    //Vector<double> w;
//    //add_compressed_column(cproblem1, 1, eq.basis().first_generator(), 1, w, jmax, strategy1, false, pmax, a, b);
//    
 
#ifdef FRAME
    //plot solution
    cout << "plotting solution" << endl;
    u_epsilon.scale(&cproblem1, -1);
    SampledMapping<2> sm1(evaluate(cproblem1.frame(), u_epsilon , true, 6));
    std::ofstream stream1("solution_ad.m");
    sm1.matlab_output(stream1);
    stream1.close();
    cout << "solution plotted" << endl;
#endif
    
#ifdef BASIS
    
    //plot solution
    cout << "plotting solution" << endl;
    u_epsilon.scale(&eq, -1);
    SampledMapping<2> sm1(evaluate(eq.basis(), u_epsilon , true, 6));
    std::ofstream stream1("solution_ad.m");
    sm1.matlab_output(stream1);
    stream1.close();
    cout << "solution plotted" << endl;
#endif
//    
//    //plot coefficients
//    cout << "plotting coefficients" << endl;
//    std::ofstream coeff_stream1;
//    coeff_stream1.open("coefficients_ad.m");
//    coeff_stream1 << "figure;" << endl;
//    plot_indices_tbasis(&eq.basis(), u_epsilon, 4, coeff_stream1,"(flipud(gray))", false, true, -8);
//    coeff_stream1 << "title('solution coefficients') " << endl;
//    coeff_stream1.close();
//    cout << "coefficients plotted" << endl;
#endif
    
#endif   
    return 0;

}


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE
#define BASIS

#define DYADIC
#define JMAX 8
#define PMAX 0
#define TWO_D

#define PRIMALORDER 3
#define DUALORDER   3

#define _WAVELETTL_USE_TBASIS 1
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0

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
#include <interval/pq_frame.h>
#include <cube/tframe.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <cube/tbasis_indexplot.h>
#include <adaptive/steepest_descent_ks.h>

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
    case 4:
      return 4*(exp(5*p[0])-1)/(exp(5)-1)*(1-(exp(5*p[0])-1)/(exp(5)-1))*4*(exp(5*p[1])-1)/(exp(5)-1)*(1-(exp(5*p[1])-1)/(exp(5)-1));
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
      //return ( (p[0] < 0.5 && p[1] < 0.5) ? 0: 2*(p[0]*(1-p[0])+p[1]*(1-p[1])) );  
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
    case 4:
      return 4*(4*5*5*exp(2*5*p[0])-5*5*(exp(5)+1)*exp(5*p[0]))/(pow(exp(5)-1,2))*4*(exp(5*p[1])-1)/(exp(5)-1)*(1-(exp(5*p[1])-1)/(exp(5)-1))+  4*(4*5*5*exp(2*5*p[1])-5*5*(exp(5)+1)*exp(5*p[1]))/(pow(exp(5)-1,2))*4*(exp(5*p[0])-1)/(exp(5)-1)*(1-(exp(5*p[0])-1)/(exp(5)-1));
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
    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
    const unsigned int dim = 2; 
    const int jmax=JMAX;
    
    
    typedef PBasis<d,dT> Basis1d;
//    typedef Basis1d::Index Index1d;
    typedef TensorBasis<Basis1d,dim> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;
    //typedef PQFrame<d,dT> Frame1d;
    //typedef TensorFrame<Frame1d,dim> Frame;
    
    
    
    
    FixedArray1D<int,2*dim> s;          //set order of boundary conditions
    s[0]=0, s[1]=0; s[2]=0, s[3]=0; 
    //cout << s << endl;
    FixedArray1D<bool, 2*dim> bc;
    bc[0]=true, bc[1]=true, bc[2]=true, bc[3]=true;
    
    Basis1d basis1d(true,true);
    basis1d.set_jmax(jmax);
    
    Basis basis(bc); 
    basis.set_jmax(jmax);
    
    mySolution<4> solution1;
    myRHS<4> rhs1;
    //Function2 function2;
#if 1
    //plot exact solution
    Point<2> p0(0,1);
    Point<2> p1(1,0);
    Grid<2> mygrid(p0,p1,100);
    SampledMapping<2> smuexact(mygrid, solution1);
    std::ofstream streamuexact("uexact.m");
    smuexact.matlab_output(streamuexact);
    streamuexact.close();
    cout << "uxact plotted" << endl;
  
    //plot right hand side
    //Point<2> p0(0,1);
    //Point<2> p1(1,0);
    //Grid<2> mygrid(p0,p1,100);
    SampledMapping<2> smrhs(mygrid, rhs1);
    
    std::ofstream streamrhs("rhs.m");
    smrhs.matlab_output(streamrhs);
    streamrhs.close();
    cout << "rhs plotted" << endl;
#endif  
    PoissonBVP<dim> poisson1(&rhs1);
    //cout << "hallo" << endl;
    TensorEquation<Basis1d,dim,Basis> eq(&poisson1, bc);
    //TensorEquation<Frame1d,dim,Frame> feq(&poisson1, bc);
    cout << eq.basis().degrees_of_freedom() << endl;
    //eq.basis().set_jmax(10);
    eq.set_jmax(jmax);                               //wichtig für cached problem
    cout << eq.basis().degrees_of_freedom() << endl;
    
#if 1
    //test support and output
    Index ind1=basis.first_generator();
    cout << ind1 << endl;
    
    
#endif
    
#if 0
    //plot stiffness matrix index version 
    int n=eq.basis().degrees_of_freedom();
    Matrix<double> L(n,n);
    int i=0;
    for (Index ind1 = eq.basis().first_generator(), itend = eq.basis().last_wavelet(2*jmax1); ind1 <= itend; ++ind1, i++){
        int j=0;
        for (Index ind2 = eq.basis().first_generator(), itend = eq.basis().last_wavelet(2*jmax1); ind2 <= itend; ++ind2, j++){
            //Index ind1=*it1;
            //Index ind2=*it2;
            L(i,j)=eq.a(ind2,ind1)/(eq.D(ind1)*eq.D(ind2));
        }     
    }
    L.matlab_output("L", "L",0);
    cout << "plot stiffness matrix done" << endl;
#endif
  
#ifdef NONADAPTIVE
    set<Index> Lambda;  //setup index set
    int zaehler=0;
    for (Index lambda=eq.basis().first_generator(), itend=eq.basis().last_wavelet(6); lambda <= itend; ++lambda)
    {
        Lambda.insert(lambda);
        //cout << lambda << endl;
        zaehler++;
    }
    cout << zaehler << endl;
    
    SparseMatrix<double> A;
    setup_stiffness_matrix(eq,Lambda,A); 
    cout << "setup stiffness matrix done" << endl;
    Vector<double> F;
    setup_righthand_side(eq, Lambda, F);
    
    Vector<double> x(Lambda.size());
    x =0;
    unsigned int iterations;
    const int maxiterations = 99;
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
    SampledMapping<2> sm1(evaluate(eq.basis(), u , true, 6));
    std::ofstream stream1("solution_na.m");
    sm1.matlab_output(stream1);
    stream1.close();
    cout << "solution plotted" << endl;
    
    //alternative plot solution
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
#ifdef ADAPTIVE
//    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > cproblem1(&eq);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > cproblem1(&eq,20.,10.);
    
    set<Index> Lambda;
  for (int i=0; i<basis.degrees_of_freedom();i++) {
    Lambda.insert(*basis.get_wavelet(i));
        cout << *basis.get_wavelet(i) << endl;
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

  stiff.matlab_output("stiff_2D_out", "stiff",1); 
  
  abort();
    
    
    
    
    
    
    
    
    
    InfiniteVector<double, Index> F_eta;
    cproblem1.RHS(1e-6, F_eta);
    const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
    double epsilon = 1e-4;
    InfiniteVector<double, Index> u_epsilon, v;
    
    //cout << eq.s_star() << endl;
    //cout << cproblem1.s_star() << endl;
    //cout << eq.space_dimension << endl;
    //cout << cproblem1.space_dimension << endl;
    //CDD1_SOLVE(cproblem1, epsilon, u_epsilon, 2*jmax1, tensor_simple);
    CDD2_SOLVE(cproblem1, nu, epsilon, u_epsilon, jmax, tensor_simple);
//    steepest_descent_ks_SOLVE(cproblem1, epsilon, u_epsilon);
    
    //APPLY(cproblem1, u_epsilon, 1e-3, v, 2*jmax1, tensor_simple);
    //APPLY_TEST(cproblem1, v, 10^-3, u_epsilon, 8, tensor_simple);
    //CompressionStrategy strategy1=tensor_simple;
    //Vector<double> w;
    //add_compressed_column(cproblem1, 1, eq.basis().first_generator(), 1, w, 2*jmax1, strategy1, false, pmax, a, b);
    
    //plot solution
    cout << "plotting solution" << endl;
    u_epsilon.scale(&eq, -1);
    SampledMapping<2> sm1(evaluate(eq.basis(), u_epsilon , true, 6));
    std::ofstream stream1("solution_ad.m");
    sm1.matlab_output(stream1);
    stream1.close();
    u_epsilon.scale(&eq,1);
    cout << "solution plotted" << endl;
    
    //plot coefficients
    cout << "plotting coefficients" << endl;
    std::ofstream coeff_stream1;
    coeff_stream1.open("coefficients_ad.m");
    coeff_stream1 << "figure;" << endl;
    plot_indices_tbasis(&eq.basis(), u_epsilon, 4, coeff_stream1,"(flipud(gray))", false, true, -8);
    coeff_stream1 << "title('solution coefficients') " << endl;
    coeff_stream1.close();
    cout << "coefficients plotted" << endl;
    
    
    //new coefficients plot
    std::ofstream coeff_stream2;
    coeff_stream2.open("coefficients_ad2.m");
    //coeff_stream2 << "figure;" << endl;
    MultiIndex<int,dim> jstart;// start=basis1.j0();
    MultiIndex<int,dim> estart;
    for (Index lambda=eq.basis().first_generator(), itend=eq.basis().last_wavelet(jmax); lambda <= itend; ++lambda){
        if(!(lambda.j()==jstart && lambda.e()==estart)){
            cout << lambda.j()[0]-1+lambda.e()[0]<<lambda.j()[1]-1+lambda.e()[1] << endl;
            jstart=lambda.j();
            estart=lambda.e();
            plot_indices_tbasis2(&eq.basis(), u_epsilon, coeff_stream2, lambda.j(), lambda.e(),"(flipud(gray))", false, true, -5);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            coeff_stream2<<"print('-djpg',sprintf('coeffs%i%i.jpg',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
    coeff_stream2.close();
    
    //plot grid
    cout << "plotting grid" << endl;
    std::ofstream grid_stream1;
    grid_stream1.open("grid_ad.m");
    grid_stream1 << "figure;" << endl;
    grid_stream1 << "box on;" << endl;
    for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it){
         //u_values_ad[i] += *it * basisframe.evaluate(0, it.index(), point)*1./cproblem1.D(it.index());
        
        //cout << ind << " : " << *it << endl;
        if(abs(*it)>1e-16){
            Index ind=it.index();
            Support supp;
            basis.support(ind, supp);
            grid_stream1 << "patch([" << (double)supp.a[0]/(1<<supp.j[0])<<"," << (double)supp.b[0]/(1<<supp.j[0]) << "," << (double)supp.b[0]/(1<<supp.j[0]) << "," << (double)supp.a[0]/(1<<supp.j[0]) << "],[" << (double)supp.a[1]/(1<<supp.j[1])<<"," << (double)supp.a[1]/(1<<supp.j[1]) << "," << (double)supp.b[1]/(1<<supp.j[1]) << "," << (double)supp.b[1]/(1<<supp.j[1]) << "],'w','edgecolor','k')" <<endl;
        }
        
    }
    grid_stream1.close();
    cout << "grid plotted" << endl;
    
    //all coefficients output
    std::ofstream stream2;
    stream2.open("coeffs_support.m");
    stream2<<"A=["<<endl;
    for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it){
         //u_values_ad[i] += *it * basisframe.evaluate(0, it.index(), point)*1./cproblem1.D(it.index());
        //cout << ind << " : " << *it << endl;
        if(abs(*it)>1e-16){
            Index ind=it.index();
            Support supp;
            basis.support(ind, supp);
            //stream2<<ind.j()[0]<<","<<ind.j()[1]<<","<<ind.e()[0]<<","<<ind.e()[1]<<","<<ind.k()[0]<<","<<ind.k()[1]<<","<<(double)supp.a[0]/(1<<supp.j[0])<<","<<(double)supp.b[0]/(1<<supp.j[0])<<","<<(double)supp.a[1]/(1<<supp.j[1])<<","<<(double)supp.b[1]/(1<<supp.j[1])<<","<<*it<<";"<< endl;
            stream2<<ind.j()[0]-1+ind.e()[0]<<","<<ind.j()[1]-1+ind.e()[1]<<","<<ind.k()[0]<<","<<ind.k()[1]<<","<<(double)supp.a[0]/(1<<supp.j[0])<<","<<(double)supp.b[0]/(1<<supp.j[0])<<","<<(double)supp.a[1]/(1<<supp.j[1])<<","<<(double)supp.b[1]/(1<<supp.j[1])<<","<<*it<<";"<< endl;
        }
        
    }
    stream2<<"];"<<endl;
    stream2.close();
    cout << "fertig" << endl;
    
#endif
    
    
    return 0;

}


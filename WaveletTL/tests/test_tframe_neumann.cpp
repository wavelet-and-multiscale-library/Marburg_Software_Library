/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#define PARALLEL 0
#undef NONADAPTIVE
#define ADAPTIVE

#undef BASIS
#define FRAME

#undef DIRICHLET
#define NEUMANN
//#ifdef NEUMANN
////#define NONZERONEUMANN
//#endif

#define DYADIC
#undef DYPLUSEN

#define JMAX 6
#define PMAX 1

#ifdef FRAME
#define _WAVELETTL_USE_TFRAME 1
#else
#define _WAVELETTL_USE_TBASIS 1
#endif
#define _DIM 2
#define TWO_D
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 2

#include <iostream>

//#include <interval/p_basis.h>
#include <interval/pq_frame.h>
//#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
//#include <cube/tbasis.h>
//#include <cube/tbasis_index.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
//#include <galerkin/tbasis_equation.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
//#include <cube/tbasis_evaluate.h>
#include <cube/tframe_evaluate.h>
#include <cube/tframe.h>
#include <cube/tframe_evaluate.h>
#include <cube/tframe_indexplot.h>
#include <cube/tframe_index.h>
#include <cube/tframe_support.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
//#include <cube/tbasis_indexplot.h>
#include <galerkin/tframe_equation.h>
#include <galerkin/cached_quarklet_tproblem.h>
//#include "TestFunctions2d.h"
#include <adaptive/steepest_descent_ks.h>
#include <adaptive/duv.h>



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
      return sin(M_PI*p[0])*sin(M_PI*p[0])*sin(M_PI*p[1])*sin(M_PI*p[1]);               
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
    case 4:
      return
           -2*M_PI*M_PI*(cos(M_PI*p[0])*cos(M_PI*p[0])-sin(M_PI*p[0])*sin(M_PI*p[0]))*sin(M_PI*p[1])*sin(M_PI*p[1])
           -2*M_PI*M_PI*(cos(M_PI*p[1])*cos(M_PI*p[1])-sin(M_PI*p[1])*sin(M_PI*p[1]))*sin(M_PI*p[0])*sin(M_PI*p[0]);
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
                    return "u(x,y)=sin^2(pi*x)*sin^2(pi*y), zero-Neumann boundary conditions";
                    break;
                   
                default:
                    return "TestProblem: N not defined.";
                    break;
            }
        }
};

template <unsigned int N>
class myBoundaryFunction
  : public Function<2,double>
{
public:
  virtual ~myBoundaryFunction() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
        if(p[0]==1) return -p[1]*(1-p[1]);
        if(p[0]==0) return -p[1]*(1-p[1]);
        if(p[1]==1) return -p[0]*(1-p[0]);
        if(p[1]==0) return -p[0]*(1-p[0]);
//        cout<<"kein randpunkt"<<endl;
        return 0;
      break;
        case 2:
      return 0;
      break;
    case 3:
      return 0;
      break;
    case 4:
      return 0;               
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

int main()
{
    cout << "Testing tframe ..." << endl;
    const int d  = 3;
    const int dT = 3;
#ifdef ADAPTIVE
    const int a=2;
    const int b=2;
#endif
    const unsigned int dim = 2; 
    const int jmax=JMAX;
    const int pmax = PMAX;


    const unsigned int testcase1=1;
    mySolution<testcase1> uexact1;
    myRHS<testcase1> rhs1;
    myBoundaryFunction<testcase1> phi1;
    
    PoissonBVP<dim> poisson1(&rhs1);
    
#if 1 //plot rhs and exact solution


    std::ofstream osrhs("rhs.m");
    std::ofstream osuexact("uexact.m");
    std::ofstream osphi("phi.m");
        Point<2> q0(0,0);
        Point<2> q1(1,1);
        Grid<2> mygrid(q0,q1,100);
        SampledMapping<2> smrhs(mygrid, rhs1); 
        SampledMapping<2> smuexact(mygrid, uexact1);
        SampledMapping<2> smphi(mygrid,phi1);
        smrhs.matlab_output(osrhs);
        smuexact.matlab_output(osuexact);
        smphi.matlab_output(osphi);
        osrhs << "surf(x,y,z)"<<endl;
        osuexact<<"surf(x,y,z)"<<endl;
        osphi<<"surf(x,y,z)"<<endl;
    osrhs << "view(30,55);"<<endl;
    osuexact << "view(30,55);"<<endl;
    osphi<<"view(30,55)"<<endl;
    osrhs.close();
    osuexact.close();
    osphi.close();
    cout << "rhs and uexact plotted" << endl;
#endif
    
    typedef PQFrame<d,dT> Frame1d;
    typedef TensorFrame<Frame1d,dim> Frame;
    typedef Frame::Index Index;
    typedef Frame::Support Support;

#ifdef DIRICHLET    
    Frame1d frame1d(true,true);   
#endif
#ifdef NEUMANN    
    Frame1d frame1d(false,false);
#endif
    
    frame1d.set_jpmax(jmax-frame1d.j0(),pmax);
    FixedArray1D<Frame1d*,2> frames;
    frames[0]=&frame1d;
    frames[1]=&frame1d;
       
    Frame frame(frames);
    frame.set_jpmax(jmax,pmax);
    cout<<"dof: "<<frame.degrees_of_freedom()<<endl;
#if 0//index tests
    cout<<frame1d.DeltaLmin(0)<<endl;
    cout<<frame1d.DeltaRmax(3,0)<<endl;
    cout<<frame1d.Deltasize(3,0)<<endl;
    cout<<frame.first_generator(frame.j0())<<endl;
    cout<<first_generator(&frame,frame.j0(),MultiIndex<int,2>())<<endl;
    cout<<frame.last_generator(frame.j0())<<endl;
    cout<<last_generator(&frame,frame.j0(),MultiIndex<int,2>())<<endl;
    cout<<frame1d.Nablamin(0)<<endl;
    cout<<frame1d.Nablamax(3,0)<<endl;
    cout<<frame.first_quarklet(frame.j0())<<endl;
    cout<<first_quarklet(&frame,frame.j0(),MultiIndex<int,2>())<<endl;
    cout<<frame.last_quarklet(frame.j0())<<endl;
    cout<<last_quarklet(&frame,frame.j0(),MultiIndex<int,2>())<<endl;
    cout<<frame.get_first_wavelet_numbers()<<endl;
    cout<<frame.get_last_wavelet_numbers()<<endl;
#endif
#ifdef DIRICHLET
    TensorFrameEquation<Frame1d,dim,Frame> eq(&poisson1, &frame, true);
#endif    
#ifdef NEUMANN
    TensorFrameEquation<Frame1d,dim,Frame> eq(&poisson1, &frame, &phi1, true);
#endif
      
#if 0 //some frame tests
    cout<<eq.frame().degrees_of_freedom()<<endl;
    cout<<eq.frame().j0()<<endl;
    cout<<eq.frame().first_generator(eq.frame().j0())<<endl;
    cout<<eq.a(eq.frame().first_generator(eq.frame().j0()),eq.frame().first_generator(eq.frame().j0()))<<endl;
    cout<<eq.f(eq.frame().first_generator(eq.frame().j0()))<<endl;
    Index ind=eq.frame().first_generator(eq.frame().j0());
    

    Support supp;
    eq.frame().support(eq.frame().first_generator(eq.frame().j0()),supp);
    cout<<supp.b[0]<<endl;
    cout<<intersect_supports(eq.frame(),0,1,supp)<<endl;
    std::list<int> inters;
    
    MultiIndex<int,2> q;
    Index ind2=first_quarklet<Frame1d,dim>(&eq.frame(), eq.frame().j0(),q);
    cout<<"ind2:"<<ind2<<endl;
//    ++ind2;
//    cout<<"ind2:"<<ind2<<endl;
//    ind2=last_quarklet<Frame1d,dim>(&eq.frame(), eq.frame().j0(),q);
//    cout<<(ind2==ind2)<<endl;
    intersecting_quarklets<Frame1d,dim>(eq.frame(),ind2,eq.frame().j0(),inters,q);
    for(list<int>::const_iterator it = inters.begin();it!=inters.end();++it){
        cout<<*it<<endl;
    }
//    cout<<supp<<endl;
//    Index ind3=frame.get_quarklet(0);
//    cout<<ind3<<endl;
#endif
 
#if 0
    //plot one function
    Index testindex=frame.get_quarklet(90);   
//    SampledMapping<dim> evalf(evaluate(frame, testindex, true, 6));
    SampledMapping<dim> evalf(frame.evaluate(testindex,6));
    
    cout << "evaluate quarklet with index " << testindex << endl;
//    evalf=frame.evaluate(testindex,6);
    
    std::ofstream osf("tframeoutput.m");
    osf << "clf;" << endl;
    osf << "axis([0 1]);" << endl;
        evalf.matlab_output(osf);
        osf << "surf(x,y,z);" << endl;
            
    osf << "view(30,55);"<<endl;
    osf.close();
#endif
    
#if 0
    //test modified f with neumann-bc
    cout<<"testing modified f with neumann-bc"<<endl;
//    Index testindex=frame.get_quarklet(0);
    cout<<"f= "<<eq.f(testindex)<<endl;
    
//    Point<2,double> x; x[0]=0;x[1]=0.5;
    Point<2,double> z(0,0.5);
    cout<<phi1.value(z)<<endl;
    cout<<eq.phi(z)<<endl;
#endif
    

       
#ifdef NONADAPTIVE
    //setup index set
    set<Index> Lambda;  
    MultiIndex<int,dim> p;p[0]=0;p[1]=0;
    Index lambda = eq.frame().first_generator(eq.frame().j0(), p);
    int zaehler=0;
    for (int l = 0; l < eq.frame().degrees_of_freedom(); l++) {
//        if(lambda<eq.frame().last_generator(p)) 
            Lambda.insert(lambda);   //using just generators for test purpose
//        cout << lambda << ", Number: " << lambda.number() << endl;
        if(lambda==eq.frame().last_quarklet(jmax, p)){
            ++p;
            lambda=eq.frame().first_generator(eq.frame().j0(), p);
        }
        else
        ++lambda;
        ++zaehler;
        
    }
    cout<<"using "<<Lambda.size()<<" dof"<<endl;
    cout << "setup stiffness matrix" << endl;
    SparseMatrix<double> A;
    setup_stiffness_matrix(eq,Lambda,A); 
    cout << "setup stiffness matrix done" << endl;
    Vector<double> F;
    setup_righthand_side(eq, Lambda, F);
//    cout<<"NONADAPTIVE f.size: "<<F.size()<<endl;
//    cout<<"l2_norm: "<<l2_norm(F)<<endl;
//    F.matlab_output("F","F");
    
//    cout<<"normA: "<<eq.norm_A()<<endl;
//    cout<<"normAinv: "<<eq.norm_Ainv()<<endl;
    
    Vector<double> x(Lambda.size());
    x =0;
    unsigned int iterations;
    const int maxiterations = 9999;
    //const double omega = 2.0 / (eq.norm_A() + 1.0/eq.norm_Ainv());
    //cout << omega << endl;
    
//    CG(A,F,x,1e-8,maxiterations,iterations); 
    Richardson(A,F,x,0.04,1e-6,maxiterations,iterations);
    cout << "iterative solution computed with " <<iterations<<" iterations"<< endl;

    {
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
    stream1<<"surf(x,y,z)"<<endl;
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
//            cout << ind.p()[0]<<ind.p()[1]<<ind.j()[0]-1+ind.e()[0]<<ind.j()[1]-1+ind.e()[1] << endl;
            jstart=ind.j();
            estart=ind.e();
            pstart=ind.p();
            plot_indices_tframe2(&eq.frame(), u,  coeff_stream2, ind.p(), ind.j(), ind.e(), "jet", false, true, -6);
            coeff_stream2<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<ind.p()[0]<<","<<ind.p()[1]<<","<<ind.j()[0]-1+ind.e()[0]<<","<<ind.j()[1]-1+ind.e()[1]<<"))"<<endl;
        }    
    }    
    coeff_stream2.close();
    cout << "coefficients plotted" << endl;
    }
 #if 0
    //all coefficients output
    u.scale(&eq, 1);
    std::ofstream os4;
    os4.open("coeffs_support.m");
    os4<<"A=["<<endl;
    for (InfiniteVector<double, Index>::const_iterator it(u.begin()); it != u.end(); ++it){
         //u_values_ad[i] += *it * basisframe.evaluate(0, it.index(), point)*1./cproblem1.D(it.index());
        //cout << ind << " : " << *it << endl;
        //if(abs(*it)>1e-16){
            Index ind=it.index();
            //Support supp;
            //frame.support(ind, supp);
            //stream2<<ind.j()[0]<<","<<ind.j()[1]<<","<<ind.e()[0]<<","<<ind.e()[1]<<","<<ind.k()[0]<<","<<ind.k()[1]<<","<<(double)supp.a[0]/(1<<supp.j[0])<<","<<(double)supp.b[0]/(1<<supp.j[0])<<","<<(double)supp.a[1]/(1<<supp.j[1])<<","<<(double)supp.b[1]/(1<<supp.j[1])<<","<<*it<<";"<< endl;
            os4<<ind.p()[0]<<","<<ind.p()[1]<<","<<ind.j()[0]-1+ind.e()[0]<<","<<ind.j()[1]-1+ind.e()[1]<<","<<ind.k()[0]<<","<<ind.k()[1]<<","<<*it<<";"<< endl;
       // }
        
    }
    os4<<"];"<<endl;
    os4.close();
#endif   
#endif

#ifdef ADAPTIVE
    CachedQuarkletTProblem<TensorFrameEquation<Frame1d,dim,Frame> > cproblem1(&eq,41,27);
    cout<<"normA: "<<cproblem1.norm_A()<<endl;
    cout<<"normAinv: "<<cproblem1.norm_Ainv()<<endl;
    
    InfiniteVector<double, int> F_eta; 
    cproblem1.RHS(1e-6, F_eta);
//    cout<<"ADAPTIVE f.size: "<<F_eta.size()<<endl;
//    cout<<"l2_norm: "<<l2_norm(F_eta)<<endl;
    double epsilon = 1e-3;
    InfiniteVector<double, Index> u_epsilon;
    InfiniteVector<double, int> u_epsilon_int;
    const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
    const unsigned int maxiter = 99;
    const double omega = 0.05;
    const double residual_stop = 1e-10;
    const double shrinkage = 0;
    
    //choose one scheme
    CDD2_QUARKLET_SOLVE(cproblem1, nu, epsilon, u_epsilon_int, jmax, tensor_simple, pmax, a, b); 
//    richardson_QUARKLET_SOLVE(cproblem1,epsilon,u_epsilon_int, maxiter, tensor_simple, a, b, shrinkage, omega, residual_stop);
//    DUV_QUARKLET_SOLVE_SD(cproblem1, nu, epsilon, u_epsilon, tensor_simple, pmax, jmax, a, b);
//    steepest_descent_ks_QUARKLET_SOLVE(cproblem1, epsilon, u_epsilon_int, tensor_simple, a, b);
            
    //    cout<<u_epsilon_int.size()<<endl;
    for (typename InfiniteVector<double,int>::const_iterator it(u_epsilon_int.begin()),
 	   itend(u_epsilon_int.end()); it != itend; ++it){
        u_epsilon.set_coefficient(*(frame.get_quarklet(it.index())), *it);
    }
//    cout<<u_epsilon.size()<<endl;
    
//    //testing APPLY
//    InfiniteVector<double, Index>  v,Av;
//    for(int i=0;i<cproblem1.frame().degrees_of_freedom();i++){
//        v.set_coefficient(cproblem1.frame().get_quarklet(i),1);
//    }
//    cout<<v.size()<<endl;
//    APPLY_QUARKLET(cproblem1, v, 1e-3, Av, jmax, tensor_simple, pmax, a,b);
//    cout<<Av.size()<<endl;
    
    {
    //plot solution
    cout << "plotting solution" << endl;
    u_epsilon.scale(&cproblem1, -1);
    SampledMapping<2> sm1(evaluate(cproblem1.frame(), u_epsilon , true, 6));
    std::ofstream stream1("solution_ad.m");
    sm1.matlab_output(stream1);
    stream1<<"surf(x,y,z)"<<endl;
    stream1.close();
    cout << "solution plotted" << endl;
    
    //new coefficients plot
    u_epsilon.scale(&cproblem1, 1);
    std::ofstream coeff_stream;
    coeff_stream.open("coefficients_ad.m");
    MultiIndex<int,dim> pstart;
    MultiIndex<int,dim> jstart;// start=basis1.j0();
    MultiIndex<int,dim> estart;
    for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it){
        Index lambda=it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart)){
            jstart=lambda.j();
            estart=lambda.e();
            plot_indices_tframe2(&frame, u_epsilon, coeff_stream, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            coeff_stream<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
    coeff_stream.close();
    cout << "coefficients plotted"<<endl;
//    cout<<u_epsilon.size()<<endl;
    }
    
#endif

   

    
    cout<<"fertig"<<endl;
    return 0;

}


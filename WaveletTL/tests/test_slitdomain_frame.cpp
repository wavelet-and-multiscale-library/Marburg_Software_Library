/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#define POISSON
#undef GRAMIAN


#undef DYADIC
#undef TRIVIAL
#define ENERGY
#undef DYPLUSEN

#define NONADAPTIVE
#undef ADAPTIVE

#ifdef ADAPTIVE
#undef SD
#undef CDD2
#define RICHARDSON
#endif

#define PARALLEL 0

#define FRAME
//#define _WAVELETTL_USE_TBASIS 1
#define _WAVELETTL_USE_TFRAME 1
#define _DIM 2
#define JMAX 7
#define PMAX 0
#define TWO_D

#define PRIMALORDER 3
#define DUALORDER   3

#include <iostream>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>
#include <slitdomain/slitdomain_frame_index.h>
#include <slitdomain/slitdomain_frame.h>
#include <slitdomain/slitdomain_frame_evaluate.h>
#include <slitdomain/slitdomain_frame_support.h>
#include <slitdomain/slitdomain_frame_indexplot.h>
#include <galerkin/slitdomain_frame_equation.h>
#include <galerkin/cached_quarklet_slitdomain_problem.h>

#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <adaptive/cdd2.h>
#include <adaptive/duv.h>
#include <adaptive/steepest_descent_ks.h>


using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

//test problems
class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    CornerSingularityRHS csrhs(Point<2>(0,0), 0.5, 2.0);
    return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])+5*csrhs.value(p);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    CornerSingularity cs(Point<2>(0,0), 0.5, 2.0);
    return sin(M_PI*p[0])*sin(M_PI*p[1])+5*cs.value(p);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};
class mySum
  : public Function<2,double>
{
public:
  virtual ~mySum() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return (1+2*M_PI*M_PI)*sin(M_PI*p[0])*sin(M_PI*p[1]);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main(){
    cout << "testing slitdomain quarklet frame" << endl;
    
    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
    const int dim = _DIM;
    const int jmax=JMAX;
    const int pmax=PMAX;
    
    //definitions
    typedef PQFrame<d,dT> Frame1d;
    typedef SlitDomainFrame<Frame1d> Frame;
    typedef Frame::Index Index;

//    typedef Frame::Support Support;

    
    //instances of 1d frames with different boundary conditions
    Frame1d frame1d(false,false);
    frame1d.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_11(true,true);
    frame1d_11.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_01(false,true);
    frame1d_01.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_10(true,false);
    frame1d_10.set_jpmax(jmax-frame1d.j0(),pmax);
    
    //call frame constructor
    Frame frame(&frame1d, &frame1d_11, &frame1d_01, &frame1d_10);
    //set jpmax
    frame.set_jpmax(jmax,pmax);
    
    //set test problem
    mySolution uexact1;
    myRHS rhs1;
    PoissonBVP<dim> poisson1(&rhs1);
#ifdef POISSON
    SlitDomainFrameEquation<Frame1d,Frame> eq(&poisson1, &frame, true);
#endif
#ifdef GRAMIAN
    LDomainFrameGramian<Frame1d,Frame> eq(&gramian1, &frame, true);
#endif
    
#if 1 //plot rhs and exact solution
    Array1D<Point<dim,int> > corners;
    corners.resize(4);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;
    corners[3][0]=0;
    corners[3][0]=0;

    std::ofstream osrhs("rhs.m");
    std::ofstream osuexact("uexact.m");
    for(int i=0;i<4;i++){
        Point<2> q0(corners[i][0],corners[i][1]);
        Point<2> q1(corners[i][0]+1,corners[i][1]+1);
        Grid<2> mygrid(q0,q1,100);
        SampledMapping<2> smrhs(mygrid, rhs1); 
        SampledMapping<2> smuexact(mygrid, uexact1);
        smrhs.matlab_output(osrhs);
        smuexact.matlab_output(osuexact);
        osrhs << "surf(x,y,z,'LineStyle','none')"<<endl;
        osuexact<<"surf(x,y,z,'LineStyle','none')"<<endl;
        osrhs << "hold on;" << endl;
        osuexact<<"hold on;"<<endl;
    }
    osrhs << "view(30,55);"<<endl;
    osrhs << "hold off;" << endl;
    osuexact << "view(30,55);"<<endl;
    osuexact<<"hold off;"<<endl;
    osrhs.close();
    osuexact.close();
    cout << "rhs and uexact plotted" << endl;
#endif  
    
#ifdef NONADAPTIVE
    //setup Index set   
    set<Index> Lambda;
    for (int i=0; i<frame.degrees_of_freedom();i++) {
        Index lambda = *(frame.get_quarklet(i));
        if(multi_degree(lambda.e())<3 && lambda.patch()!=-1){ 
            Lambda.insert(*(frame.get_quarklet(i)));
        }    
//    cout << *(frame.get_quarklet(i)) << endl;
    }
    
    //setup stiffness matrix and rhs
    SparseMatrix<double> A;
    setup_stiffness_matrix(eq,Lambda,A);
    Vector<double> F;
    setup_righthand_side(eq, Lambda, F); 
    A.compress(1e-14);
    F.compress(1e-14);
    A.matlab_output("A","A",0,A.row_dimension()-1,A.column_dimension()-1 );
    F.matlab_output("F","F");
    Vector<double> x(Lambda.size());
    
    
    //richardson iteration
    unsigned int iterations;
    const int maxiterations = 9999;
    cout << eq.norm_A() << endl;
    cout << eq.norm_Ainv() << endl;
    const double omega = 2.0 / (eq.norm_A() + 1.0/eq.norm_Ainv());
    cout << "omega: "<<omega<<endl;
//    const double omega2 = 0.04;
    Richardson(A,F,x,omega,1e-3,maxiterations,iterations);
//    CG(A,F,x,1e-3,maxiterations,iterations);
    cout << "iterations:" << iterations << endl;
//    cout<<x<<endl;
    unsigned int nontrivial = 0; 
    for(unsigned int i=0; i<x.size();i++){
        if(x[i]!=0)
        nontrivial++;
    }
    cout << "Number of nontrivial entries: " << nontrivial << endl;
    
    InfiniteVector<double,Index> u;
    unsigned int i2 = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i2){
        if(x[i2]!=0)
        u.set_coefficient(*it, x[i2]);
    }
    
    //plot solution
    //u.COARSE(1e-6,v);
    u.scale(&eq, -1);
    Array1D<SampledMapping<dim> > eval(3);
    eval=eq.frame().evaluate(u,4);
    std::ofstream os2("solution_na.m");
    os2 << "figure;" << endl;
    for(int i=0;i<4;i++){
        eval[i].matlab_output(os2);        
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "title('nonadaptive solution" << ", dof=" << nontrivial << "')" << endl;
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close();
    
    //new coefficients plot
    u.scale(&eq, 1);
    std::ofstream coeff_stream;
    coeff_stream.open("coefficients_na.m");
    //coeff_stream2 << "figure;" << endl;
    MultiIndex<int,dim> pstart;
    MultiIndex<int,dim> jstart;// start=basis1.j0();
    MultiIndex<int,dim> estart;
    for (InfiniteVector<double, Index>::const_iterator it(u.begin()); it != u.end(); ++it){
        Index lambda=it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart)){
            //cout <<lambda.p()[0]<<lambda.p()[1]<< lambda.j()[0]-1+lambda.e()[0]<<lambda.j()[1]-1+lambda.e()[1] << endl;
            jstart=lambda.j();
            estart=lambda.e();
            plot_indices_slitdomain(&frame, u, coeff_stream, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            coeff_stream<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
    coeff_stream.close();
    cout << "coefficients plotted"<<endl;
    
//    //new coefficients plot
//    std::ofstream coeff_stream("coefficients_na.m");
//    MultiIndex<int,dim> pstart;
//    MultiIndex<int,dim> jstart;// start=basis1.j0();
//    MultiIndex<int,dim> estart;
//    for (InfiniteVector<double, Index>::const_iterator it(u.begin()); it != u.end(); ++it){
//        Index lambda=it.index();
//        if(!(lambda.j()==jstart && lambda.e()==estart)){
//            //cout <<lambda.p()[0]<<lambda.p()[1]<< lambda.j()[0]-1+lambda.e()[0]<<lambda.j()[1]-1+lambda.e()[1] << endl;
//            jstart=lambda.j();
//            estart=lambda.e();
//            plot_indices_slitdomain(&frame, u, coeff_stream, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
//            //coeff_stream2 << "title('solution coefficients') " << endl;
//            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
//            coeff_stream<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
//        }
//    }
//    coeff_stream.close();
//    cout << "coefficients plotted"<<endl;
#endif
    
#ifdef ADAPTIVE
    //setup Index set   
    set<Index> Lambda;
    for (int i=0; i<frame.degrees_of_freedom();i++) {
        Index lambda = *(frame.get_quarklet(i));
        if(multi_degree(lambda.e())<3 && lambda.patch()!=-1){ 
            Lambda.insert(*(frame.get_quarklet(i)));
        }    
//    cout << *(frame.get_quarklet(i)) << endl;
    }
    
    //setup problem
    CachedQuarkletSlitDomainProblem<SlitDomainFrameEquation<Frame1d,Frame> > cproblem1(&eq);
    cout<<"normA: "<<cproblem1.norm_A()<<endl;
    cout<<"normAinv: "<<cproblem1.norm_Ainv()<<endl;
    
    
    InfiniteVector<double, Index> F_eta;
    cproblem1.RHS(1e-6, F_eta);
    const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
    cout<<"nu: "<<nu<<endl;
    
    InfiniteVector<double, Index> u_epsilon;
    InfiniteVector<double, int> u_epsilon_int;
    
    const double a=2;
    const double b=2;
    double epsilon = 1e-20;
    
#ifdef CDD2
    CDD2_QUARKLET_SOLVE(cproblem1, nu, epsilon, u_epsilon_int, jmax, tensor_simple, pmax, a, b);
#endif
#ifdef RICHARDSON
    const unsigned int maxiter = 500;
    richardson_QUARKLET_SOLVE(cproblem1,epsilon,u_epsilon_int, maxiter, tensor_simple, a, b, 0, 0.2);
#endif
    for (typename InfiniteVector<double,int>::const_iterator it(u_epsilon_int.begin()),
 	   itend(u_epsilon_int.end()); it != itend; ++it){
        u_epsilon.set_coefficient(*(frame.get_quarklet(it.index())), *it);
    }
    
    //plot solution
    //u.COARSE(1e-6,v);
    u_epsilon.scale(&cproblem1, -1);
    unsigned int i2 = 0;
    Vector<double> x(Lambda.size());
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i2){
        x[i2]=u_epsilon.get_coefficient(*it);        
    }
    x.matlab_output("x","x");
    Array1D<SampledMapping<dim> > eval(4);
    eval=cproblem1.frame().evaluate(u_epsilon,6);
    std::ofstream os2("solution_ad.m");
    os2 << "figure;"<< endl;
    for(int i=0;i<4;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z,'LineStyle','none');" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "title('slitdomain Poisson Equation: adaptive solution to test problem ("
            << "CDD2" << "), " << "pmax= " << pmax << ", jmax= " << jmax << ", d= " << d << ", dT= " << dT << "');" << endl;
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close();
    
    //new coefficients plot
    u_epsilon.scale(&cproblem1, 1);
    std::ofstream coeff_stream;
    coeff_stream.open("coefficients_ad.m");
    //coeff_stream2 << "figure;" << endl;
    MultiIndex<int,dim> pstart;
    MultiIndex<int,dim> jstart;// start=basis1.j0();
    MultiIndex<int,dim> estart;
    for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it){
        Index lambda=it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart)){
            //cout <<lambda.p()[0]<<lambda.p()[1]<< lambda.j()[0]-1+lambda.e()[0]<<lambda.j()[1]-1+lambda.e()[1] << endl;
            jstart=lambda.j();
            estart=lambda.e();
            plot_indices_slitdomain(&frame, u_epsilon, coeff_stream, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            coeff_stream<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
    coeff_stream.close();
    cout << "coefficients plotted"<<endl;
#endif
    

 
#if 0 //test index
    MultiIndex<int,dim> p;p[0]=0;p[1]=0;
    cout<<"coarsest level j="<<frame.j0()<<endl;
    Index ind1=frame.first_generator(frame.j0(),p);
    cout<<"first generator index="<<ind1<<endl;
    ++ind1;
    cout<<"second generator index="<<ind1<<endl;
    cout<<"last generator index="<<frame.last_generator(frame.j0(),p)<<endl;
    ind1=frame.first_quarklet(frame.j0(),p);
    cout<<"first quarklet index="<<ind1<<endl;
    ++ind1;
    cout<<"second quarklet index"<<ind1<<endl;
    cout<<"last quarklet index for p=0:"<<frame.last_quarklet(jmax,p)<<endl;
//    cout<<frame.frames(6,1)->Deltasize(3,0)<<endl;
//    cout<<*(frame.get_quarklet(1070))<<endl;
#endif 
    
#if 0 //test rhs
    const int resolution2=7;
    Index testindex3=frame.get_quarklet(1063);
    cout<<testindex3<<endl;
    Array1D<SampledMapping<dim> > eval3(4);
    eval3=frame.evaluate(testindex3, resolution2);
    Array1D<Point<dim,int> > corners;
    corners.resize(4);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;
    corners[3][0]=0;
    corners[3][1]=0;

    std::ofstream os3("testf.m");
    os3<<"Iexakt="<<eq.f(testindex3)<<endl;
    cout<<"f="<<eq.f(testindex3)<<endl;
    os3<<"I=0;"<<endl;
    for(int i=0;i<4;i++){
        eval3[i].matlab_output(os3);
        os3<<"zalt=z;"<<endl;
        Point<2> q0(corners[i][0],corners[i][1]);
        Point<2> q1(corners[i][0]+1,corners[i][1]+1);
        Grid<2> mygrid(q0,q1,pow(2,resolution2));
        SampledMapping<2> smf(mygrid, rhs1); 
        smf.matlab_output(os3);
        os3<<"z=z.*zalt;"<<endl;
        os3<<"xx=x(1,:);"<<endl;
        os3<<"yy=y'(1,:);"<<endl;
        os3<<"I=I+trapz(yy,trapz(xx,z,2)');"<<endl;
    }
    os3<<"I"<<endl;
    os3<<"relative_error=abs(I-Iexakt)/Iexakt"<<endl; 
    os3.close();
#endif
    
#if 0 //test bilinearform
    const int resolution=7;
    Index testindex4=eq.frame().get_quarklet(929); //1063,1064
    Index testindex5=eq.frame().get_quarklet(1057);
    cout<<testindex4<<endl;
    cout<<testindex5<<endl;
    Array1D<SampledMapping<dim> > eval1(4);
    Array1D<SampledMapping<dim> > eval2(4);
    eval1=frame.evaluate(testindex4,resolution);
    eval2=frame.evaluate(testindex5,resolution);
    std::ofstream os("testa.m");
    os<<"Iexakt="<<eq.a(testindex4,testindex5)<<endl;
    cout<<"a:"<<eq.a(testindex4,testindex5)<<endl;
    os<<"I=0;"<<endl;
    for(int i=0;i<4;i++){
        eval1[i].matlab_output(os);
        os<<"z1=z;"<<endl;
        eval2[i].matlab_output(os);
        os<<"z2=z;"<<endl;
        os<<"xx=x(1,:);"<<endl;
        os<<"yy=y'(1,:);"<<endl;
        os<<"[gx1,gy1]=gradient(z1,2^-"<<resolution<<",2^-"<<resolution<<");"<<endl;
        os<<"[gx2,gy2]=gradient(z2,2^-"<<resolution<<",2^-"<<resolution<<");"<<endl;
        os<<"L=gx1.*gx2.+gy1.*gy2;"<<endl;
        os<<"I=I+trapz(yy,trapz(xx,L,2)');"<<endl;
    }
    os<<"I"<<endl;
    os<<"relative_error=abs(I-Iexakt)/Iexakt"<<endl;
    os.close();
#endif
    
#if 0 //test methods of slitdomain_frame_equation
    Index testindex=frame.get_quarklet(1071); 
    cout << "testindex: "<<testindex<<endl;
    cout << "a(lambda_0,lambda_0) : "<<eq.a(testindex,testindex) << endl;
    cout << "f(lambda_0): "<<eq.f(testindex)<<endl;
    cout << "D(lambda_0): "<<eq.D(testindex)<<endl;
    //cout << "F_norm: "<<eq.F_norm() << endl;
    //cout << "normA: "<<eq.norm_A() <<endl;
    //cout << "normAinv: "<<eq.norm_Ainv() << endl;
#endif

    
#if 0   //test evaluate and plot one function
    Array1D<SampledMapping<dim> > evalf(4);
    Index testindex=frame.get_quarklet(1063);   
    cout << "evaluate quarklet with index " << testindex << endl;
    evalf=frame.evaluate(testindex,6);
    std::ofstream osf("slitdomainoutput.m");
    osf << "clf;" << endl;
    osf << "axis([-1 1 -1 1 0 1]);" << endl;
    for(int i=0;i<4;i++){
        evalf[i].matlab_output(osf);
        osf << "surf(x,y,z);" << endl;
        
        osf << "hold on;" << endl;
    }
    osf << "view(30,55);"<<endl;
    osf << "hold off" << endl;
    osf.close(); 
    
    //support test
    typedef Frame::Support Support;
    Index testindex2(testindex); 
    cout<<"test support of index "<<testindex2<<endl;
    Support supp;
    frame.support(testindex2, supp);
    for(int i=0;i<4;i++){
        supp.xmin[i]!=-1 ? cout<<"("<<supp.xmin[i]<<","<<supp.xmax[i]<<")/"<<(1<<supp.j[0])<<" x ("<<supp.ymin[i]<<","<<supp.ymax[i]<<")/"<<(1<<supp.j[1])<<endl : cout<<"0"<<endl;
    }
#endif   
    
  
    
    cout<<"fertig"<<endl;
    
}

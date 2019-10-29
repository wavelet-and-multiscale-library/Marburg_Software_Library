/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


// To reproduce the results of Chapter 8 Diss Keding use DYPLUSEN preconditioner, DELTA1=4, DELTA2 = 2, JMAX=9, PMAX=3



#define POISSON
#undef GRAMIAN


#define DYADIC

//define the parameters \delta_1 and \delta_2 for the H^s weights, cf. Diss Keding Formula (6.1.20) 
//and Theorem 7.12
#ifdef DYADIC
#define DELTA1 6
#define DELTA2 2
#endif

#undef TRIVIAL
#undef ENERGY
#undef DYPLUSEN

#ifdef DYPLUSEN
#define DELTA1 4
#define DELTA2 2
#endif


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
#define JMAX 6
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
#include <general_domain/domain_frame_index.h>
#include <general_domain/domain_frame.h>
#include <general_domain/domain_frame_evaluate.h>
#include <general_domain/domain_frame_support.h>
#include <general_domain/domain_frame_indexplot.h>
#include <galerkin/domain_frame_equation.h>
//#include <galerkin/domain_frame_gramian.h>
//#include <galerkin/cached_quarklet_domain_problem.h>
#include <galerkin/infinite_preconditioner.h>


#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <adaptive/cdd2.h>
#include <adaptive/duv.h>
#include <adaptive/steepest_descent_ks.h>

//#include "ldomain_solutions.h"

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;



template <unsigned int N>
class mySolution
  : public Function<2,double>
{
public:
  virtual ~mySolution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    CornerSingularity cs(Point<2>(0,0), 0.5, 1.5);
//    return sin(M_PI*p[0])*sin(M_PI*p[1])+5*cs.value(p);
    switch(N) {
    case 1:
      return sin(M_PI*p[0])*sin(M_PI*p[1]);
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
    CornerSingularityRHS csrhs(Point<2>(0,0), 0.5, 1.5);
//    return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])+5*csrhs.value(p);
    switch(N) {
    case 1:
      return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
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



int main(){    
    
    cout << "testing quarklet frame on general 2d domain" << endl;
    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
    const int dim = _DIM;
    const int jmax=JMAX;
    const int pmax=PMAX;
    
    //set corners and extensions
    Array1D<Point<dim,int> > corners;
    Array1D<FixedArray1D<int,2> > extensions;

    
//            //high rectangle
//            const int Npatches=2;
//            corners.resize(Npatches);
//            corners[0][0]=-1;
//            corners[0][1]=0;
//            corners[1][0]=-1;
//            corners[1][1]=-1;
//            
//            const int Nextensions=1;
//            extensions.resize(Nextensions);
//            extensions[0][0]=0;
//            extensions[0][1]=1;
//            cout<<"Domain: high rectangle"<<endl;

//            //L-Domain
//            const int Npatches=3;
//            corners.resize(Npatches);
//            corners[0][0]=-1;
//            corners[0][1]=0;
//            corners[1][0]=-1;
//            corners[1][1]=-1;
//            corners[2][0]=0;
//            corners[2][1]=-1;
//            
//            const int Nextensions=2;
//            extensions.resize(Nextensions);
//            extensions[0][0]=0;
//            extensions[0][1]=1;
//            extensions[1][0]=2;
//            extensions[1][1]=1;
//            cout<<"Domain: L-shaped"<<endl;


//            //T-Domain
//            const int Npatches=4;
//            corners.resize(Npatches);
//            corners[0][0]=-1;
//            corners[0][1]=0;
//            corners[1][0]=-1;
//            corners[1][1]=-1;
//            corners[2][0]=0;
//            corners[2][1]=-1;
//            corners[3][0]=-1;
//            corners[3][1]=-2;
//            
//            const int Nextensions=3;
//            extensions.resize(Nextensions);
//            extensions[0][0]=0;
//            extensions[0][1]=1;
//            extensions[1][0]=2;
//            extensions[1][1]=1;
//            extensions[2][0]=3;
//            extensions[2][1]=1;   
//            cout<<"Domain: T-shaped"<<endl;

            //+-Domain
            const int Npatches=5;
            corners.resize(Npatches);
            corners[0][0]=-1;
            corners[0][1]=0;
            corners[1][0]=-1;
            corners[1][1]=-1;
            corners[2][0]=0;
            corners[2][1]=-1;
            corners[3][0]=-1;
            corners[3][1]=-2;
            corners[4][0]=-2;
            corners[4][1]=-1;
            
            const int Nextensions=4;
            extensions.resize(Nextensions);
            extensions[0][0]=0;
            extensions[0][1]=1;
            extensions[1][0]=2;
            extensions[1][1]=1;
            extensions[2][0]=3;
            extensions[2][1]=1; 
            extensions[3][0]=4;
            extensions[3][1]=1;
            cout<<"Domain: +-shaped"<<endl;
            
//            //snake-shaped Domain
//            const int Npatches=5;
//            corners.resize(Npatches);
//            corners[0][0]=0;
//            corners[0][1]=0;
//            corners[1][0]=0;
//            corners[1][1]=1;
//            corners[2][0]=1;
//            corners[2][1]=1;
//            corners[3][0]=2;
//            corners[3][1]=1;
//            corners[4][0]=2;
//            corners[4][1]=2;
//            
//            const int Nextensions=4;
//            extensions.resize(Nextensions);
//            extensions[0][0]=0;
//            extensions[0][1]=1;
//            extensions[1][0]=2;
//            extensions[1][1]=3;
//            extensions[2][0]=4;
//            extensions[2][1]=3;
//            extensions[3][0]=2;
//            extensions[3][1]=1;
//            cout<<"snake-shaped domain"<<endl;
            
//            //big cube (not yet working)
//            const int Npatches=4;
//            corners.resize(Npatches);
//            corners[0][0]=-1;
//            corners[0][1]=-1;
//            corners[1][0]=0;
//            corners[1][1]=-1;
//            corners[2][0]=0;
//            corners[2][1]=0;
//            corners[3][0]=-1;
//            corners[1][1]=0;
//            
//            const int Nextensions=4;
//            extensions.resize(Nextensions);
//            extensions[0][0]=0;
//            extensions[0][1]=1;
//            extensions[1][0]=3;
//            extensions[1][1]=2;
//            extensions[2][0]=1;
//            extensions[2][1]=2;
//            extensions[3][0]=0;
//            extensions[3][1]=3;
 

//            //unit cube (not working)
//            const int Npatches=1;
//            corners.resize(Npatches);
//            corners[0][0]=0;
//            corners[0][1]=0;
//            
//            const int Nextensions=0;
//            extensions.resize(Nextensions);
//            cout<<"Domain: unit cube"<<endl;

  
      
    
    //setup problem
    const int N1=1;
    mySolution<N1> uexact1;
    myRHS<N1> rhs1;
    
    typedef PQFrame<d,dT> Frame1d;
    Frame1d frame1d(false,false);
    frame1d.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_11(true,true);
    frame1d_11.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_01(false,true);
    frame1d_01.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_10(true,false);
    frame1d_10.set_jpmax(jmax-frame1d.j0(),pmax);

    typedef DomainFrame<Frame1d,Npatches> Frame;
    typedef Frame::Index Index;
    Frame frame(&frame1d, &frame1d_11, &frame1d_01, &frame1d_10, corners, extensions);
    cout<<frame.num_real_patches()<<" real patches"<<endl;
    cout<<frame.num_logical_patches()<<" logical patches"<<endl;
    frame.set_jpmax(jmax,pmax);
    
    
    
#if 1 //frame tests
    for(int i=0;i<frame.num_real_patches()+frame.num_logical_patches();i++){
        cout<<"boundary conditions on patch "<<i<<": ("<<frame.frames(i,0)->get_s0()<<","<<frame.frames(i,0)->get_s1()<<")x("<<frame.frames(i,1)->get_s0()<<","<<frame.frames(i,1)->get_s1()<<")"<<endl;
    }
    cout<<"deltasize: "<<frame.Deltasize(frame.j0()[0])<<endl;
#endif
    
#if 1
    {
    //output index set
    set<Index> Lambda;
  for (int i=0; i<frame.degrees_of_freedom();i++) {

    Lambda.insert(*(frame.get_quarklet(i)));
//    cout << *(frame.get_quarklet(i)) << endl;


  }
    }
#endif
    
    PoissonBVP<dim> poisson1(&rhs1);
    DomainFrameEquation<Frame1d,Npatches,Frame> eq(&poisson1, &frame, true);
    
    
#if 1
    //plot one function
    Array1D<SampledMapping<dim> > evalf(Npatches);
//    Index testindex=frame.last_generator(frame.j0());
//    ++testindex;
    Index testindex=frame.get_quarklet(0);
//    Index testindex=frame.last_quarklet(7,frame.j0(),0); //careful with number
    cout << "evaluate quarklet with index " << testindex << endl;
    evalf=frame.evaluate(testindex,6);
    std::ofstream osf("domainoutput.m");
    osf << "clf;" << endl;
//    osf << "axis([-2 2 -2 2 0 1]);" << endl;
    for(int i=0;i<Npatches;i++){
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
    for(int i=0;i<Npatches;i++){
        supp.xmin[i]!=-1 ? cout<<"("<<supp.xmin[i]<<","<<supp.xmax[i]<<")/"<<(1<<supp.j[0])<<" x ("<<supp.ymin[i]<<","<<supp.ymax[i]<<")/"<<(1<<supp.j[1])<<endl : cout<<"0"<<endl;
    }

#endif 
#if 1 //test bilinearform
    const int resolution=6;
    Index testindex1(testindex);
    testindex2=eq.frame().get_quarklet(0);
    //    Support supp;
    if(intersect_supports(frame,testindex1,testindex2,supp)){
        cout<<"support wird geschnitten"<<endl;
    }
    else{
        cout<<"support wird nicht geschnitten"<<endl;
    }
    Array1D<SampledMapping<dim> > eval1(Npatches);
    Array1D<SampledMapping<dim> > eval2(Npatches);
    eval1=frame.evaluate(testindex1,resolution);
    eval2=frame.evaluate(testindex2,resolution);
    std::ofstream os("testa.m");
    os<<"Iexakt="<<eq.a(testindex1,testindex2)<<endl;
    os<<"I=0;"<<endl;
    for(int i=0;i<Npatches;i++){
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
    
    
        
    
#if 1 //plot rhs and exact solution
    std::ofstream osrhs("rhs.m");
    std::ofstream osuexact("uexact.m");
    for(int i=0;i<Npatches;i++){
        Point<2> q0(corners[i][0],corners[i][1]);
        Point<2> q1(corners[i][0]+1,corners[i][1]+1);
        Grid<2> mygrid(q0,q1,100);
        SampledMapping<2> smrhs(mygrid, rhs1); 
        SampledMapping<2> smuexact(mygrid, uexact1);
        smrhs.matlab_output(osrhs);
        smuexact.matlab_output(osuexact);
        osrhs << "surf(x,y,z)"<<endl;
        osuexact<<"surf(x,y,z)"<<endl;
        osrhs << "hold on;" << endl;
        osuexact<<"hold on;"<<endl;
    }
    osrhs << "view(30,55);"<<endl;
    osrhs << "hold off;" << endl;
    osrhs << "grid off;" << endl;
    osrhs << "shading('flat');" << endl;
//    osrhs << "colormap([flipud(jet);jet]);" << endl;
//    osrhs << "set(gca,'CLim', [- min(abs(get(gca,'CLim')))  min(abs(get(gca,'CLim')))]);" << endl;
    osuexact << "view(30,55);"<<endl;
    osuexact<<"hold off;"<<endl;
    osuexact << "grid off;" << endl;
    osuexact << "shading('flat');" << endl;
//    osuexact << "colormap([flipud(jet);jet]);" << endl;
//    osuexact << "set(gca,'CLim', [- min(abs(get(gca,'CLim')))  min(abs(get(gca,'CLim')))]);" << endl;
    osrhs.close();
    osuexact.close();
    cout << "rhs and uexact plotted" << endl;
#endif
    

    
#ifdef NONADAPTIVE
    //setup index set
    set<Index> Lambda;  
    MultiIndex<int,dim> p;p[0]=0;p[1]=0;
    Index lambda = eq.frame().first_generator(eq.frame().j0(), p);
    int zaehler=0;
    for (int l = 0; l < eq.frame().degrees_of_freedom(); l++) {
//        if(lambda.patch()<Npatches) 
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
    const int maxiterations = 999;
    //const double omega = 2.0 / (eq.norm_A() + 1.0/eq.norm_Ainv());
    //cout << omega << endl;
    
    CG(A,F,x,1e-8,maxiterations,iterations); 
//    Richardson(A,F,x,0.04,1e-6,maxiterations,iterations);
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
    Array1D<SampledMapping<dim> > eval(Npatches);
    eval=eq.frame().evaluate(u,6);
//    SampledMapping<2> sm1(evaluate(eq.frame(), u , true, 6));
//    std::ofstream stream1("solution_na.m");
//    sm1.matlab_output(stream1);
//    stream1<<"surf(x,y,z)"<<endl;
//    stream1.close();
    std::ofstream os2("solution_na.m");
    os2 << "figure;" << endl;
    for(int i=0;i<Npatches;i++){
        eval[i].matlab_output(os2);        
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2 << "grid off;" << endl;
    os2 << "shading('flat');" << endl;
    os2.close(); 
    cout << "solution plotted" << endl;
    
    //newer compact coefficients plot
    cout << "plotting quarklet supports" << endl;
    std::ofstream coeff_stream2;
    coeff_stream2.open("coefficients_na.m");
    plot_indices_domain2(&frame, u, coeff_stream2,1e-15);
    coeff_stream2.close();
    cout << "quarklet supports plotted"<<endl;
    
    }
   
#endif
    

    

    cout << "end of test_domain_frame" << endl;
    return 0;
}

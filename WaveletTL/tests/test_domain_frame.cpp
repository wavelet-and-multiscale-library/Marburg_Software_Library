/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


// To reproduce the results of Chapter 8 Diss Keding use DYPLUSEN preconditioner, DELTA1=4, DELTA2 = 2, JMAX=9, PMAX=3



#define POISSON
#undef GRAMIAN


#undef DYADIC

//define the parameters \delta_1 and \delta_2 for the H^s weights, cf. Diss Keding Formula (6.1.20) 
//and Theorem 7.12
#ifdef DYADIC
#define DELTA1 6
#define DELTA2 2
#endif

#undef TRIVIAL
#undef ENERGY
#define DYPLUSEN

#ifdef DYPLUSEN
#define DELTA1 4
#define DELTA2 2
#endif


#undef NONADAPTIVE
#define ADAPTIVE

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
#define PMAX 1
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
//#include <general_domain/domain_frame_indexplot.h>
//#include <galerkin/domain_frame_equation.h>
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
 

//            //unit cube
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
    
#if 0
    //output index set
    set<Index> Lambda;
  for (int i=0; i<frame.degrees_of_freedom();i++) {

    Lambda.insert(*(frame.get_quarklet(i)));
    cout << *(frame.get_quarklet(i)) << endl;


  }
#endif
    
#if 1
    //plot one function
    Array1D<SampledMapping<dim> > evalf(Npatches);
//    Index testindex=frame.last_generator(frame.j0());
//    ++testindex;
//    Index testindex=frame.get_quarklet(207);
    Index testindex=frame.last_quarklet(7,frame.j0(),0); //careful with number
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
    

    

    cout << "end of test_domain_frame" << endl;
    return 0;
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE

#define _WAVELETTL_USE_TBASIS 1
#define _WAVELETTL_USE_TFRAME 1
#define _DIM 2

#include <iostream>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>
#include <Ldomain/ldomain_frame_index.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_frame_evaluate.h>
#include <galerkin/ldomain_frame_equation.h>
#include <galerkin/cached_quarklet_ldomain_problem.h>

#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <adaptive/cdd2.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;


class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
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
    return sin(M_PI*p[0])*sin(M_PI*p[1]);
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
    
    cout << "testing L-domain quarklet frame" << endl;
    const int d  = 3;
    const int dT = 3;
    const int dim = 2;
    const int jmax = 6;
    const int pmax = 0;
    
    typedef PQFrame<d,dT> Frame1d;
    Frame1d frame1d(false,false);
    typedef LDomainFrame<Frame1d> Frame;
    typedef Frame::Index Index;
    //Frame frame(frame1d);
    Frame frame(frame1d);
    frame.set_jpmax(jmax,pmax);
    
    //CornerSingularity uexact1(Point<2>(0,0), 0.5, 1.5);
    mySolution uexact1;
    //CornerSingularityRHS rhs1(Point<2>(0,0), 0.5, 1.5); //Fnorm=0?
    //Vector<double> val(1, "1.0");
    //ConstantFunction<dim> rhs1(val);
    myRHS rhs1;
    
    PoissonBVP<dim> poisson1(&rhs1);
    LDomainFrameEquation<Frame1d,Frame> eq(&poisson1,frame, true);
    eq.set_jpmax(jmax,pmax);
    
    
    Index testindex=frame.get_quarklet(71);
#if 0 //plot rhs and exact solution
    Array1D<Point<dim,int> > corners;
    corners.resize(3);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;

    std::ofstream osrhs("rhs.m");
    std::ofstream osuexact("uexact.m");
    for(int i=0;i<3;i++){
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
    osuexact << "view(30,55);"<<endl;
    osuexact<<"hold off;"<<endl;
    osrhs.close();
    osuexact.close();
    cout << "rhs and uexact plotted" << endl;
#endif
    
#if 1
    //setup index set
    //Vector<double> a;
    //a.resize(eq.frame().degrees_of_freedom(),true);
    //setup Index set
    set<Index> Lambda;  
    MultiIndex<int,dim> p;p[0]=0;p[1]=0;
    Index lambda = eq.frame().first_generator(eq.frame().j0(), p);
    int zaehler=0;
    for (int l = 0; l < eq.frame().degrees_of_freedom(); l++) {
        //cout << lambda << " : "<<lambda.number()<< endl;
        //cout << lambda << " : " << eq.a(lambda,lambda) << endl;
        //if(multi_degree(lambda.e())<1 && lambda.patch()<4){
            Lambda.insert(lambda);
            ++zaehler; 
        //}
        if(lambda==eq.frame().last_quarklet(jmax, p)){
            ++p;
            lambda=eq.frame().first_generator(eq.frame().j0(), p);
        }
        else
        //a(zaehler)=eq.a(eq.frame().get_quarklet(126),lambda);
        ++lambda;
            
    }
    //a.matlab_output("a","a");
    
    cout << "size of Lambda: " << zaehler << endl;
#endif
        
#if 0 //test methods of ldomain_frame_equation
    cout << "testindex: "<<testindex<<endl;
    cout << "a(lambda_0,lambda_0) : "<<eq.a(testindex,testindex) << endl;
    cout << "f(lambda_0): "<<eq.f(testindex)<<endl;
    cout << "D(lambda_0): "<<eq.D(testindex)<<endl;
    //cout << "F_norm: "<<eq.F_norm() << endl;
    //cout << "normA: "<<eq.norm_A() <<endl;
    //cout << "normAinv: "<<eq.norm_Ainv() << endl;
    
    
#endif
 #ifdef NONADAPTIVE   
    //setup stiffness matrix and rhs
    SparseMatrix<double> A;
    setup_stiffness_matrix(eq,Lambda,A); 
    cout << "setup stiffness matrix done" << endl;
    Vector<double> F;
    setup_righthand_side(eq, Lambda, F); 
    A.compress(1e-14);
    F.compress(1e-14);
    A.matlab_output("A","A",0,A.row_dimension()-1,A.column_dimension()-1 );
    F.matlab_output("F","F");
    Vector<double> x(Lambda.size());
    x =0;
    
    //richardson iteration
    unsigned int iterations;
    const int maxiterations = 9999;
    const double omega = 2.0 / (eq.norm_A() + 1.0/eq.norm_Ainv());
    cout << "omega: "<<omega<<endl;
    Richardson(A,F,x,omega,1e-6,maxiterations,iterations);
    //CG(A,F,x,1e-8,maxiterations,iterations);
    cout << "iterations:" << iterations << endl;
    
    
    InfiniteVector<double,Index> u;
    unsigned int i2 = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i2){
      u.set_coefficient(*it, x[i2]);
    }
    
    //plot solution
    //u.COARSE(1e-6,v);
    u.scale(&eq, -1);
    Array1D<SampledMapping<dim> > eval(3);
    eval=eq.frame().evaluate(u,6);
    std::ofstream os2("solution_na.m");
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close(); 
 #endif   

#ifdef ADAPTIVE
    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame> > cproblem1(&eq);
    InfiniteVector<double, Index> F_eta; 
    cproblem1.RHS(1e-6, F_eta);
    const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
    double epsilon = 1e-3;
    InfiniteVector<double, Index> u_epsilon, v;
    const double a=2;
    const double b=2;
    CDD2_QUARKLET_SOLVE(cproblem1, nu, epsilon, u_epsilon, jmax, tensor_simple, pmax, a, b);
    
    //plot solution
    //u.COARSE(1e-6,v);
    u_epsilon.scale(&cproblem1, -1);
    Array1D<SampledMapping<dim> > eval(3);
    eval=cproblem1.frame().evaluate(u_epsilon,6);
    std::ofstream os2("solution_ad.m");
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close();
#endif
    
    

    

#if 0 //test bilinearform
    const int resolution=6;
    Index testindex1(testindex);
    Index testindex2=eq.frame().get_quarklet(192);
    Array1D<SampledMapping<dim> > eval1(3);
    Array1D<SampledMapping<dim> > eval2(3);
    eval1=frame.evaluate(testindex1,resolution);
    eval2=frame.evaluate(testindex2,resolution);
    std::ofstream os("testa.m");
    os<<"Iexakt="<<eq.a(testindex1,testindex2)<<endl;
    os<<"I=0;"<<endl;
    for(int i=0;i<3;i++){
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
    
#if 0 //test rhs
    const int resolution2=6;
    Index testindex3(testindex);
    Array1D<SampledMapping<dim> > eval3(3);
    eval3=frame.evaluate(testindex3, resolution2);
    Array1D<Point<dim,int> > corners;
    corners.resize(3);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;

    std::ofstream os3("testf.m");
    os3<<"Iexakt="<<eq.f(testindex)<<endl;
    os3<<"I=0;"<<endl;
    for(int i=0;i<3;i++){
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
    
    
#if 0
    //plot one function
    Array1D<SampledMapping<dim> > evalf(3);
    //Index ind=frame.get_quarklet(159);    //0-26:generatoren auf patches,
                                        //27-32:überlappende generatoren, indiziert mit p=3,4
                                        //33:überlappendes wavelet
                                        //34:nicht-überlappendes wavelet
    cout << "evaluate quarklet with index " << testindex << endl;
    evalf=frame.evaluate(testindex,6);
    std::ofstream osf("Ldomainoutput.m");
    osf << "clf;" << endl;
    osf << "axis([-1 1 -1 1 0 1]);" << endl;
    for(int i=0;i<3;i++){
        evalf[i].matlab_output(osf);
        osf << "surf(x,y,z);" << endl;
        
        osf << "hold on;" << endl;
    }
    osf << "view(30,55);"<<endl;
    osf << "hold off" << endl;
    osf.close(); 
#endif  
    
    
    
    return 0;
}
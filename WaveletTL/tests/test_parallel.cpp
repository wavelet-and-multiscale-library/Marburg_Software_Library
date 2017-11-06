/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE

#undef DYADIC
#define ENERGY
#define PARALLEL 0

#define FRAME
//#define _WAVELETTL_USE_TBASIS 1
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
#include <Ldomain/ldomain_frame_indexplot.h>
#include <galerkin/ldomain_frame_equation.h>
#include <galerkin/cached_quarklet_ldomain_problem.h>

#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <adaptive/cdd2.h>

#include "ldomain_solutions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;


class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    CornerSingularityRHS csrhs(Point<2>(0,0), 0.5, 1.5);
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
    CornerSingularity cs(Point<2>(0,0), 0.5, 1.5);
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
    
    
    
    cout << "testing L-domain quarklet frame" << endl;
    const int d  = 3;
    const int dT = 3;
    const int dim = 2;
    const int jmax=6;
    const int pmax=1;
    typedef PQFrame<d,dT> Frame1d;
    //Frame1d frame1d(false,false);
    typedef LDomainFrame<Frame1d> Frame;
    typedef Frame::Index Index;
    //typedef Frame::Support Support;
    //Frame frame(frame1d);
    //Frame frame;
    //frame.set_jpmax(jmax,pmax);
    
    //CornerSingularity uexact1(Point<2>(0,0), 0.5, 1.5);
    mySolution uexact1;
    //CornerSingularityRHS rhs1(Point<2>(0,0), 0.5, 1.5); 
    //Vector<double> val(1, "1.0");
    //ConstantFunction<dim> rhs1(val);
    myRHS rhs1;
    
    PoissonBVP<dim> poisson1(&rhs1);
    Frame1d frame1d(false,false);
    frame1d.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_11(true,true);
    frame1d_11.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_01(false,true);
    frame1d_01.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_10(true,false);
    frame1d_10.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame frame(&frame1d, &frame1d_11, &frame1d_01, &frame1d_10);
    frame.set_jpmax(jmax,pmax);
    LDomainFrameEquation<Frame1d,Frame> eq(&poisson1, &frame, false);
//    LDomainFrameEquation<Frame1d,Frame> eq(&poisson1, false);
//    eq.set_jpmax(jmax,pmax);
//    for(int i=0;i<1000;i++){
//        cout<<i<<endl;
//        eq.set_jpmax(jmax,pmax);
//    }
    
#if 0 //compare rhs
    InfiniteVector<double, Index> F_eta;
    eq.RHS(1e-6, F_eta);
//    cout<<F_eta<<endl;
    std::ofstream stream1;
#if PARALLEL==1
    stream1.open("F_etapar.m");
#else
    stream1.open("F_etaseq.m");
#endif
    for(InfiniteVector<double, Index>::const_iterator it(F_eta.begin()); it!=F_eta.end(); ++it){
        stream1<<*it<<endl;
    }
    stream1.close();  
#endif 
    
#if 0
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
     
#if 0 //runtime test for parallel setup
    SparseMatrix<double> A;
    tic=clock();
    setup_stiffness_matrix(eq, Lambda, A);  
    toc = clock();
    time = (double)(toc-tic);
    cout << "time: " << (time/CLOCKS_PER_SEC) << " s\n"<<endl; 
    Vector<double> x0(A.row_dimension(), false), yseq(A.row_dimension(), false), ypar(A.row_dimension());
    x0=1;    
#if 0    
    cout <<"testing apply"<<endl;
//    for(int i=0;i<1e4;i++){
//        A.apply<Vector<double> >(x0,y);
//    }
    A.apply<Vector<double> >(x0,ypar);
    A.apply_transposed<Vector<double> >(x0,yseq);
    ypar-=yseq;
    cout<<"l2-error="<<l2_norm(ypar)<<endl;
    cout<<"linfty-error="<<linfty_norm(ypar)<<endl;
#endif  
    cout<<"computing eigenvalues"<<endl;
    unsigned int iterations;
    double lambdamax=PowerIteration<Vector<double>,SparseMatrix<double> >(A, x0, 1e-3, 100, iterations);
    cout<<"lambdamax= "<<lambdamax<<endl;
    double lambdamin=InversePowerIteration<Vector<double>,SparseMatrix<double> >(A, x0, 1e-10,  1e-3, 100, iterations);
    cout<<"lambdamin= "<<lambdamin<<endl;
#endif
    
#if 1
    //cout << "bin hier"<<endl;
//    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame> > cproblem1(&eq,66,23);
    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame> > cproblem1(&eq);
//    cout<<"hallo"<<endl;

    cout<<"normA="<<cproblem1.norm_A()<<endl;
    cout<<"normAinv="<<cproblem1.norm_Ainv()<<endl;
    
    InfiniteVector<double, Index> F_eta;
    cproblem1.RHS(1e-6, F_eta);

#if 0 //compare rhs
    std::ofstream stream1;
#if PARALLEL==1
    stream1.open("F_etapar.m");
#else
    stream1.open("F_etaseq.m"); 
#endif
    for(InfiniteVector<double, Index>::const_iterator it(F_eta.begin()); it!=F_eta.end(); ++it){
        stream1<<*it<<endl;
    }
    stream1.close();  
#endif 
    
    
    //const double nu = cproblem1.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
    InfiniteVector<double, Index> u_epsilon, v,w,Av,Avseq,Avopar,Avipar;
    for(int i=0;i<eq.frame().degrees_of_freedom();i++){
        v.set_coefficient(eq.frame().get_quarklet(i),1);
    }
    const double a=2;
    const double b=2;
    
    
#if 1
    cout<<"comparing APPLY"<<endl;
    
    std::ofstream stream2;
    APPLY_QUARKLET(cproblem1, v, 1e-3, Av, jmax, tensor_simple, pmax, a,b);
#if PARALLEL==1
//    APPLY_QUARKLET_PARALLEL_INNER(cproblem1, v, 1e-1, Av, jmax, tensor_simple, pmax, a,b);
    stream2.open("Av_ipar.m");
#else
//    APPLY_QUARKLET_SEQUENTIAL(cproblem1, v, 1e-1, Av, jmax, tensor_simple, pmax, a,b);
    stream2.open("Av_seq.m"); 
#endif
    for(InfiniteVector<double, Index>::const_iterator it(Av.begin()); it!=Av.end(); ++it){
        stream2<<*it<<endl;
    }
    stream2.close();
#endif 
    
    
//    APPLY_QUARKLET_PARALLEL_INNER(cproblem1, v, 1e-1, Avipar, jmax, tensor_simple, pmax, a,b);
//    APPLY_QUARKLET_PARALLEL_OUTER(cproblem1, v, 1e-1, Avopar, jmax, tensor_simple, pmax, a,b);
//    APPLY_QUARKLET_SEQUENTIAL(cproblem1, v, 1e-1, Avseq, jmax, tensor_simple, pmax, a,b);
    
//    Avseq-=Avipar;
//    cout<<"l2-error="<<l2_norm(Avseq)<<endl;
//    cout<<"linfty-error="<<linfty_norm(Avseq)<<endl;
    
//     CDD2_QUARKLET_SOLVE(cproblem1, nu, epsilon, u_epsilon, jmax, tensor_simple, pmax, a, b);
#endif   
    cout << "fertig" << endl;
    return 0;
}
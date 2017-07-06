#undef BASIS
#define FRAME
#ifdef FRAME
#define _WAVELETTL_USE_TFRAME 1
#else
#define _WAVELETTL_USE_TBASIS 1
#endif

#define _DIM 2

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/function.h>
#include <utils/fixed_array1d.h>
#include <numerics/bvp.h>
#include <numerics/corner_singularity.h>
#include <geometry/sampled_mapping.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 1
#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_index.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_frame_index.h>
#include <Ldomain/ldomain_frame_support.h>
#include <Ldomain/ldomain_frame_evaluate.h>
#include <galerkin/ldomain_frame_equation.h>
#include <cube/tframe.h>
#include <galerkin/tframe_equation.h>
#include <galerkin/cached_quarklet_ldomain_problem.h>
#include <adaptive/cdd2.h>


#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <galerkin/ldomain_equation.h>
#include <galerkin/ldomain_helmholtz_equation.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


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


int main()
{
  cout << "Testing wavelet-Galerkin solution of an elliptic equation on the L-shaped domain ..." << endl;

  const int d  = 3;
  const int dT = 3;
  
  const int jmax=6;
  const int pmax=0;
  const int a = 2;
  const int b = 2;
//  typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
  
#ifdef FRAME
  typedef PQFrame<d,dT> Frame1D;
//  typedef PBasis<d,dT> Frame1D;
  Frame1D frame1D_11(true,true);
  
  Frame1D frame1D;
//  typedef Frame1D::Index Index1D;
  
  typedef LDomainFrame<Frame1D> LFrame;
  typedef LFrame::Index Index;
  typedef Index::polynomial_type polynomial_type;
  polynomial_type p;
//  typedef Index::level_type level_type;
  LFrame frame(frame1D);
//  LFrame frame2();
  
  myRHS rhs1;
    
//  PoissonBVP<dim> poisson1(&rhs1);
  
//  Vector<double> val(1, "1.0");
//  ConstantFunction<2> rhs1(val);
  PoissonBVP<2> poisson1(&rhs1);
  LDomainFrameEquation<Frame1D,LFrame> eq(&poisson1, frame, true);
  Index myindex=eq.frame().get_quarklet(192);
//  for(int i= 0;i<8;i++,++myindex){
//  cout << "RHS with quarklet: " <<  myindex << ", " << eq.f(myindex) << endl;
//  }
//  abort();
//  cout << "DeltaLmin=" << eq.frame().frame1d().DeltaLmin() << endl; 
  eq.set_jpmax(jmax, pmax);
  CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1D,LFrame> > ceq(&eq);
  InfiniteVector<double, Index> F_eta; 
  ceq.RHS(1e-6, F_eta);
////    cout << F_eta << endl;
    const double nu = ceq.norm_Ainv() * l2_norm(F_eta);   //benötigt hinreichend großes jmax
//    cout << "TEST: " << l2_norm(F_eta) << ", " << cproblem1.norm_Ainv() << endl;
  double epsilon = 1e-3;
  InfiniteVector<double, Index> u_epsilon, v, test, test2;
    
    //cout << eq.s_star() << endl;
    //cout << cproblem1.s_star() << endl;
    //cout << eq.space_dimension << endl;
    //cout << cproblem1.space_dimension << endl;
//    //CDD1_SOLVE(cproblem1, epsilon, u_epsilon, 2*jmax1, tensor_simple);
    
    
    


    CDD2_QUARKLET_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, tensor_simple, pmax, a, b);
    //plot solution
    cout << "plotting solution" << endl;
    u_epsilon.scale(&ceq, -1);
    Array1D<SampledMapping<2> > eval(3);
    eval=eq.frame().evaluate(u_epsilon,6);
    std::ofstream os2("solution_cdd2.m");
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close(); 
//    SampledMapping<2> sm1(evaluate(ceq.frame(), u_epsilon , true, 6));
//    std::ofstream stream1("solution_cdd2.m");
//    sm1.matlab_output(stream1);
//    stream1.close();
    cout << "solution plotted" << endl;
    abort();
    
//    Array1D<SampledMapping<dim> > eval(3);
//    eval=eq.frame().evaluate(u,6);
//    std::ofstream os2("solution_na.m");
//    for(int i=0;i<3;i++){
//        eval[i].matlab_output(os2);
//        os2 << "surf(x,y,z);" << endl;
//        os2 << "hold on;" << endl;
//    }  
//    os2 << "view(30,55);"<<endl;
//    os2 << "hold off" << endl;
//    os2.close(); 
  
  //  plot one function
//    polynomial_type p2(0,0);
//    level_type j2(3,3);
//    Array1D<SampledMapping<2> > eval(3);
//    Index ind=frame.last_quarklet(j2,p2); 
//    Index ind=ceq.frame().get_quarklet(192), indrun=ceq.frame().get_quarklet(70), indrunnew=ceq.frame().get_quarklet(1);   //0-26:generatoren auf patches,
//                                        //27-32:überlappende generatoren, indiziert mit p=3,4
//                                        //
//    cout << "bilinearform of all quarklets with " << ind << endl;
////    for(int i=0; i<800; i++, ++indrun){
////    
////    eval=frame.evaluate(ind,6);
////    std::ofstream os("Ldomainoutput.m");
////    os << "clf;" << endl;
////    os << "axis([-1 1 -1 1 0 1]);" << endl;
////    for(int i=0;i<3;i++){
////        eval[i].matlab_output(os);
////        os << "surf(x,y,z);" << endl;
////        
////        os << "hold on;" << endl;
////    }
////    os << "view(30,55);"<<endl;
////    os << "hold off" << endl;
////    os.close(); 
//    cout << "Bilinearform Cached/Uncached for quarklet with index " << indrun  << ",  " /*<< ceq.a(ind,indrun) << ", " */<< eq.a(ind,indrun) << endl;
//cout << "Bilinearform Cached/Uncached for quarklet with index " << indrunnew  << ",  " /*<< ceq.a(ind,indrunnew) << ", "*/ << eq.a(ind,indrunnew) << endl;
//    //    }
//    
    
    
  
  
  
//  InfiniteVector<double, Index> F_eta; 
//  ceq.RHS(1e-6, F_eta);
//////    cout << F_eta << endl;
//    const double nu = ceq.norm_Ainv() * l2_norm(F_eta);   
////    cout << "TEST: " << l2_norm(F_eta) << ", " << ceq.norm_Ainv() << endl;
//    double epsilon = 1e-3;
//    InfiniteVector<double, Index> u_epsilon;
//    
//    //cout << eq.s_star() << endl;
//    //cout << ceq.s_star() << endl;
//    //cout << eq.space_dimension << endl;
//    //cout << ceq.space_dimension << endl;
//    CDD2_QUARKLET_SOLVE(ceq, nu, epsilon, u_epsilon, jmax, tensor_simple, pmax, a, b); 
    
  
  
//  frame.set_jpmax(jmax,pmax);
//  abort();
//  eq.norm_A();
  
//  typedef PQFrame<d,dT> Frame1d;
//    typedef Basis1d::Index Index1d;
#if 1
    typedef TensorFrame<Frame1D,2> CubeFrame;
    
    typedef CubeFrame::Index CubeIndex;
    FixedArray1D<bool, 4> bc;
    bc[0]=true, bc[1]=true, bc[2]=false, bc[3]=true;
    CubeFrame cubeframe(bc);
    TensorFrameEquation<Frame1D,2,CubeFrame> cubeeq(&poisson1, bc);
    CubeIndex ind2=cubeeq.frame().get_quarklet(0), ind2run=cubeeq.frame().get_quarklet(0);
    cout << "Cubeevaluate quarklet with index: " << ind2 << endl;
    for(int i=0;i<272;i++, ++ind2run){

    cout << "Bilinearform Cube for q with index: " << ind2run << ", " << cubeeq.a(ind2,ind2run) << endl;
    }
//    cubeeq.norm_A();
#endif
    
//    abort();
//    MultiIndex<int,2> lp(0,0), lj(3,3), le(0,1), lk(1,0), mp(0,0), mj(3,3), me(0,0), mk(1,1);
//    Index lambdal(lp,lj,le,1,lk,0,&eq.frame()), mul(mp,mj,me,1,mk,0,&eq.frame());
//    
//    cout << last_generator(&ceq.frame(), lj,lp) << last_generator(&ceq.frame(), lj,lp).number() << endl;
//    cout << ceq.frame().last_generator(lj,lp) << ceq.frame().last_generator(lj,lp).number() << endl;
//    cout << first_generator(&ceq.frame(), lj,lp) << first_generator(&ceq.frame(), lj,lp).number() << endl;
//    cout << ceq.frame().first_generator(lj,lp) << ceq.frame().first_generator(lj,lp).number() << endl;
//    cout << first_quarklet(&ceq.frame(), lj,lp) << first_quarklet(&ceq.frame(), lj,lp).number() << endl;
//    cout << ceq.frame().first_quarklet(lj,lp) << ceq.frame().first_quarklet(lj,lp).number() << endl;
//    cout << last_quarklet(&ceq.frame(), lj,lp) << last_quarklet(&ceq.frame(), lj,lp).number() << endl;
//    cout << ceq.frame().last_quarklet(lj,lp) << ceq.frame().last_quarklet(lj,lp).number() << endl;
//    ceq.a(ceq.frame().get_quarklet(210),ceq.frame().get_quarklet(49));
//    
////    cout << "Nablasize: " << ceq.frame().get_Nablasize() << endl;
//    cout << ceq.frame().get_first_wavelet_numbers()[0] << endl;
//    
////    cout << "Cached: " << ceq.a(ceq.frame().get_quarklet(210),ceq.frame().get_quarklet(49)) << ", Uncached: " << eq.a(ceq.frame().get_quarklet(210),ceq.frame().get_quarklet(49)) << endl;
////    ceq.norm_A();
//    double time2=0.0, tstart2; 
//    tstart2 = clock();
//    cout << "Cached: " << endl;
//    cout << "norm_A: " << ceq.norm_A() <<endl;
//    cout << "norm_Ainv: " << ceq.norm_Ainv() <<endl;
//    time2 += (clock() - tstart2)/CLOCKS_PER_SEC;
//    cout << "time = " << time2 << " sec." << endl;
//    double time1=0.0, tstart1;
//    tstart1 = clock();
//    cout << "Uncached: " << endl;
//    cout << "norm_A: " << eq.norm_A() <<endl;
//    cout << "norm_Ainv: " << eq.norm_Ainv() <<endl;
//    time1 += (clock() - tstart1)/CLOCKS_PER_SEC;
//    cout << "time = " << time1 << " sec." << endl;
//    cout << "Cached Norms: " << ceq.norm_A() << ", " << ceq.norm_Ainv() << ", Uncached norms: " << eq.norm_A() << ", " << eq.norm_Ainv() << endl;

    abort();
    //    eq.frame().suppcache()[lambdal];
    
//  eq.a(lambdal,lambdal);
  
  InfiniteVector<double, Index> coeffs;
  Index lambda = frame.first_generator(frame.j0());
//  typedef LFrame::Support Support;
#if 0
  CubeIndex lambdacube(lp,lj,le,lk,&cubeeq.frame()), mucube(mp,mj,me,mk,&cubeeq.frame());
  typedef CubeFrame::Support CubeSupport;
  Support supp2;
  CubeSupport cubesupp;
    
  intersect_supports(cubeeq.frame(), lambdacube, mucube, cubesupp);
    //    cout << "some tests: " << eq.frame().get_jmax() << eq.frame().j0() << endl;
     cout << "wichtig evaluate quarklet with index " << lambdacube << ", " << mucube << endl;
    cout << "j[0]: " << cubesupp.j[0] << ", j[1]: " << cubesupp.j[1] << endl;
    cout << "[" << cubesupp.a[0] <<" , " <<cubesupp.b[0]<<"]x["<<cubesupp.a[1]<<" , "<<cubesupp.b[1]<<"]"<<  endl;
//    cout << "Cube: " << cubeeq.a(lambdacube,lambdacube) << ", LDomain: " << eq.a(lambdal,lambdal) << endl; 
#endif
    
//    Support supp
//    support(eq.frame(),lambdal,supp);
////    support(eq.frame(),lambdal,supp);
//    //intersect_supports(frame,frame.get_quarklet(27),frame.get_quarklet(29),supp);
//    cout << "wichtig evaluate quarklet with index " << lambdal << endl;
//    cout << "j[0]: " << supp.j[0] << ", j[1]: " << supp.j[1] << endl;
//    cout << "patch 0: [" << supp.xmin[0] <<" , " <<supp.xmax[0]<<"]x["<<supp.ymin[0]<<" , "<<supp.ymax[0]<<"]"<<  endl;
//    cout << "patch 1: [" << supp.xmin[1] <<" , " <<supp.xmax[1]<<"]x["<<supp.ymin[1]<<" , "<<supp.ymax[1]<<"]"<<  endl;
//    cout << "patch 2: [" << supp.xmin[2] <<" , " <<supp.xmax[2]<<"]x["<<supp.ymin[2]<<" , "<<supp.ymax[2]<<"]"<<  endl;
//  
// 
    
#if 0
  intersect_supports(eq.frame(), lambdal, mul, supp2);
  cout << "intersect support: " << endl;
    cout << "j[0]: " << supp2.j[0] << ", j[1]: " << supp2.j[1] << endl;
    cout << "patch 0: [" << supp2.xmin[0] <<" , " <<supp2.xmax[0]<<"]x["<<supp2.ymin[0]<<" , "<<supp2.ymax[0]<<"]"<<  endl;
    cout << "patch 1: [" << supp2.xmin[1] <<" , " <<supp2.xmax[1]<<"]x["<<supp2.ymin[1]<<" , "<<supp2.ymax[1]<<"]"<<  endl;
    cout << "patch 2: [" << supp2.xmin[2] <<" , " <<supp2.xmax[2]<<"]x["<<supp2.ymin[2]<<" , "<<supp2.ymax[2]<<"]"<<  endl;
//    eq.a(lambdal,lambdal);
    cout << "Cube: " << cubeeq.a(lambdacube,mucube) << ", LDomain: " << eq.a(lambdal,mul) << endl;
    
    eq.norm_A();
//    cubeeq.norm_A();
    cout << "Norm: " << eq.norm_A() << ", Norm_Inv: " << eq.norm_Ainv() << endl;
#endif
  
  
  
  
  
  
//  Support suppconst;
//  suppconst.j[0]=0, suppconst.j[1]=0;
//  for(int i=0;i<3;i++){
//      suppconst.xmin[i]=0, suppconst.xmax[i]=0,suppconst.ymin[i]=0,suppconst.ymax[i]=0;
//  }
//  
//  
//  set<Index> Lambda;
//  int counter = 0;
//  for (Index lambda = frame.first_generator(frame.j0());counter<=3000;) {
//      Support supp = suppconst;
//      frame.support(lambda,supp);
//      cout << "("<< supp.j[0] << "," << supp.j[1] << ", [" << supp.xmin[0]  << "," << supp.xmax[0]<< "]x["   << supp.ymin[0]<< "," << supp.ymax[0]<< "]"<< ", ["<< supp.xmin[1]<< "," << supp.xmax[1]<< "]x["<< supp.ymin[1]
//           << ","<< supp.ymax[1]<< "]"<< ", ["<< supp.xmin[2]<< ","<< supp.xmax[2]<< "]x["<< supp.ymin[2]<< ","<< supp.ymax[2]<< "])"<< endl;
//    Lambda.insert(lambda);
////    cout << lambda << ", " << lambda.number() << endl;
//    if(lambda==frame.last_quarklet(5,p)){
//            ++p;
//            lambda=frame.first_generator(frame.j0(), p, counter+1);
//        }
//            else
//                ++lambda;
//    
//    counter++;
////    if (lambda == frame.last_quarklet(frame.j0())) break;
//  }
//  cout << "Test" << endl;
//  polynomial_type p2(0,0);
//  level_type j2(2,3);
//  Support supp;
//  Index lambda2 = frame.first_generator(frame.j0(), p2);
//  Index lambda3 = frame.first_quarklet(j2,p2);
//  intersect_supports(frame,lambda2,lambda3,supp);
//  cout << "("<< supp.j[0] << "," << supp.j[1] << ", [" << supp.xmin[0]  << "," << supp.xmax[0]<< "]x["   << supp.ymin[0]<< "," << supp.ymax[0]<< "]"<< ", ["<< supp.xmin[1]<< "," << supp.xmax[1]<< "]x["<< supp.ymin[1]
//       << ","<< supp.ymax[1]<< "]"<< ", ["<< supp.xmin[2]<< ","<< supp.xmax[2]<< "]x["<< supp.ymin[2]<< ","<< supp.ymax[2]<< "])"<< endl;
//  cout << lambda2 << ", " << lambda2.number() << endl;
  
  

    
//        cout << "support of quarklet on patches:"  << endl;
//    Support supp;
//    frame.support(ind, supp);
//    //intersect_supports(frame,frame.get_quarklet(27),frame.get_quarklet(29),supp);
//    cout << "j[0]: " << supp.j[0] << ", j[1]: " << supp.j[1] << endl;
//    cout << "patch 0: [" << supp.xmin[0] <<" , " <<supp.xmax[0]<<"]x["<<supp.ymin[0]<<" , "<<supp.ymax[0]<<"]"<<  endl;
//    cout << "patch 1: [" << supp.xmin[1] <<" , " <<supp.xmax[1]<<"]x["<<supp.ymin[1]<<" , "<<supp.ymax[1]<<"]"<<  endl;
//    cout << "patch 2: [" << supp.xmin[2] <<" , " <<supp.xmax[2]<<"]x["<<supp.ymin[2]<<" , "<<supp.ymax[2]<<"]"<<  endl;
    
//    const char* filenameSolution2 = "dim1output.m";
//    Index1D lambda1D(0,4,1,15,&frame1D_11);
//    InfiniteVector<double, Index1D> u_epsilon;
//    u_epsilon[lambda1D]=1;
////    u_epsilon.scale(&ceq, -1);
//    SampledMapping<1> s2(evaluate(frame1D_11, u_epsilon, true, 7));
//    std::ofstream u_stream2(filenameSolution2);
//    s2.matlab_output(u_stream2);
//    u_stream2 << "figure;\nplot(x,y);"
//              << "title('Testausgabe');"
//              << endl;
//
//    u_stream2.close();
//    cout << "Kontrolle: mit RB: " << frame1D_11.j0() << ", ohne RB: " << frame1D.j0() << endl;
     

  
#endif
#ifdef BASIS
   typedef PBasis<d,dT> Basis1D;
  
  Basis1D basis1D;
//  basis1D.orthogonalize_boundary_wavelets(); // cf. [Barsch]

  typedef LDomainBasis<Basis1D> LBasis;
  typedef LBasis::Index Index;
  typedef LBasis::Support Support;
  Support supp;

  CornerSingularity    u_sing(Point<2>(0,0), 0.5, 1.5);
  CornerSingularityRHS f_sing(Point<2>(0,0), 0.5, 1.5);

//#if 0
//  PoissonBVP<2> poisson(&f_sing);
//#else
  EigenRHS rhs;
  PoissonBVP<2> poisson(&rhs);
//#endif
//
  LBasis basis(basis1D);
//
//#if 1
  LDomainEquation<Basis1D> eq(&poisson, basis);
//#else
//  LDomainHelmholtzEquation<Basis1D> eq(&rhs, basis, 0.0);
//#endif
//  
  InfiniteVector<double, Index> coeffs;
  eq.RHS(1e-8, coeffs);

  set<Index> Lambda;
  for (Index lambda = eq.basis().first_generator(eq.basis().j0());; ++lambda) {
    basis.support(lambda,supp);
      cout << "("
         << supp.j
         << ", ["
         << supp.xmin[0]
         << ","
         << supp.xmax[0]
         << "]x["
         << supp.ymin[0]
         << ","
         << supp.ymax[0]
         << "]"
         << ", ["
         << supp.xmin[1]
         << ","
         << supp.xmax[1]
         << "]x["
         << supp.ymin[1]
         << ","
         << supp.ymax[1]
         << "]"
         << ", ["
         << supp.xmin[2]
         << ","
         << supp.xmax[2]
         << "]x["
         << supp.ymin[2]
         << ","
         << supp.ymax[2]
         << "])"
         << endl;
      
      Lambda.insert(lambda);
    cout << lambda << endl;
    if (lambda == eq.basis().last_wavelet(eq.basis().j0()+1)) break;
  }

   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
     cout << *it << endl;

  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();

  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
//#if 0
//  std::ofstream ofs("stiff_out.m");
//  ofs << "M=";
//  print_matrix(A,ofs);
//  ofs.close();
//#endif
//
//  tend = clock();
//  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//  cout << "  ... done, time needed: " << time << " seconds" << endl;
////   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;
//
  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- right hand side: " << b << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-8, 100, iterations);
//  
//// //   cout << "- solution coefficients: " << x;
//  cout << " with residual (infinity) norm ";
//  A.apply(x, err);
//  err -= b;
//  cout << linfty_norm(err) << endl;
//
//#if 1
  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);

  u.scale(&eq, -1);
  Array1D<SampledMapping<2> > s((eq.basis()).evaluate(u, 6));
//  std::ofstream u_Lambda_stream("u_lambda.m");
//  matlab_output(u_Lambda_stream, s);
//  u_Lambda_stream.close();
//  cout << "  ... done, see file 'u_lambda.m'" << endl;
//   EigenSolution u_Lambda;
//   s.add(-1.0, SampledMapping<2>(s, u_Lambda));
//   cout << "  ... done, pointwise error: " << row_sum_norm(s.values()) << endl;
//#endif
#endif
  return 0;
}

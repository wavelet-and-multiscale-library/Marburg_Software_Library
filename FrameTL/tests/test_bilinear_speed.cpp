#define TWO_D
#define PARALLEL 0

#define QUARKLET
#undef AGGREGATED

#define DYADIC
#undef TRIVIAL
#undef ENERGY

#define _DIM 2

#ifdef QUARKLET
#define FRAME
#define _WAVELETTL_USE_TFRAME 1
#endif


#include <map>
#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/spline_basis.h>
#include <numerics/corner_singularity.h>

#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>

#include <galerkin/galerkin_utils.h>
#ifdef AGGREGATED
#include <simple_elliptic_equation.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/cached_problem.h>
#include <frame_support.h>
#endif
#include <interval/pq_frame.h>
#include <interval/pq_evaluate.h>
#ifdef QUARKLET
#include <Ldomain/ldomain_frame_index.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_frame_evaluate.h>
#include <Ldomain/ldomain_frame_indexplot.h>
#include <galerkin/ldomain_frame_equation.h>
#include <galerkin/cached_quarklet_ldomain_problem.h>
#endif




using std::cout;
using std::endl;
#ifdef AGGREGATED
using FrameTL::EvaluateFrame;
using FrameTL::SimpleEllipticEquation;
using FrameTL::AggregatedFrame;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;
#endif
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::CornerSingularityRHS;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;


using namespace std;
#ifdef AGGREGATED
//using namespace FrameTL;
#endif
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  
  cout << "Testing class SimpleEllipticEquation..." << endl;
  
  const int DIM = 2;
  const int jmax = 6;
  
#ifdef QUARKLET
  const int pmax = 0;
#endif
//  const int jmax1d=jmax-j0;

  const int d  = 3;
  const int dT = 3;

  //typedef DSBasis<2,2> Basis1D;
#ifdef AGGREGATED
  typedef PBasis<d,dT> Basis1D;
  //typedef SplineBasis<d,dT,P_construction> Basis1D;

  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  //typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

  EvaluateFrame<Basis1D,2,2> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.;
  b[1] = -1.;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.;
  A2(1,1) = 2.;
  Point<2> b2;
  b2[0] = -1.;
  b2[1] = -1.;
  AffineLinearMapping<2> affineP2(A2,b2);
  //##############################

  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;
 
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;//2

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;//2
  bound_2[2] = 1;
  bound_2[3] = 1;

  bc[1] = bound_2;

//to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;
 
  bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, jmax);
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);
#endif
#ifdef QUARKLET
    typedef PQFrame<d,dT> Frame1d;
    //Frame1d frame1d(false,false);
    typedef LDomainFrame<Frame1d> Frame2D;
    
    typedef Frame2D::Index Index;


    
#endif

  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;

  CornerSingularityRHS singRhs(origin, 0.5, 1.5);
  CornerSingularity sing2D(origin, 0.5, 1.5);
  
  PoissonBVP<DIM> poisson(&singRhs);
  //PoissonBVP<DIM> poisson(&const_fun);

  // BiharmonicBVP<DIM> biahrmonic(&singRhs);
#ifdef AGGREGATED
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson/*, 5.0048, 1.0/0.01*/);
#endif
#ifdef QUARKLET
  
    Frame1d frame1d(false,false);
    frame1d.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_11(true,true);
    frame1d_11.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_01(false,true);
    frame1d_01.set_jpmax(jmax-frame1d.j0(),pmax);
    Frame1d frame1d_10(true,false);
    frame1d_10.set_jpmax(jmax-frame1d.j0(),pmax);
    
//    typedef typename Frame1d::Index Index1D;
//    
//    for (Index1D it=frame1d_01.first_generator(3);it<=frame1d_11.last_wavelet(frame1d.get_jmax_(), frame1d_11.get_pmax_());++it){
//        cout << it << endl;
//    }
    
//    abort();
  
//    Frame2D frame(frame1d);
    Frame2D frame(&frame1d, &frame1d_11, &frame1d_01, &frame1d_10);
    frame.set_jpmax(jmax,pmax);
    LDomainFrameEquation<Frame1d,Frame2D> discrete_poisson(&poisson, &frame, false);
//    discrete_poisson.set_jpmax(jmax,pmax);
//    Frame2D frame = discrete_poisson.frame();
//    
    
    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame2D> > problem(&discrete_poisson, 1., 1.);
#endif
  

  cout.precision(12);
  
  InfiniteVector<double,Frame2D::Index> rhs;
  //Vector<double> v(225);
  
  set<Index> Lambda;
  for (int i=0; i<frame.degrees_of_freedom();i++) {
//  for (int i=0; i<817;i++) {    
#ifdef AGGREGATED
    Lambda.insert(*frame.get_wavelet(i));
    //    cout << *frame.get_wavelet(i) << endl;
#endif
#ifdef QUARKLET
//    if(((*(frame.get_quarklet(i))).k()[1]>1 && (*(frame.get_quarklet(i))).k()[1]<6) || (*(frame.get_quarklet(i))).p()[1]==0){
//    if(i>=208){
        Lambda.insert(*(frame.get_quarklet(i)));
        cout << *(frame.get_quarklet(i)) << endl;
//    }
//    }
//    cout << discrete_poisson.a(*(frame.get_quarklet(i)),*(frame.get_quarklet(i)))<<endl;
#endif

  }
//  abort();
  
  
//  for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
//       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, jmax); ++lambda) {
//    Lambda.insert(lambda);
//    cout << lambda << endl;
//  }
  
  
  
  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;
  
  
  clock_t tstart, tend;
  double time;
  tstart = clock();
//   WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff, false);
  WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
// WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff, false);
//  WaveletTL::setup_stiffness_matrix(discrete_poisson, Lambda, stiff);
  
    
  
  
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  stiff.matlab_output("stiff_2D_out", "stiff",1); 
//  cout << "Norm/InvNorm: " << problem.norm_A() << ", " << problem.norm_Ainv() << endl;
  abort();
  
 
  
  cout << "setting up full right hand side..." << endl;
  Vector<double> rh;
  WaveletTL::setup_righthand_side(problem, Lambda, rh);
  cout << rh << endl;
  

  
  Vector<double> x(Lambda.size()); x = 1;
  //double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  double lmax = 1;
//  cout << "lmax = " << lmax << endl;

  x = 1;
  //double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  
  cout << "performing iterative scheme to solve projected problem..." << endl;
  Vector<double> xk(Lambda.size()), err(Lambda.size()); xk = 0;
  
 

  //CG(stiff, rh, xk, 1.0e-6, 100, iter);
  //cout << "CG iterations needed: "  << iter << endl;
  //Richardson(stiff, rh, xk, 2. / lmax - 0.01, 0.0001, 2000, iter);
  double alpha_n = 2. / lmax - 0.001;
  
  Vector<double> resid(xk.size());
  Vector<double> help(xk.size());
  for (int i = 0; i < 5000;i++) {
    stiff.apply(xk,help);
    resid = rh - help;
    cout << "i = " << i << " " << sqrt(resid*resid) << endl;
    stiff.apply(resid,help);
    alpha_n = (resid * resid) * (1.0 / (resid * help));
    resid *= alpha_n;
    xk = xk + resid;
  }

  cout << "performing output..." << endl;
  
  
  InfiniteVector<double,Frame2D::Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i){
//    if(i>=208 && i<800)
    u.set_coefficient(*it, xk[i]);
//    if(i==799) break;
  }
  
//  cout << u << endl;
  u.scale(&discrete_poisson,-1);
  cout << "u_new: " << endl << u << endl;
  
#ifdef AGGREGATED
  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u, true, 6);//expand in primal basis
 
  
  std::ofstream ofs5("approx_solution_out.m");
  matlab_output(ofs5,U);
  ofs5.close();
#endif
  
#ifdef QUARKLET
  
 Array1D<SampledMapping<2> > eval(3);
    eval=frame.evaluate(u,6);
    std::ofstream os2("approx_solution_out.m");
    os2 << "figure\n" << endl;
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }  
    os2 << "view(30,55);"<<endl;
    os2 << "hold off;" << endl;
//    os2 << "figure;" << endl;
    os2.close();  
#endif
  
   return 0;
}

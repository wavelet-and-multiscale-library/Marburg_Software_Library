#define TWO_D
#define PARALLEL 0

#define QUARKLET
#undef AGGREGATED

#undef DYADIC
#undef TRIVIAL
#define ENERGY

#define _DIM 2
#define JMAX 8
#define PRIMALORDER 3
#define DUALORDER   3
#define TWO_D

#ifdef QUARKLET
#define _WAVELETTL_USE_TFRAME 1
#define PMAX 0
#define FRAME
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
#include <adaptive/apply.h>
#ifdef AGGREGATED
#include <simple_elliptic_equation.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/cached_problem.h>
#include <frame_support.h>
#include <steepest_descent.h>
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
#include <adaptive/steepest_descent_ks.h>

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

class myRHS
  : public Function<2,double>
{
public:
  virtual ~myRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    CornerSingularityRHS csrhs(Point<2>(0,0), 0.5, 1.5);
    return /*2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])+*/5*csrhs.value(p);
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

int main()
{
  
  cout << "Testing class Adaptive SD..." << endl;
  
  const int DIM = 2;
  const int jmax = JMAX;
  
  
  
  double epsilon = 1e-6;
  
#ifdef QUARKLET
  const int pmax = PMAX;
#endif
//  const int jmax1d=jmax-j0;

  const int d  = PRIMALORDER;
  const int dT = DUALORDER;
  
  
  
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
  A(0,0) = /*2.*/ 1.25;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = /*-1.*/ -0.25;
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
  
//  myRHS rhs1;
//    
//  PoissonBVP<DIM> poisson(&rhs1);
  
  PoissonBVP<DIM> poisson(&singRhs);
  //PoissonBVP<DIM> poisson(&const_fun);

  // BiharmonicBVP<DIM> biahrmonic(&singRhs);
#ifdef AGGREGATED
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
//  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson/*, 5.0048, 1.0/0.01*/);
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson,1.,1.);
  
  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  

  //steepest_descent_SOLVE(problem, epsilon, u_epsilon, approximations);
  
  clock_t tstart_sd, tend_sd;
  double time_sd;
  tstart_sd = clock();
    steepest_descent_SOLVE(problem, epsilon, approximations);
    tend_sd = clock();
  time_sd = (double)(tend_sd-tstart_sd)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed for SD AGGREGATED: " << time_sd << " seconds" << endl;
  
//  abort();

#endif
  #if 1  
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
    LDomainFrameEquation<Frame1d,Frame2D> discrete_poisson(&poisson, &frame, true);
//    discrete_poisson.set_jpmax(jmax,pmax);
//    Frame2D frame = discrete_poisson.frame();
//    
    
    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame2D> > problem(&discrete_poisson, 1., 1.);
    
    InfiniteVector<double, int> u_epsilon_int;
    
    clock_t tstart_sd, tend_sd;
  double time_sd;
  tstart_sd = clock();
    steepest_descent_ks_QUARKLET_SOLVE(problem, epsilon, u_epsilon_int, tensor_simple, 2, 2);
    tend_sd = clock();
  time_sd = (double)(tend_sd-tstart_sd)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed for SD QUARKLET: " << time_sd << " seconds" << endl;
  
//  abort();
  InfiniteVector<double, Index> u_epsilon;
  for (typename InfiniteVector<double,int>::const_iterator it(u_epsilon_int.begin()),
 	   itend(u_epsilon_int.end()); it != itend; ++it){
        u_epsilon.set_coefficient(*(frame.get_quarklet(it.index())), *it);
  }
  
#endif
  
  
  

  
#ifdef AGGREGATED
  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&problem,-1);



  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 6);//expand in primal basis
  

  

  std::ofstream ofs5("adaptive_solution_out.m");
  matlab_output(ofs5,U);
  ofs5 << "title('Adaptive solution to test problem (Aggregated Frame)');" << endl;
  
  ofs5.close();
  
#endif
  
#ifdef QUARKLET
  //  cout << u << endl;
    u_epsilon.scale(&discrete_poisson,-1);
//  cout << "u_new: " << endl << u_epsilon << endl;
    Array1D<SampledMapping<2> > eval(3);
    eval=frame.evaluate(u_epsilon,6);
    std::ofstream os2("adaptive_solution_out.m");
    os2 << "figure\n" << endl;
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    } 
    os2 << "title('Adaptive solution to test problem (Quarklet Frame), jmax=" << JMAX << ", pmax=" << PMAX << "');" << endl;
    os2 << "view(30,55);"<<endl;
    os2 << "hold off;" << endl;
//    os2 << "figure;" << endl;
    os2.close(); 
    
    //new coefficients plot
    u_epsilon.scale(&problem, 1);
    std::ofstream coeff_stream;
    coeff_stream.open("adaptive_coefficients.m");
    //coeff_stream2 << "figure;" << endl;
    MultiIndex<int,DIM> pstart;
    MultiIndex<int,DIM> jstart;// start=basis1.j0();
    MultiIndex<int,DIM> estart;
    for (InfiniteVector<double, Index>::const_iterator it(u_epsilon.begin()); it != u_epsilon.end(); ++it){
        Index lambda=it.index();
        if(!(lambda.j()==jstart && lambda.e()==estart)){
            //cout <<lambda.p()[0]<<lambda.p()[1]<< lambda.j()[0]-1+lambda.e()[0]<<lambda.j()[1]-1+lambda.e()[1] << endl;
            jstart=lambda.j();
            estart=lambda.e();
            plot_indices(&frame, u_epsilon, coeff_stream, lambda.p(), lambda.j(), lambda.e(),"(jet)", false, true, -6);
            //coeff_stream2 << "title('solution coefficients') " << endl;
            //coeff_stream2 << "title(sprintf('coefficients on level (%i,%i)',"<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"));"<<endl;
            coeff_stream<<"print('-djpg',sprintf('coeffs%i%i%i%i.jpg',"<<lambda.p()[0]<<","<<lambda.p()[1]<<","<<lambda.j()[0]-1+lambda.e()[0]<<","<<lambda.j()[1]-1+lambda.e()[1]<<"))"<<endl;
        }
    }
    coeff_stream.close();
    cout << "coefficients plotted"<<endl;
#endif
#endif
  cout << "done plotting adaptive solution" << endl;
   return 0;
}

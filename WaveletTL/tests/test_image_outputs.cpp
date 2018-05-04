#include <iostream>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1



#define DYADIC
#undef ENERGY
#undef TRIVIAL

//#undef SD
//#define RICHARDSON
//#undef CONGRAD
//
//#undef SHRINKAGE


#define JMAX 8
#define PMAX 2
#define _DIM 1

#define PRIMALORDER 3
#define DUALORDER   3
#include <fstream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>


#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/i_q_index.h>
#include <interval/i_q_indexplot.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#include <interval/pq_frame.h>

#undef BASIS
#define FRAME
#undef RFRAME
#undef CDF


#include <Rd/cdf_basis.h>
#include <Rd/quarklet_frame.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/TestProblem.h>
#include <galerkin/TestFunctions.h>
#ifdef BASIS
#include <interval/p_expansion.h>
#include <galerkin/cached_problem.h>
#endif
#ifdef FRAME
#endif



using namespace std;
using namespace WaveletTL;







int main()
{
    
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

    const unsigned int testcase=1;
    TestProblem<testcase> T;    
      

      
  const int d  = PRIMALORDER;
  const int dT = DUALORDER; // be sure to use a continuous dual here, otherwise the RHS test will fail
  
  const int jmax = JMAX;
  const int pmax = PMAX;
    
    

  
  
#ifdef BASIS
  typedef PBasis<d,dT> Basis;
  Basis basis(1,1);
  typedef Basis::Index Index;
  const char* basis_type = "Primbs basis";
  basis.set_jmax(jmax);
#endif
  
#ifdef CDF
  typedef CDFBasis<d,dT> Basis;
  Basis basis;
  typedef Basis::Index Index;
#endif
  
#ifdef RFRAME
  typedef QuarkletFrame<d,dT> Basis;
  Basis basis;
  typedef Basis::Index Index;
#endif  
  
#ifdef FRAME
  typedef PQFrame<d,dT> Basis;
  bool dirichlet_left=1;
  bool dirchlet_right=1;
  Basis basis(dirichlet_left, dirchlet_right, true);
  typedef Basis::Index Index;
//  const char* basis_type = "Primbs quarklet frame";
  basis.set_jpmax(jmax,pmax);
  bool primal = true;
  
#endif
  cout << "setup equation.." << endl;
  SturmEquation<Basis> eq(T, basis);
  cout << "end setup equation!" << endl;  

  cout << "M_j0: " << endl << basis.get_Mj0() << endl;
  cout << "M_j1: " << endl << basis.get_Mj1() << endl;
  cout << "tildeM_j0: " << endl << basis.get_Mj0T() << endl; 
  
  cout << "- writing M_j0 to the file refinementmatrix.m ..." << endl;
   std::ofstream Astream("Image_outputs/refinementmatrix.m");
   Astream << "A=";
   print_matrix(basis.get_Mj0(), Astream);
   Astream << ";" << endl;
   Astream.close();
   cout << "  ... done!" << endl;
//  abort();
 

#ifdef RFRAME
  const int p = 1;
  const int j = 0;
  const int e = 1;
  const int k = 0;
  Index lambda(p,j,e,k);
  char plotnameInd[128];
  const char* filenameInd = "Image_outputs/indexplot.m";
  if(e==1)
    sprintf(plotnameInd, "%s%d%s%d%s%d%s", "rquarklet_m_", d, "_mT_" , dT, "_p_" , p, ".png");
  else
    sprintf(plotnameInd, "%s%d%s%d%s", "rquark_m_", d, "_p_" , p, ".png");
  
  std::ofstream ind_stream1(filenameInd);
  const int leftx = -2;
  const int rightx= 3;
  Grid<1> grid(leftx, rightx, 1024);
  Array1D<double> values;
  basis.evaluate(0, lambda, grid.points(), values);
  SampledMapping<1>(grid,values).matlab_output(ind_stream1);
  ind_stream1 << "figure;\nplot(x,y,'-');\nset(gca,'fontsize',20)\nxlim([" <<leftx << " " << rightx << "])"
          ";\nylim([-inf inf]);\nprint -dpng " << plotnameInd  << endl;
  ind_stream1.close();
#endif
  
  
  
#ifdef CDF
  const int j = 0;
  const int e = 1;
  const int k = 0;
  Index lambda(j,e,k);
  char plotnameInd[128];
  const char* filenameInd = "Image_outputs/indexplot.m";
  if(e==1)
    sprintf(plotnameInd, "%s%d%s%d%s", "cdf_motherwavelet_m_", d, "_mT_" , dT, ".png");
  else
    sprintf(plotnameInd, "%s%d%s", "cdf_generator_m_", d, ".png");
  
  std::ofstream ind_stream1(filenameInd);
  const int leftx = -4;
  const int rightx= 5;
  Grid<1> grid(leftx, rightx, 1024);
  Array1D<double> values;
  basis.evaluate(0, lambda, grid.points(), values);
  SampledMapping<1>(grid,values).matlab_output(ind_stream1);
  ind_stream1 << "figure;\nplot(x,y,'-');\nset(gca,'fontsize',20)\nxlim([" <<leftx << " " << rightx << "]);;\nprint -dpng " << plotnameInd  << endl;
  ind_stream1.close();
#endif
  
#ifdef FRAME
  const int p = 1;
  const int j = basis.j0();
  const int e = 1;
  const int k = 5;
  
  InfiniteVector<double,Index> indplot;
//  const char* filenameInd = "Image_outputs/indexplot.m";  
//  std::ofstream ind_stream1(filenameInd);
//  ind_stream1 << "figure;" << endl;   
  cout << basis.Nablamin() << ", " << basis.Nablamax(j,p) << endl;
          
//  for(int k = basis.Nablamin();k<=basis.Nablamax(j,p);k++){
  
      
    Index lambda(p,j,e,k,&basis);
    char filenameInd[128], plotnameInd[128];
    int number=lambda.number();
    sprintf(filenameInd, "%s%d%s%d%s%d%s%d%s%d%s", "Image_outputs/indexplot_bc_", dirichlet_left, "_case_" , primal , "_m_", d, "_mT_" , dT, "_number_" , number , ".m");
    sprintf(plotnameInd, "%s%d%s%d%s%d%s%d%s%d%s", "quarklet_bc_", dirichlet_left, "_case_" , primal , "_m_", d, "_mT_" , dT, "_number_" , number , ".png");
    std::ofstream ind_stream1 (filenameInd);
      
    

    /* plot solution*/ 
//    InfiniteVector<double,Index> indplot;
//    onst char* filenameInd = "Image_outputs/indexplot.m";  
  ////  u.scale(&eq, -1);
  //  //u.set_coefficient(Index(1,3,0,1,&basis), u[Index(1,3,0,1,&basis)]+1);
    indplot.clear();
    indplot.set_coefficient(lambda, 1);
  ////  basis.reconstruct(w,4,testv);
  ////  cout << "reconstruction: " << testv << endl;
  ////  cout << "scaled solution: " << endl << u << endl;
    SampledMapping<1> indsamp(evaluate(basis, indplot, primal, 8));
  //  SampledMapping<1> s(evaluate(eq.basis(), Index(1,3,0,1,&basis), true, 7));
    indsamp.matlab_output(ind_stream1);
    ind_stream1 << "figure;\nplot(x,y,'-');\nset(gca,'fontsize',20)\nprint -dpng " << plotnameInd << endl;
    ind_stream1.close();
    
    
//    if(k<=basis.DeltaLmin()+1 || k>=basis.DeltaRmax(j,p)-1)
//        ind_stream1 << "\nplot(x,y,'--');\nhold on\nclear x;\nclear y;" << endl;
//    else
//        ind_stream1 << "\nplot(x,y,'-');\nhold on\nclear x;\nclear y;" << endl;            
    cout << "Index " << lambda << " plotted" << endl;
//  }
//  ind_stream1 << "\nhold off" << endl; 
//  ind_stream1.close();
#endif 
  
#if 0
  int mysize=PMAX;
  Vector<double> Squarenorm;
  Squarenorm.resize(mysize);
  for(int p=0; p<mysize; p++){
      Index lambdag(p,basis.j0(),0,3,&basis);
      Squarenorm[p]=eq.a(lambdag,lambdag);
  }
  
  cout << "- writing Squarenorm to the file squarenorm.m ..." << endl;
   std::ofstream sqstream("squarenorm.m");
   sqstream << "sqnorm=";
   print_vector(Squarenorm, sqstream);
   sqstream << ";" << endl;
   sqstream.close();
   cout << "  ... done!" << endl;

  
//  abort();
 


  
  //optional outputs
//  Matrix<double> evecs;
//  Vector<double> evals;
//  SymmEigenvalues(A,evals,evecs);
//  cout<< "Eigenwerte A: " << evals << endl;
//  cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;
//  cout << "- right hand side: " << b << endl << endl;
//  cout << "  point values of Galerkin solution: " << u_values << endl;
//  cout << "  point values of exact solution: " << uexact_values << endl;
//  cout << " pointwise error (exact solution): " << u_values-uexact_values << endl;
//  cout << "solution: " << endl << u << endl;
//  cout << "coarsed solution: " << endl << v << endl;
  cout << "Number of iterations: " << iterations << endl;

  
  //plots
  
  
  
  /* plot solution coefficients */
  // output for Frame has to be renewed
  
#ifdef FRAME

    
  
 /* plot solution*/ 
  char filenameSolution1[128];
  sprintf(filenameSolution1, "%s%d%s", "quarklet_with_coordinates.m");
//  const char* filenameSolution1 = "sturm_bvp_solution.m";  
  u.scale(&eq, -1);
  //u.set_coefficient(Index(1,3,0,1,&basis), u[Index(1,3,0,1,&basis)]+1);
//  w.set_coefficient(Index(0,3,1,3,&basis), 1);
//  basis.reconstruct(w,4,testv);
//  cout << "reconstruction: " << testv << endl;
//  cout << "scaled solution: " << endl << u << endl;
  SampledMapping<1> s(evaluate(eq.basis(), u, true, 7));
//  SampledMapping<1> s(evaluate(eq.basis(), Index(1,3,0,1,&basis), true, 7));
  std::ofstream u_stream1(filenameSolution1);
  s.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('Sturm bvp: solution to test problem ("
            << basis_type << "), " << "pmax= " << pmax << ", jmax= " << jmax << ", d= " << d << "');" << endl;
  u_stream1.close();
  
  
  
  
#endif
  
#endif


  return 0;
}

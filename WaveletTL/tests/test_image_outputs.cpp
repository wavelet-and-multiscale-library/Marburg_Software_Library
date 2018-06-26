#include <iostream>
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1


#define JMAX 8
#define PMAX 0
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
#ifdef BASIS
#include <interval/p_expansion.h>
#endif




using namespace std;
using namespace WaveletTL;







int main()
{
    
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

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
  bool dirichlet_left=0;
  bool dirchlet_right=0;
  Basis basis(dirichlet_left, dirchlet_right, true);
  
  typedef Basis::Index Index;
  basis.set_jpmax(jmax,pmax);
  bool primal = true;
 
#endif 

#ifdef RFRAME
  const int p = 2;
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
  const int p = 2;
  const int j = basis.j0();
  const int e = 0;
  const int k = 5;
  
  InfiniteVector<double,Index> indplot;
  cout << basis.Nablamin() << ", " << basis.Nablamax(j,p) << endl;
          
      
    Index lambda(p,j,e,k,&basis);
    char filenameInd[128], plotnameInd[128];
    int number=lambda.number();
    sprintf(filenameInd, "%s%d%s%d%s%d%s%d%s%d%s", "Image_outputs/indexplot_bc_", dirichlet_left, "_case_" , primal , "_m_", d, "_mT_" , dT, "_number_" , number , ".m");
    sprintf(plotnameInd, "%s%d%s%d%s%d%s%d%s%d%s", "quarklet_bc_", dirichlet_left, "_case_" , primal , "_m_", d, "_mT_" , dT, "_number_" , number , ".png");
    std::ofstream ind_stream1 (filenameInd);
      
    

    /* plot solution*/ 
   indplot.clear();
    indplot.set_coefficient(lambda, 1);
    SampledMapping<1> indsamp(evaluate(basis, indplot, primal, 8));
    indsamp.matlab_output(ind_stream1);
    ind_stream1 << "figure;\nplot(x,y,'-');\nset(gca,'fontsize',20)\nprint -dpng " << plotnameInd << endl;
    ind_stream1.close();
    
    
    cout << "Index " << lambda << " plotted" << endl;

#endif 
  


  return 0;
}

// more output for the cached problemm (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// Code in CDD1.h uses APPLY or APPLY_CAched, depending on this makro
#define _WAVELETTL_USE_TBASIS 1
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0
#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
//#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>

#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>

#include <geometry/sampled_mapping.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tbasis_indexplot.h>

/*
 * for comparisn with the isotropic case
 */
//#include <cube/cube_basis.h>
//#include <galerkin/cube_equation.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>

// incompatible with tbasis_indexplot:
//#include <cube/cube_indexplot.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};

/*
  Some test problems for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: u(x,y) = x(1-x)y(1-y), -Delta u(x,y) = 2(x(1-x)+y(1-y))
  2: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y)
  3: u(x,y) = x(1-x)^2y^2(1-y), -Delta u(x,y)= 4*(1-x)*y^2*(1-y)-2*x*y^2*(1-y)-2*x*(1-x)^2*(1-y)+4*x*(1-x)^2*y
*/
template <unsigned int N>
class TestRHS
  : public Function<2,double>
{
public:
  virtual ~TestRHS() {};
  double value(const Point<2>& p, const unsigned int component = 0) const {
    switch(N) {
    case 2:
      return
	(200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
	* exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
      break;
    case 3:
      return
	4*(1-p[0])*p[1]*p[1]*(1-p[1])
	- 2*p[0]*p[1]*p[1]*(1-p[1])
	- 2*p[0]*(1-p[0])*(1-p[0])*(1-p[1])
	+ 4*p[0]*(1-p[0])*(1-p[0])*p[1];
      break;
    case 1:
    default:
      return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

template<class T> struct index_cmp2 {
index_cmp2(const T arr) : arr(arr) {}
bool operator()(const size_t a, const size_t b) const
{ return arr[a] < arr[b]; }
const T arr;
};


int main()
{
  const int d  = 3;
  const int dT = 3;
  const int dim = 2;
  const int offset = 0; // Radius for the index set Lambda
  
  //typedef DSBasis<d,dT> Basis1d;
  typedef PBasis<d,dT> Basis1d;
  typedef TensorBasis<Basis1d,dim> TBasis;
  typedef TBasis::Index Index;
  typedef Index::level_type index_lt;
  
  //typedef TensorEquation<Basis1d,dim,TBasis>::Index EqIndex;
  //typedef CachedTProblem<TensorEquation<Basis1d,dim,TBasis> >::Index CachIndex;

  // uncoment for constant RHS ...
  //ConstantFunction<dim> constant_rhs(Vector<double>(1, "2.0"));
  //PoissonBVP<dim> poisson(&constant_rhs);

  // uncomment in 2 space dimensions
  TestRHS<2> rhs;
  PoissonBVP<dim> poisson(&rhs);

  FixedArray1D<bool,(2*dim)> bc;
  if (dim==1)
  {
      bc[0] = bc[1] = true;
  }
  else
  {
      bc[0] = bc[1] = bc[2] = bc[3] = true;
  }

  clock_t tstart, tend;
  double time;

  cout << "main:: setting up equation" << endl;
  TensorEquation<Basis1d,dim,TBasis> eq(&poisson, bc);
  eq.set_jmax(multi_degree(eq.basis().j0())+offset); //calls setup_full_collection

  cout << "main:: setting up problem" << endl;
  CachedTProblem<TensorEquation<Basis1d,dim,TBasis> > ctproblem(&eq);

  tstart = clock();
  
  cout << "Compute index Set up to the wavelet";
  
  
  Index tempindex(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));
  cout << tempindex<< " number = "<<tempindex.number()<<endl;  
  
  set<Index> Lambda;
  for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == itend) break;
  }

#if 0
  cout << "testing the index set" <<endl;
  cout << "j0 = " << (ctproblem.basis().j0()) << " multidegree = " << (multi_degree(ctproblem.basis().j0())) << endl;
  set<Index> Lambda2;
   for (Index lambda ( eq.basis().first_generator() ), itend(eq.basis().last_wavelet(multi_degree(eq.basis().j0())+offset));; ++lambda) {
    Lambda2.insert(lambda);
    if (lambda == itend) break;
  }

  Index itend1(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset)), itend2(eq.basis().last_wavelet(multi_degree(eq.basis().j0())+offset));
  if (itend1 != itend2)
  {
      cout << itend1 << " , " << itend2;
      cout << " End of index sets is not equal!" << endl;
  }

  for (set<Index>::iterator it1(Lambda.begin()),it2(Lambda2.begin()), itend(Lambda.end()); it1!=itend;it1++, it2++)
  {
      //cout << (*it1) << " , " << (*it2);
      if (*it1 != *it2)
      {
          cout << (*it1) << " , " << (*it2) << " warning: should be equal"<<endl;
      }
      //cout << endl;
  }
#else
  cout << "not testing the index set" << endl;
#endif

#if 0
  cout << "testing D ..."<<endl;;
  for (set<Index>::iterator it(Lambda.begin()), itend(Lambda.end()); it != itend; it++)
  {
      if (abs(ctproblem.D(*it)-eq.D(*it)) != 0)
          cout << "  Problem for lambda = "<<(*it)<<endl;
  }
  cout << " done"<<endl;
#else
  cout << "skipping testing of D"<<endl;
#endif

  double ctmax(0),eqmax(0);
  Index it1max(0,&eq.basis()),it2max(0,&eq.basis());

#if 0
  cout << "testing symmetry of cached a ..."<<endl;;

  double cttemp1,cttemp2;
  for (set<Index>::iterator it1(Lambda.begin()), itend(Lambda.end()); it1 != itend; it1++)
  {
      //cout << "  doing row it1 = " << (*it1).number() <<endl;
      for (set<Index>::iterator it2(Lambda.begin()); it2 != itend; it2++)
      {
          cttemp1 = ctproblem.a(*it1,*it2); //, eqtemp1(eq.a(*it1,*it2));
          cttemp2 = ctproblem.a(*it2,*it1); //, eqtemp2(eq.a(*it2,*it1));
          if (abs(cttemp1-cttemp2) != 0)
          {
              if (abs(cttemp1-cttemp2)>abs(ctmax-eqmax))
              {
                  ctmax=cttemp1;eqmax=cttemp2;
                  it1max=*it1;it2max=*it2;
              }
              cout << "  Problem for (lambda1,lambda2) =  ("<< (*it1).number() << " = "<<(*it1)<<" , "<< (*it2).number() << " = "<<(*it2)<<" )" <<endl;
              cout << "  difference = "<< (abs(cttemp1-cttemp2))<< " ct.a(i,j) = " <<cttemp1<<" ct.a(j,i) = "<<cttemp2<<endl;
          }
      }
  }
  if (abs(ctmax)+abs(eqmax)> 0)
  {

      cout << "  largest difference for (lambda1,lambda2) =  ("<<(it1max)<<" , "<<(it2max)<<" )"<<endl;
      cout << "  difference = "<< (abs(ctmax-eqmax))<< " ct.a(i,j) = " <<ctmax<<" ct.a(j,i) = "<<eqmax<<endl;
  }
  cout << " done"<<endl;
#else
  cout << "skipping symmetry test for cached a" << endl;
#endif

#if 0
  cout << "comparing cached and uncached a ..."<<endl;;
  double cttemp3, eqtemp3;
  for (set<Index>::iterator it1(Lambda.begin()), itend(Lambda.end()); it1 != itend; it1++)
  {
      for (set<Index>::iterator it2(Lambda.begin()); it2 != itend; it2++)
      {
          cttemp3 = ctproblem.a(*it1,*it2);
          eqtemp3 = eq.a(*it1,*it2);
          if (abs(cttemp3-eqtemp3) != 0)
          {
              if (abs(cttemp3-eqtemp3)>abs(ctmax-eqmax))
              {
                  ctmax=cttemp3;eqmax=eqtemp3;
                  it1max=*it1;it2max=*it2;
              }
              // uncomment if errors above tolerance are present
              //cout << "  Problem for (lambda1,lambda2) =  ("<< (*it1).number() << " = "<<(*it1)<<" , "<< (*it2).number() << " = "<<(*it2)<<" )" <<endl;
              //cout << "  difference = "<< (abs(cttemp-eqtemp))<< " ct.a = " <<cttemp<<" eq.a = "<<eqtemp<<endl;
          }
      }
  }
  if (abs(ctmax)+abs(eqmax)> 1e-16)
  {

      cout << "  largest difference for (lambda1,lambda2) =  ("<<(it1max)<<" , "<<(it2max)<<" )"<<endl;
      cout << "  difference = "<< (abs(ctmax-eqmax))<< " ct.a = " <<ctmax<<" eq.a = "<<eqmax<<endl;
  }
  else
  {
      cout << "  differences below tolerance of 1e-16" << endl;
  }
  cout << " done"<<endl;
#else
  cout << "skipping cache test for a" << endl;
#endif

#if 0
  // useful if test of cached and uncached a found an error
  cout << "testing cached and uncached a for specific indices" << endl;
  Index ind1(15,&ctproblem.basis());
  Index ind2(14,&ctproblem.basis());
  cout << "ind1 = " << ind1 << " ind2 = " << ind2 << endl;
  ctproblem.a(ind1,ind2); // for breakpoints
  eq.a(ind1,ind2); // for breakpoints
  cout << "ct.a(i,j) = " << ctproblem.a(ind1,ind2) << " ct.a(j,i) = " << ctproblem.a(ind2,ind1) << endl;
  cout << "eq.a(i,j) = " << eq.a(ind1,ind2) << " eq.a(j,i) = " << eq.a(ind2,ind1) << endl;
#else
  cout << "not testing a for specific indices" << endl;
#endif

#if 0
  cout << "testing calls to the cache ..."<<endl;
  index_lt it1_level,first_level;
  int it1_blocknumber;
  double ctemp1,ctemp2,ctemp3;
  for (set<Index>::iterator it1(Lambda.begin()), itend(Lambda.end()); it1 != itend; it1++)
  {
      for (set<Index>::iterator it2(Lambda.begin()); it2 != itend; it2++)
      {
          //cout << "lambda1 = " << (*it1) << " lambda2 = " << (*it2)<<endl;
          
          it1_level = (*it1).j();
          first_level = ctproblem.basis().j0();
          for (int k=0;k<dim;k++)
          {
              it1_level[k] = it1_level[k]-first_level[k];
          }
          it1_blocknumber = it1_level.number();

          //ctproblem.entries_cache[(*it2).number()][it1_blocknumber][(*it1).number()];

          ctemp1 = ctproblem.a(*it1,*it2);
          ctemp2 = ctproblem.a(*it1,*it2);
          ctemp3 = ctproblem.entries_cache[(*it2).number()][it1_blocknumber][(*it1).number()];
          if ((abs(ctemp1-ctemp2) != 0) || (abs(ctemp1-ctemp3) != 0) || (abs(ctemp2-ctemp3) != 0))
          {
              //if (abs(cttemp-eqtemp)>abs(ctmax-eqmax))
              //{
              //    ctmax=cttemp;eqmax=eqtemp;
              //    it1max=*it;it2max=*it2;
              //}
              cout << "  Problem for (lambda1,lambda2) =  ("<<(*it1)<<" , "<<(*it2)<<" )"<<endl;
              cout << "  temp1 = " <<ctemp1<<" temp2 = "<<ctemp2<<" temp3 = "<<ctemp3<<endl;
          }
      }
  }
  //cout << "  largest difference for (lambda1,lambda2) =  ("<<(it1max)<<" , "<<(it2max)<<" )"<<endl;
  //cout << "  difference = "<< (abs(ctmax-eqmax))<< " ct.a = " <<ctmax<<" eq.a = "<<eqmax<<endl;
  cout << " done"<<endl;
#else
  cout << "skipping testing of cache calls"<<endl;
#endif


#if 0
  cout << "testing norm_A, norm_Ainv ..."<<endl;
  
  cout << "compute eq.norm_A()" << endl;
  double eqnormA (eq.norm_A());  
  cout << "eq.norm_A() = " << eqnormA << endl;
  //if (abs(ctproblem.norm_A()-3.40341)>1e-6)
  if (abs(ctproblem.norm_A() - eqnormA) > 1e-6)
  {
      //cout << "  HIGH ERROR! difference = "<< (abs(ctproblem.norm_A() - 3.40341))<< " ct.norm_A = " <<ctproblem.norm_A()<<" eq.norm_A = 3.40341" <<endl;
      cout << "  HIGH ERROR! difference = "<< (abs(ctproblem.norm_A() - eqnormA))<< " ct.norm_A = " <<ctproblem.norm_A()<<" eq.norm_A = " << eqnormA << endl;
  }
  else
  {
      cout << "  difference below tolerance of 1e-6. ct.norm_A = "<< ctproblem.norm_A()<<endl;
  }
  cout << "compute eq.norm_Ainv()." << endl;
  double eqnormAinv (eq.norm_Ainv());
  cout << "eq.norm_Ainv() = " << eqnormAinv << endl;
  //if (abs(ctproblem.norm_A()-3.40341)>1e-6)
  if (abs(ctproblem.norm_Ainv() - eqnormAinv) > 1e-6)
  {
      //cout << "  HIGH ERROR! difference = "<< (abs(ctproblem.norm_A() - 3.40341))<< " ct.norm_A = " <<ctproblem.norm_A()<<" eq.norm_A = 3.40341" <<endl;
      cout << "  HIGH ERROR! difference = "<< (abs(ctproblem.norm_Ainv() - eqnormAinv))<< " ct.norm_Ainv = " <<ctproblem.norm_Ainv()<<" eq.norm_Ainv = " << eqnormAinv << endl;
  }
  else
  {
      cout << "  difference below tolerance of 1e-6. ct.norm_Ainv = "<< ctproblem.norm_Ainv()<<endl;
  }

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
  cout << "skipping testing of norm_A"<<endl;
#endif

#if 0
  // computation of normA with a Lanczos iteration with tolerance tol produces only results up to tolerance tol !!
  cout << "computing normA in the main routine because of high differences in the above test" << endl;
  // For the cached equation
  double normA1;
  set<Index> Lambda_NormA1;
  for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda)
  {
    Lambda_NormA1.insert(lambda);
    if (lambda == itend) break;
  }
  SparseMatrix<double> A_Lambda1;
  setup_stiffness_matrix(ctproblem, Lambda_NormA1, A_Lambda1);
  double help1;
  unsigned int iterations1;
  LanczosIteration(A_Lambda1, 1e-6, help1, normA1, 200, iterations1);
  

  // For the uncached equation:
  double normA2;
  std::set<Index> Lambda_NormA2;
  for (Index lambda = eq.basis().first_generator(), itend = eq.basis().last_wavelet(multi_degree(eq.basis().j0())+offset);; ++lambda) 
  {
      Lambda_NormA2.insert(lambda);
      if (lambda == itend) break;
  }
  SparseMatrix<double> A_Lambda2;
  setup_stiffness_matrix(eq, Lambda_NormA2, A_Lambda2);
  double help2;
  unsigned int iterations2;
  LanczosIteration(A_Lambda2, 1e-8, help2, normA2, 200, iterations2);

  //cout << normA1;
  //cout << normA2;
  cout << "Differences between cached and uncached evaluation:" << endl;
  cout << "  difference = "<< (abs(normA1-normA2))<< " normA1 = " <<normA1<<" normA2 = "<<normA2<<endl;
  cout << "Differences between call of norm_A() and computation here:" << endl;
  cout << "  difference cached = "<< (abs(ctproblem.norm_A()-normA1))<< "  difference uncached = "<< (abs(eq.norm_A()-normA2)) << endl;
  cout << "Calling norm_A_test ..." << endl;
  eq.norm_A_test(A_Lambda2,normA2);
#endif

#if 0
  cout << "for the statistic: " << endl;
  cout << "compute ctproblem.normA" << endl;
  ctproblem.norm_A();
  cout << "... done. computing ctproblem.normAinv"<<endl;
  ctproblem.norm_Ainv();
  cout << "... done." << endl;
  cout << "normA = " << ctproblem.norm_A() << " normAinv = " << ctproblem.norm_Ainv() << " kondition = " << ctproblem.norm_A()*ctproblem.norm_Ainv() << endl;
#endif

#if 0 // I forgot to implement&call set_jmax for thee underlying equation!
  // non the less: could be moved to test basis:
  cout << "first_gen =("<< ctproblem.basis().first_generator().number() << ")="<< ctproblem.basis().first_generator() << endl;
  Index temp_index(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));
  cout << "last_wav  =("<< temp_index.number() << ")="<< temp_index << endl;

  ctproblem.basis().get_wavelet(30);
  ctproblem.basis().get_wavelet(31);

  for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda)
  {

      //if ( (ctproblem.basis().get_wavelet(lambda.number()))->number() != lambda.number())
      {
          cout << " get_wavelet = " <<  (ctproblem.basis().get_wavelet(lambda.number()))->number() << " lambda.number = " << lambda.number() << endl;
      }
      if (lambda == itend) break;
  }
#endif

#if 0
  cout << "testing apply"<<endl;

  set<Index> Lambda_apply;
  set<int> window;
  for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda)
  {
      window.insert(lambda.number());
      Lambda_apply.insert(lambda);
      if (lambda == itend) break;
  }
  SparseMatrix<double> A_apply;
  //A_apply.resize(Lambda_apply.size(), Lambda_apply.size());
  cout << "  setting up stiffness matrix ..." << endl;
  setup_stiffness_matrix(ctproblem, Lambda_apply, A_apply);
  cout << "  done" << endl;
  Vector<double> x(Lambda_apply.size()), res_apply1(Lambda_apply.size()), res_apply2(Lambda_apply.size());
  Vector<double> apply_diff;
  double temp_apply1(0), temp_apply2(0);
  ctmax = 0;
  eqmax = 0;
  for (set<Index>::iterator it1(Lambda_apply.begin()), itend(Lambda_apply.end());it1!=itend;it1++)
  {
      x[(*it1).number()] = 1;
      A_apply.apply(x,res_apply1);
      ctproblem.apply(window,x,res_apply2);
      apply_diff = res_apply1;
      apply_diff.subtract(res_apply2);
      temp_apply1 = linfty_norm(apply_diff);
      temp_apply2 = l2_norm(apply_diff);
      if (temp_apply1 > ctmax)
      {
          ctmax = temp_apply1;
          it1max = (*it1);
      }
      if (temp_apply2 > eqmax)
      {
          eqmax = temp_apply2;
          it2max = (*it1);
      }
      if (abs(temp_apply1)+abs(temp_apply2)>1e-15)
      {
          cout << "index (" << (*it1).number() << ") = " << (*it1) << " difference of Ax: linfty =" << temp_apply1 << " l2 norm = " << temp_apply2 << endl;
      }
      /*
      if (temp_apply1 != temp_apply2)
      {
          cout << apply_diff << endl;
      }
       */
      x[(*it1).number()] = 0;
      
  }
  if (abs(ctmax)+abs(eqmax)>1e-16)
  {
      cout << "largest errors in apply:" <<endl;
      cout << "linfty = "<< ctmax << " at index ("<< (it1max.number()) <<") = " << it1max << endl;
      cout << "l2     = "<< eqmax << " at index ("<< (it2max.number()) <<") = " << it2max << endl;
  }
  else
  {
      cout << " all errors below tolerance of 1e-16" << endl;
  }
  cout << " done"<<endl;
#else
  cout << "skipping test of apply"<<endl;
#endif

#if 0 // use only if apply test is active
  cout << "testing apply for specific indices" << endl;
  x.resize(Lambda_apply.size());
  it1max = new Index(14,&eq.basis());
  x[it1max.number()]=1;
  A_apply.apply(x,res_apply1);
  ctproblem.apply(window,x,res_apply2);
  Vector<double> res_apply3(Lambda_apply.size());
  ctproblem.apply(window,x,res_apply3);
  if (res_apply2 != res_apply3)
  {
      cout << " Error: repeated call to apply does not return same value" << endl;
  }
  apply_diff = res_apply1;
  apply_diff.subtract(res_apply2);
  temp_apply1 = linfty_norm(apply_diff);
  temp_apply2 = l2_norm(apply_diff);
  if (temp_apply1+temp_apply2 > 0)
  {
      cout << " Different results for apply! linfty = " << temp_apply1 << " l2 = " << temp_apply2 << endl;
  }
  // compute row by hand:
  double help1(0),help2(0);
  int i =(*Lambda_apply.begin()).number();
  for (set<Index>::iterator it2(Lambda_apply.begin()), itend(Lambda_apply.end());it2!=itend;it2++)
  {
      help1 += ctproblem.a(it1max,(*it2))/ctproblem.D(it1max)/ctproblem.D(*it2)*x[i];
      help2 += A_apply.get_entry(it1max.number(),i)*x[i];
      i++;
  }
  if (abs(help1-help2)>0)
  {
      cout << " defferent row sum for row = (" << it1max << ")=" << it1max << " diff = " << abs(help1-help2) << " ct_sum = " << help1 << " matrix_sum = " << help2 << endl;
  }
  
#endif

#if 0
  cout << "testing add_ball"<<endl;
  set<Index> Lambda_aball;
  //set<int> window;
  for (Index lambda ( ctproblem.basis().first_generator() ), itend(ctproblem.basis().last_wavelet(multi_degree(ctproblem.basis().j0())+offset));; ++lambda)
  {
      //window.insert(lambda.number());
      Lambda_aball.insert(lambda);
      if (lambda == itend) break;
  }
  SparseMatrix<double> A_aball;
  cout << "  setting up stiffness matrix ..." << endl;
  setup_stiffness_matrix(ctproblem, Lambda_aball, A_aball); // uses the cache.
  //setup_stiffness_matrix(eq, Lambda_aball, A_aball);
  cout << "  done" << endl;
  Vector<double> w1(Lambda_aball.size()), w2(Lambda_aball.size()), x_aball(Lambda_aball.size()), y_aball(Lambda_aball.size());
  double temp_linfty(0),temp_l2(0),max_linfty(0),max_l2(0);
  Index max_linfty_pos, max_l2_pos;
  index_lt temp_level;
  int radius, diff;
  bool noerror(true);
  set<Index>::iterator itend(Lambda_aball.end()),itendtemp(Lambda_aball.begin());;
  ++itendtemp;
  //for (set<Index>::iterator it(Lambda_aball.begin()),itend(Lambda_aball.end());it!=itend;++it)
  //for (set<Index>::iterator it(Lambda_aball.begin());it!=itendtemp;++it)
  for (set<Index>::iterator it(Lambda_aball.begin());it!=itend;++it)
  {
      w1.resize(Lambda_aball.size());
      //w2.resize(Lambda_aball.size());
      radius=20;
      //cout << "radius/2 = " << (radius/2) << endl;
      ctproblem.add_ball(*it,w1,radius,1,multi_degree(eq.basis().j0())+offset);
      x_aball[(*it).number()]=1;
      A_aball.apply(x_aball,w2);
      x_aball[(*it).number()]=0;
      y_aball = w1;
      y_aball.subtract(w2);
      temp_linfty = linfty_norm(y_aball);
      temp_l2 = l2_norm(y_aball);
      if (temp_linfty > max_linfty)
      {
          max_linfty = temp_linfty;
          max_linfty_pos = *it;
      }
      if (temp_l2 > max_l2)
      {
          max_l2 = temp_l2;
          max_l2_pos = *it;
      }
      //cout << " iterNr " << (*it).number() << "="<< *it << endl;
      for (set<Index>::iterator it2(Lambda_aball.begin()); it2!=itend;++it2)
      {
          //cout << " (*it2).j() = " << (*it2).j() << endl;
          temp_level = (*it).j();
          //cout << " temp_level = " << temp_level << endl;
          for (int i=0;i<dim;i++)
          {
              temp_level[i]=abs(temp_level[i]-(*it2).j()[i]);
          }
          diff = multi_degree(temp_level);
          if (diff <= radius)
          {
              if (y_aball[(*it2).number()] > 1e-15)
              {
                  noerror=false;
                  
                  cout << "  error! it="<< (*it).number() << "=" << *it << "y_aball[" << (*it2).number() << "=" << *it2;
                  cout << "]="<< y_aball[(*it2).number()] << " w1[]=" << w1[(*it2).number()];
                  cout << " w2[]=" << w2[(*it2).number()] << " A[][]=" << A_aball.get_entry((*it2).number(),(*it).number()) << endl;
              }
          }
          else
          {
              cout << "  shouldn't happen!!" << endl;
          }
      }
      //if (temp_linfty + temp_l2 > 1e-15)
      if (false)
      //if ( ((*it).number() >1184) && (*it).number()<1441)
      {

          cout << "row w1 w2 A_aball column = " << (*it).number() << " = " << (*it) << endl;
          for (set<Index>::iterator it2(Lambda_aball.begin());it2!=itend;++it2)
          {
              cout << " ("<<(*it2).number()<<")="<<(*it2)<<" " << w1[(*it2).number()] << " "<< w2[(*it2).number()] << " " << A_aball.get_entry((*it2).number(),(*it).number()) << endl;
          }
      }
  }
  //if (max(max_linfty,max_l2) > 1e-15)
  if (!noerror)
  {
      cout << "  large errors. linfty = " << max_linfty << " at (" << max_linfty_pos.number() << ")=" << max_linfty_pos << " l2 = " << max_l2 << " at (" << max_l2_pos.number() << ")=" << max_l2_pos << endl;
  }
  else
  {
      cout << "  error below tolerance of 1e-15" << endl;
  }
  //int ro(124),co(125),bl(3);
  //cout << " a specific entry of the cache: " << ctproblem.entries_cache[co][bl][ro] << " with precond: " << ctproblem.entries_cache[co][bl][ro]/ctproblem.D(ctproblem.basis().get_wavelet(co))/ctproblem.D(ctproblem.basis().get_wavelet(ro))<< endl;
#else
  cout << "skipping test of add_ball" << endl;
#endif

#if 0
  cout << "testing tiny tools twotothejhalf " << endl;
  for (unsigned int i=0;i<100;i++)
  {
      cout << "i = " << i << " twotothejhalf, pow = " << twotothejhalf(i) << " " << pow(2,i/2.0) << endl;
  }
#endif

#if 0
  cout << "testing tiny_tools intpower" << endl;
  for (unsigned int dime = 1; dime < 4; dime++)
  {
      for (unsigned int i=0; i< 100; i++)
      {
          if (intpower(i,dime) != pow(i,dime))
          {
              cout << "intpower failed! intpower("<<i<<", "<<dime<<") = " << intpower(i,dime) << " pow("<<i<<", "<<dime<<") = " << pow(i,dime) << endl;
          }
      }
  }
#endif
  
# if 0
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD2_SOLVE ..." << endl;
  tstart = clock();
  cout << "main:: setting up RHS" << endl;
  InfiniteVector<double, Index> F_eta;
  ctproblem.RHS(1e-6, F_eta);
  cout << "l2norm(f) = " << l2_norm(F_eta) << endl;
  cout << "main:: calling norm_Ainv" << endl;
  const double nAinv = ctproblem.norm_Ainv();
  //cout << "main:: calling norm_A" << endl;
  //const double nA = ctproblem.norm_A();
  //cout << "main:: calling l2_norm" << endl;
  //const double l2normF = l2_norm(F_eta);
  //cout << "l2norm(f_eta) = " << l2normF << " normA = "<< nA << " normAinv = " << nAinv << " kond = "<<(nA*nAinv)<<endl;
  double nu = ctproblem.norm_Ainv() * l2_norm(F_eta);
  cout << "nu = " << nu << endl;

  cout << "degrees_of_freedom = " << ctproblem.basis().degrees_of_freedom() << endl;

  //tend = clock();
  //time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  //cout << "  ... done, time needed: " << time << " seconds" << endl;

  //tstart = clock();
  
  InfiniteVector<double, Index> u_epsilon;
  cout << "calling CDD2_SOLVE_TENSOR" << endl;
  unsigned int maxlevel = ctproblem.basis().get_jmax();

  //nu = min(nu,2.0);
  //nu = 3000;
  CDD2_SOLVE(ctproblem, nu, 1e-1, u_epsilon,maxlevel);
  cout << "  ... done. Saving output in file 'u_adaptCDD2.m'" << endl;
  cout << "main:: u_epsilon.size() = " << u_epsilon.size() << endl;
  u_epsilon.scale(&ctproblem, -1);
  cout << "main:: u_epsilon.size() = " << u_epsilon.size() << endl;
  SampledMapping<dim> s(evaluate(eq.basis(), u_epsilon, true, d+dim+offset+1));
  std::ofstream u_stream("u_adaptCDD2.m");
  s.matlab_output(u_stream);
  u_stream.close();
//  cout << "  ... done" << endl;
  cout << "  ... done. Coarsening with tol = 1e-6. Saving output in file 'u_adaptCDD2c.m'" << endl;

  InfiniteVector<double, Index> u_epsilon_coarse;
  u_epsilon.COARSE(1.0e-6, u_epsilon_coarse);
  cout << "main:: u_epsilon_coarse.size() = " << u_epsilon_coarse.size() << endl;
  SampledMapping<dim> s2(evaluate(eq.basis(), u_epsilon_coarse, true, d+dim+offset+1));
  std::ofstream u_stream2("u_adaptCDD2c.m");
  s2.matlab_output(u_stream2);
  u_stream2.close();
  cout << "  ... done" << endl;

  if (dim == 2)
  {
      std::ofstream plotstream;
      plotstream.open("coefficient_plotCDD2.m");
      cout << "Writing active coefficients to coefficient_plotCDD2.m" << endl;
      plot_indices(&ctproblem.basis(), u_epsilon, offset, plotstream, "jet", true, true);
      plotstream.close();
  }

  //cout << "u_epsilon = " << u_epsilon << endl;

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
  cout << "skipping CDD2 test"<<endl;
#endif
  
# if 1
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD1_SOLVE ..." << endl;
  tstart = clock();
  cout << "main:: setting up RHS" << endl;
  InfiniteVector<double, Index> F_eta3;
  ctproblem.RHS(1e-6, F_eta3);
  cout << "l2norm(f) = " << l2_norm(F_eta3) << endl;
  //cout << "main:: calling norm_Ainv" << endl;
  //const double nAinv = ctproblem.norm_Ainv();
  //cout << "main:: calling norm_A" << endl;
  //const double nA = ctproblem.norm_A();
  //cout << "main:: calling l2_norm" << endl;
  //const double l2normF = l2_norm(F_eta);
  //cout << "l2norm(f_eta) = " << l2normF << " normA = "<< nA << " normAinv = " << nAinv << " kond = "<<(nA*nAinv)<<endl;
  //const double nu = ctproblem.norm_Ainv() * l2_norm(F_eta3);
  //cout << "nu = " << nu << endl;

  cout << "degrees_of_freedom = " << ctproblem.basis().degrees_of_freedom() << endl;
  InfiniteVector<double, Index> u_epsilon3;
  cout << "calling CDD1_SOLVE" << endl;
  unsigned int maxlevel2 = ctproblem.basis().get_jmax();

  CDD1_SOLVE(ctproblem, 1e-2, u_epsilon3,maxlevel2,tensor_simple);
  cout << "  ... done. u_epsilon.size() = " << u_epsilon3.size() << endl;

  u_epsilon3.scale(&ctproblem, -1);
  //cout << " u_epsilon = " << endl << u_epsilon3 << endl;

  //std::ofstream plotstream;
  //plotstream.open("u_adapt_coefficient_plot.m");
  //plot_indices(&ctproblem.basis(), u_epsilon3, maxlevel2, plotstream, "jet", true, true);
  //plotstream.close();
  if (dim == 2)
  {
      std::ofstream plotstream2;
      plotstream2.open("coefficient_plotCDD1.m");
      cout << "Writing active coefficients to coefficient_plotCDD1.m" << endl;
      plot_indices(&ctproblem.basis(), u_epsilon3, offset, plotstream2, "jet", true, true);
      plotstream2.close();
  }

  cout << "saving computed solution to u_adaptCDD1.m" << endl;
  SampledMapping<dim> s3(evaluate(eq.basis(), u_epsilon3, true, d+dim+offset+1));
  std::ofstream u_stream3("u_adaptCDD1.m");
  s3.matlab_output(u_stream3);
  u_stream3.close();

  cout << "  ... done. Coarsening u_epsilon with tol = 1e-6. " << endl;

  InfiniteVector<double, Index> u_epsilon3_coarse;
  u_epsilon3.COARSE(1.0e-6, u_epsilon3_coarse);
  cout << "main:: u_epsilon_coarse.size() = " << u_epsilon3_coarse.size() << endl;
  cout << "Saving output in file 'u_adaptCDD1c.m'" << endl;
  SampledMapping<dim> s4(evaluate(eq.basis(), u_epsilon3_coarse, true, d+dim+offset+1));
  std::ofstream u_stream4("u_adaptCDD1c.m");
  s4.matlab_output(u_stream4);
  u_stream4.close();
  cout << "  ... done" << endl;

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
  cout << "skipping CDD1 test"<<endl;
#endif

  /*
  //Some sorting
  Array1D<unsigned int> test(3), b2(3);
  test[0]=7;test[1]=2;test[2]=3;
  std::vector<int> a;
  a.push_back(7); a.push_back(2); a.push_back(3);
  vector<size_t> b;
  for (unsigned i = 0; i < a.size(); ++i)
  {
      b.push_back(i);
      b2[i]=i;
  }
  cout << "test = [ ";
  for (unsigned int i = 0 ; i< a.size(); i++)
      cout << test[i] << " ";
  cout << "]" << endl;
  cout << "b = [ ";
  for (unsigned int i = 0 ; i< b.size(); i++)
      cout << b[i] << " ";
  cout << "]" << endl;
  // b = [0, 1, 2]
  //sort(b.begin(), b.end(), index_cmp<vector<int>&>(a));
  sort(b.begin(), b.end(), index_cmp2<Array1D<unsigned int>&>(test));
  // b = [0, 2, 1]
  cout << "test = [ ";
  for (unsigned int i = 0 ; i< a.size(); i++)
      cout << test[i] << " ";
  cout << "]" << endl;
  cout << "b = [ ";
  for (unsigned int i = 0 ; i< b.size(); i++)
      cout << b[i] << " ";
  cout << "]" << endl;
  */
  return 0;
}
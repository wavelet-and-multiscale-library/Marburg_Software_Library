/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE

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
    
    cout << "testing L-domain quarklet norm of the laplacian" << endl;
    const int d  = 3;
    const int dT = 3;
    const int dim = 2;
    const int jmax=8;
    const int pmax=2;
    const int offsetj=0;
    const int offsetp=0;
    
    typedef PQFrame<d,dT> Frame1d;
    typedef LDomainFrame<Frame1d> Frame;
    myRHS rhs1;
    
    PoissonBVP<dim> poisson1(&rhs1);
    LDomainFrameEquation<Frame1d,Frame> eq(&poisson1, false);
    eq.set_jpmax(jmax,pmax,false);
    CachedQuarkletLDomainProblem<LDomainFrameEquation<Frame1d,Frame> > ceq(&eq);
    
    ceq.normtest(offsetj,offsetp);
      
    
    
    return 0;
}
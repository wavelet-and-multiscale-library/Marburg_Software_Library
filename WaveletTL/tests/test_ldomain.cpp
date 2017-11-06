/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#undef NONADAPTIVE
#define ADAPTIVE

#define _WAVELETTL_USE_TBASIS 1
#define _WAVELETTL_CDD1_VERBOSITY 1
#define _APPLY_TENSOR_DEBUGMODE 0
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0

#include <iostream>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <cube/tbasis_index.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_lproblem.h>
#include <galerkin/cached_problem.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <cube/tbasis_evaluate.h>
#include <interval/pq_frame.h>
#include <cube/tframe.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <cube/tbasis_indexplot.h>
#include <numerics/corner_singularity.h>

#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_evaluate.h>
#include <Ldomain/ldomain_evaluate.h>
#include <galerkin/ldomain_equation.h>

#include "ldomain_solutions.h"

//#include "TestFunctions2d.h"

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;


int main(){
    cout << "testing ldomain wavelet basis" << endl;
    const int d  = 2;
    const int dT = 2;
    const int dim = 2;
    const int jmax = 4;
    
    typedef PBasis<d,dT> Basis1d;
    typedef Basis1d::Index Index1d;
    typedef TensorBasis<Basis1d,dim> TBasis;
    typedef TBasis::Index TIndex;
    typedef LDomainBasis<Basis1d> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;
    
    Basis1d basis1d;
    basis1d.set_jmax(jmax);
    Basis basis(basis1d);
    basis.set_jmax(jmax);
    
    //setup index set
    set<Index> Lambda;  
    int zaehler=0;
    for (Index lambda=basis.first_generator(basis.j0()), itend=basis.last_wavelet(jmax); lambda <= itend; ++lambda)
    {
        Lambda.insert(lambda);
        //cout << lambda << " : " << lambda.number() << endl;
        zaehler++;
    }
    cout << zaehler << endl;
    
    //plot one function
    Array1D<SampledMapping<dim> > eval(3);
    Index ind=basis.get_wavelet(27);    //0-26:generatoren auf patches,
                                        //27-32:überlappende generatoren, indiziert mit p=3,4
                                        //33:überlappendes wavelet
                                        //34:nicht-überlappendes wavelet
    cout << "evaluate wavelet with index " << ind << endl;
    eval=basis.evaluate(ind,6);
    std::ofstream os("Ldomainoutput.m");
    os << "clf;" << endl;
    os << "axis([-1 1 -1 1 0 1]);" << endl;
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os);
        os << "surf(x,y,z);" << endl;
        
        os << "hold on;" << endl;
    }
    os << "view(30,55);"<<endl;
    os << "hold off" << endl;
    os.close(); 
     
    cout << "support of wavelet on patches:"  << endl;
    Support supp;
    basis.support(ind, supp);
    //intersect_supports(basis,basis.get_wavelet(27),basis.get_wavelet(29),supp);
    cout << "patch 0: [" << supp.xmin[0] <<" , " <<supp.xmax[0]<<"]x["<<supp.ymin[0]<<" , "<<supp.ymax[0]<<"]"<<  endl;
    cout << "patch 1: [" << supp.xmin[1] <<" , " <<supp.xmax[1]<<"]x["<<supp.ymin[1]<<" , "<<supp.ymax[1]<<"]"<<  endl;
    cout << "patch 2: [" << supp.xmin[2] <<" , " <<supp.xmax[2]<<"]x["<<supp.ymin[2]<<" , "<<supp.ymax[2]<<"]"<<  endl;
    
    
    //plot solution and right-hand-side
    //PolyRHS rhs1;
    //PolySolution solution1;
    CornerSingularity    solution1(Point<2>(0,0), 0.5, 1.5);
    CornerSingularityRHS rhs1(Point<2>(0,0), 0.5, 1.5);
    Array1D<Point<dim,int> > corners;
    corners.resize(3);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;
    std::ofstream osrhs("rhs.m");
    std::ofstream ossolution("uexakt.m");
    for(int i=0;i<3;i++){
        Point<2> q0(corners[i][0],corners[i][1]);
        Point<2> q1(corners[i][0]+1,corners[i][1]+1);
        Grid<2> mygrid(q0,q1,100);
        SampledMapping<2> smuexakt(mygrid, solution1);
        SampledMapping<2> smrhs(mygrid, rhs1); 
        smuexakt.matlab_output(ossolution);
        smrhs.matlab_output(osrhs);
        ossolution << "surf(x,y,z)"<<endl;
        osrhs << "surf(x,y,z)"<<endl;
        ossolution << "hold on;" << endl;
        osrhs << "hold on;" << endl;
    }
    ossolution << "hold off;" << endl;
    osrhs << "hold off;" << endl;
    ossolution.close();
    osrhs.close();
    cout << "solution and rhs plotted" << endl;
#if 0   
    //quick test of cdd1
    PoissonBVP<dim> poisson1(&rhs1);
    LDomainEquation<Basis1d> eq(&poisson1, basis);
    CachedLProblem<LDomainEquation<Basis1d> > cproblem1(&eq);
    //cout << basis.first_generator() << endl;
    //cout << cproblem1.norm_A() << endl;
    //cout << cproblem1.norm_Ainv() << endl;
    
    
    InfiniteVector<double, Index> u_epsilon;
    CDD1_SOLVE(cproblem1, 1e-3, u_epsilon, u_epsilon, jmax, tensor_simple);
    //u_epsilon.scale(-1);
    Array1D<SampledMapping<dim> > eval2(3);
    eval2=basis.evaluate(u_epsilon, 6);
    std::ofstream os2("solutioncdd1.m");
    for(int i=0;i<3;i++){
        eval2[i].matlab_output(os2);
        os2 << "surf(x,y,z);" << endl;
        os2 << "hold on;" << endl;
    }
    //os2 << "view(30,55);"<<endl;
    os2 << "hold off" << endl;
    os2.close();
#endif  
    cout << "fertig" << endl;
    return 0;
}
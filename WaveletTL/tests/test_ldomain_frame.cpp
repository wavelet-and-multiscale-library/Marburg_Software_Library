/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#define NONADAPTIVE
#undef ADAPTIVE

#define _WAVELETTL_USE_TBASIS 1

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
//#include <galerkin/cached_lproblem.h>

#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <adaptive/cdd2.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main(){
    
    cout << "testing L-domain quarklet frame" << endl;
    const int d  = 2;
    const int dT = 2;
    const int dim = 2;
    const int jmax = 6;
    const int pmax = 2;
    
    typedef PQFrame<d,dT> Frame1d;
    Frame1d frame1d(false,false);
    typedef LDomainFrame<Frame1d> Frame;
    typedef Frame::Index Index;
    //Frame frame(frame1d);
    Frame frame;
    cout << frame.j0() << endl;
    frame.set_jpmax(jmax,pmax);
    //cout << frame.frame1d().j0() << endl;
    
    //test LDomainEquation
    //CornerSingularityRHS rhs1(Point<2>(0,0), 0.5, 1.5);
    Vector<double> val(1, "1.0");
    ConstantFunction<dim> rhs1(val);
#if 1
    Array1D<Point<dim,int> > corners;
    corners.resize(3);
    corners[0][0]=-1;
    corners[0][1]=0;
    corners[1][0]=-1;
    corners[1][1]=-1;
    corners[2][0]=0;
    corners[2][1]=-1;
    std::ofstream osrhs("rhs.m");
    for(int i=0;i<3;i++){
        Point<2> q0(corners[i][0],corners[i][1]);
        Point<2> q1(corners[i][0]+1,corners[i][1]+1);
        Grid<2> mygrid(q0,q1,100);
        SampledMapping<2> smrhs(mygrid, rhs1); 
        smrhs.matlab_output(osrhs);
        osrhs << "surf(x,y,z)"<<endl;
        osrhs << "hold on;" << endl;
    }
    osrhs << "hold off;" << endl;
    osrhs.close();
    cout << "rhs plotted" << endl;
#endif
    PoissonBVP<dim> poisson1(&rhs1);
    LDomainFrameEquation<Frame1d,Frame> eq(&poisson1, true);
    eq.set_jpmax(jmax,pmax);
    cout << "normA: "<<eq.norm_A() <<endl;
    cout << "normAinv: "<<eq.norm_Ainv() << endl;
    //cout << "F_norm: "<<eq.F_norm() << endl;
    
    Index testindex=frame.get_quarklet(161);
    cout << "testindex: "<<testindex<<endl;
    cout << "a(lambda_0,lambda_0) : "<<eq.a(testindex,testindex) << endl;
    cout << "f(lambda_0): "<<eq.f(testindex)<<endl;
    cout << "D(lambda_0): "<<eq.D(testindex)<<endl;
    
    
    
#if 1
    //plot one function
    Array1D<SampledMapping<dim> > eval(3);
    //Index ind=frame.get_quarklet(159);    //0-26:generatoren auf patches,
                                        //27-32:überlappende generatoren, indiziert mit p=3,4
                                        //33:überlappendes wavelet
                                        //34:nicht-überlappendes wavelet
    cout << "evaluate quarklet with index " << testindex << endl;
    eval=frame.evaluate(testindex,6);
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
#endif  
    
    
    
    return 0;
}
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#define PARALLEL 0
#undef NONADAPTIVE
#define ADAPTIVE

#undef BASIS
#define FRAME

#define CUBE
#undef LDOMAIN

#ifdef FRAME
#define _WAVELETTL_USE_TFRAME 1
#else
#define _WAVELETTL_USE_TBASIS 1
#endif
#define _DIM 2
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 2

#include <iostream>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <cube/tbasis_index.h>
#include <utils/multiindex.h>
#include <numerics/bvp.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tframe_evaluate.h>
#include <interval/pq_frame.h>
#include <cube/tframe.h>
#include <cube/tframe_evaluate.h>
#include <cube/tframe_indexplot.h>
#include <Ldomain/ldomain_frame_index.h>
#include <Ldomain/ldomain_frame.h>
#include <Ldomain/ldomain_frame_evaluate.h>
#include <Ldomain/ldomain_frame_indexplot.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <cube/tbasis_indexplot.h>
#include <galerkin/tframe_equation.h>
#include <galerkin/cached_quarklet_tproblem.h>
//#include "TestFunctions2d.h"




using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;



int main()
{
    
    const int d  = 3;
    const int dT = 3;
    const unsigned int dim = 2; 
    const int jmax=6;
    const int pmax = 2;
    
#ifdef LDOMAIN
    cout << "Testing ldomainframe image outputs" << endl;
    
    
   

    
#ifdef FRAME
    
    typedef PQFrame<d,dT> Frame1d;
    typedef LDomainFrame<Frame1d> Frame;
    typedef Frame::Index Index;
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
#endif  
    
    InfiniteVector<double,Index> u;
    Index lambda=frame.get_quarklet(3997);
    cout << lambda << endl;
    u.set_coefficient(lambda, 1);
//    cout << u << endl;
    
    cout << "plotting quarklet on ldomain" << endl;
    
    
    Array1D<SampledMapping<dim> > eval(3);
    eval=frame.evaluate(u,8);
    std::ofstream os2("Image_outputs/quarklet_ldomain.m");
    os2 << "figure;"<< endl;
//    os2 <<  "cm=colormap('jet'); cm=[flipud(cm);cm]; colormap(cm); " << endl;
    for(int i=0;i<3;i++){
        eval[i].matlab_output(os2);        
           os2 << "surf(x,y,z,'LineStyle','none');"
           << "hold on;" << endl;
    }  
    os2 << "hold off" << endl;
    os2.close();
    
//    SampledMapping<2> sm1(evaluate(frame, u , true, 8));
//    std::ofstream stream1("Image_outputs/quarklet_ldomain.m");
//    sm1.matlab_output(stream1);
//    stream1 << "cm=colormap('jet'); cm=[flipud(cm);cm]; colormap(cm); "
////            << "surf(x,y,z,'LineStyle','none');"
//            << "h=surf(x,y,z);"
//            << "set(h, 'LineStyle','none');"
//             << "\ncolormap(cm);"
//              << "view(30,55);"
////              << "matlab2tikz('myquarklet.tex', 'standalone', true);"
//              <<endl;
//        
//    
//              
////    stream1   << "\nprint('-djpg', 'quarklet_cube.png');" << endl;
////              << endl;
//    stream1.close();
//    
//    cout << "quarklet plotted" << endl;
    
    
#endif        
    
    
#ifdef CUBE
    
    cout << "Testing tframe image outputs" << endl;
//        FixedArray1D<int,2*dim> s;          //set order of boundary conditions
//    s[0]=0, s[1]=0; s[2]=0, s[3]=0; 
//    //cout << s << endl;
    FixedArray1D<bool, 2*dim> bc;
    bc[0]=true, bc[1]=true, bc[2]=true, bc[3]=true;
    
    
   

    
#ifdef FRAME
    
    typedef PQFrame<d,dT> Frame1d;
//    typedef Basis1d::Index Index1d;
    typedef TensorFrame<Frame1d,dim> Frame;
    typedef Frame::Index Index;
//    TensorFrameEquation<Frame1d,dim,Frame> eq(&poisson1, bc);

    
    Frame frame(bc);
    frame.set_jpmax(jmax, pmax);
#endif  
    
    InfiniteVector<double,Index> u;
    Index lambda=frame.get_quarklet(15);
    cout << lambda << endl;
    u.set_coefficient(lambda, 1);
    cout << u << endl;
    
    cout << "plotting quarklet on cube" << endl;
    SampledMapping<2> sm1(evaluate(frame, u , true, 6));
    std::ofstream stream1("Image_outputs/quarklet_cube.m");
    sm1.matlab_output(stream1);
    stream1 //<< "figure;\nsurf(x,y,z,'LineStyle','none');"
            << "figure;\ns=surf(x,y,z);"            
            << "\nview(30,55);" 
//            << "\ngrid off;"
//            << "\nshading interp;"
//            << "print('-djpg',sprintf('jetztaber.jpg'));"

//              << "matlab2tikz('myquarklet.tex', 'standalone', true);"
              << endl;
        
    
              
//    stream1   << "\nprint('-djpg', 'quarklet_cube.png');" << endl;
//              << endl;
    stream1.close();
    
    cout << "quarklet plotted" << endl;
    
    
#endif    
    

    return 0;

}


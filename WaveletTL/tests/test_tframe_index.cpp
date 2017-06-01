#include <iostream>

#include <interval/pq_frame.h>
#include <utils/fixed_array1d.h>
#include <cube/tframe.h>
#include <cube/tframe_index.h>
#include <utils/multiindex.h>

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
    const int d  = 2;
    const int dT = 2;
#if 1
    const unsigned int dim = 2; const int levelrange(2); const int pmax(2);
//    const unsigned int dim = 2; const int levelrange(1);
    typedef PQFrame<d,dT> Frame1d;
//    typedef DSFrame<d,dT> Frame1d;
    typedef TensorFrame<Frame1d,dim> Frame;
//    typedef Frame::Index Index;
//    typedef TensorQIndex<Frame1d,dim>::polynomial_type polynomial_type;
//    polynomial_type p;
//    for(int i=0; i<100; i++){
//        ++p;
//        cout << p << endl;
//    }
//    
//    abort();
//    typedef TensorQIndex<Frame1d,dim>::level_type level_type;
//    typedef TensorQIndex<Frame1d,dim>::type_type type_type;
//    typedef TensorQIndex<Frame1d,dim>::translation_type translation_type;
    
    Frame frameH;
    FixedArray1D<int,2*dim> s; //,st;
    for (unsigned int i(0);i<(2*dim);i++) //s[0]=s[1]=s[2]=s[3]=s[4]=s[5]=0;
    {
        s[i]=0;
    }
    // st[0]=st[1]=st[2]=st[3]=0;
    Frame frameS(s);
    FixedArray1D<bool,2*dim> bc;
    for (unsigned int i=0;i<2*dim;i++)
    {
        bc[i]=false;
    }
    //bc[0]=bc[1]=bc[2]=bc[3]=bc[4]=bc[5]=true;
    Frame frameBC(bc);
    Frame1d fra00(false,false); Frame1d fra10(true,false); Frame1d fra01(false,true); Frame1d fra11(true,true);
    fra00.set_jpmax(5,0); fra10.set_jpmax(5,0); fra01.set_jpmax(5,0); fra11.set_jpmax(5,0);
    FixedArray1D<Frame1d*,dim> framesArray;
    for(unsigned int i=0;i<dim;i++)
    {
        framesArray[i] = &fra00;
    }
    


    // Run tests with this frame:
    //Frame frame;
    //Frame frame(s);
    //Frame frame(bc);
    Frame frame(framesArray);
    
    frame.set_jpmax(multi_degree(frame.j0())+levelrange, pmax);
    
    
    
    Frame frameFra(framesArray); frameFra.set_jpmax(multi_degree(frameFra.j0())+levelrange, pmax);
    framesArray[0] = &fra10;
    Frame frameFra2(framesArray);frameFra2.set_jpmax(multi_degree(frameFra2.j0())+levelrange, pmax);
    framesArray[0] = &fra01;
    Frame frameFra3(framesArray);frameFra3.set_jpmax(multi_degree(frameFra3.j0())+levelrange, pmax);
    framesArray[dim-1] = &fra10;
    Frame frameFra4(framesArray);frameFra4.set_jpmax(multi_degree(frameFra4.j0())+levelrange, pmax);
    framesArray[1] = &fra01;
    Frame frameFra5(framesArray);frameFra5.set_jpmax(multi_degree(frameFra5.j0())+levelrange, pmax);
    framesArray[0] = &fra11; framesArray[dim-1] = &fra01;
    Frame frameFra6(framesArray);frameFra6.set_jpmax(multi_degree(frameFra6.j0())+levelrange, pmax);

    FixedArray1D<Frame*,6> TFrameArray;
    TFrameArray[0] = &frameFra;
    TFrameArray[1] = &frameFra2;
    TFrameArray[2] = &frameFra3;
    TFrameArray[3] = &frameFra4;
    TFrameArray[4] = &frameFra5;
    TFrameArray[5] = &frameFra6;

    
    cout << "Testing tframe_index" << endl;
    for (unsigned int i=0; (int)i<TFrameArray[0]->degrees_of_freedom(); ++i)
    {
        cout << "i = " << i << "; lambba(i) = " << *TFrameArray[0]->get_wavelet(i) << endl;
    }
    abort();
#endif


    return 0;
}

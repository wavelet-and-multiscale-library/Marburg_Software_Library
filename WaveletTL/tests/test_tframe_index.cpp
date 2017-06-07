#include <iostream>

#include <interval/pq_frame.h>
#include <utils/fixed_array1d.h>
#include <cube/tframe.h>
#include <cube/tframe_index.h>
#include <cube/tframe_support.h>
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
    const unsigned int DIM = 2;
//    const unsigned int dim = 2; const int levelrange(1);
    typedef PQFrame<d,dT> Frame1D;
//    typedef DSFrame<d,dT> Frame1D;
    typedef TensorFrame<Frame1D,dim> Frame;
    typedef Frame::Index Index;
//    typedef TensorQIndex<Frame1D,dim>::polynomial_type polynomial_type;
//    polynomial_type p;
//    for(int i=0; i<100; i++){
//        ++p;
//        cout << p << endl;
//    }
//    
//    abort();
//    typedef TensorQIndex<Frame1D,dim>::level_type level_type;
//    typedef TensorQIndex<Frame1D,dim>::type_type type_type;
//    typedef TensorQIndex<Frame1D,dim>::translation_type translation_type;
    
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
    Frame1D fra00(false,false); Frame1D fra10(true,false); Frame1D fra01(false,true); Frame1D fra11(true,true);
    fra00.set_jpmax(5,pmax); fra10.set_jpmax(5,pmax); fra01.set_jpmax(5,pmax); fra11.set_jpmax(5,pmax);
    FixedArray1D<Frame1D*,dim> framesArray;
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
    
    #if 1
    
//    typedef Frame1D::Index Index1D;
        typedef TensorQIndex<Frame1D,DIM,Frame>::level_type level_type;
        typedef TensorQIndex<Frame1D,DIM,Frame>::polynomial_type polynomial_type;
    typedef std::list<Index> SupportList;
    for (unsigned int b = 0; b < TFrameArray.size(); b++)
    {
        Index lambda(TFrameArray[b]->get_quarklet(2000));

        level_type lambda_j(TFrameArray[b]->j0());
        polynomial_type p(1,1);

        SupportList intwav1, intwav2, intwav3, intel1, intel2;

        Frame::Support supp_lambda;
        TFrameArray[b]->support(lambda, supp_lambda);
        cout << "intersecting generators" << endl;
        cout << "lambda :" << endl;
        cout << "    N = " << lambda.number() << " == " << lambda << "; support = 2^{-(" << supp_lambda.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp_lambda.j[i];
            }
            cout << ")}[" << supp_lambda.a[0] << "," << supp_lambda.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp_lambda.a[i] << "," << supp_lambda.b[i] ;
            }
            cout << "]" << endl;


        Index lambda1(TFrameArray[b]->get_quarklet(81));
        Frame::Support supp_lambda1;
        TFrameArray[b]->support(lambda1, supp_lambda1);
        cout << "lambda1 :" << endl;
        cout << "    N = " << lambda1.number() << " == " << lambda1 << "; support = 2^{-(" << supp_lambda1.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp_lambda1.j[i];
            }
            cout << ")}[" << supp_lambda1.a[0] << "," << supp_lambda1.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp_lambda1.a[i] << "," << supp_lambda1.b[i] ;
            }
            cout << "]" << endl;

        Frame::Support temp_supp;
        bool temp_bool;
        temp_bool = intersect_supports(*TFrameArray[b], lambda, lambda1, temp_supp);
        if (!temp_bool)
        {
            cout << "supports do not intersect. output of intersect_supports:" << endl;
        }
        cout << "    support intersection = 2^{-(" << temp_supp.j[0];
        for (unsigned int  i=1; i<DIM; ++i)
        {
            cout << ", " << temp_supp.j[i];
        }
        cout << ")}[" << temp_supp.a[0] << "," << temp_supp.b[0];
        for (unsigned int  i=1; i<DIM; ++i)
        {
            cout << "]x[" << temp_supp.a[i] << "," << temp_supp.b[i] ;
        }
        cout << "]" << endl;
    
    // compute all intersecting_generators for lambda and print the support and support intersection
        cout << "intersecting_generators for lambda = " << lambda << "; .number() = " << lambda.number() << "; on level = " << lambda_j << 
                ", and polynomial degree = " << p << endl;
        intersecting_quarklets<Frame1D,DIM>(*TFrameArray[b], lambda, lambda_j, true, intwav1, p);
        cout << "Number of intersecting generators is " << intwav1.size() << endl;
        for (SupportList::const_iterator it(intwav1.begin()); it != intwav1.end(); ++it)
        {
            Frame::Support supp, supp_it;
            TFrameArray[b]->support(*it, supp_it);
            assert (intersect_supports(*TFrameArray[b], lambda, *it, supp) == true);
            cout << "    N = " << (*it).number() << " == " << (*it) << "; support = 2^{-(" << supp_it.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp_it.j[i];
            }
            cout << ")}[" << supp_it.a[0] << "," << supp_it.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp_it.a[i] << "," << supp_it.b[i] ;
            }
            cout << "]; support intersection = 2^{-(" << supp.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp.j[i];
            }
            cout << ")}[" << supp.a[0] << "," << supp.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp.a[i] << "," << supp.b[i] ;
            }
            cout << "]" << endl;

        }
        // compute all intersecting_wavelets for lambda and print the support and support intersection
        cout << "intersecting wavelets on level j0 = " << TFrameArray[b]->j0() << endl;
        intersecting_quarklets<Frame1D,DIM>(*TFrameArray[b], lambda, lambda_j, false, intwav2, p);
        cout << "Number of intersecting wavelets is " << intwav2.size() << endl;
        for (SupportList::const_iterator it(intwav2.begin()); it != intwav2.end(); ++it)
        {
            Frame::Support supp, supp_it;
            TFrameArray[b]->support(*it, supp_it);
            intersect_supports(*TFrameArray[b], lambda, *it, supp);
            cout << "    N = " << (*it).number() << " == " << (*it) << "; support = 2^{-(" << supp_it.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp_it.j[i];
            }
            cout << ")}[" << supp_it.a[0] << "," << supp_it.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp_it.a[i] << "," << supp_it.b[i] ;
            }
            cout << "]; support intersection = 2^{-(" << supp.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp.j[i];
            }
            cout << ")}[" << supp.a[0] << "," << supp.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp.a[i] << "," << supp.b[i] ;
            }
            cout << "]" << endl;

        }

        lambda_j[DIM-1]++;

        cout << "intersecting wavelets on level " << lambda_j << endl;
        intersecting_quarklets<Frame1D,DIM>(*TFrameArray[b], lambda, lambda_j, false, intwav3, p);
        cout << "Number of intersectiong wavelets is " << intwav3.size() << endl;
        for (SupportList::const_iterator it(intwav3.begin()); it != intwav3.end(); ++it)
        {
            Frame::Support supp, supp_it;
            /*
            if ((*it).number() == 361)
            {
                Frame::Support temp_supp, temp_inters;
                basis.support(*it, temp_supp);
                intersect_supports(basis, lambda, *it, temp_inters);

                assert (intersect_supports(basis, lambda, *it, temp_inters) == true);

                cout << "    N = " << (*it).number() << " == " << (*it) << "; support = 2^{-(" << temp_supp.j[0];
                for (unsigned int  i=1; i<DIM; ++i)
                {
                    cout << ", " << temp_supp.j[i];
                }
                cout << ")}[" << temp_supp.a[0] << "," << temp_supp.b[0];
                for (unsigned int  i=1; i<DIM; ++i)
                {
                    cout << "]x[" << temp_supp.a[i] << "," << temp_supp.b[i] ;
                }
                cout << "]; support intersection = 2^{-(" << temp_inters.j[0];
                for (unsigned int  i=1; i<DIM; ++i)
                {
                    cout << ", " << temp_inters.j[i];
                }
                cout << ")}[" << temp_inters.a[0] << "," << temp_inters.b[0];
                for (unsigned int  i=1; i<DIM; ++i)
                {
                    cout << "]x[" << temp_inters.a[i] << "," << temp_inters.b[i] ;
                }
                cout << "]" << endl;
            }
            */
            TFrameArray[b]->support(*it, supp_it);
            intersect_supports(*TFrameArray[b], lambda, *it, supp);
            cout << "    N = " << (*it).number() << " == " << (*it) << "; support = 2^{-(" << supp_it.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp_it.j[i];
            }
            cout << ")}[" << supp_it.a[0] << "," << supp_it.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp_it.a[i] << "," << supp_it.b[i] ;
            }
            cout << "]; support intersection = 2^{-(" << supp.j[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << ", " << supp.j[i];
            }
            cout << ")}[" << supp.a[0] << "," << supp.b[0];
            for (unsigned int  i=1; i<DIM; ++i)
            {
                cout << "]x[" << supp.a[i] << "," << supp.b[i] ;
            }
            cout << "]" << endl;

        }
    }
    
#endif
    
#if 0
    
    int tempA, tempB;
    //Index1D temp_mu(4,1,6,TFrameArray[0]->frames()[1]);
    Index1D first(TFrameArray[0]->frames()[1]->get_wavelet(0));
    //typedef Frame1D::Support Support1D;
    //Support1D supp1D;
    TFrameArray[0]->frames()[1]->support(first, tempA, tempB);
    cout << "wavelet = " << first << "; 1D support = (" << tempA << ", " << tempB << ")" << endl;
            
            
    Index1D temp_mu(TFrameArray[0]->frames()[1]->get_wavelet(20));
    cout << "temp_mu = " << temp_mu << endl;
    
    get_intersecting_quarklets_on_level(*(TFrameArray[0]->frames()[1]),
                    temp_mu,
            4,true,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
    get_intersecting_quarklets_on_level(*(TFrameArray[0]->frames()[1]),
                    temp_mu,
            4,false,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
     
    cout << "evaluate wavelet mu = " << temp_mu << endl;
    for (unsigned int i=0; i< 100; ++i)
    {
        cout << "f(" << (0.01*i) << ") = " << TFrameArray[0]->frames()[1]->evaluate(0, temp_mu, 0.01*i) << endl;
    }

    cout << "min = " << tempA << "; max = " << tempB << endl;
    for (unsigned int i = 0; i< 400; ++i)
    {
        Index1D temp_index(TFrameArray[0]->frames()[1]->get_wavelet(i));
        TFrameArray[0]->frames()[1]->support(temp_index, tempA, tempB);
        cout << "N = " << i << "; lam = " << temp_index << "; k1 = " << tempA << "; k2 = " << tempB << endl;
    }
    abort();
#endif

    
    cout << "Testing tframe_index" << endl;
    cout << "DoF: " << TFrameArray[0]->degrees_of_freedom() << endl;
    for (unsigned int i=0; (int)i<TFrameArray[0]->degrees_of_freedom(); ++i)
    {
        cout << "i = " << i << "; lambba(i) = " << *TFrameArray[0]->get_quarklet(i) << endl;
    }
    abort();
#endif


    return 0;
}

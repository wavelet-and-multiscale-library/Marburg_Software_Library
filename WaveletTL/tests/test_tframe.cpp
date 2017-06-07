#include <iostream>


#include <utils/fixed_array1d.h>
#include <cube/tframe.h>
#include <utils/multiindex.h>
#include <cube/tframe_evaluate.h>
#include <interval/pq_frame.h>
#include <interval/pq_evaluate.h>
#include <cube/tframe_indexplot.h>
#include <cube/tframe_support.h>
#include <cube/tframe_index.h>
using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
    cout << "Testing tframe." << endl;

    typedef PQFrame<2,2> Frame1d;
//    typedef Frame1d::Index IIndex;
    const unsigned int DIM = 3;
    typedef TensorFrame<Frame1d,DIM> Frame;
//    typedef Frame::Support Support;
//    typedef Frame::Index Index;
//    typedef Index::level_type level_type;
//    typedef Index::polynomial_type polynomial_type;
//    typedef Index::type_type type_type;
//    typedef Index::translation_type translation_type;


    const unsigned int tempj(9), tempp(0);
    

    cout << "Testing constructors" << endl;
    cout << "default constructor (with homogeneous b.c.)" << endl;
    Frame frameH;
    cout << "done" << endl;

    //
    cout << "frameH ...";
    frameH.set_jpmax(tempj, tempp);
    cout << "done" << endl;

    cout << "Constructor with specified boundary condition orders"<<endl;
    cout << "Initialize FixedArray<int,4> with b.c.'s" << endl;

    FixedArray1D<int,2*DIM> s; // ,st; // BC for the dual frame can only be specified for dsframe, not for pframe
    for (unsigned int i=0; i<2*DIM; i++)
    {
        s[i]=0;
        //st[i]=0;
    }
    //s[0] = 1;
    //s[2*DIM-1]=1;

    cout << "Call constructor with s" << endl;
    Frame frameS(s);
    cout << "done" << endl;

    //
    cout << "frameS ...";
    frameS.set_jpmax(tempj, tempp);
    cout << "done" << endl;
    
    cout << "Constructor with specified Dirichlet boundary conditions for the primal functions" << endl;
    cout << "Initialize FixedArray<bool,4> with b.c.'s" << endl;
    FixedArray1D<bool,2*DIM> bc;
    for (unsigned int i=0; i<2*DIM; i++)
    {
        bc[i]=false;
    }
    //bc[0] = true;
    //bc[2*DIM-1] = true;
    
    cout << "Call constructor with bc" << endl;
    Frame frameBC(bc);
    cout << "done" << endl;

    //
    cout << "frameBC ...";
    frameBC.set_jpmax(tempj, tempp);
    cout << "done" << endl;
    
    cout << "Constructor with precomputed instances of the 1D frames" << endl;
    cout << "Initialize 1d frames" << endl;

    //Frame1d fra0();
    Frame1d fra00(false,false);
    Frame1d fra10(true,false);
    Frame1d fra01(false,true);
    Frame1d fra11(true,true);
    fra00.set_jpmax(5, tempp);
    fra10.set_jpmax(5, tempp);
    fra01.set_jpmax(5, tempp);
    fra11.set_jpmax(5, tempp);
    

    FixedArray1D<Frame1d*,DIM> framesArray;
    for (unsigned int i=0; i<DIM; i++)
    {
        framesArray[i] = &fra00;
    }
    //cout << "test: " << frameS.j0() << endl;
    //framesArray[0] = &fra10;
    //cout << "test: " << frameS.j0() << endl;
    //framesArray[DIM-1] = &fra01;
    
    cout << "Call constructor with array of 1d frames" << endl;
    Frame frameBas(framesArray);
    cout << "done" << endl;
    
    //
    cout << "frameBas ...";
    frameBas.set_jpmax(tempj, tempp);
    cout << "done" << endl;
    

    // Use this frame for testing:
    //Frame frame;
    //Frame frame(s);
    //Frame frame(s,st);
    //Frame frame(bc);
    Frame frame(framesArray);

    // test equality of the constructed frames
    //FixedArray1D<Frame*,4> tensorframes;
    //tensorframes[0] = &frameH; tensorframes[1] = &frameS; tensorframes[2] = &frameBC; tensorframes[3] = &frameBas;
    FixedArray1D<Frame*,3> tensorframes;
    tensorframes[0] = &frameS; tensorframes[1] = &frameBC; tensorframes[2] = &frameBas;
    
    unsigned int maxlevelrange(4); // levels up to \|j0\|+maxlevelrange are tested
    unsigned int jpmax(multi_degree(frame.j0())+maxlevelrange);
    
    // Run tests with this frame:

    frame.set_jpmax(multi_degree(frame.j0())+maxlevelrange);
    frameBas.set_jpmax(multi_degree(frameBas.j0())+maxlevelrange);
    
    framesArray[0] = &fra10;
    Frame frameBas2(framesArray);
    frameBas2.set_jpmax(multi_degree(frameBas2.j0())+maxlevelrange);
    framesArray[0] = &fra01;
    Frame frameBas3(framesArray);frameBas3.set_jpmax(multi_degree(frameBas3.j0())+maxlevelrange);
    framesArray[DIM-1] = &fra10;
    Frame frameBas4(framesArray);frameBas4.set_jpmax(multi_degree(frameBas4.j0())+maxlevelrange);
    framesArray[DIM-1] = &fra01;
    Frame frameBas5(framesArray);frameBas5.set_jpmax(multi_degree(frameBas5.j0())+maxlevelrange);
    framesArray[0] = &fra11; framesArray[DIM-1] = &fra01;
    Frame frameBas6(framesArray);frameBas6.set_jpmax(multi_degree(frameBas6.j0())+maxlevelrange);
    
    FixedArray1D<Frame*,6> TFrameArray;
    TFrameArray[0] = &frameBas;
    TFrameArray[1] = &frameBas2;
    TFrameArray[2] = &frameBas3;
    TFrameArray[3] = &frameBas4;
    TFrameArray[4] = &frameBas5;
    TFrameArray[5] = &frameBas6;
    
    
    cout << "Testing set_jpmax" << endl;
    cout << "frameH ...";
    frameH.set_jpmax(jpmax);
    cout << "done" << endl;
    cout << "frameS ...";
    frameS.set_jpmax(jpmax);
    cout << "done" << endl;
    cout << "frameBC ...";
    frameBC.set_jpmax(jpmax);
    cout << "done" << endl;
    cout << "frameBas ...";
    frameBas.set_jpmax(jpmax);
    cout << "done" << endl;
   
    for (unsigned int i=0; i< tensorframes.size(); ++i)
    {
        cout << "jpmax= " << jpmax << "; " << tensorframes[i]->last_quarklet(tempj).number() << "; " << tensorframes[i]->last_quarklet(jpmax) << endl;
    }
       
    cout << "" << endl;
    
#if 0
    for (unsigned int b=0; b<3; ++b)
    {
        cout << "List differences of tensorframe " << b << " and " << (b+1) << endl;
        if (tensorframes[b]->degrees_of_freedom() 
                != tensorframes[b+1]->degrees_of_freedom())
        {
            cout << " DOF[" << b << "] = " << tensorframes[b]->degrees_of_freedom() << "; "
                 << " DOF[" << (b+1) << "] = " << tensorframes[b+1]->degrees_of_freedom() << endl;
        }
        
        if (tensorframes[b]->j0()
                != tensorframes[b+1]->j0())
        {
            cout << " j0[" << b << "] = " << tensorframes[b]->j0() << "; "
                 << " j0[" << (b+1) << "] = " << tensorframes[b+1]->j0() << endl;
        }
        
        if (tensorframes[b]->get_jpmax()
                != tensorframes[b+1]->get_jpmax())
        {
            cout << " get_jpmax[b=" << b << "] = " << tensorframes[b]->get_jpmax() << "; "
                 << " get_jpmax[b=" << (b+1) << "] = " << tensorframes[b+1]->get_jpmax() << endl;
        }
        
        for (unsigned int d=0; d<DIM; d++)
        {
            cout << "List differences of tensorframe " << b << " and " << (b+1) << " in direction " << d << endl;
            if (tensorframes[b]->frames()[d]->DeltaLmin() 
                    != tensorframes[b+1]->frames()[d]->DeltaLmin()  )
            {
                cout << " DeltaLmin()[b=" << b     << "][d=" << d << "] = " << tensorframes[b]->frames()[d]->DeltaLmin()   << "; "
                     << " DeltaLmin()[b=" << (b+1) << "][d=" << d << "] = " << tensorframes[b+1]->frames()[d]->DeltaLmin() << endl;
            }
            if (tensorframes[b]->frames()[d]->DeltaLmax() 
                    != tensorframes[b+1]->frames()[d]->DeltaLmax()  )
            {
                cout << " DeltaLmax()[b=" << b     << "][d=" << d << "] = " << tensorframes[b]->frames()[d]->DeltaLmax()   << "; "
                     << " DeltaLmax()[b=" << (b+1) << "][d=" << d << "] = " << tensorframes[b+1]->frames()[d]->DeltaLmax() << endl;
            }
            
            
            if (tensorframes[b]->frames()[d]->Nablamin()
                    != tensorframes[b+1]->frames()[d]->Nablamin()  )
            {
                cout << " Nablamin()[b=" << b     << "][d=" << d << "] = " << tensorframes[b]->frames()[d]->Nablamin()   << "; "
                     << " Nablamin()[b=" << (b+1) << "][d=" << d << "] = " << tensorframes[b+1]->frames()[d]->Nablamin() << endl;
            }
            
            
            for (unsigned int j= multi_degree(tensorframes[0]->j0()); j< (multi_degree(tensorframes[0]->j0()) + maxlevelrange); ++j)
            {
                cout << "List differences of tensorframe " << b << " and " << (b+1) << " in direction " << d << " on level " << j << endl;
                if (tensorframes[b]->frames()[d]->Deltasize(j)
                        != tensorframes[b+1]->frames()[d]->Deltasize(j)  )
                {
                    cout << " Deltasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorframes[b]->frames()[d]->Deltasize(j)  << "; "
                         << " Deltasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorframes[b+1]->frames()[d]->Deltasize(j) << endl;
                }
                
                if (tensorframes[b]->frames()[d]->Nablamax(j)
                        != tensorframes[b+1]->frames()[d]->Nablamax(j)  )
                {
                    cout << " Nablamax(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorframes[b]->frames()[d]->Nablamax(j)  << "; "
                         << " Nablamax(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorframes[b+1]->frames()[d]->Nablamax(j) << endl;
                }
                
                if (tensorframes[b]->frames()[d]->Nablasize(j)
                        != tensorframes[b+1]->frames()[d]->Nablasize(j)  )
                {
                    cout << " Nablasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorframes[b]->frames()[d]->Nablasize(j)  << "; "
                         << " Nablasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorframes[b+1]->frames()[d]->Nablasize(j) << endl;
                }
                if (tensorframes[b]->frames()[d]->Deltasize(j)
                        != tensorframes[b+1]->frames()[d]->Deltasize(j)  )
                {
                    cout << " Deltasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorframes[b]->frames()[d]->Deltasize(j)  << "; "
                         << " Deltasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorframes[b+1]->frames()[d]->Deltasize(j) << endl;
                }
            }
        }
    }
    
#endif

    
    /*        
    cout << "Testing j0()" << endl
         << "frameH " << frameH.j0() << endl
         << "frameS " << frameS.j0() << endl
//             << "frameSST "<< frameSST.j0() << endl
         << "frameBC "<< frameBC.j0() << endl
         << "frameBas "<< frameBas.j0() << endl;
*/
    

#if 0
    for (unsigned int b=0; b<4; ++b)
    {
    cout << "b = " << b 
            << "; last_quarklet = " << tensorframes[b]->last_quarklet(tensorframes[b]->get_jpmax()) 
            << "; .number() = " << tensorframes[b]->last_quarklet(tensorframes[b]->get_jpmax()).number() 
            << "; dof = " << tensorframes[b]->degrees_of_freedom() << endl;
    }
#endif

#if 0
    // see test_tframe_support.cpp
    cout << "Testing support" << endl;
    // void support(const Index& lambda, Support& supp) const;
    Support supp;

    //cout << "DeltaLmax()  " << frame.frames()[0]->DeltaLmax() << " " << frame.frames()[1]->DeltaLmax() << endl;
    //cout << frameH.frames()[0]->DeltaLmax() << endl;
    //cout << "frame.first_q_generator();"  << frame.first_q_generator()<< endl;
    //frameH.first_q_generator();

    int maxnumber(10);
    for (Index it(frame.first_q_generator());it < frame.last_quarklet(jpmax);++it)
    {
        frame.support(it,supp);
        //cout << "supp = " << supp << endl;
        cout << "    N = " << it.number() << " == " << it << "; support = 2^{-j=(" << it.j()[0];
        for (unsigned int  i=1; i<DIM; ++i)
        {
            cout << ", " << it.j()[i];
        }
        cout << ") -e=(" << it.e()[0];
        for (unsigned int  i=1; i<DIM; ++i)
        {
            cout << ", " << it.e()[i];
        }
        cout << ")}[" << supp.a[0] << "," << supp.b[0];
        for (unsigned int  i=1; i<DIM; ++i)
        {
            cout << "]x[" << supp.a[i] << "," << supp.b[i] ;
        }
        cout << "]" << endl;

        if (it.number() == maxnumber) break;
    }
#endif


#if 0

    // static double primal_regularity() { return IFRAME::primal_regularity(); }
    // static unsigned int primal_polynomial_degree() { return IFRAME::primal_polynomial_degree(); }
    // static unsigned int primal_vanishing_moments() { return IFRAME::primal_vanishing_moments(); }
    cout << "Testing primal_ functions"<<endl;
    cout << "primal_regularity "<< frame.primal_regularity() << endl
         << "primal_polynomial_degree " << frame.primal_polynomial_degree() << endl
         << "primal_vanishing_moments " << frame.primal_vanishing_moments() << endl;
#endif
#if 0
    cout << "Testing first/last routines" << endl;
    Index tempi1, tempi2;
    tempi1 = frame.first_q_generator();
    tempi2 = first_q_generator<Frame1d,DIM,Frame> (&frame);
    assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ) );

    tempi1 = frame.last_generator();
    tempi2 = last_q_generator<Frame1d,DIM,Frame> (&frame);
    assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

    level_type j0(frame.j0());
    level_type currentlevel(j0);
    int range(0), maxrange(2); // number of level steps we have climbed so far.

    while (range <= maxrange)
    {
        tempi1 = frame.first_quarklet(currentlevel);
        tempi2 = first_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

        tempi1 = frame.last_quarklet(currentlevel);
        tempi2 = last_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

        // increase the level j,  meaning: "currentlevel++"
        //cout << currentlevel << endl;
        for (int i(DIM-1); i >= 0; i--)
        {
            if (i != 0)
            {
                if (currentlevel[i] != j0[i])
                {
                    // increase left neighbor
                    currentlevel[i-1]=currentlevel[i-1]+1;
                    int temp = currentlevel[i]-j0[i];
                    currentlevel[i]=j0[i];
                    currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                    break;
                }
            } else // i == 0. "big loop" arrived at the last index. We have to increase range!
            {
                range = range +1;
                if (DIM == 1)
                {
                    currentlevel[0]=currentlevel[0]+1;
                }
                else
                {
                    //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                    currentlevel[DIM-1]=j0[DIM-1]+range;
                    currentlevel[0]=j0[0];
                }
            }
        } // end of "curentlevel++"
    } // end of while (range < maxrange)

    cout << "Testing first_quarklet" << endl;
    unsigned int maxrange2(2);
    assert (multi_degree(frame.j0()) + maxrange2 <= frame.get_jmax());
    
    for (unsigned int lev = multi_degree(frame.j0()); lev < multi_degree(frame.j0()) + maxrange2; lev++)
    {
        tempi1 = frame.first_quarklet(lev);
        tempi2 = first_quarklet<Frame1d,DIM,Frame>(&frame,lev);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

        tempi1 = frame.last_quarklet(lev);
        tempi2 = last_quarklet<Frame1d,DIM,Frame>(&frame,lev);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));
    }
    
#endif    
#if 0
    
    cout << "Comparing speed of first/last routines" << endl
            << " - in tframe.h and tframe_index.h" << endl;
    Index tempi3, tempi4;
    int range2(0), maxrange3(12), repetitions(1000);
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; rep++)
    {
        //cout << "testing tframe.h; rep = " << rep << "; range2 = ";
        range2=0;
        currentlevel=j0;
        while (range2 <= maxrange3)
        {
            //cout << range2 << "; ";
            tempi3 = frame.first_quarklet(currentlevel);
            tempi4 = frame.last_quarklet(currentlevel);
            //tempi3 = first_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);
            //tempi4 = last_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);

            // increase the level j,  meaning: "currentlevel++"
            //cout << currentlevel << endl;
            for (int i(DIM-1); i >= 0; i--)
            {
                if (i != 0)
                {
                    if (currentlevel[i] != j0[i])
                    {
                        // increase left neighbor
                        currentlevel[i-1]=currentlevel[i-1]+1;
                        int temp = currentlevel[i]-j0[i];
                        currentlevel[i]=j0[i];
                        currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                        break;
                    }
                } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                {
                    range2 = range2 +1;
                    if (DIM == 1)
                    {
                        currentlevel[0]=currentlevel[0]+1;
                    }
                    else
                    {
                        //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                        currentlevel[DIM-1]=j0[DIM-1]+range2;
                        currentlevel[0]=j0[0];
                    }
                }
            } // end of "curentlevel++"
        } // end of while (range < maxrange)
        //cout << endl;
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; rep++)
    {
        range2=0;
        currentlevel=j0;
        while (range2 <= maxrange3)
        {
            tempi3 = first_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);
            tempi4 = last_quarklet<Frame1d,DIM,Frame>(&frame, currentlevel);

            // increase the level j,  meaning: "currentlevel++"
            //cout << currentlevel << endl;
            for (int i(DIM-1); i >= 0; i--)
            {
                if (i != 0)
                {
                    if (currentlevel[i] != j0[i])
                    {
                        // increase left neighbor
                        currentlevel[i-1]=currentlevel[i-1]+1;
                        int temp = currentlevel[i]-j0[i];
                        currentlevel[i]=j0[i];
                        currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                        break;
                    }
                } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                {
                    range2 = range2 +1;
                    if (DIM == 1)
                    {
                        currentlevel[0]=currentlevel[0]+1;
                    }
                    else
                    {
                        //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                        currentlevel[DIM-1]=j0[DIM-1]+range2;
                        currentlevel[0]=j0[0];
                    }
                }
            } // end of "curentlevel++"
        } // end of while (range < maxrange)
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; levelrange = " << maxrange3 << endl
            << "tframe.h: " << time1 << "sec; tframe_index.h: " << time2 << "sec" << endl;
#endif
    
    



#if 0
    // this test works only for 1 tensorindex at a time and uses matlab
    // therefore it is discouraged
    // evaluate is rather simple anyways!
    cout << "Test evaluate by creating matlab output" << endl
         << "1. by hand in each coordinate direction" << endl
         << "2. directly for a tensorindex" << endl;
    cout << "Create a sampling of a 1D function given by a grid" << endl;


    InfiniteVector<double, IIndex> coeff1d;
    IIndex index1d(first_q_generator(frame.frames()[0], frame.frames()[0]->j0()));
    coeff1d.set_coefficient(index1d,1.0);
    int dil=8;
    Array1D<double> grid1d;
    Array1D<double> funval1d;
    grid1d.resize((1<<dil)+1);
    for(unsigned int k=0;k<grid1d.size();k++)
        grid1d[k]=k*(1.0/(1<<dil));
    Grid<1>gitter(grid1d);
    //evaluate(fra11,0,index1d,grid1d,funval1d);
    //evaluate(*(framesArray[0]),0,index1d,grid1d,funval1d);
    evaluate(*(frame.frames()[0]),0,index1d,grid1d,funval1d);
    SampledMapping<1> my_quarklet1d (gitter, funval1d);
    ofstream fileA("my_quarklet1d.m");
    my_quarklet1d.matlab_output(fileA);
    fileA.close();

    cout << "Fix an tensorindex" << endl;
    cout << "compute by hand 1D samples corresponding to each 1d quarklet in the tensorquarklet" << endl;
    cout << "store them as matlab output" << endl;

    int res(7);

    level_type my_j;
    type_type my_e;
    translation_type my_k;
    polynomial_type my_p;
    
    my_j = frame.j0();
    //my_j[0]=frame.frames()[0]->j0(); // == 4
    //my_j[1]=frame.frames()[1]->j0();
    for (unsigned int i=0; i<DIM; i++)
    {
        assert(my_e[i] == 0);
        my_k[i] = frame.frames()[i]->DeltaLmin();
    }
    my_k[DIM-1] += 1;
    //my_e[0]=0;
    //my_e[1]=0;
    //my_k[0]=frame.frames()[0]->DeltaLmin()+0; // DeltaLmin ==2, Deltasize(4) == 12
    //my_k[1]=frame.frames()[1]->DeltaLmin()+1;

    //TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORFRAME* frame);
    Index myindex (my_p,my_j,my_e,my_k,&frame);


    // simpler:
    cout << "N = " << myindex.number() << endl;
    Index myindex2(myindex.number(),&frame);
    
    assert(my_p == myindex2.p());
    assert(my_j == myindex2.j());
    assert(my_e == myindex2.e());
    assert(my_k == myindex2.k());
    my_p = myindex.p();
    my_j = myindex.j();
    my_e = myindex.e();
    my_k = myindex.k();


    cout << "j0 " <<frame.frames()[0]->j0() << endl;
    cout << "deltasize " << frame.frames()[0]->Deltasize(frame.frames()[0]->j0())<<endl;
    cout << "deltaLmin " << frame.frames()[0]->DeltaLmin() << endl;
    cout << "firstgen " << first_q_generator(frame.frames()[0],frame.frames()[0]->j0()) << endl;
    cout << Frame1d::Index(my_p[0],my_j[0],my_e[0],my_k[0],frame.frames()[0]) << endl;

    SampledMapping<1> my_quarklet1d00(evaluate(*(frame.frames()[0]), Frame1d::Index(my_p[0],my_j[0],my_e[0],my_k[0],frame.frames()[0]), false,res ));
//      SampledMapping<1> my_quarklet1d00(evaluate(*(frame.frames()[0]), first_q_generator(frame.frames()[0],4), false,res ));
    SampledMapping<1> my_quarklet1d01(evaluate(*(frame.frames()[0]), Frame1d::Index(my_p[0],my_j[0],my_e[0],my_k[0],frame.frames()[0]), true ,res ));
    SampledMapping<1> my_quarklet1d10(evaluate(*(frame.frames()[1]), Frame1d::Index(my_p[1],my_j[1],my_e[1],my_k[1],frame.frames()[1]), false,res ));
    SampledMapping<1> my_quarklet1d11(evaluate(*(frame.frames()[1]), Frame1d::Index(my_p[1],my_j[1],my_e[1],my_k[1],frame.frames()[1]), true ,res ));

    ofstream file1d00("my_quarklet1d0false.m");
    ofstream file1d01("my_quarklet1d0true.m");
    ofstream file1d10("my_quarklet1d1false.m");
    ofstream file1d11("my_quarklet1d1true.m");

    my_quarklet1d00.matlab_output(file1d00);
    my_quarklet1d01.matlab_output(file1d01);
    my_quarklet1d10.matlab_output(file1d10);
    my_quarklet1d11.matlab_output(file1d11);

    file1d00.close();
    file1d01.close();
    file1d10.close();
    file1d11.close();

    if (DIM > 2)
    {
        SampledMapping<1> my_quarklet1d20(evaluate(*(frame.frames()[1]), Frame1d::Index(my_p[1],my_j[1],my_e[1],my_k[1],frame.frames()[1]), false,res ));
        SampledMapping<1> my_quarklet1d21(evaluate(*(frame.frames()[1]), Frame1d::Index(my_p[1],my_j[1],my_e[1],my_k[1],frame.frames()[1]), true ,res ));
        ofstream file1d20("my_quarklet1d2false.m");
        ofstream file1d21("my_quarklet1d2true.m");
        my_quarklet1d20.matlab_output(file1d20);
        my_quarklet1d21.matlab_output(file1d21);
        file1d20.close();
        file1d21.close();
    }

    cout << "Create matlab output via the tensorquarklet method" << endl;

    SampledMapping<DIM> my_quarklet0(evaluate(frame, myindex,false,res));
    SampledMapping<DIM> my_quarklet1(evaluate(frame, myindex,true,res));

    ofstream file0("my_quarkletfalse.m");
    ofstream file1("my_quarklettrue.m");
//    my_quarklet0.matlab_output(file0);
//    my_quarklet1.matlab_output(file1);
    file0.close();
    file1.close();
#endif
    
#if 0
    cout << "Check code for precomputation of first/last quarklets" << endl;
    
    MultiIndex<int,DIM> offset_it, level_it;
    Index temp_ind1, temp_ind2;
    for (unsigned int b=0; b< TFrameArray.size(); ++b)
    {
        cout << "Frame b = " << b << endl;
        for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
        while(multi_degree(offset_it) <= maxlevelrange)
        {
            for (unsigned int i=0; i<DIM; ++i)
                level_it[i] = offset_it[i] + TFrameArray[b]->j0()[i];
            //cout << "level = " << level_it << endl;
            temp_ind1 = TFrameArray[b]->first_quarklet( level_it );
            temp_ind2 = TFrameArray[b]->first_quarklet2( level_it );
            //cout << "j = " << level_it << "; first_quarklet" << temp_ind1 << "; first_quarklet2 = " << temp_ind2 << endl;
            assert ( (temp_ind1 == temp_ind2) && (temp_ind1.number() == temp_ind2.number()) );
            temp_ind1 = TFrameArray[b]->last_quarklet( level_it );
            temp_ind2 = TFrameArray[b]->last_quarklet2( level_it );
            //cout << "j = " << level_it << "; last_quarklet" << temp_ind1 << "; last_quarklet2 = " << temp_ind2 << endl;
            assert ( (temp_ind1 == temp_ind2) && (temp_ind1.number() == temp_ind2.number()) );
            ++offset_it;
        }
    }
            
    int repetitions(100);
    clock_t tstart, tend;
    double time1(0), time2(0);
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< TFrameArray.size(); ++b)
        {
            for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
            while(multi_degree(offset_it) <= maxlevelrange)
            {
                for (unsigned int i=0; i<DIM; ++i)
                {
                    level_it[i] = offset_it[i] + TFrameArray[b]->j0()[i];
                }
                temp_ind1 = TFrameArray[b]->first_quarklet( level_it );
                temp_ind2 = TFrameArray[b]->last_quarklet( level_it );
                ++offset_it;
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< TFrameArray.size(); ++b)
        {
            for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
            while(multi_degree(offset_it) <= maxlevelrange)
            {
                for (unsigned int i=0; i<DIM; ++i)
                {
                    level_it[i] = offset_it[i] + TFrameArray[b]->j0()[i];
                }
                temp_ind1 = TFrameArray[b]->first_quarklet2( level_it );
                temp_ind2 = TFrameArray[b]->last_quarklet2( level_it );
                ++offset_it;
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; old code " << time1 << "sec; new code " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
    
#endif
    return 0;
}
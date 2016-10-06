#include <iostream>


#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <utils/multiindex.h>
#include <cube/tbasis_evaluate.h>
#include <interval/p_basis.h>
#include <interval/p_evaluate.h>
#include <interval/ds_basis.h>
#include <interval/ds_evaluate.h>
using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
    cout << "Testing tbasis." << endl;

    //typedef DSBasis<2,2> Basis1d;
    typedef PBasis<2,2> Basis1d;
//    typedef Basis1d::Index IIndex;
    const unsigned int DIM = 3;
    typedef TensorBasis<Basis1d,DIM> Basis;
//    typedef Basis::Support Support;
    typedef Basis::Index Index;
    typedef Index::level_type level_type;
//    typedef Index::type_type type_type;
//    typedef Index::translation_type translation_type;


    const unsigned int tempj(9);
    

    cout << "Testing constructors" << endl;
    cout << "default constructor (with homogeneous b.c.)" << endl;
    Basis basisH;
    cout << "done" << endl;

    //
    cout << "basisH ...";
    basisH.set_jmax(tempj);
    cout << "done" << endl;

    cout << "Constructor with specified boundary condition orders"<<endl;
    cout << "Initialize FixedArray<int,4> with b.c.'s" << endl;

    FixedArray1D<int,2*DIM> s; // ,st; // BC for the dual basis can only be specified for dsbasis, not for pbasis
    for (unsigned int i=0; i<2*DIM; i++)
    {
        s[i]=0;
        //st[i]=0;
    }
    //s[0] = 1;
    //s[2*DIM-1]=1;

    cout << "Call constructor with s" << endl;
    Basis basisS(s);
    cout << "done" << endl;

    //
    cout << "basisS ...";
    basisS.set_jmax(tempj);
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
    Basis basisBC(bc);
    cout << "done" << endl;

    //
    cout << "basisBC ...";
    basisBC.set_jmax(tempj);
    cout << "done" << endl;
    
    cout << "Constructor with precomputed instances of the 1D bases" << endl;
    cout << "Initialize 1d bases" << endl;

    //Basis1d bas0();
    Basis1d bas00(false,false);
    Basis1d bas10(true,false);
    Basis1d bas01(false,true);
    Basis1d bas11(true,true);
    bas00.set_jmax(5);
    bas10.set_jmax(5);
    bas01.set_jmax(5);
    bas11.set_jmax(5);
    

    FixedArray1D<Basis1d*,DIM> basesArray;
    for (unsigned int i=0; i<DIM; i++)
    {
        basesArray[i] = &bas00;
    }
    //cout << "test: " << basisS.j0() << endl;
    //basesArray[0] = &bas10;
    //cout << "test: " << basisS.j0() << endl;
    //basesArray[DIM-1] = &bas01;
    
    cout << "Call constructor with array of 1d bases" << endl;
    Basis basisBas(basesArray);
    cout << "done" << endl;
    
    //
    cout << "basisBas ...";
    basisBas.set_jmax(tempj);
    cout << "done" << endl;
    

    // Use this basis for testing:
    //Basis basis;
    //Basis basis(s);
    //Basis basis(s,st);
    //Basis basis(bc);
    Basis basis(basesArray);

    // test equality of the constructed bases
    //FixedArray1D<Basis*,4> tensorbases;
    //tensorbases[0] = &basisH; tensorbases[1] = &basisS; tensorbases[2] = &basisBC; tensorbases[3] = &basisBas;
    FixedArray1D<Basis*,3> tensorbases;
    tensorbases[0] = &basisS; tensorbases[1] = &basisBC; tensorbases[2] = &basisBas;
    
    unsigned int maxlevelrange(4); // levels up to \|j0\|+maxlevelrange are tested
    unsigned int jmax(multi_degree(basis.j0())+maxlevelrange);
    
    // Run tests with this basis:

    basis.set_jmax(multi_degree(basis.j0())+maxlevelrange);
    basisBas.set_jmax(multi_degree(basisBas.j0())+maxlevelrange);
    
    basesArray[0] = &bas10;
    Basis basisBas2(basesArray);
    basisBas2.set_jmax(multi_degree(basisBas2.j0())+maxlevelrange);
    basesArray[0] = &bas01;
    Basis basisBas3(basesArray);basisBas3.set_jmax(multi_degree(basisBas3.j0())+maxlevelrange);
    basesArray[DIM-1] = &bas10;
    Basis basisBas4(basesArray);basisBas4.set_jmax(multi_degree(basisBas4.j0())+maxlevelrange);
    basesArray[DIM-1] = &bas01;
    Basis basisBas5(basesArray);basisBas5.set_jmax(multi_degree(basisBas5.j0())+maxlevelrange);
    basesArray[0] = &bas11; basesArray[DIM-1] = &bas01;
    Basis basisBas6(basesArray);basisBas6.set_jmax(multi_degree(basisBas6.j0())+maxlevelrange);
    
    FixedArray1D<Basis*,6> TBasisArray;
    TBasisArray[0] = &basisBas;
    TBasisArray[1] = &basisBas2;
    TBasisArray[2] = &basisBas3;
    TBasisArray[3] = &basisBas4;
    TBasisArray[4] = &basisBas5;
    TBasisArray[5] = &basisBas6;
    
    
    cout << "Testing set_jmax" << endl;
    cout << "basisH ...";
    basisH.set_jmax(jmax);
    cout << "done" << endl;
    cout << "basisS ...";
    basisS.set_jmax(jmax);
    cout << "done" << endl;
    cout << "basisBC ...";
    basisBC.set_jmax(jmax);
    cout << "done" << endl;
    cout << "basisBas ...";
    basisBas.set_jmax(jmax);
    cout << "done" << endl;
   
    for (unsigned int i=0; i< tensorbases.size(); ++i)
    {
        cout << "jmax= " << jmax << "; " << tensorbases[i]->last_wavelet(tempj).number() << "; " << tensorbases[i]->last_wavelet(jmax) << endl;
    }
       
    cout << "" << endl;
    
#if 0
    for (unsigned int b=0; b<3; ++b)
    {
        cout << "List differences of tensorbasis " << b << " and " << (b+1) << endl;
        if (tensorbases[b]->degrees_of_freedom() 
                != tensorbases[b+1]->degrees_of_freedom())
        {
            cout << " DOF[" << b << "] = " << tensorbases[b]->degrees_of_freedom() << "; "
                 << " DOF[" << (b+1) << "] = " << tensorbases[b+1]->degrees_of_freedom() << endl;
        }
        
        if (tensorbases[b]->j0()
                != tensorbases[b+1]->j0())
        {
            cout << " j0[" << b << "] = " << tensorbases[b]->j0() << "; "
                 << " j0[" << (b+1) << "] = " << tensorbases[b+1]->j0() << endl;
        }
        
        if (tensorbases[b]->get_jmax()
                != tensorbases[b+1]->get_jmax())
        {
            cout << " get_jmax[b=" << b << "] = " << tensorbases[b]->get_jmax() << "; "
                 << " get_jmax[b=" << (b+1) << "] = " << tensorbases[b+1]->get_jmax() << endl;
        }
        
        for (unsigned int d=0; d<DIM; d++)
        {
            cout << "List differences of tensorbasis " << b << " and " << (b+1) << " in direction " << d << endl;
            if (tensorbases[b]->bases()[d]->DeltaLmin() 
                    != tensorbases[b+1]->bases()[d]->DeltaLmin()  )
            {
                cout << " DeltaLmin()[b=" << b     << "][d=" << d << "] = " << tensorbases[b]->bases()[d]->DeltaLmin()   << "; "
                     << " DeltaLmin()[b=" << (b+1) << "][d=" << d << "] = " << tensorbases[b+1]->bases()[d]->DeltaLmin() << endl;
            }
            if (tensorbases[b]->bases()[d]->DeltaLmax() 
                    != tensorbases[b+1]->bases()[d]->DeltaLmax()  )
            {
                cout << " DeltaLmax()[b=" << b     << "][d=" << d << "] = " << tensorbases[b]->bases()[d]->DeltaLmax()   << "; "
                     << " DeltaLmax()[b=" << (b+1) << "][d=" << d << "] = " << tensorbases[b+1]->bases()[d]->DeltaLmax() << endl;
            }
            
            
            if (tensorbases[b]->bases()[d]->Nablamin()
                    != tensorbases[b+1]->bases()[d]->Nablamin()  )
            {
                cout << " Nablamin()[b=" << b     << "][d=" << d << "] = " << tensorbases[b]->bases()[d]->Nablamin()   << "; "
                     << " Nablamin()[b=" << (b+1) << "][d=" << d << "] = " << tensorbases[b+1]->bases()[d]->Nablamin() << endl;
            }
            
            
            for (unsigned int j= multi_degree(tensorbases[0]->j0()); j< (multi_degree(tensorbases[0]->j0()) + maxlevelrange); ++j)
            {
                cout << "List differences of tensorbasis " << b << " and " << (b+1) << " in direction " << d << " on level " << j << endl;
                if (tensorbases[b]->bases()[d]->Deltasize(j)
                        != tensorbases[b+1]->bases()[d]->Deltasize(j)  )
                {
                    cout << " Deltasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorbases[b]->bases()[d]->Deltasize(j)  << "; "
                         << " Deltasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorbases[b+1]->bases()[d]->Deltasize(j) << endl;
                }
                
                if (tensorbases[b]->bases()[d]->Nablamax(j)
                        != tensorbases[b+1]->bases()[d]->Nablamax(j)  )
                {
                    cout << " Nablamax(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorbases[b]->bases()[d]->Nablamax(j)  << "; "
                         << " Nablamax(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorbases[b+1]->bases()[d]->Nablamax(j) << endl;
                }
                
                if (tensorbases[b]->bases()[d]->Nablasize(j)
                        != tensorbases[b+1]->bases()[d]->Nablasize(j)  )
                {
                    cout << " Nablasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorbases[b]->bases()[d]->Nablasize(j)  << "; "
                         << " Nablasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorbases[b+1]->bases()[d]->Nablasize(j) << endl;
                }
                if (tensorbases[b]->bases()[d]->Deltasize(j)
                        != tensorbases[b+1]->bases()[d]->Deltasize(j)  )
                {
                    cout << " Deltasize(j)[b=" << b     << "][d=" << d << "][j=" << j << "] = " << tensorbases[b]->bases()[d]->Deltasize(j)  << "; "
                         << " Deltasize(j)[b=" << (b+1) << "][d=" << d << "][j=" << j << "] = " << tensorbases[b+1]->bases()[d]->Deltasize(j) << endl;
                }
            }
        }
    }
    
#endif

    
    /*        
    cout << "Testing j0()" << endl
         << "basisH " << basisH.j0() << endl
         << "basisS " << basisS.j0() << endl
//             << "basisSST "<< basisSST.j0() << endl
         << "basisBC "<< basisBC.j0() << endl
         << "basisBas "<< basisBas.j0() << endl;
*/
    

#if 0
    for (unsigned int b=0; b<4; ++b)
    {
    cout << "b = " << b 
            << "; last_wavelet = " << tensorbases[b]->last_wavelet(tensorbases[b]->get_jmax()) 
            << "; .number() = " << tensorbases[b]->last_wavelet(tensorbases[b]->get_jmax()).number() 
            << "; dof = " << tensorbases[b]->degrees_of_freedom() << endl;
    }
#endif

#if 0
    // see test_tbasis_support.cpp
    cout << "Testing support" << endl;
    // void support(const Index& lambda, Support& supp) const;
    Support supp;

    //cout << "DeltaLmax()  " << basis.bases()[0]->DeltaLmax() << " " << basis.bases()[1]->DeltaLmax() << endl;
    //cout << basisH.bases()[0]->DeltaLmax() << endl;
    //cout << "basis.first_generator();"  << basis.first_generator()<< endl;
    //basisH.first_generator();

    int maxnumber(10);
    for (Index it(basis.first_generator());it < basis.last_wavelet(jmax);++it)
    {
        basis.support(it,supp);
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

    // static double primal_regularity() { return IBASIS::primal_regularity(); }
    // static unsigned int primal_polynomial_degree() { return IBASIS::primal_polynomial_degree(); }
    // static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }
    cout << "Testing primal_ functions"<<endl;
    cout << "primal_regularity "<< basis.primal_regularity() << endl
         << "primal_polynomial_degree " << basis.primal_polynomial_degree() << endl
         << "primal_vanishing_moments " << basis.primal_vanishing_moments() << endl;
#endif

    cout << "Testing first/last routines" << endl;
    Index tempi1, tempi2;
    tempi1 = basis.first_generator();
    tempi2 = first_generator<Basis1d,DIM,Basis> (&basis);
    assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ) );

    tempi1 = basis.last_generator();
    tempi2 = last_generator<Basis1d,DIM,Basis> (&basis);
    assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

    level_type j0(basis.j0());
    level_type currentlevel(j0);
    int range(0), maxrange(2); // number of level steps we have climbed so far.

    while (range <= maxrange)
    {
        tempi1 = basis.first_wavelet(currentlevel);
        tempi2 = first_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

        tempi1 = basis.last_wavelet(currentlevel);
        tempi2 = last_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);
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

    cout << "Testing first_wavelet" << endl;
    unsigned int maxrange2(2);
    assert (multi_degree(basis.j0()) + maxrange2 <= basis.get_jmax());
    
    for (unsigned int lev = multi_degree(basis.j0()); lev < multi_degree(basis.j0()) + maxrange2; lev++)
    {
        tempi1 = basis.first_wavelet(lev);
        tempi2 = first_wavelet<Basis1d,DIM,Basis>(&basis,lev);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));

        tempi1 = basis.last_wavelet(lev);
        tempi2 = last_wavelet<Basis1d,DIM,Basis>(&basis,lev);
        assert( (tempi1 == tempi2) && (tempi1.number() == tempi2.number() ));
    }
    
    
#if 0
    
    cout << "Comparing speed of first/last routines" << endl
            << " - in tbasis.h and tbasis_index.h" << endl;
    Index tempi3, tempi4;
    int range2(0), maxrange3(12), repetitions(1000);
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; rep++)
    {
        //cout << "testing tbasis.h; rep = " << rep << "; range2 = ";
        range2=0;
        currentlevel=j0;
        while (range2 <= maxrange3)
        {
            //cout << range2 << "; ";
            tempi3 = basis.first_wavelet(currentlevel);
            tempi4 = basis.last_wavelet(currentlevel);
            //tempi3 = first_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);
            //tempi4 = last_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);

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
            tempi3 = first_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);
            tempi4 = last_wavelet<Basis1d,DIM,Basis>(&basis, currentlevel);

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
            << "tbasis.h: " << time1 << "sec; tbasis_index.h: " << time2 << "sec" << endl;
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
    IIndex index1d(first_generator(basis.bases()[0], basis.bases()[0]->j0()));
    coeff1d.set_coefficient(index1d,1.0);
    int dil=8;
    Array1D<double> grid1d;
    Array1D<double> funval1d;
    grid1d.resize((1<<dil)+1);
    for(unsigned int k=0;k<grid1d.size();k++)
        grid1d[k]=k*(1.0/(1<<dil));
    Grid<1>gitter(grid1d);
    //evaluate(bas11,0,index1d,grid1d,funval1d);
    //evaluate(*(basesArray[0]),0,index1d,grid1d,funval1d);
    evaluate(*(basis.bases()[0]),0,index1d,grid1d,funval1d);
    SampledMapping<1> my_wavelet1d (gitter, funval1d);
    ofstream fileA("my_wavelet1d.m");
    my_wavelet1d.matlab_output(fileA);
    fileA.close();

    cout << "Fix an tensorindex" << endl;
    cout << "compute by hand 1D samples corresponding to each 1d wavelet in the tensorwavelet" << endl;
    cout << "store them as matlab output" << endl;

    int res(7);

    level_type my_j;
    type_type my_e;
    translation_type my_k;

    my_j = basis.j0();
    //my_j[0]=basis.bases()[0]->j0(); // == 4
    //my_j[1]=basis.bases()[1]->j0();
    for (unsigned int i=0; i<DIM; i++)
    {
        assert(my_e[i] == 0);
        my_k[i] = basis.bases()[i]->DeltaLmin();
    }
    my_k[DIM-1] += 1;
    //my_e[0]=0;
    //my_e[1]=0;
    //my_k[0]=basis.bases()[0]->DeltaLmin()+0; // DeltaLmin ==2, Deltasize(4) == 12
    //my_k[1]=basis.bases()[1]->DeltaLmin()+1;

    //TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
    Index myindex (my_j,my_e,my_k,&basis);


    // simpler:
    cout << "N = " << myindex.number() << endl;
    Index myindex2(myindex.number(),&basis);
    assert(my_j == myindex2.j());
    assert(my_e == myindex2.e());
    assert(my_k == myindex2.k());
    my_j = myindex.j();
    my_e = myindex.e();
    my_k = myindex.k();


    cout << "j0 " <<basis.bases()[0]->j0() << endl;
    cout << "deltasize " << basis.bases()[0]->Deltasize(basis.bases()[0]->j0())<<endl;
    cout << "deltaLmin " << basis.bases()[0]->DeltaLmin() << endl;
    cout << "firstgen " << first_generator(basis.bases()[0],basis.bases()[0]->j0()) << endl;
    cout << Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]) << endl;

    SampledMapping<1> my_wavelet1d00(evaluate(*(basis.bases()[0]), Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]), false,res ));
//      SampledMapping<1> my_wavelet1d00(evaluate(*(basis.bases()[0]), first_generator(basis.bases()[0],4), false,res ));
    SampledMapping<1> my_wavelet1d01(evaluate(*(basis.bases()[0]), Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]), true ,res ));
    SampledMapping<1> my_wavelet1d10(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), false,res ));
    SampledMapping<1> my_wavelet1d11(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), true ,res ));

    ofstream file1d00("my_wavelet1d0false.m");
    ofstream file1d01("my_wavelet1d0true.m");
    ofstream file1d10("my_wavelet1d1false.m");
    ofstream file1d11("my_wavelet1d1true.m");

    my_wavelet1d00.matlab_output(file1d00);
    my_wavelet1d01.matlab_output(file1d01);
    my_wavelet1d10.matlab_output(file1d10);
    my_wavelet1d11.matlab_output(file1d11);

    file1d00.close();
    file1d01.close();
    file1d10.close();
    file1d11.close();

    if (DIM > 2)
    {
        SampledMapping<1> my_wavelet1d20(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), false,res ));
        SampledMapping<1> my_wavelet1d21(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), true ,res ));
        ofstream file1d20("my_wavelet1d2false.m");
        ofstream file1d21("my_wavelet1d2true.m");
        my_wavelet1d20.matlab_output(file1d20);
        my_wavelet1d21.matlab_output(file1d21);
        file1d20.close();
        file1d21.close();
    }

    cout << "Create matlab output via the tensorwavelet method" << endl;

    SampledMapping<DIM> my_wavelet0(evaluate(basis, myindex,false,res));
    SampledMapping<DIM> my_wavelet1(evaluate(basis, myindex,true,res));

    ofstream file0("my_waveletfalse.m");
    ofstream file1("my_wavelettrue.m");
    my_wavelet0.matlab_output(file0);
    my_wavelet1.matlab_output(file1);
    file0.close();
    file1.close();
#endif
    
#if 0
    cout << "Check code for precomputation of first/last wavelets" << endl;
    
    MultiIndex<int,DIM> offset_it, level_it;
    Index temp_ind1, temp_ind2;
    for (unsigned int b=0; b< TBasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
        while(multi_degree(offset_it) <= maxlevelrange)
        {
            for (unsigned int i=0; i<DIM; ++i)
                level_it[i] = offset_it[i] + TBasisArray[b]->j0()[i];
            //cout << "level = " << level_it << endl;
            temp_ind1 = TBasisArray[b]->first_wavelet( level_it );
            temp_ind2 = TBasisArray[b]->first_wavelet2( level_it );
            //cout << "j = " << level_it << "; first_wavelet" << temp_ind1 << "; first_wavelet2 = " << temp_ind2 << endl;
            assert ( (temp_ind1 == temp_ind2) && (temp_ind1.number() == temp_ind2.number()) );
            temp_ind1 = TBasisArray[b]->last_wavelet( level_it );
            temp_ind2 = TBasisArray[b]->last_wavelet2( level_it );
            //cout << "j = " << level_it << "; last_wavelet" << temp_ind1 << "; last_wavelet2 = " << temp_ind2 << endl;
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
        for (unsigned int b=0; b< TBasisArray.size(); ++b)
        {
            for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
            while(multi_degree(offset_it) <= maxlevelrange)
            {
                for (unsigned int i=0; i<DIM; ++i)
                {
                    level_it[i] = offset_it[i] + TBasisArray[b]->j0()[i];
                }
                temp_ind1 = TBasisArray[b]->first_wavelet( level_it );
                temp_ind2 = TBasisArray[b]->last_wavelet( level_it );
                ++offset_it;
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< TBasisArray.size(); ++b)
        {
            for (unsigned int i=0; i<DIM; ++i)
            {
                offset_it[i] = 0;
            }
            while(multi_degree(offset_it) <= maxlevelrange)
            {
                for (unsigned int i=0; i<DIM; ++i)
                {
                    level_it[i] = offset_it[i] + TBasisArray[b]->j0()[i];
                }
                temp_ind1 = TBasisArray[b]->first_wavelet2( level_it );
                temp_ind2 = TBasisArray[b]->last_wavelet2( level_it );
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
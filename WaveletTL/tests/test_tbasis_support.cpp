#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/jl_basis.h>
#include <interval/jl_support.h>
#include <interval/jl_evaluate.h>

#include <interval/i_index.h>
#include <cube/tbasis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
    cout << "Testing support calculations for wavelet bases on the cube..." << endl;
    const int d  = 3;
    const int dT = 3;
    //const unsigned int DIM = 1; const int levelrange(4);
    const unsigned int DIM = 2; const int levelrange(4);
    //const unsigned int DIM = 3; const int levelrange(2);
    typedef PBasis<d,dT> Basis1D;
    //typedef DSBasis<d,dT> Basis1D;
//   typedef JLBasis Basis1D;
    
    typedef TensorBasis<Basis1D,DIM> Basis;
    typedef Basis::Index Index;
//    typedef TensorIndex<Basis1D,DIM>::level_type level_type;
//    typedef TensorIndex<Basis1D,DIM>::type_type type_type;
    typedef TensorIndex<Basis1D,DIM>::translation_type translation_type;
    typedef Basis1D::Index Index1D;
    
    FixedArray1D<bool,2*DIM> bc;
    for (unsigned int i=0; i< 2*DIM; ++i)
    {
        bc[i] = false;
    }
    //bc[0] = true;
    //bc[1] = true;
    //bc[2*DIM-2] = true;
    //bc[2*DIM-1] = true;
    //Basis basis(bc);
    //basis.set_jmax(multi_degree(basis.j0())+levelrange);
    
    Basis1D bas00(false,false); Basis1D bas10(true,false); Basis1D bas01(false,true); Basis1D bas11(true,true);
    bas00.set_jmax(5); bas10.set_jmax(5); bas01.set_jmax(5); bas11.set_jmax(5);
    FixedArray1D<Basis1D*,DIM> basesArray;
    for(unsigned int i=0;i<DIM;i++)
    {
        basesArray[i] = &bas00;
    }
    Basis basisBas1(basesArray);basisBas1.set_jmax(multi_degree(basisBas1.j0())+levelrange);
    basesArray[0] = &bas10;
    Basis basisBas2(basesArray);basisBas2.set_jmax(multi_degree(basisBas2.j0())+levelrange);
    basesArray[0] = &bas01;
    Basis basisBas3(basesArray);basisBas3.set_jmax(multi_degree(basisBas3.j0())+levelrange);
    basesArray[DIM-1] = &bas10;
    Basis basisBas4(basesArray);basisBas4.set_jmax(multi_degree(basisBas4.j0())+levelrange);
    if (DIM > 1) basesArray[1] = &bas01;
    Basis basisBas5(basesArray);basisBas5.set_jmax(multi_degree(basisBas5.j0())+levelrange);
    basesArray[0] = &bas11; basesArray[DIM-1] = &bas01;
    Basis basisBas6(basesArray);basisBas6.set_jmax(multi_degree(basisBas6.j0())+levelrange);

    FixedArray1D<Basis*,6> TBasisArray;
    TBasisArray[0] = &basisBas1;    TBasisArray[1] = &basisBas2;    TBasisArray[2] = &basisBas3;    TBasisArray[3] = &basisBas4;    TBasisArray[4] = &basisBas5;    TBasisArray[5] = &basisBas6;
    
    //FixedArray1D<Basis*,1> TBasisArray;
    //TBasisArray[0] = &basisBas3;
    

#if 0
    /*
     * get_intersecting_wavelets_on_level does not work for DSBasis!
     * The following code is copied to the dsbasis_testfile
     */
    typedef Basis1D::Index Index1D;
    int tempA, tempB;
    //Index1D temp_mu(4,1,6,TBasisArray[0]->bases()[1]);
    Index1D first(TBasisArray[0]->bases()[1]->get_wavelet(0));
    //typedef Basis1D::Support Support1D;
    //Support1D supp1d;
    TBasisArray[0]->bases()[1]->support(first, tempA, tempB);
    cout << "wavelet = " << first << "; 1d support = (" << tempA << ", " << tempB << ")" << endl;
            
            
    Index1D temp_mu(TBasisArray[0]->bases()[1]->get_wavelet(20));
    cout << "temp_mu = " << temp_mu << endl;
    
    get_intersecting_wavelets_on_level(*(TBasisArray[0]->bases()[1]),
                    temp_mu,
            4,true,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
    get_intersecting_wavelets_on_level(*(TBasisArray[0]->bases()[1]),
                    temp_mu,
            4,false,tempA,tempB);
    cout << "temp_mu -> get_intersecting_wavelets (4, true) (min, max) = (" << tempA << ", " << tempB << ")" << endl;
     
    cout << "evaluate wavelet mu = " << temp_mu << endl;
    for (unsigned int i=0; i< 100; ++i)
    {
        cout << "f(" << (0.01*i) << ") = " << TBasisArray[0]->bases()[1]->evaluate(0, temp_mu, 0.01*i) << endl;
    }

    cout << "min = " << tempA << "; max = " << tempB << endl;
    for (unsigned int i = 0; i< 40; ++i)
    {
        Index1D temp_index(TBasisArray[0]->bases()[1]->get_wavelet(i));
        TBasisArray[0]->bases()[1]->support(temp_index, tempA, tempB);
        cout << "N = " << i << "; lam = " << temp_index << "; k1 = " << tempA << "; k2 = " << tempB << endl;
    }
    abort();
#endif
    
#if 1
    cout << "- testing calculation of supports:" << endl;
    Basis::Support supp;
    //int maxnumber(203);
    int minnumber(0);
    int maxnumber(1000);
    //for (unsigned int b = 0; b < TBasisArray.size(); b++)
    for (unsigned int b = 5; b < 6; b++)
    {
        cout << "Basis Nr b = " << b << endl;
        for (Index lambda(TBasisArray[b]->get_wavelet(minnumber));; ++lambda) 
        {
            support<Basis1D,DIM>(*TBasisArray[b], lambda, supp);
            //support(basis, lambda, supp);
            cout << "    N = " << lambda.number() << " == " << lambda << "; support = 2^{-(" << supp.j[0];
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

            //if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0()+1)) break;
            if ((int)lambda.number() == maxnumber) break;
        }
    }

    double aa(-0.1), bb(1.1); 
    int steps (100);
    double h = (bb-aa) / steps;
    Index temp_ind;
    int basnum(5), waveletnum(248); // test a suspicious wavelet
    {
        temp_ind = TBasisArray[basnum]->get_wavelet(waveletnum);
        cout << "Testing support of" << endl; // lambda = " << temp_ind << "; .number() = " << temp_ind.number() <<  endl;
        support<Basis1D,DIM>(*TBasisArray[basnum], temp_ind, supp);
        //support(basis, lambda, supp);
        cout << "    N = " << temp_ind.number() << " == " << temp_ind << "; support = 2^{-(" << supp.j[0];
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
        
        Index1D temp_ind1D(temp_ind.j()[1], temp_ind.e()[1], temp_ind.k()[1], TBasisArray[basnum]->bases()[1]);
        int k11,k12,k21,k22;
        TBasisArray[basnum]->bases()[1]->support(temp_ind.j()[1], temp_ind.e()[1], temp_ind.k()[1],k21,k22);
        TBasisArray[basnum]->bases()[1]->support(temp_ind1D,k11,k12);
        assert ((k11 == k21) && (k12 == k22));
        cout << "temp_ind1D. has support j = " << temp_ind1D.j() << "; e = " << temp_ind1D.e() << "; k1 = " << k11 << "; k2 = " << k12 << endl;
        
        for (unsigned int i=0; (int)i<steps; ++i)
        {
            double val1, val2;
            double x(aa+i*h);
            val1 = TBasisArray[basnum]->bases()[1]->evaluate(0,
                    temp_ind1D.j(), temp_ind1D.e(), temp_ind1D.k(),
                    x);
            val2 = TBasisArray[basnum]->bases()[1]->evaluate(0,
                    temp_ind1D,
                    x);
            assert (val1 == val2);
            if ((val1 != 0) || (val2 != 0))
            {
                cout << "i=" << i << "; x = " << x << "; val1 = " << val1 << "; val2 = " << val2 <<  endl;
            }
        }
        
        for (unsigned int i=52; (int)i<steps; ++i)
        {
            Point<DIM> point(aa+i*h, aa+i*h);
            cout << "x = " << point << endl;
            //Point<DIM> point(aa+i*h, 1-aa-i*h);
            
            double value1(1.0), value2(1.0);

            cout << "j[0] = " << temp_ind.j()[0] << "; e[0] = " << temp_ind.e()[0] << "; k[0] = " << temp_ind.k()[0] << "; x[0] = " << point[0] << endl;
            value1 *= TBasisArray[basnum]->bases()[0]->evaluate(0,
                    temp_ind.j()[0], temp_ind.e()[0], temp_ind.k()[0],
                    point[0]);
            cout << "j[1] = " << temp_ind.j()[1] << "; e[1] = " << temp_ind.e()[1] << "; k[1] = " << temp_ind.k()[1] << "; x[1] = " << point[1] << endl;
            value2 *= TBasisArray[basnum]->bases()[1]->evaluate(0,
                    temp_ind.j()[1], temp_ind.e()[1], temp_ind.k()[1],
                    point[1]);
            double value3;
            value3 *= TBasisArray[basnum]->bases()[1]->evaluate(0,
                    Index1D(temp_ind.j()[1], temp_ind.e()[1], temp_ind.k()[1], TBasisArray[basnum]->bases()[1]),
                    point[1]);
            InfiniteVector<double,int> gcoeffs;
            double r(0);
            TBasisArray[basnum]->bases()[1]->reconstruct_1(temp_ind.j()[1], temp_ind.e()[1], temp_ind.k()[1], temp_ind.j()[1]+1, gcoeffs); 
            const int Lmin(TBasisArray[basnum]->bases()[1]->DeltaLmin());
            for (InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                    it != gcoeffs.end(); ++it)
            {
                // gcoeffs contains only coeffs related to generators on level j+1
                // j_ = j+1; e_ = 0; k_ = DeltaLmin() + num
                r += *it * evaluate(*(TBasisArray[basnum]->bases()[1]), 0, temp_ind.j()[1]+1,0,Lmin+it.index(), point[1]);
            }
      
            double temp_d = TBasisArray[basnum]->evaluate(0, temp_ind, point);
            double temp_d2 = TBasisArray[basnum]->evaluate(1, temp_ind, point);

             if  ( (abs(temp_d) != 0) || (abs(temp_d2) != 0) )
                cout << "i=" << i << "; x = " << point << "; evaluate(0,x) = " << temp_d << "; evaluate(1,x) = " << temp_d2 <<  endl;
        }
    }
     
    cout << "" << endl;
#endif

//    typedef TensorIndex<Basis1D,DIM,Basis>::level_type level_type;
//    typedef std::list<Index> SupportList;
    
#if 0
    for (unsigned int b = 0; b < TBasisArray.size(); b++)
    {
        Index lambda(6000,TBasisArray[b]);

        level_type lambda_j(TBasisArray[b]->j0());

        SupportList intwav1, intwav2, intwav3, intel1, intel2;

        Basis::Support supp_lambda;
        TBasisArray[b]->support(lambda, supp_lambda);
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


        Index lambda1(81,TBasisArray[b]);
        Basis::Support supp_lambda1;
        TBasisArray[b]->support(lambda1, supp_lambda1);
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

        Basis::Support temp_supp;
        bool temp_bool;
        temp_bool = intersect_supports(*TBasisArray[b], lambda, lambda1, temp_supp);
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
        cout << "intersecting_generators for lambda = " << lambda << "; .number() = " << lambda.number() << "; on level = " << lambda_j << endl;
        intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda, lambda_j, true, intwav1);
        cout << "Number of intersecting generators is " << intwav1.size() << endl;
        for (SupportList::const_iterator it(intwav1.begin()); it != intwav1.end(); ++it)
        {
            Basis::Support supp, supp_it;
            TBasisArray[b]->support(*it, supp_it);
            assert (intersect_supports(*TBasisArray[b], lambda, *it, supp) == true);
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
        cout << "intersecting wavelets on level j0 = " << TBasisArray[b]->j0() << endl;
        intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda, lambda_j, false, intwav2);
        cout << "Number of intersecting wavelets is " << intwav2.size() << endl;
        for (SupportList::const_iterator it(intwav2.begin()); it != intwav2.end(); ++it)
        {
            Basis::Support supp, supp_it;
            TBasisArray[b]->support(*it, supp_it);
            intersect_supports(*TBasisArray[b], lambda, *it, supp);
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
        intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda, lambda_j, false, intwav3);
        cout << "Number of intersectiong wavelets is " << intwav3.size() << endl;
        for (SupportList::const_iterator it(intwav3.begin()); it != intwav3.end(); ++it)
        {
            Basis::Support supp, supp_it;
            /*
            if ((*it).number() == 361)
            {
                Basis::Support temp_supp, temp_inters;
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
            TBasisArray[b]->support(*it, supp_it);
            intersect_supports(*TBasisArray[b], lambda, *it, supp);
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

        // throws assertion
        // intersecting_wavelets<Basis1D,DIM>(basis, lambda, lambda_j, true, nus);

        lambda_j = TBasisArray[b]->j0();

        cout << " Testing intersecting_elements on level j0 = " << TBasisArray[b]->j0() << lambda_j <<endl;
        intersecting_elements<Basis1D,DIM>(*TBasisArray[b], lambda, lambda_j, intel1);
        cout << "Number of intersectiong elements is " << intel1.size() << endl;
        for (SupportList::const_iterator it(intel1.begin()); it != intel1.end(); ++it)
        {
            Basis::Support supp, supp_it;
            TBasisArray[b]->support(*it, supp_it);
            intersect_supports(*TBasisArray[b], lambda, *it, supp);
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
        cout << " Testing intersecting_elements on level "<< lambda_j <<endl;
        intersecting_elements<Basis1D,DIM>(*TBasisArray[b], lambda, lambda_j, intel2);
        cout << "Number of intersectiong elements is " << intel2.size() << endl;
        for (SupportList::const_iterator it(intel2.begin()); it != intel2.end(); ++it)
        {
            Basis::Support supp, supp_it;
            TBasisArray[b]->support(*it, supp_it);
            intersect_supports(*TBasisArray[b], lambda, *it, supp);
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

        cout << "intwav1.size() = " << intwav1.size() << endl;
        for (SupportList::const_iterator it(intwav1.begin()); it != intwav1.end(); ++it)
        {
            cout << "N = " << (*it).number() << " == " << *it << endl;
        }
        cout << "intwav2.size() = " << intwav2.size() << endl;
        for (SupportList::const_iterator it(intwav2.begin()); it != intwav2.end(); ++it)
        {
            cout << "N = " << (*it).number() << " == " << *it << endl;
        }
        cout << "intwav3.size() = " << intwav3.size() << endl;
        for (SupportList::const_iterator it(intwav3.begin()); it != intwav3.end(); ++it)
        {
            cout << "N = " << (*it).number() << " == " << *it << endl;
        }
        cout << "intel1.size()  = " << intel1.size()  << endl;
        for (SupportList::const_iterator it(intel1.begin()); it != intel1.end(); ++it)
        {
            cout << "N = " << (*it).number() << " == " << *it << endl;
        }
        cout << "intel2.size()  = " << intel2.size()  << endl;
        for (SupportList::const_iterator it(intel2.begin()); it != intel2.end(); ++it)
        {
            cout << "N = " << (*it).number() << " == " << *it << endl;
        }
    }
#endif

#if 0
    // intersecting_wavelets and intersecting_elements produce the same sets
    // (on level j0: of course intersecting _wavs and _gens together)
    // in this test not the sets are compared, but the number of elements in the sets
    // proof:
            
    cout << "Testing intersecting_wavelets and intersecting_elements" << endl;
    for (unsigned int b = 0; b < TBasisArray.size(); b++)
    {
        cout << "Basis Nr b = " << b << endl;
        //Index lambda_run(basis.first_generator());
        level_type j0(TBasisArray[b]->j0());
        level_type currentlevel(j0);
        int range(0), minlambdanum(6000), maxlambdanum(6005), count(0); // number of level steps we have climbed so far.
        TBasisArray[b]->set_jmax(multi_degree(TBasisArray[b]->j0())+levelrange);
        assert(maxlambdanum < TBasisArray[b]->degrees_of_freedom());
        Index lambda_run(TBasisArray[b]->get_wavelet(minlambdanum));
        SupportList nus;
        bool isordered = true;
        SupportList::const_iterator it2;

        cout << "levelrange = " << levelrange << "; minlambdanum = " << minlambdanum << "; maxlambdanum = " << maxlambdanum << endl;

        while (lambda_run.number() <= maxlambdanum)
        {
            // compute intersecting_wavelets and intersecting_elements for all levels up to \|j0\|+maxrange
            //cout << "lambda_run = " << lambda_run << "; currentlevel = ";
            //cout.flush();
            while (range <= levelrange)
            {
                //cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                //cout  << currentlevel << "; ";
                //cout.flush();
                cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                count = 0;
                if (range == 0)
                {
                    nus.clear();
                    assert (nus.size() == 0 );
                    cout << "main: intersecting_gens:" << endl;
                    //intersecting_wavelets2<Basis1D,DIM>(basis, lambda_run, currentlevel, true, nus);
                    //cout << "intersecting_generators2 found " << nus.size() << " intersections" << endl;

                    intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run, currentlevel, true, nus);
                    cout << "intersecting_generators found " << nus.size() << " intersections" << endl;
                    count += nus.size();
                    //abort();
              /*      
               
                    cout << "nus =" << endl;
                    for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it)
                    {
                        cout << " N = " << (*it).number() << "; " << (*it) << endl;
                    }
              */      
            //   /*     
                    isordered = true;
                    it2 = nus.begin();
                    ++it2;
                    for (SupportList::const_iterator it(nus.begin()); it2 != nus.end(); ++it, ++it2)
                    {
                        isordered = (isordered && ((*it) < (*it2)) );
                    }
                    if (!isordered)
                    {
                        cout << "intersecting_generators is NOT ordered" << endl;
                        cout << "N = " << lambda_run.number() << "; lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                        abort();
                    }
         //   */        

                }
                nus.clear();
                assert (nus.size() == 0);
                cout << "main: intersecting_wavs:" << endl;
                //intersecting_wavelets2<Basis1D,DIM>(basis, lambda_run, currentlevel, false, nus);
                //cout << "intersecting_wavelets2 found " << nus.size() << " intersections" << endl;
                intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run, currentlevel, false, nus);
                cout << "intersecting_wavelets found " << nus.size() << " intersections" << endl;

                count += nus.size();
              /*      
                //cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                cout << "nus =" << endl;
                for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it)
                {
                    cout << " N = " << (*it).number() << "; " << (*it) << endl;
                }
              */  

       //   /*      
                isordered = true;
                it2 = nus.begin();
                ++it2;
                for (SupportList::const_iterator it(nus.begin()); it2 != nus.end(); ++it, ++it2)
                {
                    isordered = (isordered && ((*it) < (*it2)) );
                }
                if (!isordered)
                {
                    cout << "intersecting_wavelets is NOT ordered" << endl;
                    cout << "N = " << lambda_run.number() << "; lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                    abort();
                }
      //     */            

                //cout << "#nus = " << nus.size() << " ";

                nus.clear();
                assert (nus.size() == 0);
                
                // intersecting_elements<Basis1D,DIM>(*TBasisArray[b], lambda_run, currentlevel, nus);
                count -= nus.size();
                cout << "intersecting_elements found " << nus.size() << " intersections" << endl;
            /*
                cout << "nus =" << endl;
                for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it)
                {
                    cout << " N = " << (*it).number() << "; " << (*it) << endl;
                }
            */
            /*
                isordered = true;
                it2 = nus.begin();
                ++it2;
                for (SupportList::const_iterator it(nus.begin()); it2 != nus.end(); ++it, ++it2)
                {
                    isordered = (isordered && ((*it) < (*it2)) );
                }
                if (!isordered)
                {
                    cout << "intersecting_elements is NOT ordered" << endl;
                    cout << "N = " << lambda_run.number() << "; lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                    abort();
                }
            */              

                //cout << "#nus = " << nus.size() << " ";
                if (count != 0)
                {
                    cout << "Problem detected" <<  endl;
                    cout <<  "intersecting_wavelets and intersecting_elements produce different sets" << endl;
                    cout << "lambda_run = " << lambda_run << "; .number() = " << lambda_run.number() << "; currentlevel = " << currentlevel << endl;
                    abort();
                }
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
            ++lambda_run;
            currentlevel = j0;
            range = 0;
            cout << endl;
            cout.flush();
        }
        //cout << "halt" << endl;
    }
#endif
    
#if 0
    // methods are equally fast! intersect_elements is marginally slower
    cout << "Comparing speed of intersect_wavelets and intersect_elements" << endl;
    level_type j02;
    level_type currentlevel2;
    SupportList nus2;
    int range2(0), minlambdanum2(7995), maxlambdanum2(8000), count2(0), repetitions(5);
    clock_t tstart, tend;
    double time1(0), time2(0);
    Index lambda_run2;
    
    tstart = clock();
    
    for (unsigned int b=0; b<TBasisArray.size(); ++b)
    {
        j02 = TBasisArray[b]->j0();
        currentlevel2 = j02;
        assert(maxlambdanum2 < TBasisArray[b]->degrees_of_freedom());

        for (int rep(0); rep < repetitions; rep++)
        {
            lambda_run2 = TBasisArray[b]->get_wavelet(minlambdanum2);
            while (lambda_run2.number() <= maxlambdanum2)
            {
                count2 = 0;
                intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run2, currentlevel2, true, nus2);
                count2 += nus2.size();
                nus2.clear();
                while (range2 <= levelrange)
                {
                    intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run2, currentlevel2, false, nus2);
                    count2 += nus2.size();
                    nus2.clear();

                    //intersecting_elements<Basis1D,DIM>(basis, lambda_run, currentlevel, nus);
                    //nus.clear();

                    // increase the level j,  meaning: "currentlevel++"
                    //cout << currentlevel << endl;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        if (i != 0)
                        {
                            if (currentlevel2[i] != j02[i])
                            {
                                // increase left neighbor
                                currentlevel2[i-1]=currentlevel2[i-1]+1;
                                int temp = currentlevel2[i]-j02[i];
                                currentlevel2[i]=j02[i];
                                currentlevel2[DIM-1]=j02[DIM-1]+temp-1;
                                break;
                            }
                        } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                        {
                            range2 = range2 +1;
                            if (DIM == 1)
                            {
                                currentlevel2[0]=currentlevel2[0]+1;
                            }
                            else
                            {
                                //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                currentlevel2[DIM-1]=j02[DIM-1]+range2;
                                currentlevel2[0]=j02[0];
                            }
                        }
                    } // end of "curentlevel++"
                } // end of while (range < maxrange)
                ++lambda_run2;
                currentlevel2 = j02;
                range2 = 0;
            }
        }
    }
    
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    for (unsigned int b=0; b<TBasisArray.size(); ++b)
    {
        for (int rep(0); rep < repetitions; rep++)
        {
            lambda_run2 = TBasisArray[b]->get_wavelet(minlambdanum2);
            while (lambda_run2.number() <= maxlambdanum2)
            {
                count2 = 0;
                while (range2 <= levelrange)
                {
                    //intersecting_elements<Basis1D,DIM>(*TBasisArray[b], lambda_run2, currentlevel2, nus2);
                    count2 += nus2.size();
                    nus2.clear();

                    // increase the level j,  meaning: "currentlevel++"
                    //cout << currentlevel << endl;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        if (i != 0)
                        {
                            if (currentlevel2[i] != j02[i])
                            {
                                // increase left neighbor
                                currentlevel2[i-1]=currentlevel2[i-1]+1;
                                int temp = currentlevel2[i]-j02[i];
                                currentlevel2[i]=j02[i];
                                currentlevel2[DIM-1]=j02[DIM-1]+temp-1;
                                break;
                            }
                        } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                        {
                            range2 = range2 +1;
                            if (DIM == 1)
                            {
                                currentlevel2[0]=currentlevel2[0]+1;
                            }
                            else
                            {
                                //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                currentlevel2[DIM-1]=j02[DIM-1]+range2;
                                currentlevel2[0]=j02[0];
                            }
                        }
                    } // end of "curentlevel++"
                } // end of while (range < maxrange)
                ++lambda_run2;
                currentlevel2 = j02;
                range2 = 0;
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; levelrange = " << levelrange << "; minlambdanum = " << minlambdanum2 << "; maxlambdanum = " << maxlambdanum2 << endl
            << "intersecting_wavelets " << time1 << "sec; intersecting_elements " << time2 << "sec" << endl;
    
#endif
    
    
    
    
    
    
    
    
    
    
    
#if 0
    unsigned int von(143),bis(143);
    cout << "- compute all intersecting generators"<<endl;
    cout << "for lambda.number() = "<<von<<" until number = "<<bis<<endl;
    cout << "with gen/wav on levels with a total_degree <= ||j0||+x" << endl;
    for (unsigned int b=0; b<TBasisArray.size(); ++b)
    {
            
    for (Index lambda = first_generator<Basis1D,2,Basis>(&basis);; ++lambda)
    {
        if (lambda.number() < von) continue;
        cout << "  * for lambda=" << lambda << ":" << endl;
        typedef std::list<Index> SupportList;
        SupportList nus;
        typedef TensorIndex<Basis1D,2,Basis>::level_type level_type;
        level_type j(basis.j0());
        intersecting_wavelets<Basis1D,DIM>(basis, lambda, j, true, nus);
        /*
        // generators:
        for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
            Basis::Support supp, supp_it;
            basis.support(*it, supp_it);
            intersect_supports(basis, lambda, *it, supp);
            
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
         * */
        //compute intersecting wavelets and increase j until j is too high
        level_type level = basis.j0();
        while (multi_degree(level) <= multi_degree(basis.j0())+1) // +2)
        {
            intersecting_wavelets<Basis1D,DIM>(basis, lambda, level, false, nus);
            for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) 
            {
                Basis::Support supp, supp_it;
                basis.support(*it, supp_it);
                assert (intersect_supports(basis, lambda, *it, supp) == true);
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
            // level++
            for (int i(DIM-1); i >= 0; i--)
            {
                if (i != 0) // increase sublevel
                {
                    if (level[i] != basis.j0()[i])
                    {
                        // increase left neighbor
                        level[i-1]=level[i-1]+1;
                        int temp = level[i]-basis.j0()[i];
                        level[i]=basis.j0()[i];
                        level[DIM-1]=basis.j0()[DIM-1]+temp-1;
                        break;
                    }
                } else // i == 0. "big loop" arrived at the last index. We have to increase the level! (will only happen 2 times in test_basis_support)
                {
                    if (DIM == 1)
                    {
                        level[0]=level[0]+1;
                    }
                    else
                    {
                        level[DIM-1]=basis.j0()[DIM-1]+level[0]-basis.j0()[0]+1;
                        level[0]=basis.j0()[0];
                    }
                    break; // unnoetig, da i==0 gilt.
                }
            }
        }   // end of while
        if (lambda.number() > bis-1) break;
        if (lambda == last_wavelet<Basis1D,DIM,Basis>(&basis, basis.j0())) break;
    }
#endif
    
#if 0
    Index lambda2(51,&basis);
    cout << "Testing intersect_singular_support with lambda = " << lambda2 << endl;

    //for (Index it = first_generator<Basis1D,DIM,Basis>(&basis);; ++lambda)
    for (Index it = first_generator<Basis1D,DIM,Basis>(&basis);it <= last_wavelet<Basis1D,2,Basis>(&basis,basis.j0()); ++it)
    {
        if (intersect_singular_support<Basis1D,DIM>(basis,lambda2,it))
            cout << "intersection with nu = " << it << endl;
    }
#endif
    
#if 0
    cout << "Comparing speed of old and new code of intersecting_wavelets" << endl;
    
    int range2(0), minlambdanum2(7995), maxlambdanum2(8000), count2(0), repetitions(5);
    level_type currentlevel2;
    SupportList nus2;
    Index lambda_run2;
    
    
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    tstart = clock();
    for (unsigned int b=0; b<TBasisArray.size(); ++b)
    {
        level_type j02(TBasisArray[b]->j0());
        assert(maxlambdanum2 < TBasisArray[b]->degrees_of_freedom());
        for (int rep(0); rep < repetitions; rep++)
        {
            lambda_run2 = TBasisArray[b]->get_wavelet(minlambdanum2);
            while (lambda_run2.number() <= maxlambdanum2)
            {
                currentlevel2 = j02;
                range2 = 0;
                count2 = 0;
                intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run2, currentlevel2, true, nus2);
                count2 += nus2.size();
                nus2.clear();
                while (range2 <= levelrange)
                {
                    intersecting_wavelets<Basis1D,DIM>(*TBasisArray[b], lambda_run2, currentlevel2, false, nus2);
                    count2 += nus2.size();
                    nus2.clear();

                    //intersecting_elements<Basis1D,DIM>(basis, lambda_run, currentlevel, nus);
                    //nus.clear();

                    // increase the level j,  meaning: "currentlevel++"
                    //cout << currentlevel << endl;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        if (i != 0)
                        {
                            if (currentlevel2[i] != j02[i])
                            {
                                // increase left neighbor
                                currentlevel2[i-1]=currentlevel2[i-1]+1;
                                int temp = currentlevel2[i]-j02[i];
                                currentlevel2[i]=j02[i];
                                currentlevel2[DIM-1]=j02[DIM-1]+temp-1;
                                break;
                            }
                        } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                        {
                            range2 = range2 +1;
                            if (DIM == 1)
                            {
                                currentlevel2[0]=currentlevel2[0]+1;
                            }
                            else
                            {
                                //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                currentlevel2[DIM-1]=j02[DIM-1]+range2;
                                currentlevel2[0]=j02[0];
                            }
                        }
                    } // end of "curentlevel++"
                } // end of while (range < maxrange)
                ++lambda_run2;
                
            }
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; levelrange = " << levelrange << "; minlambdanum = " << minlambdanum2 << "; maxlambdanum = " << maxlambdanum2 << endl
            << "intersecting_wavelets " << time1 << "sec; intersecting_elements " << time2 << "sec" << endl;
    
#endif
    cout << endl;
}

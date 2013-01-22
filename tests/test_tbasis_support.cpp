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
  //  typedef DSBasis<3,5> Basis1D;
  typedef PBasis<2,2> Basis1D;
//   typedef JLBasis Basis1D;

    const unsigned int DIM = 2;
    typedef TensorBasis<Basis1D,DIM> Basis;
    typedef Basis::Index Index;

    FixedArray1D<bool,2*DIM> bc;
    for (unsigned int i=0; i< 2*DIM; ++i)
    {
        bc[i] = false;
    }
    //bc[0] = true;
    bc[1] = true;
    //bc[2*DIM-2] = true;
    //bc[2*DIM-1] = true;
    
    Basis basis(bc);

#if 1
    cout << "- testing calculation of supports:" << endl;
    Basis::Support supp;
    int maxnumber(203);
    for (Index lambda(first_generator<Basis1D,DIM,Basis>(&basis));; ++lambda) 
    {
        support<Basis1D,DIM>(basis, lambda, supp);
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
        if (lambda.number() == maxnumber) break;
    }
    

#endif

    typedef TensorIndex<Basis1D,DIM,Basis>::level_type level_type;
    typedef std::list<Index> SupportList;
    
#if 0
    Index lambda(202,&basis);
    
    level_type lambda_j(basis.j0());
    
    SupportList intwav1, intwav2, intwav3, intel1, intel2;
    
    Basis::Support supp_lambda;
    basis.support(lambda, supp_lambda);
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
        
    
    Index lambda1(81,&basis);
    Basis::Support supp_lambda1;
    basis.support(lambda1, supp_lambda1);
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
    temp_bool = intersect_supports(basis, lambda, lambda1, temp_supp);
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
        
        
    intersecting_wavelets<Basis1D,DIM>(basis, lambda, lambda_j, true, intwav1);
    cout << "Number of intersectiong generators is " << intwav1.size() << endl;
    for (SupportList::const_iterator it(intwav1.begin()); it != intwav1.end(); ++it)
    {
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
    cout << "intersecting wavelets on level j0 = " << basis.j0() << endl;
    intersecting_wavelets<Basis1D,DIM>(basis, lambda, lambda_j, false, intwav2);
    cout << "Number of intersectiong wavelets is " << intwav2.size() << endl;
    for (SupportList::const_iterator it(intwav2.begin()); it != intwav2.end(); ++it)
    {
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
    
    lambda_j[DIM-1]++;
    
    cout << "intersecting wavelets on level " << lambda_j << endl;
    intersecting_wavelets<Basis1D,DIM>(basis, lambda, lambda_j, false, intwav3);
    cout << "Number of intersectiong wavelets is " << intwav3.size() << endl;
    for (SupportList::const_iterator it(intwav3.begin()); it != intwav3.end(); ++it)
    {
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
    
    // throws assertion
    // intersecting_wavelets<Basis1D,DIM>(basis, lambda, lambda_j, true, nus);
    
    lambda_j = basis.j0();
    
    cout << " Testing intersecting_elements on level j0 = " << basis.j0() << lambda_j <<endl;
    intersecting_elements<Basis1D,DIM>(basis, lambda, lambda_j, intel1);
    cout << "Number of intersectiong elements is " << intel1.size() << endl;
    for (SupportList::const_iterator it(intel1.begin()); it != intel1.end(); ++it)
    {
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
    
    lambda_j[DIM-1]++;
    cout << " Testing intersecting_elements on level "<< lambda_j <<endl;
    intersecting_elements<Basis1D,DIM>(basis, lambda, lambda_j, intel2);
    cout << "Number of intersectiong elements is " << intel2.size() << endl;
    for (SupportList::const_iterator it(intel2.begin()); it != intel2.end(); ++it)
    {
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
        
#endif

#if 1
    // intersecting_wavelets and intersecting_elements produce the same sets
    // (on level j0: of course intersecting _wavs and _gens together)
    // in this test not the sets are compared, but the number of elements in the sets
    // proof:
            
    cout << "Testing intersect_wavelets and intersect_elements" << endl;
    //Index lambda_run(basis.first_generator());
    level_type j0(basis.j0());
    level_type currentlevel(j0);
    int range(0), maxrange(2), minlambdanum(1456), maxlambdanum(1466), count(0); // number of level steps we have climbed so far.
    basis.set_jmax(multi_degree(basis.j0())+maxrange);
    assert(maxlambdanum < basis.degrees_of_freedom());
    Index lambda_run(basis.get_wavelet(minlambdanum));
    SupportList nus;
    bool isordered = true;
    SupportList::const_iterator it2;

    cout << "levelrange = " << maxrange << "; minlambdanum = " << minlambdanum << "; maxlambdanum = " << maxlambdanum << endl;
    
    while (lambda_run.number() <= maxlambdanum)
    {
        // compute intersecting_wavelets and intersecting_elements for all levels up to \|j0\|+maxrange
        cout << "lambda_run = " << lambda_run << "; currentlevel = ";
        cout.flush();
        while (range <= maxrange)
        {
            //cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
            cout  << currentlevel << "; ";
            cout.flush();
            count = 0;
            if (range == 0)
            {
                //cout << nus.size() << endl;
                intersecting_wavelets<Basis1D,DIM>(basis, lambda_run, currentlevel, true, nus);
                //cout << nus.size() << endl;
          /*      
                //cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
                
                //cout << "nus =" << endl;
                //for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it)
                //{
                //    cout << " N = " << (*it).number() << "; " << (*it) << endl;
                //}
                
                
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
                }
        */        
                count += nus.size();
                //cout << "#nus = " << nus.size() << " ";
                nus.clear();
            }
            intersecting_wavelets<Basis1D,DIM>(basis, lambda_run, currentlevel, false, nus);
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
                cout << "intersecting_wavelets is NOT ordered" << endl;
                cout << "N = " << lambda_run.number() << "; lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
            }
    */            
            count += nus.size();
            //cout << "#nus = " << nus.size() << " ";
            nus.clear();
            intersecting_elements<Basis1D,DIM>(basis, lambda_run, currentlevel, nus);
            
        ///*
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
        //*/              
            count -= nus.size();
            //cout << "#nus = " << nus.size() << " ";
            nus.clear();
            if (count != 0)
            {
                cout << "Problem detected" <<  endl;
                cout <<  "intersecting_wavelets and intersecting_elements produce different sets" << endl;
                cout << "lambda_run = " << lambda_run << "; currentlevel = " << currentlevel << endl;
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
#endif
    
    
    

#if 0
    // methods are equally fast! intersect_elements is marginally slower
    cout << "Comparing speed of intersect_wavelets and intersect_elements" << endl;
    level_type j02(basis.j0());
    level_type currentlevel2(j02);
    SupportList nus2;
    Index lambda_run2(basis.first_generator());
    int range2(0), maxrange2(4), minlambdanum2(14567), maxlambdanum2(14567), count2(0), repetitions(5);
    
    basis.set_jmax(multi_degree(basis.j0())+maxrange2 );
    assert(maxlambdanum < basis.degrees_of_freedom());
    
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    tstart = clock();
    for (int rep(0); rep < repetitions; rep++)
    {
        lambda_run2 = basis.get_wavelet(minlambdanum);
        while (lambda_run2.number() <= maxlambdanum2)
        {
            count2 = 0;
            intersecting_wavelets<Basis1D,DIM>(basis, lambda_run2, currentlevel2, true, nus2);
            count2 += nus2.size();
            nus2.clear();
            while (range2 <= maxrange2)
            {
                intersecting_wavelets<Basis1D,DIM>(basis, lambda_run2, currentlevel2, false, nus2);
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
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    for (int rep(0); rep < repetitions; rep++)
    {
        lambda_run2 = basis.get_wavelet(minlambdanum);
        while (lambda_run2.number() <= maxlambdanum2)
        {
            count2 = 0;
            while (range2 <= maxrange2)
            {
                intersecting_elements<Basis1D,DIM>(basis, lambda_run2, currentlevel2, nus2);
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
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; levelrange = " << maxrange2 << "; minlambdanum = " << minlambdanum2 << "; maxlambdanum = " << maxlambdanum2 << endl
            << "intersecting_wavelets " << time1 << "sec; intersecting_elements " << time2 << "sec" << endl;
    
#endif
    
    
    
    
    
    
    
    
    
    
    
#if 0
    unsigned int von(143),bis(143);
    cout << "- compute all intersecting generators"<<endl;
    cout << "for lambda.number() = "<<von<<" until number = "<<bis<<endl;
    cout << "with gen/wav on levels with a total_degree <= ||j0||+x" << endl;
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
            for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
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
    cout << endl;
}

#include <iostream>

#include <interval/p_basis.h>
#include <interval/p_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
    cout << "Testing support calculation of [P] generators and wavelets..." << endl;
#define d 3
#define dT 3
    //const int d = 3;
    //const int dT = 3;

    typedef PBasis<d,dT> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;

//   Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
//    Basis basis(1, 2); // complementary b.c.'s

        
    Array1D<Basis*> BasisArray;
    Basis basis(0,0); // 1st order complementary b.c.'s at x=0 and x=1
    Basis basis01(0,1);
    Basis basis10(1,0);
    Basis basis11(1,1);
#if d > 2
    Basis basis02(0,2);
    Basis basis20(2,0);
    Basis basis12(1,2);
    Basis basis21(2,1);
    Basis basis22(2,2);
    BasisArray.resize(9);
    BasisArray[0] = &basis;
    BasisArray[1] = &basis01;
    BasisArray[2] = &basis10;
    BasisArray[3] = &basis11;
    BasisArray[4] = &basis02;
    BasisArray[5] = &basis20;
    BasisArray[6] = &basis12;
    BasisArray[7] = &basis21;
    BasisArray[8] = &basis22;
#else // d=2
    BasisArray.resize(4);
    BasisArray[0] = &basis;
    BasisArray[1] = &basis01;
    BasisArray[2] = &basis10;
    BasisArray[3] = &basis11;
#endif    
    
    const unsigned int maxlevelrange(7);
    
    for (unsigned int b=0; b<BasisArray.size(); ++b)
    {
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+maxlevelrange);
    }
    
#if 0
    cout << "validate the supports computed by support()" << endl;
    int k1, k2, l1,l2, j,e,k, steps(100);
    
    Index temp_ind;
    double x1, x2, val1, val2, val3, val4, h(1.0/steps), temp_d;
    for (unsigned int b=0; b<BasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        for (unsigned int i=0; i< BasisArray[b]->degrees_of_freedom(); ++i)
        //for (unsigned int i=0; i< 20; ++i)
        {
            temp_ind = BasisArray[b]->get_wavelet(i);
            j = temp_ind.j();
            e = temp_ind.e();
            k = temp_ind.k();
            BasisArray[b]->support(temp_ind,k1,k2);
            //cout << "N = " << i << "; temp_ind = " << temp_ind << "; k1 k2 = " << k1 << " " << k2 << endl;
            /*
            BasisArray[b]->support(j,e,k,l1,l2);
            if (!((k1 == l1) && (k2 == l2)))
            {
                cout << "Support computation strange: lambda = " << temp_ind << endl;
                BasisArray[b]->support(temp_ind,k1,k2);
                BasisArray[b]->support(j,e,k,l1,l2);
            }
            assert ((k1 == l1) && (k2 == l2));
             */
            // support = 2^{-(j+e)} [k1,k2]
            // check meaning of k1
            x1 = ldexp (k1, -(j+e));
            
            /*
            int s0(BasisArray[b]->get_s0());
            int s1(BasisArray[b]->get_s1());
            int sT0(BasisArray[b]->get_sT0());
            int sT1(BasisArray[b]->get_sT1());
            
            int r0(temp_ind.basis()->get_s0());
            int r1(temp_ind.basis()->get_s1());
            int rT0(temp_ind.basis()->get_sT0());
            int rT1(temp_ind.basis()->get_sT1());
            
            cout << "alive!" << endl;
            
            int deltasize1(BasisArray[b]->Deltasize(j));
            int deltasize2(temp_ind.basis()->Deltasize(j));
            int deltasize3(BasisArray[b]->Deltasize(j+1));
            int deltasize4(temp_ind.basis()->Deltasize(j+1));
            
            int temp_int2(BasisArray[b]->get_Mj0().row_dimension() );
            
            cout << "alive!" << endl;
            //const SparseMatrix<double>& temp_Mj0 = BasisArray[b]->get_Mj0();
            //bool temp_bool = temp_Mj0.empty();
            
            InfiniteVector<double,int> gcoeffs;
            BasisArray[b]->reconstruct_1(j, e, k, j+1, gcoeffs);
      
            InfiniteVector<double,Index> gcoeffs2;
            BasisArray[b]->reconstruct_1(temp_ind,j+1, gcoeffs2);
            */
            val1 = BasisArray[b]->evaluate(0, temp_ind, x1);
            x2 = ldexp (k2, -(j+e));
            val2 = BasisArray[b]->evaluate(0, temp_ind, x2);
            if (k < (d+dT)/2-1) // left boundary wavelet
            {
                if (val1 == 0)
                {
                    // this is not a bug:
                    //cout << "N = " << i << "; lam = " << temp_ind << " is left boundary wavelet and strange at left boundary: k_1 = " << k1 << "; val1 = " << val1 << endl;
                }
                if (k1 != 0)
                {
                    // bug! (not every time ... does not check for generators)
                    cout << "N = " << i << "; lam = " << temp_ind << " is no boundary wavelet but should be!" << endl;
                }
                if (val2 != 0)
                {
                    if (x2 != 1)
                    {
                        // bug!
                        cout << "N = " << i << "; lam = " << temp_ind << " is left boundary wavelet and nonzero at right boundary!" << endl;
                    }
                    else
                    {
                        // this is not a bug, but indicates a very big support
                        cout << "N = " << i << "; lam = " << temp_ind << " is left boundary wavelet with supp [0,1]! strange!" << endl;
                    }
                }
            }
            else
            {
                if ((1<<j)-k <= (d+dT)/2-1) // right boundary wavelet
                {
                    if  (val2 == 0)
                    {
                        // this is not a bug:
                        //cout << "N = " << i << "; lam = " << temp_ind << " is right boundary wavelet and strange at right boundary: k_2 = " << k2 << "; val2 = " << val2 << endl;
                    }
                    if (k2 != (1<<(j+e)))
                    {
                        // bug! (not every time ... does not check for generators)
                        cout << "N = " << i << "; lam = " << temp_ind << " is no boundary wavelet but should be!" << endl;
                    }
                    if (val1 != 0)
                    {
                        if (x1 != 0)
                        {
                            // bug!
                            cout << "N = " << i << "; lam = " << temp_ind << " is right boundary wavelet and nonzero at left boundary!" << endl;
                        }
                        else
                        {
                            // this is not a bug, but indicates a very big support
                            cout << "N = " << i << "; lam = " << temp_ind << " is right boundary wavelet with supp [0,1]! strange!" << endl;
                        }
                    }
                }
                else
                {
                    if ((val1 != 0) || (val2 != 0))
                    {
                        // bug!
                        cout << "N = " << i << "; lam = " << temp_ind << " is inner wavelet and nonzero at the support boundaries!" << endl;
                    }
                }

            }
            val1=0;
            val2=0;
            val3=0;
            val4=0;
            
            for (unsigned int j=1; j<(steps/2); ++j)
            {
                temp_d = ldexp(h*j, (-j+e) );
                // check left of left support boundary (val1 == 0)
                val1 += abs(BasisArray[b]->evaluate(0, temp_ind, x1 - temp_d));
                // check right of right support boundary (val2 == 0)
                val2 += abs(BasisArray[b]->evaluate(0, temp_ind, x2 + temp_d));
                // check right of left support boundary (val3 != 0)
                val3 += abs(BasisArray[b]->evaluate(0, temp_ind, x1 + temp_d));
                // check left of right support boundary (val4 != 0)
                val4 += abs(BasisArray[b]->evaluate(0, temp_ind, x2 - temp_d));
            }
            if (val1 != 0)
            {
                cout << "N = " << i << "; lam = " << temp_ind << " is nonzero left of k1! => k1 too big (underestimate!)" << endl;
            }
            if (val2 != 0)
            {
                cout << "N = " << i << "; lam = " << temp_ind << " is nonzero right of k2! => k2 too small (underestimate!)" << endl;
            }
            if (val3 == 0)
            {
                cout << "N = " << i << "; lam = " << temp_ind << " is zero right of k1! => k1 is too small (overestimate)" << endl;
            }
            if (val4 == 0)
            {
                cout << "N = " << i << "; lam = " << temp_ind << " is zero left of k2! => k2 is too big (overestimate)" << endl;
            }
        }
    }
#endif
#if 0
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- computing the supports of all generators and wavelets on level j=" << level << ":" << endl;
    
    Index lambda(first_generator(&basis, level));
    for (;; ++lambda) {
      int k1, k2;
      support(basis, lambda, k1, k2);
      cout << "  lambda=" << lambda << ", supp(psi_lambda)=2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << k1
	   << ","
	   << k2
	   << "]"
	   << endl;
      
      if (lambda == last_wavelet(&basis, level)) break;
    }
  }
#endif
  
#if 0
  cout << "- calculating some support intersections:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2
	   << "]"
	   << endl;

      cout << "support intersection with first generator on level j0: ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif

#if 0
  cout << "- compute all intersecting wavelets:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      cout << "  * for lambda=" << lambda << ":" << endl;
      typedef std::list<std::pair<Index, Basis::Support> > SupportList;
      SupportList nus;
      intersecting_wavelets(basis, lambda, basis.j0(), true, nus);
      for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	cout << "    nu=" << it->first 
	     << " with support intersection "
	     << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
      }
      for (int level = basis.j0(); level <= basis.j0()+2; level++) {
	intersecting_wavelets(basis, lambda, level, false, nus);
	for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	  cout << "    nu=" << it->first 
	       << " with support intersection "
	       << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
	}
      }
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif  

#if 0
  cout << "- checking intersection of singular supports:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0()+2);; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has the support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2307681640 0 0

	   << "]"
	   << endl;
      
      Support supp_0;
      support(basis, first_generator(&basis, basis.j0()), supp_0.k1, supp_0.k2);
      cout << "* first generator on level j0 has the support 2^{-"
	   << basis.j0()
	   << "}["
	   << supp_0.k1
	   << ","
	   << supp_0.k2
	   << "]"
	   << endl;

      cout << "* support intersection with first generator on level j0         : ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      cout << "* singular support intersection with first generator on level j0: ";
      inter = intersect_singular_support(basis, lambda, first_generator(&basis, basis.j0()), supp.j, supp.k1, supp.k2);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      if (lambda == last_wavelet(&basis, basis.j0()+2)) break;
    }
#endif

#if 0
  
    cout << "Check code modification for intersect_singular_support(bas, jlam,elam,klam, jmu,emu,kmu, j,k1,k2)" << endl;
    
    Index temp_ind1, temp_ind2;
    int j1,j2, k11,k12,k21,k22;
    bool temp_b1, temp_b2;
    
    cout << "validate new code" << endl;
    for (unsigned int b=0; b< BasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        
        for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
        {
            temp_ind1 = BasisArray[b]->get_wavelet(i);
            BasisArray[b]->support(temp_ind1, k11,k21);
            BasisArray[b]->support(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), k12,k22);
            
            if ((k11 != k12) || (k21 != k22) )
            {
                cout << "i = " << i << "; temp_ind1 = " << temp_ind1 << "; k11 k21 = " << k11 << " " << k21 << "; k12 k22 = " << k12 << " " << k22 << endl;
                BasisArray[b]->support(temp_ind1, k11,k21);
                BasisArray[b]->support(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), k12,k22);
            }
            temp_b1 = intersect_singular_support(*BasisArray[b], temp_ind1, temp_ind1, j1,k11, k21);
            if (!temp_b1)
            {
                cout << "i = " << i << "; temp_ind1 = " << temp_ind1 << "; k11 k21 = " << k11 << " " << k21
                        << "; self_intersecting_support: temp_b1 = " << temp_b1 << "; j1 k11 k21 = " << j1 << " " << k11 << " " << k21 << endl;
                temp_b1 = intersect_singular_support(*BasisArray[b], temp_ind1, temp_ind1, j1,k11, k21);
            }
        }
        
        for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
        {
            temp_ind1 = BasisArray[b]->get_wavelet(i);
            //BasisArray[b]->support(temp_ind1, k11,k12);
            //temp_b1 = intersect_singular_support(*BasisArray[b],temp_ind1,  temp_ind1.j()+temp_ind1.e(), k11, k12, j1, k12, k22);
            
            for (unsigned int j=0; j< BasisArray[b]->degrees_of_freedom(); ++j)
            {
                k11 = BasisArray[b]->degrees_of_freedom();
                temp_ind2 = BasisArray[b]->get_wavelet(j);
                BasisArray[b]->support(temp_ind1, k11,k21);
                BasisArray[b]->support(temp_ind2, k12,k22);
                support(*BasisArray[b], temp_ind2.j(), temp_ind2.e(), temp_ind2.k(), k12, k22);
                BasisArray[b]->support(temp_ind2.j(), temp_ind2.e(), temp_ind2.k(), k12, k22);
                
                //cout << "i j = " << i << " " << j << "; temp_ind1 = " << temp_ind1 << "; temp_ind2 = " << temp_ind2 << endl;
                //cout << "k11 k21 = " << k11 << " " << k21 << "; k12 k22 = " << k12 << " " << k22 << endl;
                
                temp_b1 = intersect_singular_support(*BasisArray[b], temp_ind1, temp_ind2, j1,k11, k21);
                temp_b2 = intersect_singular_support(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), temp_ind2.j(), temp_ind2.e(), temp_ind2.k(), j2,k12,k22);
                if (temp_b1 || temp_b2)
                {
                    if ((((temp_b1 != temp_b2) || (j1 != j2))
                        || (k11 != k12))
                        || (k21 != k22))
                    {
                        cout << "i = " << i << "; lam = " << temp_ind1 
                                << "; j = " << j << "; mu = " << temp_ind2
                                << "; temp_b1 = " << temp_b1 << "; temp_b2 = " << temp_b2
                                << "; j1 k11 k21 = " << j1 << " " << k11 << " " << k21 
                                << "; j2 k12 k22 = " << j1 << " " << k12 << " " << k22 << endl;
                    }
                }
            }
        }
    }

    int repetitions(0);
    clock_t tstart, tend;
    double time1(0), time2(0);
    cout << "Comparing speed of old and new code for intersect_singular_support" << endl;
    cout << "running old code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            //cout << "Basis b = " << b << endl;
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int j = 0; j< BasisArray[b]->degrees_of_freedom(); ++j)
                {
                    temp_ind2 = BasisArray[b]->get_wavelet(j);
                    temp_b1 = intersect_singular_support(*BasisArray[b], temp_ind1, temp_ind2, j1,k11, k21);
                    //temp_b2 = intersect_singular_support(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), temp_ind2.j(), temp_ind2.e(), temp_ind2.k(), j2,k12,k22);
                }
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "running new code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            //cout << "Basis b = " << b << endl;
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int j = 0; j< BasisArray[b]->degrees_of_freedom(); ++j)
                {
                    temp_ind2 = BasisArray[b]->get_wavelet(j);
                    //temp_b1 = intersect_singular_support(*BasisArray[b], temp_ind1, temp_ind2, j1,k11, k21);
                    temp_b2 = intersect_singular_support(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), temp_ind2.j(), temp_ind2.e(), temp_ind2.k(), j2,k12,k22);
                }
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; old code " << time1 << "sec; new code " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
    
    
    
#endif
    
#if 0
    cout << "Check new code for get_intersecting_wavelets_on_level (bas, jlam,elam,klam, j, gen?, k1,k2)" << endl;
    cout << "It is only checked whether the new code computes smaller or equal supports than the old code (not the correctness)" << endl;
    Index temp_ind1;
    int k11,k12,k21,k22;
    bool gen;
    
    cout << "validate new code" << endl;
    for (unsigned int b=0; b< BasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        
        for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
        {
            temp_ind1 = BasisArray[b]->get_wavelet(i);
            for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
            {
                gen = true;
                get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                
                // d=2,b=3 => bass11: here one would expect the same output for old and new code:
                // d=3,b=8 => bass22: here one would expect the same output for old and new code:
                // in any other case one would expect less or equal intersection with the new code (since there is no more overestimation)
                if (( (((b == 3) && (d==2)  ) && ((k11 != k12) || (k21 != k22)))
                    || ((k11 > k12) || (k21 < k22)) )
                    || (((b == 8) && (d==3)  ) && ((k11 != k12) || (k21 != k22))) )
                {
                    cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; k11 k21 k12 k22 = " << k11 << " " << k21 << " " << k12 << " " << k22 << endl;
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                    abort();
                }
                
                gen = false;
                get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                // d=2,b=3 => bass11: here one would expect the same output for old and new code:
                // d=3,b=8 => bass22: here one would expect the same output for old and new code:
                // in any other case one would expect less or equal intersection with the new code (since there is no more overestimation)
                
                if (((b == 3) && (d==2)  ) && ((k11 != k12) || (k21 != k22)))
                {
                    cout << "A" << endl;
                }
                if ((k11 > k12) || (k21 < k22))
                
                {
                    cout << "B" << endl;;
                }
                if (((b == 8) && (d==3)  ) && ((k11 != k12) || (k21 != k22)))
                {
                    cout << "C" << endl;
                }
                if (( (((b == 3) && (d==2)  ) && ((k11 != k12) || (k21 != k22)))
                    || ((k11 > k12) || (k21 < k22)) )
                    || (((b == 8) && (d==3)  ) && ((k11 != k12) || (k21 != k22))) )
                {
                    cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; wavelets; level = " << level << "; k11 k21 k12 k22 = " << k11 << " " << k21 << " " << k12 << " " << k22 << endl;
                    //Index temp_ind2 = BasisArray[b]->last_wavelet(level);
                    //cout << "last wavelet = " << temp_ind2 << "; .number() = " << temp_ind2.number() << "; .Nablamax(" << level << ") = " << BasisArray[b]->Nablamax(level) << endl;
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                    abort();
                }
            }
        }
    }

    int repetitions(200);
    clock_t tstart, tend;
    double time1(0), time2(0);
    cout << "Comparing speed of old and new code for intersect_singular_support" << endl;
    cout << "running old code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
                {
                    gen = true;
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    //get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                    gen = false;
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    //get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                }
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "running new code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
                {
                    gen = true;
                    //get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                    gen = false;
                    //get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1, level, gen, k11,k21);
                    get_intersecting_wavelets_on_level(*BasisArray[b], temp_ind1.j(), temp_ind1.e(), temp_ind1.k(), level, gen, k12,k22);
                }
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; old code " << time1 << "sec; new code " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
#endif

    
#if 0
    cout << "Test a suspicious wavelet" << endl;
    for (unsigned int b=0; b< BasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        Index temp_ind1(4,0,1, BasisArray[b]);
        Index temp_ind2(4,1,6, BasisArray[b]);
        int k11, k21, k12, k22;
        BasisArray[b]->support(temp_ind1,k11,k21);
        BasisArray[b]->support(temp_ind2,k12,k22);
        cout << "temp_ind1 = " << temp_ind1 << "; k11 k21 = " << k11 << " " << k21 << endl;
        cout << "temp_ind2 = " << temp_ind2 << "; k12 k22 = " << k12 << " " << k22 << endl;
    }
#endif
    
    
#if 0
    cout << "Check new code for intersecting_wavelets (both versions)" << endl;
    list<Index> intersecting1, intersecting2;
    
    list<pair<Index, Support> > intersectingpairs1, intersectingpairs2;
    
    Index temp_ind1;
    bool gen;
    
    cout << "validate new code" << endl;
    for (unsigned int b=0; b< BasisArray.size(); ++b)
    {
        cout << "Basis b = " << b << endl;
        
        for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
        {
            temp_ind1 = BasisArray[b]->get_wavelet(i);
            for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
            {
                gen = true;
                intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                
                intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                // d=2,b=3 => bass11: here one would expect the same output for old and new code:
                // d=3,b=8 => bass22: here one would expect the same output for old and new code:
                // in any other case one would expect less or equal intersection with the new code (since there is no more overestimation)
                //cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; k11 k21 k12 k22 = " << k11 << " " << k21 << " " << k12 << " " << k22 << endl;
                if (((b == 3) && (d==2)) || ((b == 8) && (d==3)))
                {
                    assert(intersecting1.size() == intersecting2.size());
                    for (list<Index>::const_iterator it1(intersecting1.begin()), itend(intersecting1.end()), it2(intersecting2.begin()); it1 != itend; ++it1)
                    {
                        assert(*it1 == *it2);
                        ++it2;
                    }
                    
                    assert(intersecting1.size() == intersectingpairs1.size());
                    assert(intersectingpairs1.size() == intersectingpairs2.size());
                    for (list<pair<Index,Support> >::const_iterator it1(intersectingpairs1.begin()), itend(intersectingpairs1.end()), it2(intersectingpairs2.begin()); it1 != itend; ++it1)
                    {
                        assert ((*it1).first == (*it2).first);
                        assert ((*it1).second.j == (*it2).second.j);
                        assert ((*it1).second.k1 == (*it2).second.k1);
                        assert ((*it1).second.k2 == (*it2).second.k2);
                        ++it2;
                    }
                }
                else
                {
                    assert(intersecting1.size() >= intersecting2.size());
                    bool temp_b(true);
                    for (list<Index>::const_iterator it1(intersecting1.begin()), itend(intersecting1.end()), it2(intersecting2.begin()); it1 != itend; ++it1)
                    {
                        //cout << "conflicting value: *it1 = " << (*it1) << "; *it2 = " << (*it2) << endl;
                        if ( (*it1 != *it2) && (temp_b || (it2 == intersecting2.end())))
                            continue;
                        else
                        {
                            temp_b = false;

                            if (*it1 != *it2)
                            {
                                cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; intersecting1.size() = " << intersecting1.size() << "; intersecting2.size() = " << intersecting2.size() << endl;
                                cout << "conflicting value: *it1 = " << (*it1) << "; *it2 = " << (*it2) << endl;
                            }
                            if (*it1 != *it2)
                            {
                                cout << "intersecting1 = " << endl;
                                for (list<Index>::const_iterator itA(intersecting1.begin()), itendA(intersecting1.end()); itA != itendA; ++itA)
                                {
                                    cout << (*itA) << " ";
                                }
                                cout << endl;
                                cout << "intersecting2 = " << endl;
                                for (list<Index>::const_iterator itB(intersecting2.begin()), itendB(intersecting2.end()); itB != itendB; ++itB)
                                {
                                    cout << (*itB) << " ";
                                }
                                cout << endl;
                            }
                            assert (*it1 == *it2);
                            ++it2;
                        }
                    }
                    
                    assert(intersecting1.size() == intersectingpairs1.size());
                    assert(intersecting2.size() == intersectingpairs2.size());
                    assert(intersectingpairs1.size() >= intersectingpairs2.size());
                    temp_b = true;
                    bool temp_b2;
                    for (list<pair<Index,Support> >::const_iterator it1(intersectingpairs1.begin()), itend(intersectingpairs1.end()), it2(intersectingpairs2.begin()); it1 != itend; ++it1)
                    {
                        temp_b2 = (*it1).first == (*it2).first;
                        temp_b2 = temp_b2 && (*it1).second.j == (*it2).second.j;
                        temp_b2 = temp_b2 && (*it1).second.k1 == (*it2).second.k1;
                        temp_b2 = temp_b2 && (*it1).second.k2 == (*it2).second.k2;
                        
                        if ( (!temp_b2) && (temp_b || (it2 == intersectingpairs2.end())))
                            continue;
                        else
                        {
                            temp_b = false;

                            if (!temp_b2)
                            {
                                cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; intersectingpairs1.size() = " << intersectingpairs1.size() << "; intersectingpairs2.size() = " << intersectingpairs2.size() << endl;
                                cout << "conflicting value: *it1 = " << (*it1).first << "; " << (*it1).second.j << " " << (*it1).second.k1 << " " << (*it1).second.k2 << "; ";
                                cout << "*it2 = " << (*it2).first << "; [" << (*it2).second.j << "; " << (*it2).second.k1 << "; " << (*it2).second.k2 << "]; ";
                                        
                            }
                            if (!temp_b2)
                            {
                                cout << "intersecting1 = " << endl;
                                for (list<pair<Index,Support> >::const_iterator itA(intersectingpairs1.begin()), itendA(intersectingpairs1.end()); itA != itendA; ++itA)
                                {
                                    cout << (*itA).first << ";[" << (*itA).second.j << "; " << (*itA).second.k1 << "; " << (*itA).second.k2 << "];  ";
                                }
                                cout << endl;
                                cout << "intersecting2 = " << endl;
                                for (list<pair<Index, Support> >::const_iterator itB(intersectingpairs2.begin()), itendB(intersectingpairs2.end()); itB != itendB; ++itB)
                                {
                                    cout << (*itB).first << ";[" << (*itB).second.j << "; " << (*itB).second.k1 << "; " << (*itB).second.k2 << "];  ";
                                }
                                cout << endl;
                            }
                            assert (temp_b2);
                            ++it2;
                        }
                    }
                }
                
                gen = false;
                intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                
                intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                
                // d=2,b=3 => bass11: here one would expect the same output for old and new code:
                // d=3,b=8 => bass22: here one would expect the same output for old and new code:
                // in any other case one would expect less or equal intersection with the new code (since there is no more overestimation)
                //cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; k11 k21 k12 k22 = " << k11 << " " << k21 << " " << k12 << " " << k22 << endl;
                if (((b == 3) && (d==2)) || ((b == 8) && (d==3)))
                {
                    assert(intersecting1.size() == intersecting2.size());
                    for (list<Index>::const_iterator it1(intersecting1.begin()), itend(intersecting1.end()), it2(intersecting2.begin()); it1 != itend; ++it1)
                    {
                        assert(*it1 == *it2);
                        ++it2;
                    }
                    
                    assert(intersecting1.size() == intersectingpairs1.size());
                    assert(intersectingpairs1.size() == intersectingpairs2.size());
                    for (list<pair<Index,Support> >::const_iterator it1(intersectingpairs1.begin()), itend(intersectingpairs1.end()), it2(intersectingpairs2.begin()); it1 != itend; ++it1)
                    {
                        assert ((*it1).first == (*it2).first);
                        assert ((*it1).second.j == (*it2).second.j);
                        assert ((*it1).second.k1 == (*it2).second.k1);
                        assert ((*it1).second.k2 == (*it2).second.k2);
                        ++it2;
                    }
                }
                else
                {
                    assert(intersecting1.size() >= intersecting2.size());
                    bool temp_b(true);
                    for (list<Index>::const_iterator it1(intersecting1.begin()), itend(intersecting1.end()), it2(intersecting2.begin()); it1 != itend; ++it1)
                    {
                        //cout << "conflicting value: *it1 = " << (*it1) << "; *it2 = " << (*it2) << endl;
                        if ( (*it1 != *it2) && (temp_b || (it2 == intersecting2.end())))
                            continue;
                        else
                        {
                            temp_b = false;

                            if (*it1 != *it2)
                            {
                                cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; intersecting1.size() = " << intersecting1.size() << "; intersecting2.size() = " << intersecting2.size() << endl;
                                cout << "conflicting value: *it1 = " << (*it1) << "; *it2 = " << (*it2) << endl;
                            }
                            if (*it1 != *it2)
                            {
                                cout << "intersecting1 = " << endl;
                                for (list<Index>::const_iterator itA(intersecting1.begin()), itendA(intersecting1.end()); itA != itendA; ++itA)
                                {
                                    cout << (*itA) << " ";
                                }
                                cout << endl;
                                cout << "intersecting2 = " << endl;
                                for (list<Index>::const_iterator itB(intersecting2.begin()), itendB(intersecting2.end()); itB != itendB; ++itB)
                                {
                                    cout << (*itB) << " ";
                                }
                                cout << endl;
                            }
                            assert (*it1 == *it2);
                            ++it2;
                        }
                    }
                    
                    assert(intersecting1.size() == intersectingpairs1.size());
                    assert(intersecting2.size() == intersectingpairs2.size());
                    assert(intersectingpairs1.size() >= intersectingpairs2.size());
                    temp_b = true;
                    bool temp_b2;
                    for (list<pair<Index,Support> >::const_iterator it1(intersectingpairs1.begin()), itend(intersectingpairs1.end()), it2(intersectingpairs2.begin()); it1 != itend; ++it1)
                    {
                        temp_b2 = (*it1).first == (*it2).first;
                        temp_b2 = temp_b2 && (*it1).second.j == (*it2).second.j;
                        temp_b2 = temp_b2 && (*it1).second.k1 == (*it2).second.k1;
                        temp_b2 = temp_b2 && (*it1).second.k2 == (*it2).second.k2;
                        
                        if ( (!temp_b2) && (temp_b || (it2 == intersectingpairs2.end())))
                            continue;
                        else
                        {
                            temp_b = false;

                            if (!temp_b2)
                            {
                                cout << " i = " << i << "; temp_ind1 = " << temp_ind1 << "; generators; intersectingpairs1.size() = " << intersectingpairs1.size() << "; intersectingpairs2.size() = " << intersectingpairs2.size() << endl;
                                cout << "conflicting value: *it1 = " << (*it1).first << "; " << (*it1).second.j << " " << (*it1).second.k1 << " " << (*it1).second.k2 << "; ";
                                cout << "*it2 = " << (*it2).first << "; [" << (*it2).second.j << "; " << (*it2).second.k1 << "; " << (*it2).second.k2 << "]; ";
                                        
                            }
                            if (!temp_b2)
                            {
                                cout << "intersecting1 = " << endl;
                                for (list<pair<Index,Support> >::const_iterator itA(intersectingpairs1.begin()), itendA(intersectingpairs1.end()); itA != itendA; ++itA)
                                {
                                    cout << (*itA).first << ";[" << (*itA).second.j << "; " << (*itA).second.k1 << "; " << (*itA).second.k2 << "];  ";
                                }
                                cout << endl;
                                cout << "intersecting2 = " << endl;
                                for (list<pair<Index, Support> >::const_iterator itB(intersectingpairs2.begin()), itendB(intersectingpairs2.end()); itB != itendB; ++itB)
                                {
                                    cout << (*itB).first << ";[" << (*itB).second.j << "; " << (*itB).second.k1 << "; " << (*itB).second.k2 << "];  ";
                                }
                                cout << endl;
                            }
                            assert (temp_b2);
                            ++it2;
                        }
                    }
                }
            }
        }
    }

    int repetitions(10);
    clock_t tstart, tend;
    double time1(0), time2(0);
    cout << "Comparing speed of old and new code for intersect_singular_support" << endl;
    cout << "running old code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
                {
                    gen = true;
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                    intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                    gen = false;
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                    intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                }
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "running new code" << endl;
    tstart = clock();
    for (unsigned int r = 0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b< BasisArray.size(); ++b)
        {
            for (unsigned int i = 0; i< BasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = BasisArray[b]->get_wavelet(i);
                for (unsigned int level=BasisArray[b]->j0(); level <= BasisArray[b]->j0()+maxlevelrange; ++level)
                {
                    gen = true;
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                    intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                    gen = false;
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersecting1);
                    //intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersecting2);
                    //intersecting_wavelets(*BasisArray[b], temp_ind1, level, gen, intersectingpairs1);
                    intersecting_wavelets2(*BasisArray[b], temp_ind1, level, gen, intersectingpairs2);
                }
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; old code " << time1 << "sec; new code " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
#endif
    return 0;
}

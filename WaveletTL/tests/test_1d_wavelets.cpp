/*
 * This file prints informations about one dimensional wavelet bases to cout.
 * E.g. number of generators/wavelets on a level, first/last wavelet on a level, which wavelets do not vanish at boundary ...
 */

#define _TEST_WITH_PBASIS 1
#define _D 3
#define _DT 3
#define _BCLEFT 0
#define _BCRIGHT 1

#include <iostream>
#if _TEST_WITH_PBASIS == 1

#include <interval/p_basis.h>
#include <interval/p_support.h>

#else

#include <interval/ds_basis.h>
#include <interval/ds_support.h>

#endif


using namespace std;
using namespace WaveletTL;

using namespace std;

/*
 * 
 */
int main()
{

    cout << "Checking criteria for a wavelet to be a true boundary wavelet" << endl;

    const int d = _D;
    const int dT = _DT;

#if _TEST_WITH_PBASIS == 1
    typedef PBasis<d,dT> Basis;
#else
    typedef DSBasis<d,dT> Basis;
#endif
    typedef Basis::Index Index;
    typedef Basis::Support Support;

    //Basis basis;
    Basis basis(_BCLEFT, _BCRIGHT);
    //Basis basis( false, false);

    if (_TEST_WITH_PBASIS == 1)
    {
        cout << "testing PBasis" << endl;
    }
    else
    {
        cout << "testing DSBasis" << endl;
    }
    cout << "d = " << d << " dT = " << dT << " s0 = " << basis.get_s0() << " left_BC = " << _BCLEFT << " s1 = " << basis.get_s1() << " right BC = " << _BCRIGHT << " j0 = " << basis.j0() << endl;

    cout << "values of Delta and Nabla functions"<< endl;
    cout << "DeltaLmin() " << basis.DeltaLmin() << endl;
    cout << "DeltaLmax() " << basis.DeltaLmax() << endl;
    cout << "Delta0min() " << basis.Delta0min() << endl;
    cout << "Delta0max " << basis.Delta0max(basis.j0()) << " " << basis.Delta0max(basis.j0()+1) << " " << basis.Delta0max(basis.j0()+2) << endl;
    cout << "DeltaRmin " << basis.DeltaRmin(basis.j0()) << " " << basis.DeltaRmin(basis.j0()+1) << " " << basis.DeltaRmin(basis.j0()+2) << endl;
    cout << "DeltaRmax " << basis.DeltaRmax(basis.j0()) << " " << basis.DeltaRmax(basis.j0()+1) << " " << basis.DeltaRmax(basis.j0()+2) << endl;

    cout << "DeltaLTmin() " << basis.DeltaLTmin() << endl;
    cout << "DeltaLTmax() " << basis.DeltaLTmax() << endl;
    cout << "DeltaRTmin " << basis.DeltaRTmin(basis.j0()) << " " << basis.DeltaRTmin(basis.j0()+1) << " " << basis.DeltaRTmin(basis.j0()+2) << endl;
    cout << "DeltaRTmax " << basis.DeltaRTmax(basis.j0()) << " " << basis.DeltaRTmax(basis.j0()+1 ) << " " << basis.DeltaRTmax(basis.j0()+2) << endl;
    cout << "Deltasize " << basis.Deltasize(basis.j0()) << " " << basis.Deltasize(basis.j0()+1) << " " << basis.Deltasize(basis.j0()+2) << endl;
    cout << "Nablamin() " << basis.Nablamin() << endl;
    cout << "Nablamax " << basis.Nablamax(basis.j0()) << " " << basis.Nablamax(basis.j0()+1) << " " << basis.Nablamax(basis.j0()+2) << endl;
    cout << "Nablasize " << basis.Nablasize(basis.j0()) << " " << basis.Nablasize(basis.j0()+1) << " " << basis.Nablasize(basis.j0()+2) << endl;


#if 1


    for (int level = basis.j0(); level <= basis.j0()+1; level++)
    {
        cout << "locate boundary generators on level j=" << level << endl;
        Index lambda(first_generator(&basis, level));
        Support supp;
        for (;; ++lambda) 
        {
            support(basis, lambda, supp.k1, supp.k2);
            cout << lambda << " " << basis.evaluate(0,lambda,0.0)
                 //<< " " << basis.evaluate(0,lambda,0.0000001)
                 //<< " " << basis.evaluate(0,lambda,0.9999999)
                 << " " << basis.evaluate(0,lambda,1.0)
                 << " supp = 2^{-" << lambda.j()+lambda.e() << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
            if (lambda == last_generator(&basis, level)) break;
        }
        cout << "locate boundary wavelets on level j=" << level << endl;
        lambda = first_wavelet(&basis, level);
        for (;; ++lambda)
        {
            support(basis, lambda, supp.k1, supp.k2);
            cout << lambda << " " << basis.evaluate(0,lambda,0.0)
                 //<< " " << basis.evaluate(0,lambda,0.0000001)
                 //<< " " << basis.evaluate(0,lambda,0.9999999)
                 << " " << basis.evaluate(0,lambda,1.0)
                 << " supp = 2^{-" << lambda.j()+lambda.e() << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
            if (lambda == last_wavelet(&basis, level)) break;
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
    return 0;
}


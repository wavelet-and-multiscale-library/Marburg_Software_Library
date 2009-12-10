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
    typedef DSBasis<2,2> Basis1D;
//  typedef PBasis<2,2> Basis1D;
//   typedef JLBasis Basis1D;

    const unsigned int DIM = 2;
    typedef TensorBasis<Basis1D,DIM> Basis;
    typedef Basis::Index Index;
#if 1
    Basis basis;
#else
    FixedArray1D<bool,4> bc;
    bc[0] = bc[1] = true;
    bc[2] = bc[3] = true;
    Basis basis(bc);
#endif

#if 1
    cout << "- testing calculation of supports:" << endl;
    Basis::Support supp;
    for (Index lambda(first_generator<Basis1D,2,Basis>(&basis));; ++lambda) {
        support<Basis1D,2>(basis, lambda, supp);
        //support(basis, lambda, supp);
        cout << lambda << " has support 2^{-(" << supp.j[0] << ", " << supp.j[1] << ")}"
             << "[" << supp.a[0] << "," << supp.b[0]
             << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
             << endl;
        //if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0()+1)) break;
        if (lambda.number() == 20) break;
    }
#endif
#if 0
    unsigned int von(143),bis(143);
    cout << "- compute all intersecting wavelets"<<endl;
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
        intersecting_wavelets<Basis1D,2>(basis, lambda, j, true, nus);
        /*
        // generators:
        for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
            Basis::Support supp;
            intersect_supports(basis, lambda, *it, supp);
            cout << "    nu=" << *it
                    << " with support intersection "
                    << "2^{-{" << supp.j[0] << ", " << supp.j[1] << "}}"
                    << "[" << supp.a[0] << "," << supp.b[0]
                    << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
                    << endl;
        }
         * */
        //compute intersecting wavelets and increase j until j is too high
        level_type level = basis.j0();
        while (multi_degree(level) <= multi_degree(basis.j0())+1) // +2)
        {
            intersecting_wavelets<Basis1D,2>(basis, lambda, level, false, nus);
            for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
                Basis::Support supp;
                intersect_supports(basis, lambda, *it, supp);
                cout << "    nu=" << *it
                        << " with support intersection "
                        << "2^{-{" << supp.j[0] << ", " << supp.j[1] << "}}"
                        << "[" << supp.a[0] << "," << supp.b[0]
                        << "]x[" << supp.a[1] << "," << supp.b[1] << "]"
                        << endl;
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
                        level[i]=level[i]+1;
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
        if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0())) break;
    }
#endif
#if 1
    Index lambda(51,&basis);
    cout << "Testing intersect_singular_support with lambda = " << lambda << endl;

    //for (Index it = first_generator<Basis1D,2,Basis>(&basis);; ++lambda)
    for (Index it = first_generator<Basis1D,2,Basis>(&basis);it <= last_wavelet<Basis1D,2,Basis>(&basis,basis.j0()); ++it)
    {
        if (intersect_singular_support<Basis1D,2>(basis,lambda,it))
            cout << "intersection with nu = " << it << endl;
    }
#endif
}

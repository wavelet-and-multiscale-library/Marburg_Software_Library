/*
 * This file is intended to test qtbasis
 * 
 * Several aspects are tested. Which tests are active is determined by MACROS
 * 
 * check out test_cached_qtproblem for useful domain decompositions
 */

#include <iostream>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>

#include <interval/p_basis.h>
#include <interval/p_evaluate.h>
#include <interval/ds_basis.h>
#include <interval/ds_evaluate.h>

#include <cube/tbasis.h>

#include <general_domain/qtbasis.h>
#include <general_domain/qtbasis_index.h>
using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{

    cout << "Testing qtbasis." << endl;

#define d  3
#define dT 3
    const int num_of_patches = 1;
    const int DIM = 2; const int maxleveloffset = 2;
    //const int DIM = 3; const int maxleveloffset = 3;
    /* implemented Decompositions:
     * 0 :: 1 patch
     * 1 :: 2 patches, homogeneous BC
     */
    const int decomposition = 0;
    

    
    //typedef DSBasis<3,5> Basis1d;
    typedef PBasis<d,dT> Basis1d;
    typedef QTBasis<Basis1d,DIM> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;
    typedef TensorBasis<Basis1d,DIM> TBasis;
    typedef TBasis::Index TIndex;
    typedef TBasis::Support TSupport;

    // construct the basis
    Array1D<Point<DIM,int> > corners;
    Array1D<FixedArray1D<int,2*DIM> > neighbours;
    Array1D<FixedArray1D<bool,2*DIM> > bc_bool;
    Array1D<FixedArray1D<int,2*DIM> > bc_int;

    corners.resize(num_of_patches);
    neighbours.resize(num_of_patches);
    bc_bool.resize(num_of_patches);
    bc_int.resize(num_of_patches);

    for (unsigned int p = 0; p < num_of_patches; ++p)
    {
        for (unsigned int i = 0; i < DIM; ++i)
        {
            // initialize with no ...
            neighbours[p][2*i] = neighbours[p][2*i+1] =-1;; // ... left and right neighbours
            bc_bool[p][2*i] = bc_bool[p][2*i+1] = true;
            bc_int[p][2*i] = bc_int[p][2*i+1] = 1;
        }
    }

    switch (decomposition) {
        case 0: // 1 domain
            assert (num_of_patches == 1);
            break;
        case 1: // 2 domains, homogeneous BC, left patch is extended to the right
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[1][0] = 1;
            }
            else if (DIM == 2)
            {
                corners[1][0] = 1;
                corners[1][1] = 0;
            }
            neighbours[0][1] = 1; // patch 0, x-right-neighbour = patch 1
            neighbours[1][0] = 0; // patch 1, x-left-neighbour = patch 0
            bc_bool[0][1] = false;
            bc_int[0][1] = 0;
            break;
        case 2: // 2 domains, homogeneous BC, patch 1 is extended to patch 0 (left in the y direction for DIM=2)
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (DIM == 2)
            {
                corners[0][0] = 3;
                corners[0][1] = 4;
                corners[1][0] = 3;
                corners[1][1] = 5;
            }
            neighbours[0][3] = 1; // patch 0 has a neighbour in the south
            neighbours[1][2] = 0; // patch 1 has a neighbour in the north
            bc_bool[1][2] = false;
            bc_int[1][2] = 0;
            break;
        case 3: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 1
            assert (num_of_patches == 3);
            assert (DIM == 2);
            corners[0][0] = 0;
            corners[0][1] = 1;
            corners[1][0] = 0;
            corners[1][1] = 0;
            corners[2][0] = 1;
            corners[2][1] = 0;
            
            neighbours[0][2] = 1; // patch 0 has a neighbour in the south
            neighbours[1][3] = 0; // patch 1 has a neighbour in the north
            neighbours[1][1] = 2; // patch 1 has a neighbour in the east
            neighbours[2][0] = 1; // patch 2 has a neighbour in the west
            bc_bool[0][2] = false;
            bc_bool[2][0] = false;
            bc_int[0][2] = 0;
            bc_int[2][0] = 0;
            break;
        case 4: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 0
            assert (num_of_patches == 3);
            assert (DIM == 2);
            corners[0][0] = 2;
            corners[0][1] = 2;
            corners[1][0] = 2;
            corners[1][1] = 1;
            corners[2][0] = 3;
            corners[2][1] = 2;
            
            neighbours[0][2] = 1; // patch 0 has a neighbour in the south
            neighbours[0][1] = 2; // patch 0 has a neighbour in the east
            neighbours[1][3] = 0; // patch 1 has a neighbour in the north
            neighbours[2][0] = 0; // patch 2 has a neighbour in the west
            bc_bool[0][2] = false;
            bc_bool[2][0] = false;
            bc_int[0][2] = 0;
            bc_int[2][0] = 0;
            break;
        default:
            cout << "main:: error! no decomposition specified!" << endl;
            abort();
            break;
    }
    
    cout << "constructor called with parameters:" << endl;
    for (unsigned int i=0; i< num_of_patches; ++i)
    {
        cout << "Patch number "<<i<<":"<<endl;
        cout << "corners["<<i<<"]="<<corners[i]<<endl;
        cout << "neighbours["<<i<<"]= ";
        for (unsigned int j=0; j<DIM-1;++j)
        {
            cout << "{"<< neighbours[i][2*j] << ", " << neighbours[i][2*j+1] <<"}x";
        }
        cout << "{"<< neighbours[i][2*(DIM-1)] << ", " << neighbours[i][2*(DIM-1)+1] <<"}" << endl;

        cout << "bc_bool["<<i<<"]= ";
        for (unsigned int j=0; j<DIM-1;++j)
        {
            cout << "{"<< bc_bool[i][2*j] << ", " << bc_bool[i][2*j+1] <<"}x";
        }
        cout << "{"<< bc_bool[i][2*(DIM-1)] << ", " << bc_bool[i][2*(DIM-1)+1] <<"}" << endl;

        cout << "bc_int["<<i<<"]= ";
        for (unsigned int j=0; j<DIM-1;++j)
        {
            cout << "{"<< bc_int[i][2*j] << ", " << bc_int[i][2*j+1] <<"}x";
        }
        cout << "{"<< bc_int[i][2*(DIM-1)] << ", " << bc_int[i][2*(DIM-1)+1] <<"}" << endl;
    }
    


#if 0
    cout << "test series 1: test equality of default, bc_bool, and bc_int constructor (on the cube)" << endl;
    cout << "test series: check whether default constructor == constructor with bool BC == constructor with int BC" << endl;
    FixedArray1D<Basis*,3> QTBasisArray;
    
    cout << "default constructor (with homogeneous b.c.)" << endl;
    Basis basis_def;
    cout << "done" << endl;

    cout << "Call constructor with bc_bool" << endl;
    Basis basis_bool(corners,neighbours,bc_bool);
    cout << "done" << endl;

    cout << "Call constructor with bc_int" << endl;
    Basis basis_int(corners,neighbours,bc_int);
    cout << "done" << endl;
    
    QTBasisArray[0] = &basis_def;
    QTBasisArray[1] = &basis_bool;
    QTBasisArray[2] = &basis_int;
#endif
    
#if 1
    cout << "test series 2: check whether basis with bool BC corresponds with tbasis (for all possible BC, on the cube)" << endl;
    assert (DIM == 2);
    FixedArray1D<Basis*,16> QTBasisArray;
    FixedArray1D<TBasis*,16> TBasisArray;
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = false;
    bc_bool[0][2] = false;
    bc_bool[0][3] = false;
    Basis qtbasis0000(corners, neighbours,bc_bool);
    TBasis tbasis0000(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = false;
    bc_bool[0][2] = false;
    bc_bool[0][3] = true;
    Basis qtbasis0001(corners, neighbours,bc_bool);
    TBasis tbasis0001(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = false;
    bc_bool[0][2] = true;
    bc_bool[0][3] = false;
    Basis qtbasis0010(corners, neighbours,bc_bool);
    TBasis tbasis0010(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = true;
    bc_bool[0][2] = false;
    bc_bool[0][3] = false;
    Basis qtbasis0100(corners, neighbours,bc_bool);
    TBasis tbasis0100(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = false;
    bc_bool[0][2] = false;
    bc_bool[0][3] = false;
    Basis qtbasis1000(corners, neighbours,bc_bool);
    TBasis tbasis1000(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = false;
    bc_bool[0][2] = true;
    bc_bool[0][3] = true;
    Basis qtbasis0011(corners, neighbours,bc_bool);
    TBasis tbasis0011(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = true;
    bc_bool[0][2] = false;
    bc_bool[0][3] = true;
    Basis qtbasis0101(corners, neighbours,bc_bool);
    TBasis tbasis0101(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = true;
    bc_bool[0][2] = true;
    bc_bool[0][3] = false;
    Basis qtbasis0110(corners, neighbours,bc_bool);
    TBasis tbasis0110(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = false;
    bc_bool[0][2] = false;
    bc_bool[0][3] = true;
    Basis qtbasis1001(corners, neighbours,bc_bool);
    TBasis tbasis1001(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = false;
    bc_bool[0][2] = true;
    bc_bool[0][3] = false;
    Basis qtbasis1010(corners, neighbours,bc_bool);
    TBasis tbasis1010(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = true;
    bc_bool[0][2] = false;
    bc_bool[0][3] = false;
    Basis qtbasis1100(corners, neighbours,bc_bool);
    TBasis tbasis1100(bc_bool[0]);
    
    bc_bool[0][0] = false;
    bc_bool[0][1] = true;
    bc_bool[0][2] = true;
    bc_bool[0][3] = true;
    Basis qtbasis0111(corners, neighbours,bc_bool);
    TBasis tbasis0111(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = false;
    bc_bool[0][2] = true;
    bc_bool[0][3] = true;
    Basis qtbasis1011(corners, neighbours,bc_bool);
    TBasis tbasis1011(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = true;
    bc_bool[0][2] = false;
    bc_bool[0][3] = true;
    Basis qtbasis1101(corners, neighbours,bc_bool);
    TBasis tbasis1101(bc_bool[0]);
    
    bc_bool[0][0] = true;
    bc_bool[0][1] = true;
    bc_bool[0][2] = true;
    bc_bool[0][3] = false;
    Basis qtbasis1110(corners, neighbours,bc_bool);
    TBasis tbasis1110(bc_bool[0]);
    
    bc_bool[0][0] = bc_bool[0][1] = bc_bool[0][2] = bc_bool[0][3] = true;
    Basis qtbasis1111(corners, neighbours,bc_bool);
    TBasis tbasis1111(bc_bool[0]);
    
    
    QTBasisArray[0]  = &qtbasis0000;
    QTBasisArray[1]  = &qtbasis0001;
    QTBasisArray[2]  = &qtbasis0010;
    QTBasisArray[3]  = &qtbasis0100;
    QTBasisArray[4]  = &qtbasis1000;
    QTBasisArray[5]  = &qtbasis0011;
    QTBasisArray[6]  = &qtbasis0101;
    QTBasisArray[7]  = &qtbasis0110;
    QTBasisArray[8]  = &qtbasis1001;
    QTBasisArray[9]  = &qtbasis1010;
    QTBasisArray[10] = &qtbasis1100;
    QTBasisArray[11] = &qtbasis0111;
    QTBasisArray[12] = &qtbasis1011;
    QTBasisArray[13] = &qtbasis1101;
    QTBasisArray[14] = &qtbasis1110;
    QTBasisArray[15] = &qtbasis1111;
    
    TBasisArray[0]  = &tbasis0000;
    TBasisArray[1]  = &tbasis0001;
    TBasisArray[2]  = &tbasis0010;
    TBasisArray[3]  = &tbasis0100;
    TBasisArray[4]  = &tbasis1000;
    TBasisArray[5]  = &tbasis0011;
    TBasisArray[6]  = &tbasis0101;
    TBasisArray[7]  = &tbasis0110;
    TBasisArray[8]  = &tbasis1001;
    TBasisArray[9]  = &tbasis1010;
    TBasisArray[10] = &tbasis1100;
    TBasisArray[11] = &tbasis0111;
    TBasisArray[12] = &tbasis1011;
    TBasisArray[13] = &tbasis1101;
    TBasisArray[14] = &tbasis1110;
    TBasisArray[15] = &tbasis1111;
    
    cout << "main :: setting jmax for the tbases" << endl;
    for (unsigned int b=0; b<TBasisArray.size(); ++b)
    {
        TBasisArray[b]->set_jmax(multi_degree(TBasisArray[b]->j0())+maxleveloffset);
    }
    cout << "done" << endl;
#endif
    
    cout << "main :: setting jmax" << endl;
    for (unsigned int b = 0; b < QTBasisArray.size(); ++b)
    {
        QTBasisArray[b] ->set_jmax( multi_degree(QTBasisArray[b]->j0()[0])+maxleveloffset);
    }
    
    cout << "degrees_of_freedom = " << endl;
    for (unsigned int b = 0; b < QTBasisArray.size(); ++b)
    {
        cout << "b = " << b << "; DOF = " << QTBasisArray[b]->degrees_of_freedom() << endl;
    }
    
    cout << "values of Delta and Nabla functions"<< endl;
    for (unsigned int b = 0; b < QTBasisArray.size(); ++b)
    {
        for (unsigned int p = 0; p < num_of_patches; ++p)
        {
            cout << "QTBasis = " << b << "; Patch = " << p << endl;
            cout << "DeltaLmin() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaLmin() << " ";}
            cout << endl;

            cout << "DeltaLmax() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaLmax() << " ";}
            cout << endl;

            cout << "Delta0min() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Delta0min() << " ";}
            cout << endl;

            cout << "Delta0max(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Delta0max(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "DeltaRmin(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaRmin(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "DeltaRmax(j0)";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaRmax(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "DeltaLTmin() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaLTmin() << " ";}
            cout << endl;

            cout << "DeltaLTmax() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaLTmax() << " ";}
            cout << endl;

            cout << "DeltaRTmin(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaRTmin(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "DeltaRTmax(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->DeltaRTmax(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "Deltasize(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Deltasize(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "Nablamin() ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Nablamin() << " ";}
            cout << endl;

            cout << "Nablamax(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Nablamax(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

            cout << "Nablasize(j0) ";
            for (unsigned int i=0; i<DIM; ++i) { cout << QTBasisArray[b]->bases()[p][i]->Nablasize(QTBasisArray[b]->bases()[p][i]->j0()) << " ";}
            cout << endl;

        }
    }

#if 0
    
    cout << "Test series 1" << endl;
    
    if (QTBasisArray.size() > 1)
    {
        cout << "Check equality of all bases in QTBasisArray" << endl;
        cout << "Intended to compare _bool and _int BC constructors (maybe also with default constructor)" << endl;
        for (unsigned int b=0; b < QTBasisArray.size()-1; ++b)
        {
            assert (QTBasisArray[b]->get_nop() == QTBasisArray[b+1]->get_nop());
            for (unsigned int p=0; p<QTBasisArray[b]->get_nop(); ++p)
            {
                assert (QTBasisArray[b]->j0()[p] == QTBasisArray[b+1]->j0()[p]);
            }
            assert (QTBasisArray[b]->get_jmax() == QTBasisArray[b+1]->get_jmax());
            assert (QTBasisArray[b]->degrees_of_freedom() == QTBasisArray[b+1]->degrees_of_freedom());
            assert (QTBasisArray[b]->get_numberoflevels() == QTBasisArray[b+1]->get_numberoflevels());
            for (unsigned int levelnum = 0; levelnum < QTBasisArray[0]->get_numberoflevels(); ++levelnum)
            {
                assert (QTBasisArray[b]->first_wavelet(levelnum) == QTBasisArray[b+1]->first_wavelet(levelnum));
                assert (QTBasisArray[b]->last_wavelet(levelnum) == QTBasisArray[b+1]->last_wavelet(levelnum));
                MultiIndex<int,DIM> temp_mi1, temp_mi2;
                unsigned int temp_i1, temp_i2;
                QTBasisArray[b]->get_level(levelnum,temp_mi1,temp_i1);
                QTBasisArray[b+1]->get_level(levelnum,temp_mi2,temp_i2);
                assert ((temp_mi1 == temp_mi2) && (temp_i1 == temp_i2));
            }
            assert (QTBasisArray[b]->first_generator() == QTBasisArray[b+1]->first_generator() );
            for (unsigned int lambdanum =0; lambdanum< QTBasisArray[b]->degrees_of_freedom();++lambdanum)
            {
                Index temp_ind1(QTBasisArray[b]->get_wavelet(lambdanum));
                Index temp_ind2(QTBasisArray[b+1]->get_wavelet(lambdanum));
                assert (temp_ind1 == temp_ind2);
                cout << "lambdanum = " << lambdanum << "; lambda = " << temp_ind1 << endl;
                assert (QTBasisArray[b]->get_levelnum(temp_ind1.j(),temp_ind1.p()) == QTBasisArray[b+1]->get_levelnum(temp_ind2.j(),temp_ind2.p()));
                Support temp_supp1, temp_supp2;
                QTBasisArray[b]->support_local(temp_ind1,temp_supp1);
                QTBasisArray[b+1]->support_local(temp_ind2,temp_supp2);
                for (unsigned int i=0; i<DIM;++i)
                {
                    assert (( (temp_supp1.j[i] == temp_supp2.j[i]) && (temp_supp1.a[i] == temp_supp2.a[i])) && (temp_supp1.b[i] == temp_supp2.b[i]));
                }
                QTBasisArray[b]->support_full(temp_ind1,temp_supp1);
                QTBasisArray[b+1]->support_full(temp_ind2,temp_supp2);
                for (unsigned int i=0; i<DIM;++i)
                {
                    assert (( (temp_supp1.j[i] == temp_supp2.j[i]) && (temp_supp1.a[i] == temp_supp2.a[i])) && (temp_supp1.b[i] == temp_supp2.b[i]));
                }
                
                
                //unsigned int steps = 30;
                unsigned int steps_per_patch(4);
                
                
                int index[DIM]; 
                for (unsigned int i = 0; i < DIM; i++)
                {
                    index[i] = 0;
                }
                Point<DIM> x;
                FixedArray1D<double,DIM> h;
                for (unsigned int i=0; i<DIM;i++)
                {
                    h[i]=ldexp(1.0, -temp_supp1.j[i]); // granularity for the quadrature
                }
                FixedArray1D<Array1D<double>,DIM> points;
        
                for (unsigned int i = 0; i < DIM; i++) 
                {
                    points[i].resize((temp_supp1.b[i]-temp_supp1.a[i])*steps_per_patch);
                    for (int patch = temp_supp1.a[i]; patch < temp_supp1.b[i]; patch++)
                    {
                        for (int n = 0; n < steps_per_patch; n++) 
                        {
                            points[i][(patch-temp_supp1.a[i])*steps_per_patch+n]
                                    = h[i]*(patch+n/(double)steps_per_patch);
                        }
                    }
                }
                double temp_d1, temp_d2;
                while (true)
                {
                    for (unsigned int i=0; i<DIM; ++i)
                    { 
                        //cout << " i = " << i << "; index[" << i << "] = " << index[i] << "; points[" << i <<"].size() = " << points[i].size() << ";" << endl;
                        x[i] = points[i][index[i]];
                    }
                    //cout << "x = (" << x[0] << ", " << x[1] << ")" << endl;
                    temp_d1 = QTBasisArray[b]->evaluate_check(0,lambdanum, x);
                    temp_d2 = QTBasisArray[b+1]->evaluate_check(0,lambdanum, x);
                    assert (temp_d1 == temp_d2); // will probably fail because of inexact evaluation
                    temp_d1 = QTBasisArray[b]->evaluate_check(1,lambdanum, x);
                    temp_d2 = QTBasisArray[b+1]->evaluate_check(1,lambdanum, x);
                    assert (temp_d1 == temp_d2);
                    temp_d1 = QTBasisArray[b]->evaluate_simple(0,lambdanum, x);
                    temp_d2 = QTBasisArray[b+1]->evaluate_simple(0,lambdanum, x);
                    assert (temp_d1 == temp_d2); // will probably fail because of inexact evaluation
                    temp_d1 = QTBasisArray[b]->evaluate_simple(1,lambdanum, x);
                    temp_d2 = QTBasisArray[b+1]->evaluate_simple(1,lambdanum, x);
                    assert (temp_d1 == temp_d2);
                    
                    
                    // "++index"
                    
                    bool exit = false;
                    for (unsigned int i = 0; i < DIM; i++) {
                        if (index[i] == (temp_supp1.b[i]-temp_supp1.a[i])*steps_per_patch -1) {
                            index[i] = 0;
                            exit = (i == DIM-1);
                        } else {
                            index[i]++;
                            break;
                        }
                    }
                    if (exit) break;
                }
            }
            assert (QTBasisArray[b]->primal_regularity() == QTBasisArray[b]->primal_regularity());
            assert (QTBasisArray[b]->primal_polynomial_degree() == QTBasisArray[b]->primal_polynomial_degree());
            assert (QTBasisArray[b]->primal_vanishing_moments() == QTBasisArray[b]->primal_vanishing_moments());
        }
        
        cout << "integrate and expand are not tested" << endl;
    }
#endif
    
#if 1
    cout << "Test series 2" << endl;
    cout << "Check equality of the bases in QTBasisArray and TBasisArray" << endl;
    cout << "Also checked: evaluate_check == evaluate_simple for QTBasis" << endl;
    cout << "It is assumed that QTBasisArray contains bases on the unit cube which correspond to the basis in TBasisArray" << endl;

    cout << "********************************" << endl;
    cout << "Test 1" << endl;
    cout << "********************************" << endl;
    cout << "(this will take some time)" << endl;
    assert (QTBasisArray.size() == TBasisArray.size());
    for (unsigned int b=0; b < QTBasisArray.size(); ++b)
    {
        cout << "--------------------------------" << endl;
        cout << "compare basis number b = " << b << endl;
        cout << "--------------------------------" << endl;
        assert (QTBasisArray[b]->get_nop() == 1);
        assert (QTBasisArray[b]->j0()[0] == TBasisArray[b]->j0());
        assert (QTBasisArray[b]->get_jmax() == TBasisArray[b]->get_jmax());
        assert (QTBasisArray[b]->degrees_of_freedom() == TBasisArray[b]->degrees_of_freedom());
        cout << "lambdanum = ";
        cout.flush();
        for (unsigned int lambdanum =0; lambdanum< QTBasisArray[b]->degrees_of_freedom();++lambdanum)
        {
            Index temp_ind1(QTBasisArray[b]->get_wavelet(lambdanum));
            TIndex temp_tind1(TBasisArray[b]->get_wavelet(lambdanum));
            //cout << "--------------------------------" << endl;
            //cout << "lambdanum = " << lambdanum << "; lambda = " << temp_ind1 << endl;
            //cout << "--------------------------------" << endl;
            cout << lambdanum << " ";
            cout.flush();
            assert (temp_ind1.j() == temp_tind1.j());
            assert (temp_ind1.e() == temp_tind1.e());
            assert (temp_ind1.k() == temp_tind1.k());
            Support temp_supp1;
            TSupport temp_tsupp1;
            QTBasisArray[b]->support_local(temp_ind1,temp_supp1);
            TBasisArray[b]->support(temp_tind1,temp_tsupp1);
            for (unsigned int i=0; i<DIM;++i)
            {
                assert (( (temp_supp1.j[i] == temp_tsupp1.j[i]) && (temp_supp1.a[i] == temp_tsupp1.a[i])) && (temp_supp1.b[i] == temp_tsupp1.b[i]));
            }
            //unsigned int steps = 30;
            unsigned int steps_per_patch(4);
            int index[DIM]; 
            for (unsigned int i = 0; i < DIM; i++)
            {
                index[i] = 0;
            }
            Point<DIM> x;
            FixedArray1D<double,DIM> h;
            for (unsigned int i=0; i<DIM;i++)
            {
                h[i]=ldexp(1.0, -temp_supp1.j[i]); // granularity for the quadrature
            }
            FixedArray1D<Array1D<double>,DIM> points;

            for (unsigned int i = 0; i < DIM; i++) 
            {
                points[i].resize((temp_supp1.b[i]-temp_supp1.a[i])*steps_per_patch);
                for (int patch = temp_supp1.a[i]; patch < temp_supp1.b[i]; patch++)
                {
                    for (int n = 0; n < steps_per_patch; n++) 
                    {
                        points[i][(patch-temp_supp1.a[i])*steps_per_patch+n]
                                = h[i]*(patch+n/(double)steps_per_patch);
                    }
                }
            }
            double temp_d1, temp_d2;
            while (true)
            {
                for (unsigned int i=0; i<DIM; ++i)
                { 
                    //cout << " i = " << i << "; index[" << i << "] = " << index[i] << "; points[" << i <<"].size() = " << points[i].size() << ";" << endl;
                    x[i] = points[i][index[i]];
                }
                //cout << "x = (" << x[0] << ", " << x[1] << ")" << endl;
                temp_d1 = QTBasisArray[b]->evaluate_check(0,lambdanum, x);
                temp_d2 = QTBasisArray[b]->evaluate_simple(0,lambdanum, x);
                assert (temp_d1 == temp_d2);
                temp_d2 = TBasisArray[b]->evaluate(0,temp_tind1, x);
                assert (temp_d1 == temp_d2); // will probably fail because of inexact evaluation
                temp_d1 = QTBasisArray[b]->evaluate_check(1,lambdanum, x);
                temp_d2 = QTBasisArray[b]->evaluate_simple(1,lambdanum, x);
                assert (temp_d1 == temp_d2);
                temp_d2 = TBasisArray[b]->evaluate(1,temp_tind1, x);
                assert (temp_d1 == temp_d2);
                // "++index"
                bool exit = false;
                for (unsigned int i = 0; i < DIM; i++) {
                    if (index[i] == (temp_supp1.b[i]-temp_supp1.a[i])*steps_per_patch -1) {
                        index[i] = 0;
                        exit = (i == DIM-1);
                    } else {
                        index[i]++;
                        break;
                    }
                }
                if (exit) break;
            }
        }
        cout << endl;
        cout.flush();
    }
    cout << "********************************" << endl;
    cout << "Test 2" << endl;
    cout << "********************************" << endl;
    for (unsigned int b=0; b < QTBasisArray.size(); ++b)
    {
        cout << "--------------------------------" << endl;
        cout << "compare basis number b = " << b << endl;
        cout << "--------------------------------" << endl;
        MultiIndex<int,DIM> temp_mi1;
        unsigned int temp_i;
        Index temp_ind1;
        TIndex temp_tind1;
        QTBasisArray[b]->get_level(QTBasisArray[b]->get_numberoflevels()-1, temp_mi1, temp_i);
        assert (multi_degree(temp_mi1) == TBasisArray[b]->get_jmax());
        for (unsigned int levelnum = 0; levelnum < QTBasisArray[b]->get_numberoflevels(); ++levelnum)
        {
            QTBasisArray[b]->get_level(levelnum, temp_mi1, temp_i);
            temp_ind1 = QTBasisArray[b]->first_wavelet(levelnum);
            if (levelnum == 0)
            {
                temp_tind1 = TBasisArray[b]->first_generator();
            }
            else
            {
                temp_tind1 = TBasisArray[b]->first_wavelet(temp_mi1);
            }
            //cout << "levelnum = " << levelnum << "; temp_ind1 = " << temp_ind1 << "; temp_tind1 = " << temp_tind1 << endl;
            assert (((temp_ind1.j() == temp_tind1.j())
                &&
                (temp_ind1.e() == temp_tind1.e()))
                &&
                (temp_ind1.k() == temp_tind1.k()));
            temp_ind1 = QTBasisArray[b]->last_wavelet(levelnum);
            temp_tind1 = TBasisArray[b]->last_wavelet(temp_mi1);
            assert (((temp_ind1.j() == temp_tind1.j())
                &&
                (temp_ind1.e() == temp_tind1.e()))
                &&
                (temp_ind1.k() == temp_tind1.k()));
        }
        temp_ind1 = QTBasisArray[b]->first_generator();
        temp_tind1 = TBasisArray[b]->first_generator(TBasisArray[b]->j0());
        assert (((temp_ind1.j() == temp_tind1.j())
                &&
                (temp_ind1.e() == temp_tind1.e()))
                &&
                (temp_ind1.k() == temp_tind1.k()));
        assert (QTBasisArray[b]->primal_regularity() == TBasisArray[b]->primal_regularity());
        assert (QTBasisArray[b]->primal_polynomial_degree() == TBasisArray[b]->primal_polynomial_degree());
        assert (QTBasisArray[b]->primal_vanishing_moments() == TBasisArray[b]->primal_vanishing_moments());
    }
    cout << "Not compared:" << endl;
    cout << "QTBasis and TBasis: integrate, expand" << endl;
    cout << "TBasis: last_generator, first_wavelet(int), last_wavelet(int)" << endl;
        
#endif
    
   
#if 0
    cout << "Plot results of setup_full_collection and precompute_firstlast_wavelets" << endl;
    
    cout << "basis.first_generator() = " << basis.first_generator() << "; .number() = " << basis.first_generator().number() << endl;
    for (unsigned int l = 0; l<basis.get_numberoflevels();++l)
    {
        cout << "basis.first_wavelet(" << l << ") = " << basis.first_wavelet(l) << "; .number() = " << basis.first_wavelet(l).number() << "; "; //endl;
        cout << "*basis.get_wavelet(" << basis.first_wavelet(l).number() << ") = " << *basis.get_wavelet(basis.first_wavelet(l).number()) << endl;
        cout << "basis.last_wavelet(" << l << ") = " << basis.last_wavelet(l) << "; .number() = " << basis.last_wavelet(l).number() << "; "; //endl;
        cout << "*basis.get_wavelet(" << basis.last_wavelet(l).number() << ") = " << *basis.get_wavelet(basis.last_wavelet(l).number()) << endl;
    }
    
#endif
#if 0
    cout << "Compare computed first/last wavelets with the results given by ++" << endl;
    Index temp_ind;
    MultiIndex<int,DIM> temp_mi1, temp_mi2, temp_mi3;
    unsigned int temp_i1, temp_i2;
    for (unsigned int l = 0; l<basis.get_numberoflevels();++l)
    {
        basis.get_level(l,temp_mi1,temp_i1);
        if (l<basis.get_numberoflevels()-1)
        {
            basis.get_level(l+1,temp_mi2,temp_i2);
            assert ( (temp_mi1 < temp_mi2) || ((temp_mi1 == temp_mi2) && (temp_i1<temp_i2)) ); // assert that we increase (j,p) in each step
        }
        assert (l == basis.get_levelnum(temp_mi1,temp_i1));
        //cout << "l = " << l << "; (j,p) = (" << temp_mi1 << ", " << temp_i1 << ")" << endl;
        
        temp_ind = basis.first_wavelet(l);
        //cout << "l = " << l << "; first_wavelet["<<l<<"] = " << temp_ind << "; get_wavelet(ind.number) = " << *basis.get_wavelet(temp_ind.number())<< endl;
        assert (temp_ind == *basis.get_wavelet(temp_ind.number()));
        temp_ind = basis.last_wavelet(l);
        //cout << "l = " << l << "; last_wavelet["<<l<<"] = " << temp_ind << "; get_wavelet(ind.number) = " << *basis.get_wavelet(temp_ind.number())<< endl;
        assert (temp_ind == *basis.get_wavelet(temp_ind.number()));
    }

    
    if (num_of_patches == 1)
    {
        for (unsigned int i=0;i<DIM;++i)
            temp_mi1[i]=0;
        cout << "Comparing QTBasis and TBasis for the special case of 1 patch" << endl;
        typedef TensorBasis<Basis1d,DIM> TBasis;
        typedef TensorBasis<Basis1d,DIM>::Index TIndex;
        TBasis tbasis(bc_bool[0]);
        tbasis.set_jmax(multi_degree(basis.j0()[0])+maxleveloffset);
        TIndex temp_tind;
        for (unsigned int l = 0; l<basis.get_numberoflevels(); ++l)
        {
            temp_mi2 = basis.j0()[0];
            for (unsigned int i=0; i<DIM; ++i)
                temp_mi2[i] += temp_mi1[i];
            basis.get_level(l,temp_mi3,temp_i1);
            assert (temp_mi2 == temp_mi3);
            assert (temp_i1 == 0);
            
            assert (l == temp_mi1.number());
            
            temp_ind = basis.first_wavelet(l);
            temp_tind = tbasis.first_wavelet(temp_mi3);
            //cout << "temp_ind = " << temp_ind << "; temp_tind = " << temp_tind << endl;
            assert(temp_ind.j() == temp_tind.j());
            if (l > 0) 
            {
                //on the first level TBasis really gives the first wavelet. However QTBasis "only" gives the first basis function, i.e., the first generator!
                assert(temp_ind.e() == temp_tind.e()); 
                assert(temp_ind.k() == temp_tind.k());
                assert(temp_ind.number() == temp_tind.number());
            }
            
            temp_ind = basis.last_wavelet(l);
            temp_tind = tbasis.last_wavelet(temp_mi3);
            assert(temp_ind.j() == temp_tind.j());
            assert(temp_ind.e() == temp_tind.e());
            assert(temp_ind.k() == temp_tind.k());
            assert(temp_ind.number() == temp_tind.number());
            
            ++temp_mi1;
        }
    }
#endif
    




    /*
     special case DIM = 1 nicht vergessen!
     */
    return 0;
}
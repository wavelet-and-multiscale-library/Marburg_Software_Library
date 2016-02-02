/*
 * This file is intended to test cached_qtproblem
 * Obviously it uses qtbasis. This test file is more advanced than test_qtbasis.
 * 
 * Several aspects are tested. Which tests are active is determined by MACROS
 * _TEST_CDD1 includes solution of elliptic PDE and MATLAB code to plot output
 * 
 * This file contains several useful domain decompositions
 */
#define _VALIDITY_OF_QT_INDICES 0
#define _ROUTINES_FROM_QTBASIS 0
#define _COMPARE_QT_a_AND_T_a 0
#define _SYMMETRY_OF_CACHED_A 0
#define _CALLS_TO_CQTPROBLEM_f 0
#define _TEST_NORM_A_AND_NORM_A_INV 0
#define _CHECK_CACHED_A_BY_HAND 0
#define _TEST_ADDBALL 0
#define _TEST_CDD2 0
#define _COMPARE_QT_RHS_AND_T_RHS 0
#define _TEST_CDD1 1


/* implemented Decompositions:
 * 0 :: 1 patch, 2d or 3d
 * 1,2,7,8 :: 2 patches, homogeneous BC
 * 1: extension to the right
 * 2: extension to the south
 * 7: extension to the north
 * 8: extension to the left
 * 3,4 :: L-shaped domain, homogeneous BC
 * 5:: L-shaped. corner domain is extended in 2 directions. inhomogeneous BC at reentrant corner. L-shaped supports are not implemented in compute_sum_over_patches !!
 * 6: big cube decomposed into 4 subcubes. Homogeneous BC
 * 9: Donut; 8 domains
 * 10: 3 domains in a line. middle is extended left and right
 * 11: Slit domain (0,2)^2 \ {1}x(1,2)  (4 sub-domains)
 * 12: Big square in 3d (8 sub-squares)
 * 13:  (0,2)^3 \ (1,2)^3 = Fichera corner domain (7 sub-squares)
 */
#define _DECOMP 0 // 2D: 3, 11, 3D: 13
//#define _NOP 1
#define _MAXOFFSET 1 // maxlevel == 1-norm of minlevel pf patch 0 + _MAXOFFSET
#define _RHS_NUMBER 1
#define _DIMENSION 2
#define _PRIMAL_ORDER_D  2
#define _DUAL_ORDER_DT 2
#define _COMPUTE_RHS 1 
#define _ONEDIMHAARCOUNT 1
//#define _RHS_FILENAME "f_prec_unso_rhs5_dec0_off1_3_3"

#define _APPLY_TENSOR_DEBUGMODE 0
//CLEANUP
#define M_ 0

// display j0, first/last wavelets, and \N -> levels
#define _DISPLAY_INFOS_ABOUT_BASIS 0
        
// more output for the cached problem (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// switch between isotropic and anisotropic Wavelets (in cdd1.h)
#define _WAVELETTL_USE_TBASIS 1
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0


#include <iostream>
#include <interval/p_basis.h>
#include <general_domain/qtbasis.h>
#include <time.h>
#include <galerkin/cached_qtproblem.h>
#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <cube/tbasis_evaluate.h>
#include <galerkin/cached_tproblem.h>
#include <galerkin/galerkin_utils.h>
#include <adaptive/apply.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <io/vector_io.h>

#if 0
#include <map>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <interval/i_index.h>
//#include <interval/ds_basis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <geometry/sampled_mapping.h>
#include <cube/tbasis_indexplot.h>

/*
 * for comparison with the isotropic case
 */
//#include <cube/cube_basis.h>
//#include <galerkin/cube_equation.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>

// incompatible with tbasis_indexplot:
//#include <cube/cube_indexplot.h>

#endif


using namespace std;
using namespace WaveletTL;


template <unsigned int DIM, unsigned int N>
class TestRHS
  : public Function<DIM,double>
{
public:
    virtual ~TestRHS() {};
    double value(const Point<DIM>& p, const unsigned int component = 0) const 
    {
        switch(N) 
        {
            case 0:
                return 0;
                break;
            case 1:
                return 1.0;
                break;
            case 2:
                return
                    (200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
                    * exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
                break;
            case 3:
                return
                    4*(1-p[0])*p[1]*p[1]*(1-p[1])
                    - 2*p[0]*p[1]*p[1]*(1-p[1])
                    - 2*p[0]*(1-p[0])*(1-p[0])*(1-p[1])
                    + 4*p[0]*(1-p[0])*(1-p[0])*p[1];
                break;
            case 4:
                return 10*sin(2*M_PI*p[0])*sin(6*M_PI*p[1]);
                break;
            case 5:
                return 3.0;
                break;
            case 6:
                return ((p[0] <= 1) && (p[1] <= 1)) ? 3.0:0;
                break;
            case 7:            
                return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
                break;
            default:
                return 1.0;
                break;
        }
    }
    void vector_value(const Point<DIM>& p, Vector<double>& values) const {
        values[0] = value(p);
    }
};

int main()
{
    cout << "Testing cached_qtproblem" << endl;

//    Vector<double> temp_v1, temp_v2;
//    Vector<int> temp_v3, temp_v4;
//    
//    temp_v1.resize(4);
//    temp_v2.resize(8);
//    temp_v3.resize(8);
//    
//    cout << "1: " << temp_v1 << endl;
//    cout << "2: " << temp_v2 << endl;
//    cout << "3: " << temp_v3 << endl;
//    
//    for (unsigned int i=0; i<temp_v1.size(); ++i)
//    {
//        temp_v1[i] = 1;
//    }
//    for (unsigned int i=0; i<temp_v2.size(); ++i)
//    {
//        temp_v2[i] = 2;
//    }
//    for (unsigned int i=0; i<temp_v3.size(); ++i)
//    {
//        temp_v3[i] = 3;
//    }
//    
//
//    cout << "1: " << temp_v1 << endl;
//    cout << "2: " << temp_v2 << endl;
//    cout << "3: " << temp_v3 << endl;
//    
//    cout << "store temp_v1 to disk" << endl;
//    std::ofstream outputstream("dummy1.file", std::ofstream::binary);
//    writeVectorToFile(temp_v1,outputstream);
//    outputstream.close();
//  
//    for (unsigned int i=0; i<temp_v1.size(); ++i)
//    {
//        temp_v1[i] = 4;
//    }
//    cout << "load old temp_v1 from disk" << endl;
//    std::ifstream inputstream ("dummy1.file", std::ifstream::binary);
//    readVectorFromFile(temp_v1,inputstream);
//    cout << temp_v1 << endl;
//    
//    cout << "store temp_v2 to disk" << endl;
//    std::ofstream outputstream2("dummy2.file", std::ofstream::binary);
//    writeVectorToFile(temp_v2,outputstream2);
//    outputstream2.close();
//    
//    
//    cout << "load temp_v2 from disk and store it in temp_v1" << endl;
//    inputstream.open("dummy2.file", std::ifstream::binary);
//    readVectorFromFile(temp_v1,inputstream);
//    cout << temp_v1 << endl;
//    
//    cout << "store temp_v3 to disk" << endl;
//    std::ofstream outputstream3("dummy3.file", std::ofstream::binary);
//    writeVectorToFile(temp_v3,outputstream3);
//    outputstream3.close();
//    
//    cout << "load temp_v3 from disk and store it in temp_v4" << endl;
//    inputstream.open("dummy3.file", std::ifstream::binary);
//    cout << "pre loading: temp_v4 = " << temp_v4 << endl;
//    readVectorFromFile(temp_v4,inputstream);
//    cout << "post loading: temp_v4 = " << temp_v4 << endl;
//    
//    cout << "load temp_v3 from disk and store it in temp_v2" << endl;
//    inputstream.open("dummy3.file", std::ifstream::binary);
//    cout << "pre loading: temp_v2 = " << temp_v2 << endl;
//    readVectorFromFile(temp_v2,inputstream);
//    cout << "post loading: temp_v2 = " << temp_v2 << endl;
//    assert(false);
    
    // specify the tests that will be performed:
    const int decomposition (_DECOMP);
    int num_of_patches; // (_NOP);
    const int maxleveloffset (_MAXOFFSET);    
    const int DIM (_DIMENSION);
    
    //typedef DSBasis<3,5> Basis1d;
    typedef PBasis<_PRIMAL_ORDER_D,_DUAL_ORDER_DT> Basis1d;
    typedef Basis1d::Index Index1d;
    typedef QTBasis<Basis1d,DIM> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;
    //typedef TensorBasis<Basis1d,DIM> TBasis;
    //typedef TBasis::Index TIndex;
    //typedef TBasis::Support TSupport;

    // construct the basis
    Array1D<Point<DIM,int> > corners;
    Array1D<FixedArray1D<int,2*DIM> > neighbours;
    Array1D<FixedArray1D<bool,2*DIM> > bc_bool;
    Array1D<FixedArray1D<int,2*DIM> > bc_int;


    switch (decomposition) {
        case 0: // 1 domain
            num_of_patches = 1;
            break;
        case 1: // 2 domains, homogeneous BC, left patch is extended to the right
            num_of_patches = 2;
            break;
        case 2: // 2 domains, homogeneous BC, patch 1 is extended southwards to patch 0 (left in the y direction for DIM=2)
            num_of_patches = 2;
            break;
        case 7: // 2 domains, homogeneous BC, patch 0 is extended northwards to patch 1 (right in the y direction for DIM=2)
            num_of_patches = 2;
            break;
        case 8: // 2 domains, homogeneous BC, patch 1 is extended westwards to patch 0 (left in the x direction for DIM=2)
            num_of_patches = 2;
            break;
        case 3: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 1
            num_of_patches = 3;
            assert (DIM == 2);
            break;
        case 4: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 0
            num_of_patches = 3;
            assert (DIM == 2);
            break;
        case 5: // L-shaped domain, homogeneous BC everywhere except around the reentrant corner
            // patch 0 is extended to right to patch 1 and to the north to patch 2
            num_of_patches = 3;
            assert (DIM == 2);
            break;
        case 6: // Big square, decomposed into 4 sub-squares
            // patch 0 is extended to right to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            num_of_patches = 4;
            assert (DIM == 2);
            break;
        case 9: // Donut! Decomposed into 8 sub-squares
            // 6 5 4
            // 7 x 3
            // 0 1 2
            // 
            // patches 1,3,5,7 are extended to their neighbours
            // homogeneous BC wherever possible
            num_of_patches = 8;
            assert (DIM == 2);
            break;
        case 10: // Line
            // 0 1 2
            // 
            // patches 1 is extended left and right
            // homogeneous BC wherever possible
            num_of_patches = 3;
            assert (DIM == 2);
            break;
        case 11: // Slit domain (0,2)^2 \ {1}x(1,2)  (4 sub-domains)
            // patch 0 is extended to right to patch 1
            // patch 1 is extended nowhere
            // patch 2 is extended to the south to patch 1. But is is not extended to the patch left of it (number 3)
            // patch 3 is extended to the south to patch 0. But is is not extended to the patch right of it (number 2)
            num_of_patches = 4;
            assert (DIM == 2);
            break;
        case 12: // Big square in 3d (8 sub-squares)
            // patches 0,1,2,3 and 4,5,6,7 are extended like in decomposition 6 to form bases for (0,2)^2x (0,1) and (0,2)^2x(1,2), respectively.
            // patches 0,1,2,3 are extended up in z direction.
            // E.g. :
            // patch 0 is extended to the east to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            num_of_patches = 8;
            assert (DIM == 3);
            break;
        case 13: // (0,2)^3 \ (1,2)^3 = Fichera corner domain (7 sub-squares)
            // patches 0,1,2,3 are extended like in decomposition 6 to form bases for (0,2)^2x (0,1) and (0,2)^2x(1,2), respectively.
            // patches 4,5,6 are extended like in decomposition 3
            // patches 4,5,6 are extended down in z direction.
            num_of_patches = 7;
            break;
        default:
            cout << "main:: error! no decomposition specified!" << endl;
            abort();
            break;
    }
    
    cout << "Decomposition = " << decomposition << "; maxleveloffset = " << maxleveloffset << "; Number of patches = " << num_of_patches << endl;
    
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
            
            bc_bool[0][0] = false;
            bc_bool[0][1] = false;
            bc_bool[0][2] = false;
            bc_bool[0][3] = false;
            
//            bc_int[0][1] = 0;
//            bc_int[0][3] = 0;
            break;
        case 1: // 2 domains, homogeneous BC, left patch is extended to the right
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[1][0] = 1;
            }
            else if (DIM == 2)
            {
                corners[0][0] = 0;
                corners[0][1] = 0;
                corners[1][0] = 1;
                corners[1][1] = 0;
            }
            neighbours[0][1] = 1; // patch 0, x-right-neighbour = patch 1
            neighbours[1][0] = 0; // patch 1, x-left-neighbour = patch 0
            bc_bool[0][1] = false;
            bc_int[0][1] = 0;
            break;
        case 2: // 2 domains, homogeneous BC, patch 1 is extended southwards to patch 0 (left in the y direction for DIM=2)
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (DIM == 2)
            {
                corners[0][0] = 0;
                corners[0][1] = 0;
                corners[1][0] = 0;
                corners[1][1] = 1;
            }
            neighbours[0][3] = 1; // patch 0 has a neighbour in the south
            neighbours[1][2] = 0; // patch 1 has a neighbour in the north
            bc_bool[1][2] = false;
            bc_int[1][2] = 0;
            break;
        case 7: // 2 domains, homogeneous BC, patch 0 is extended northwards to patch 1 (right in the y direction for DIM=2)
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (DIM == 2)
            {
                corners[0][0] = 0;
                corners[0][1] = 0;
                corners[1][0] = 0;
                corners[1][1] = 1;
            }
            neighbours[0][3] = 1; // patch 0 has a neighbour in the south
            neighbours[1][2] = 0; // patch 1 has a neighbour in the north
            bc_bool[0][3] = false;
            bc_int[0][3] = 0;
            break;
        case 8: // 2 domains, homogeneous BC, patch 1 is extended westwards to patch 0 (left in the x direction for DIM=2)
            assert (num_of_patches == 2);
            if (DIM == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (DIM == 2)
            {
                corners[0][0] = 0;
                corners[0][1] = 0;
                corners[1][0] = 1;
                corners[1][1] = 0;
            }
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            bc_bool[1][0] = false;
            bc_int[1][0] = 0;
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
            corners[0][0] = 0;
            corners[0][1] = 1;
            corners[1][0] = 0;
            corners[1][1] = 0;
            corners[2][0] = 1;
            corners[2][1] = 1;
            
            neighbours[0][2] = 1; // patch 0 has a neighbour in the south
            neighbours[0][1] = 2; // patch 0 has a neighbour in the east
            neighbours[1][3] = 0; // patch 1 has a neighbour in the north
            neighbours[2][0] = 0; // patch 2 has a neighbour in the west
            bc_bool[0][2] = false;
            bc_bool[2][0] = false;
            bc_int[0][2] = 0;
            bc_int[2][0] = 0;
            break;
        case 5: // L-shaped domain, homogeneous BC everywhere except around the reentrant corner
            // patch 1 is extended to right to patch 2 and to the north to patch 0
            assert (num_of_patches == 3);
            assert (DIM == 2);
            corners[0][0] = 19;
            corners[0][1] = 3;
            corners[1][0] = 19;
            corners[1][1] = 2;
            corners[2][0] = 20;
            corners[2][1] = 2;
            
            neighbours[0][2] = 0; // patch 0 has a neighbour in the south
            neighbours[1][1] = 1; // patch 1 has a neighbour in the east
            neighbours[1][3] = 2; // patch 2 has a neighbour in the north
            neighbours[2][0] = 0; // patch 2 has a neighbour in the west
            bc_bool[0][1] = false;
            bc_bool[1][1] = false;
            bc_bool[1][3] = false;
            bc_bool[2][3] = false;
            bc_int[0][1] = 0;
            bc_int[1][1] = 0;
            bc_int[1][3] = 0;
            bc_int[2][3] = 0;
            
            
            break;
        case 6: // Big square, decomposed into 4 sub-squares
            // patch 0 is extended to right to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            assert (num_of_patches == 4);
            assert (DIM == 2);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[2][0] = 1;
            corners[2][1] = 1;
            corners[3][0] = 0;
            corners[3][1] = 1;
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[0][3] = 3; // patch 0 has a neighbour in the north
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][3] = 2; // patch 1 has a neighbour in the north
            neighbours[2][0] = 3; // patch 2 has a neighbour in the west
            neighbours[2][2] = 1; // patch 2 has a neighbour in the south
            neighbours[3][1] = 2; // patch 3 has a neighbour in the east
            neighbours[3][2] = 0; // patch 3 has a neighbour in the south
            
            bc_bool[0][1] = false;
            bc_bool[0][3] = false;
            bc_bool[1][3] = false;
            bc_bool[2][0] = false;
            bc_int[0][1] = 0;
            bc_int[0][3] = 0;
            bc_int[1][3] = 0;
            bc_int[2][0] = 0;
            break;
        case 9: // Donut! Decomposed into 8 sub-squares
            // 6 5 4
            // 7 x 3
            // 0 1 2
            // 
            // patches 1,3,5,7 are extended to their neighbours
            // homogeneous BC wherever possible
            assert (num_of_patches == 8);
            assert (DIM == 2);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[2][0] = 2;
            corners[2][1] = 0;
            corners[3][0] = 2;
            corners[3][1] = 1;
            corners[4][0] = 2;
            corners[4][1] = 2;
            corners[5][0] = 1;
            corners[5][1] = 2;
            corners[6][0] = 0;
            corners[6][1] = 2; // 6 5 4
            corners[7][0] = 0; // 7 x 3
            corners[7][1] = 1; // 0 1 2
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[0][3] = 7; // patch 0 has a neighbour in the north
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][1] = 2; // patch 1 has a neighbour in the east
            neighbours[2][0] = 1; // patch 2 has a neighbour in the west
            neighbours[2][3] = 3; // patch 2 has a neighbour in the north
            neighbours[3][2] = 2; // patch 3 has a neighbour in the south
            neighbours[3][3] = 4; // patch 3 has a neighbour in the north
            neighbours[4][0] = 5; // patch 4 has a neighbour in the west
            neighbours[4][2] = 3; // patch 4 has a neighbour in the south
            neighbours[5][0] = 6; // patch 5 has a neighbour in the west
            neighbours[5][1] = 4; // patch 5 has a neighbour in the east
            neighbours[6][1] = 5; // patch 6 has a neighbour in the east
            neighbours[6][2] = 7; // patch 6 has a neighbour in the south
            neighbours[7][2] = 0; // patch 7 has a neighbour in the south
            neighbours[7][3] = 6; // patch 7 has a neighbour in the north
            
            bc_bool[1][0] = false; // 6 5 4
            bc_bool[1][1] = false; // 7 x 3
            bc_bool[3][2] = false; // 0 1 2
            bc_bool[3][3] = false; 
            bc_bool[5][0] = false; 
            bc_bool[5][1] = false;
            bc_bool[7][2] = false;
            bc_bool[7][3] = false;
            
            bc_int[1][0] = 0; // 6 5 4
            bc_int[1][1] = 0; // 7 x 3
            bc_int[3][2] = 0; // 0 1 2
            bc_int[3][3] = 0; 
            bc_int[5][0] = 0; 
            bc_int[5][1] = 0;
            bc_int[7][2] = 0;
            bc_int[7][3] = 0;
            break;
        case 10: // Line
            // 0 1 2
            // 
            // patches 1 is extended left and right
            // homogeneous BC wherever possible
            assert (num_of_patches == 3);
            assert (DIM == 2);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[2][0] = 2;
            corners[2][1] = 0;
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][1] = 2; // patch 1 has a neighbour in the east
            neighbours[2][0] = 1; // patch 2 has a neighbour in the west
            
            bc_bool[1][0] = false;
            bc_bool[1][1] = false;
            
            bc_int[1][0] = 0;
            bc_int[1][1] = 0;
            break;
        case 11: // Slit domain (0,2)^2 \ {1}x(1,2)  (4 sub-domains)
            // patch 0 is extended to right to patch 1
            // patch 1 is extended nowhere
            // patch 2 is extended to the south to patch 1. But is is not extended to the patch left of it (number 3)
            // patch 3 is extended to the south to patch 0. But is is not extended to the patch right of it (number 2)
            assert (num_of_patches == 4);
            assert (DIM == 2);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[2][0] = 1;
            corners[2][1] = 1;
            corners[3][0] = 0;
            corners[3][1] = 1;
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[0][3] = 3; // patch 0 has a neighbour in the north
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][3] = 2; // patch 1 has a neighbour in the north
            neighbours[2][2] = 1; // patch 2 has a neighbour in the south
            neighbours[3][2] = 0; // patch 3 has a neighbour in the south
            
            bc_bool[0][1] = false;
            bc_bool[2][2] = false;
            bc_bool[3][2] = false;
            bc_int[0][1] = 0;
            bc_int[2][2] = 0;
            bc_int[3][2] = 0;
            break;
        case 12: // Big square in 3d (8 sub-squares)
            // patches 0,1,2,3 and 4,5,6,7 are extended like in decomposition 6 to form bases for (0,2)^2x (0,1) and (0,2)^2x(1,2), respectively.
            // patches 0,1,2,3 are extended up in z direction.
            // E.g. :
            // patch 0 is extended to the east to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            assert (num_of_patches == 8);
            assert (DIM == 3);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[0][2] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[1][2] = 0;
            corners[2][0] = 1;
            corners[2][1] = 1;
            corners[2][2] = 0;
            corners[3][0] = 0;
            corners[3][1] = 1;
            corners[3][2] = 0;
            corners[4][0] = 0;
            corners[4][1] = 0;
            corners[4][2] = 1;
            corners[5][0] = 1;
            corners[5][1] = 0;
            corners[5][2] = 1;
            corners[6][0] = 1;
            corners[6][1] = 1;
            corners[6][2] = 1;
            corners[7][0] = 0;
            corners[7][1] = 1;
            corners[7][2] = 1;
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[0][3] = 3; // patch 0 has a neighbour in the north
            neighbours[0][5] = 4; // patch 0 has a neighbour above
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][3] = 2; // patch 1 has a neighbour in the north
            neighbours[1][5] = 5; // patch 1 has a neighbour above
            neighbours[2][0] = 3; // patch 2 has a neighbour in the west
            neighbours[2][2] = 1; // patch 2 has a neighbour in the south
            neighbours[2][5] = 6; // patch 2 has a neighbour above
            neighbours[3][1] = 2; // patch 3 has a neighbour in the east
            neighbours[3][2] = 0; // patch 3 has a neighbour in the south
            neighbours[3][5] = 7; // patch 3 has a neighbour above
            
            neighbours[4][1] = 5; // patch 4 has a neighbour in the east
            neighbours[4][3] = 7; // patch 4 has a neighbour in the north
            neighbours[4][4] = 0; // patch 4 has a neighbour below
            neighbours[5][0] = 4; // patch 5 has a neighbour in the west
            neighbours[5][3] = 6; // patch 5 has a neighbour in the north
            neighbours[5][4] = 1; // patch 5 has a neighbour below
            neighbours[6][0] = 7; // patch 6 has a neighbour in the west
            neighbours[6][2] = 5; // patch 6 has a neighbour in the south
            neighbours[6][4] = 2; // patch 6 has a neighbour below
            neighbours[7][1] = 6; // patch 7 has a neighbour in the east
            neighbours[7][2] = 4; // patch 7 has a neighbour in the south
            neighbours[7][4] = 3; // patch 7 has a neighbour below
            
            bc_bool[0][1] = false;
            bc_bool[0][3] = false;
            bc_bool[0][5] = false;
            bc_bool[1][3] = false;
            bc_bool[1][5] = false;
            bc_bool[2][0] = false;
            bc_bool[2][5] = false;
            bc_bool[3][5] = false;
            bc_bool[4][1] = false;
            bc_bool[4][3] = false;
            bc_bool[5][3] = false;
            bc_bool[6][0] = false;
            
            bc_int[0][1] = 0;
            bc_int[0][3] = 0;
            bc_int[0][5] = 0;
            bc_int[1][3] = 0;
            bc_int[1][5] = 0;
            bc_int[2][0] = 0;
            bc_int[2][5] = 0;
            bc_int[3][5] = 0;
            bc_int[4][1] = 0;
            bc_int[4][3] = 0;
            bc_int[5][3] = 0;
            bc_int[6][0] = 0;
            break;
        case 13: // (0,2)^3 \ (1,2)^3 = Fichera corner domain (7 sub-squares)
            // patches 0,1,2,3 are extended like in decomposition 6 to form bases for (0,2)^2x (0,1) and (0,2)^2x(1,2), respectively.
            // patches 4,5,6 are extended like in decomposition 3
            // patches 4,5,6 are extended down in z direction.
            assert (num_of_patches == 7);
            assert (DIM == 3);
            corners[0][0] = 0;
            corners[0][1] = 0;
            corners[0][2] = 0;
            corners[1][0] = 1;
            corners[1][1] = 0;
            corners[1][2] = 0;
            corners[2][0] = 1;
            corners[2][1] = 1;
            corners[2][2] = 0;
            corners[3][0] = 0;
            corners[3][1] = 1;
            corners[3][2] = 0;
            corners[4][0] = 0;
            corners[4][1] = 0;
            corners[4][2] = 1;
            corners[5][0] = 1;
            corners[5][1] = 0;
            corners[5][2] = 1;
            corners[6][0] = 0;
            corners[6][1] = 1;
            corners[6][2] = 1;
            
            neighbours[0][1] = 1; // patch 0 has a neighbour in the east
            neighbours[0][3] = 3; // patch 0 has a neighbour in the north
            neighbours[0][5] = 4; // patch 0 has a neighbour above
            neighbours[1][0] = 0; // patch 1 has a neighbour in the west
            neighbours[1][3] = 2; // patch 1 has a neighbour in the north
            neighbours[1][5] = 5; // patch 1 has a neighbour above
            neighbours[2][0] = 3; // patch 2 has a neighbour in the west
            neighbours[2][2] = 1; // patch 2 has a neighbour in the south
            // patch 2 has no neighbor above
            neighbours[3][1] = 2; // patch 3 has a neighbour in the east
            neighbours[3][2] = 0; // patch 3 has a neighbour in the south
            neighbours[3][5] = 6; // patch 3 has a neighbour above
            neighbours[4][1] = 5; // patch 4 has a neighbour in the east
            neighbours[4][3] = 6; // patch 4 has a neighbour in the north
            neighbours[4][4] = 0; // patch 4 has a neighbour below
            neighbours[5][0] = 4; // patch 5 has a neighbour in the west
            neighbours[5][4] = 1; // patch 5 has a neighbour below
            neighbours[6][2] = 4; // patch 6 has a neighbour in the south
            neighbours[6][4] = 3; // patch 6 has a neighbour below
            
            bc_bool[0][1] = false;
            bc_bool[0][3] = false;
            bc_bool[1][3] = false;
            bc_bool[2][0] = false;
            bc_bool[4][4] = false;
            bc_bool[5][0] = false;
            bc_bool[5][4] = false;
            bc_bool[6][2] = false;
            bc_bool[6][4] = false;
            
            bc_int[0][1] = 0;
            bc_int[0][3] = 0;
            bc_int[1][3] = 0;
            bc_int[2][0] = 0;
            bc_int[4][4] = 0;
            bc_int[5][0] = 0;
            bc_int[5][4] = 0;
            bc_int[6][2] = 0;
            bc_int[6][4] = 0;
            break;
        default:
            cout << "main:: error! no decomposition specified!" << endl;
            abort();
            break;
    }
    
    Basis qtbasis(corners, neighbours,bc_bool);
    cout << multi_degree(qtbasis.j0()[0]) << " " << maxleveloffset << endl;
    qtbasis.set_jmax(multi_degree(qtbasis.j0()[0])+maxleveloffset);
    
#if _DISPLAY_INFOS_ABOUT_BASIS
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

    for (unsigned int i=0; i<DIM; ++i)
    {
        for (unsigned int p=0; p<num_of_patches; p++)
        {
            cout << " qtbasis.j0_[" << p<< "][" << i << "] = " << qtbasis.j0_[p][i] << endl;
        }
    }
    cout << "first_wavelets = " << qtbasis.first_wavelets_ << endl;
    cout << "last_wavelets = " << qtbasis.last_wavelets_ << endl;
    cout << "level_to_num = " << qtbasis.level_to_num_ << endl;
    cout << "num_to_level = " << endl;
    for (unsigned int n=0; n< qtbasis.num_to_level_.size(); ++n)
    {
        cout << n << " : " << "(" << qtbasis.num_to_level_[n].first << ", " << qtbasis.num_to_level_[n].second << ")" << endl;
    }
#endif
    
    //setup_full_collection();
    //cout << "full_collection_ = " << qtbasis.full_collection_ << endl;

    // test methods that were written for cached_qtproblem, but were added added to qtbasis.h /.cpp
    Index temp_nu, temp_mu;
    int mui_basisnum, nui_basisnum, kmin, kmax, temp_i1, temp_i2, nui_k1, nui_k2;
    int kmingen_i, kmaxgen_i, kminwav_i, kmaxwav_i;
    bool temp_b, exit, reflected, nu_and_mu_intersect, gen_intersection_i, wav_intersection_i;
    
    MultiIndex<unsigned int, DIM> intinfo, intinfo2;
    MultiIndex<int, DIM> kmingen, kmaxgen,kminwav, kmaxwav;
    MultiIndex<bool, DIM> gen_intersection, wav_intersection;
    unsigned int geometry_type;
    int  centerpatchnumber;
    MultiIndex<bool, DIM> orientation;
    
    clock_t tstart, tend;
    double time;
 
    
#if _VALIDITY_OF_QT_INDICES
    cout << "testing vailidity of qtbasis indices" << endl;
    for (int n=0; n< qtbasis.degrees_of_freedom(); n++)
    {
        temp_nu= qtbasis.get_wavelet(n);
        //cout << "n = " << n << "; temp_nu = " << temp_nu << endl;
        for (unsigned int i=0; i< DIM;++i)
        {
            if (temp_nu.e()[i] == 0)
            {
                assert(temp_nu.k()[i] >= qtbasis.get_bases(temp_nu.p(),i)->DeltaLmin());
                assert(temp_nu.k()[i] <= qtbasis.get_bases(temp_nu.p(),i)->DeltaRmax(temp_nu.j()[i]));
            }
            else
            {
                assert(temp_nu.k()[i] >= qtbasis.get_bases(temp_nu.p(),i)->Nablamin());
                assert(temp_nu.k()[i] <= qtbasis.get_bases(temp_nu.p(),i)->Nablamax(temp_nu.j()[i]));
            }
        }
    }
    cout << "done" << endl;
#else
    cout << "skipping initial index test" << endl;
#endif
    
/*   // a test for the pbasis

    //Basis1d bas01(0,1);
    
    //Index1d lambda1d(3,0,0,&bas01);
    cout << qtbasis.get_bases_infact()[0]->Deltasize(3) << endl;
    cout << qtbasis.get_bases_infact()[2]->Deltasize(3) << endl;
    //cout << qtbasis.get_bases_infact()[2]->Deltasize(3) << endl;
    //cout << qtbasis.get_bases_infact()[3]->Deltasize(3) << endl;
    
    qtbasis.bases_infact_[0]->set_jmax(5);
    Index1d temp_ind(qtbasis.get_bases_infact()[0]->first_generator(3));
    for (unsigned int n=0; n<20;++n)
    {
        cout << "n = " << n << "; temp_ind  = " << temp_ind 
                << "; .number() = " << temp_ind.number() << endl;
        cout << "get_wavelet(" << n << ")   = " << *qtbasis.get_bases_infact()[0]->get_wavelet(n)
                << "; .number() = " << qtbasis.get_bases_infact()[0]->get_wavelet(n)->number() << endl;
        Index1d temp_ind2(temp_ind.j(), temp_ind.e(), temp_ind.k(), qtbasis.get_bases_infact()[0]);
        cout << "Index1d(" << temp_ind.j() << ", " << temp_ind.e() << ", " << temp_ind.k() << ") = " << temp_ind2 
                << "; .number() = " << temp_ind2.number() << endl << endl;
        ++temp_ind;
    }
    cout << "done" << endl;
    //abort();
  */  
    
#if _ROUTINES_FROM_QTBASIS
#define _n -1
#define _m 278
    cout << "begin testing routines from qtbasis" << endl;
    cout << "n =";
    for (int n=0; n< qtbasis.degrees_of_freedom(); n++)
    {
        //cout << "n = " << n << endl;
        cout << " " << n;
        temp_nu= qtbasis.get_wavelet(n);
        for (int m = 0; m < qtbasis.degrees_of_freedom(); m++)
        {
            nu_and_mu_intersect = true;
            temp_mu= qtbasis.get_wavelet(m);
            // cout << "m = " << m << endl;
            //cout << "temp_nu = " << temp_nu << endl;
            //cout << "temp_mu = " << temp_mu << endl;
            //cout << "n = " << n << ";m = " << m << "; temp_nu = " << temp_nu << "; temp_mu = " << temp_mu << endl;
            //if (false)
            if ((n == _n) && (m == _m))
            {
                cout << "temp_nu = " << temp_nu << endl;
                cout << "temp_mu = " << temp_mu << endl;
            }
            temp_i2 = qtbasis.get_levelnum(temp_mu.j(), temp_mu.p());
            
            nu_and_mu_intersect = qtbasis.intersecting_wavelets(n, 
                    temp_i2,
                    kmingen,
                    kmaxgen,
                    kminwav,
                    kmaxwav,
                    gen_intersection,
                    wav_intersection,
                    intinfo2);
            // cache ok?
            MultiIndex<unsigned int, DIM> intinfo3;
            MultiIndex<int, DIM> kmingen2, kmaxgen2,kminwav2, kmaxwav2;
            MultiIndex<bool, DIM> gen_intersection2, wav_intersection2;
            temp_b = qtbasis.intersecting_wavelets(n, 
                    temp_i2,
                    kmingen2,
                    kmaxgen2,
                    kminwav2,
                    kmaxwav2,
                    gen_intersection2,
                    wav_intersection2,
                    intinfo3);
            if ((n == _n) && (m == _m))
            {
                cout << "kmingen = " << kmingen << endl;
                cout << "kmaxgen = " << kmaxgen << endl;
                cout << "kminwav = " << kminwav << endl;
                cout << "kmaxwav = " << kmaxwav << endl;
                cout << "gen_intersection = " << gen_intersection << endl;
                cout << "wav_intersection = " << wav_intersection << endl;
                cout << "intinfo2 = " << intinfo2 << endl;
                cout << "output = " << nu_and_mu_intersect << endl;
                
                cout << "kmingen2 = " << kmingen2 << endl;
                cout << "kmaxgen2 = " << kmaxgen2 << endl;
                cout << "kminwav2 = " << kminwav2 << endl;
                cout << "kmaxwav2 = " << kmaxwav2 << endl;
                cout << "gen_intersection2 = " << gen_intersection2 << endl;
                cout << "wav_intersection2 = " << wav_intersection2 << endl;
                cout << "intinfo3 = " << intinfo3 << endl;
                cout << "output " << temp_b << endl;
                
            }
            assert (temp_b == nu_and_mu_intersect);
            if (temp_b)
            {
                assert (gen_intersection == gen_intersection2);
                assert (wav_intersection == wav_intersection2);
                for (unsigned int j=0; j<DIM; ++j)
                {
                    if (gen_intersection[j])
                    {
                        assert (wav_intersection[j]);
                        assert (kmingen[j] == kmingen2[j]);
                        assert (kmaxgen[j] == kmaxgen2[j]);
                    }
                    if (wav_intersection[j])
                    {
                        assert (kminwav[j] == kminwav2[j]);
                        assert (kmaxwav[j] == kmaxwav2[j]);
                    }
                }
                assert (intinfo2 == intinfo3);
            }
            
            
            temp_mu= qtbasis.get_wavelet(m);
            
            temp_b = qtbasis.get_LMR_info(n, temp_mu.j(), temp_mu.p(), intinfo);

            assert (!nu_and_mu_intersect || temp_b ); // temp_b == false => temp_b2 == false
            if (!temp_b)
            {
                continue;
            }
            
            for (unsigned int i = 0; i<DIM; i++)
            {
                assert (intinfo[i] == intinfo2[i]);
            }
            
            for (unsigned int i = 0; i<DIM; ++i)
            {
                nui_basisnum = (((qtbasis.get_bc()[temp_nu.p()][2*i])?0:2) + ((qtbasis.get_bc()[temp_nu.p()][2*i+1])?0:1));
                mui_basisnum = (((qtbasis.get_bc()[temp_mu.p()][2*i])?0:2) + ((qtbasis.get_bc()[temp_mu.p()][2*i+1])?0:1));
                
                //muss vor aufruf von get_cached_onedim_intersections geprüft werden (und nur falls ergebnis von get_LMR_info == true)
                
                //if (false)
                //if (n==272)
                if ((n == _n) && (m == _m))
                {
                    cout << "onedimtest" << endl;
                }
                qtbasis.get_onedim_intersections(intinfo[i],
                        temp_nu.j()[i],
                        temp_nu.e()[i],
                        temp_nu.k()[i],
                        nui_basisnum,
                        (temp_mu.j()[i] == qtbasis.get_j0(temp_mu.p(),i)), // mintype_mui
                        temp_mu.j()[i],
                        mui_basisnum,
                        kmingen_i,
                        kmaxgen_i,
                        kminwav_i,
                        kmaxwav_i,
                        gen_intersection_i,
                        wav_intersection_i);
                temp_b = temp_b && (wav_intersection_i || gen_intersection_i);
                
                int kmingen_i2, kmaxgen_i2, kminwav_i2, kmaxwav_i2;
                bool gen_intersection_i2, wav_intersection_i2;
                
                qtbasis.get_onedim_intersections(intinfo[i],
                        temp_nu.j()[i],
                        temp_nu.e()[i],
                        temp_nu.k()[i],
                        nui_basisnum,
                        (temp_mu.j()[i] == qtbasis.get_j0(temp_mu.p(),i)), // mintype_mui
                        temp_mu.j()[i],
                        mui_basisnum,
                        kmingen_i2,
                        kmaxgen_i2,
                        kminwav_i2,
                        kmaxwav_i2,
                        gen_intersection_i2,
                        wav_intersection_i2);
                
                assert (gen_intersection_i == gen_intersection_i2);
                assert (wav_intersection_i == wav_intersection_i2);
                if (gen_intersection_i)
                {
                    assert (wav_intersection_i);
                    assert (kmingen_i == kmingen_i2);
                    assert (kmaxgen_i == kmaxgen_i2);
                }
                if (wav_intersection_i)
                {
                    assert (kminwav_i == kminwav_i2);
                    assert (kmaxwav_i == kmaxwav_i2);
                }

                if (!temp_b)
                {
                    assert (!nu_and_mu_intersect);
                    //nu_and_mu_intersect = false;
                    break;
                }
                
                //if (false)
                //if (n==272)
                if ((n == _n) && (m == _m))
                {
                    cout << "intinfo[i] = " << intinfo[i] << endl;
                    cout << "temp_nu.j()[i] = " << temp_nu.j()[i] << endl;
                    cout << "temp_nu.e()[i] = " << temp_nu.e()[i] << endl;
                    cout << "temp_nu.k()[i] = " << temp_nu.k()[i] << endl;
                    cout << "nui_basisnum = " << nui_basisnum << endl;
                    cout << "mintype_mui = (temp_mu.j()[i] == qtbasis.get_j0(temp_mu.p(),i)) = " << (temp_mu.j()[i] == qtbasis.get_j0(temp_mu.p(),i)) << endl;
                    cout << "temp_mu.j()[i] = " << temp_mu.j()[i] << endl;
                    cout << "mui_basisnum = " << mui_basisnum << endl;
                    cout << "kmingen_i = " << kmingen_i << endl;
                    cout << "kmaxgen_i = " << kmaxgen_i << endl;
                    cout << "kminwav_i = " << kminwav_i << endl;
                    cout << "kmaxwav_i = " << kmaxwav_i << endl;
                    cout << "gen_intersection_i = " << gen_intersection_i << endl;
                    cout << "wav_intersection_i = " << wav_intersection_i << endl;
                    cout << "kmingen[i] = " << kmingen[i] << endl;
                    cout << "kmaxgen[i] = " << kmaxgen[i] << endl;
                    cout << "kminwav[i] = " << kminwav[i] << endl;
                    cout << "kmaxwav[i] = " << kmaxwav[i] << endl;
                    cout << "gen_intersection[i] = " << gen_intersection[i] << endl;
                    cout << "wav_intersection[i] = " << wav_intersection[i] << endl;
                    qtbasis.get_bases(0,0)->support(temp_nu.j()[i], temp_nu.e()[i], temp_nu.k()[i], temp_i1, temp_i2);
                    cout << "support nu_i = 2^{-(temp_nu.j()[i]+ temp_nu.e()[i])}[" << temp_i1 << ", " << temp_i2 << "]" << endl;
                    Support temp_supp;
                    qtbasis.support_local(temp_nu, temp_supp);
                    cout << "support temp_nu: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                            "]; a = [" << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                            "]; b = [" << temp_supp.b[0] << ", " << temp_supp.b[1] << "]" << endl;
                    
                    Index temp_ind(qtbasis.get_wavelet(m));
                    qtbasis.support_local(temp_ind, temp_supp);
                    cout << "support temp_mu: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                            "]; a = [" << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                            "]; b = [" << temp_supp.b[0] << ", " << temp_supp.b[1] << "]" << endl;
                    /*
                    temp_ind = qtbasis.get_wavelet(m+1);
                    qtbasis.support_local(temp_ind, temp_supp);
                    cout << "support temp_mu+1: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                            "]; a = " << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                            "]; b = " << temp_supp.b[0] << ", " << temp_supp.b[1] << endl;
                    
                    temp_ind = qtbasis.get_wavelet(m+2);
                    qtbasis.support_local(temp_ind, temp_supp);
                    cout << "support temp_mu+2: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                            "]; a = " << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                            "]; b = " << temp_supp.b[0] << ", " << temp_supp.b[1] << endl;
                    
                    temp_ind = qtbasis.get_wavelet(m+3);
                    qtbasis.support_local(temp_ind, temp_supp);
                    cout << "support temp_mu+3: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                            "]; a = " << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                            "]; b = " << temp_supp.b[0] << ", " << temp_supp.b[1] << endl;
                    */
                     
                }
                
                if (temp_mu.j()[i] == qtbasis.get_j0(temp_mu.p(),i))
                {
                    assert (gen_intersection[i] == gen_intersection_i);
                    if (gen_intersection_i)
                    {
                        assert (kmingen[i] == kmingen_i);
                        assert (kmaxgen[i] == kmaxgen_i);
                    }
                }
                assert (wav_intersection[i] == wav_intersection_i);
                if (wav_intersection_i)
                {
                    assert (kminwav[i] == kminwav_i);
                    assert (kmaxwav[i] == kmaxwav_i);
                }
                
                if (i == DIM-1 )
                {
                    assert (temp_b == nu_and_mu_intersect);
                }
                else
                {
                    assert (temp_b >= nu_and_mu_intersect);
                }
            }
            
            for (unsigned int i=0; i<DIM; ++i)
            {
                nui_basisnum = (((qtbasis.get_bc()[temp_nu.p()][2*i])?0:2) + ((qtbasis.get_bc()[temp_nu.p()][2*i+1])?0:1));
                mui_basisnum = (((qtbasis.get_bc()[temp_mu.p()][2*i])?0:2) + ((qtbasis.get_bc()[temp_mu.p()][2*i+1])?0:1));
                
                // primitive test whether nu is reflected:
                //if (false)
                if ( (n == _n) && (m == _m))
                {
                    //IBASIS::Index lambda(lambda_j,lambda_e,lambda_k,bases_infact_[lambda_basisnum]);
                    //IBASIS::Index lambda(3,1,0,bases_infact_[0]); -> NR = 8 ???
                    cout << "temp_nu = " << temp_nu << endl;
                    cout << "temp_mu = " << temp_mu << endl;
                }
                reflected = (temp_mu.p() != temp_nu.p()) && (qtbasis.get_corners(temp_nu.p(),i) != qtbasis.get_corners(temp_mu.p(),i) );
                if ( (temp_nu.j()[i] == 3) && ((temp_nu.e()[i] == 0) && ( (temp_nu.k()[i] == 8) && (nui_basisnum == 0)  )))
                {
                    cout << "error!!!" << endl;
                    cout << "n = " << n << "; temp_nu = " << temp_nu << endl;
                    abort();
                }
                    
                qtbasis.get_cached_onedim_intersections(nui_basisnum,
                        temp_nu.j()[i],
                        temp_nu.e()[i],
                        temp_nu.k()[i],
                        reflected,
                        mui_basisnum,
                        temp_mu.j()[i],
                        (temp_mu.e()[i] == 0),
                        kmin,
                        kmax);
                
                qtbasis.get_cached_onedim_intersections(nui_basisnum,
                        temp_nu.j()[i],
                        temp_nu.e()[i],
                        temp_nu.k()[i],
                        reflected,
                        mui_basisnum,
                        temp_mu.j()[i],
                        (temp_mu.e()[i] == 0),
                        temp_i1,
                        temp_i2);

                assert (kmin == temp_i1);
                assert (kmax == temp_i2);
                //if (false)
                if ( (n == _n) && (m == _m))
                {
                    cout << "nui_basisnum = " << nui_basisnum << endl;
                    cout << "temp_nu.j()[i] = " << temp_nu.j()[i] << endl;
                    cout << "temp_nu.e()[i] = " << temp_nu.e()[i] << endl;
                    cout << "temp_nu.k()[i] = " << temp_nu.k()[i] << endl;
                    cout << "reflected = " << reflected << endl;
                    cout << "mui_basisnum = " << mui_basisnum << endl;
                    cout << "temp_mu.j()[i] = " << temp_mu.j()[i] << endl;
                    cout << "(temp_mu.e()[i] == 0) = " << (temp_mu.e()[i] == 0) << endl;
                    cout << "kmin = " << kmin << endl;
                    cout << "kmax = " << kmax << endl;
                }
                
                qtbasis.get_bases_infact()[nui_basisnum]->support(temp_nu.j()[i],temp_nu.e()[i],temp_nu.k()[i],nui_k1,nui_k2);
                if (reflected)
                {
                    temp_i1 = 1<< (temp_nu.j()[i] + temp_nu.e()[i]);
                    temp_i2 = temp_i1-nui_k2;
                    nui_k2 = temp_i1-nui_k1;
                    nui_k1 = temp_i2;
                }
                get_intersecting_wavelets_on_level(*qtbasis.get_bases_infact()[mui_basisnum],
                        temp_nu.j()[i],
                        temp_nu.e()[i],
                        nui_k1,
                        nui_k2,
                        temp_mu.j()[i], 
                        (temp_mu.e()[i] == 0), 
                        temp_i1, 
                        temp_i2);
                //if (false)
                if ( (n == _n) && (m == _m))
                {
                    cout << "temp_i1 = " << temp_i1 << endl;
                    cout << "temp_i2 = " << temp_i2 << endl;
                }
                assert (kmin == temp_i1);
                assert (kmax == temp_i2);
            }
            
            if (!nu_and_mu_intersect)
                continue;
            
            // up to now  nu_and_mu_intersect meant "there exist a mu that intersects"
            // now we test if the current mu does! I.e., if it lies in the intersecting range computed above
            
            if ( (n == _n) && (m == _m))
            {
                Support temp_supp;
                qtbasis.support_local(temp_nu, temp_supp);
                cout << "support temp_nu: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                        "]; a = [" << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                        "]; b = [" << temp_supp.b[0] << ", " << temp_supp.b[1] << "]" << endl;

                Index temp_ind(qtbasis.get_wavelet(m));
                qtbasis.support_local(temp_ind, temp_supp);
                cout << "support temp_mu: j = [" << temp_supp.j[0] << ", " << temp_supp.j[1] << 
                        "]; a = [" << temp_supp.a[0] << ", " << temp_supp.a[1] <<
                        "]; b = [" << temp_supp.b[0] << ", " << temp_supp.b[1] << "]" << endl;
                cout << kmingen << endl;
                cout << kmaxgen << endl;
                cout << kminwav << endl;
                cout << kmaxwav << endl;
            }
            for (unsigned int i=0; i<DIM; i++)
            {
                if (temp_mu.e()[i] == 0)
                {
                    nu_and_mu_intersect = nu_and_mu_intersect && ( gen_intersection[i] && ( (temp_mu.k()[i] <= kmaxgen[i]) && (temp_mu.k()[i] >= kmingen[i]) ) );
                }
                else
                {
                    nu_and_mu_intersect = nu_and_mu_intersect && ( wav_intersection[i] && ( (temp_mu.k()[i] <= kmaxwav[i]) && (temp_mu.k()[i] >= kminwav[i]) ) );
                }
            }
            
            if (!nu_and_mu_intersect)
                continue;
            
            // include information about mu into LMR info
            // -> store update in intinfo2
            for (unsigned int i=0; i<DIM; ++i)
            {
                // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                switch (intinfo[i])
                {
                    case 0:
                        // LL: check whether mu is a left boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw() ) ) ) )
                        {
                            // mu is not extended left
                            intinfo2[i] = 4; // MM
                        }
                        
                        break;
                    case 2:
                        // LR: check whether mu is a right boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                          ||
                          ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                        {
                            intinfo2[i] = 1; // LM
                        }
                        break;
                    case 3:
                        // ML : check whether mu is a left boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw()) ) ) )
                        {
                            // mu is not extended left. There is no intersection of mu and nu!
                            if (temp_mu.e()[i] == 0)
                            {
                                assert(temp_mu.k()[i] > kmaxgen[i]);
                            }
                            else
                            {
                                assert(temp_mu.k()[i] < kmaxwav[i]);
                            }
                        }
                        break;
                    case 5:
                        // MR: check whether mu is a right boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                          ||
                          ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                        {
                            // mu is not extended right. There is no intersection of mu and nu!
                            if (temp_mu.e()[i] == 0)
                            {
                                assert(temp_mu.k()[i] < kmingen[i]);
                            }
                            else
                            {
                                assert(temp_mu.k()[i] < kminwav[i]);
                            }
                        }
                        break;
                    case 6:
                        // RL : check whether mu is a left boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw()) ) ) )
                        {
                            // mu is not extended left. There is no intersection of mu and nu!
                            intinfo2[i] = 7; // RM
                        }
                        break;
                    case 8:
                        // RR: check whether mu is a right boundary gen/wav
                        if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                          ||
                          ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                        {
                            // nu is not extended right
                            intinfo2[i] = 4; // MM
                        }
                        break;
                    case 1:
                    case 4:
                    case 7:
                        break;
                    default:
                        abort();
                        break;
                }
            }
            
            //if (false)
            if ((n == _n) && (m == _m))
            {
                cout << "temp_nu = " << temp_nu << endl;
                cout << "temp_mu = " << temp_mu << endl;
                cout << "intinfo = " << intinfo << endl;
                cout << "intinfo2 = " << intinfo2 << endl;
            }
            qtbasis.get_intersection_geometry(temp_nu.p(), 
                    intinfo2, 
                    geometry_type, 
                    centerpatchnumber, 
                    orientation);
            //if (false)
            //if (n == 273)
            if ((n == _n) && (m == _m))
            {
                cout << "geometry_type = " << geometry_type << endl;
                cout << "centerpatchnumber = " << centerpatchnumber << endl;
                cout << "orientation = " << orientation << endl;
            }

            /*
            cout << intinfo << endl;
            cout << geometry_type << endl;
            cout << centerpatchnumber << endl;
            cout << orientation << endl;
            */
// TODO : interpret (check for correctness of ) LMR Info
            
            // intinfo2
            if (false)
            //if (n == 277)
            {
                cout << "test" << endl;
            }
            
            FixedArray1D<bool,2> mu_lbw, mu_rbw, nu_lbw, nu_rbw;
            
            switch (decomposition) 
            {
                case 0: // 1 domain
                    assert (centerpatchnumber == 0);
                    assert (geometry_type == 0);
                    assert (orientation[0] == orientation[1]);
                    assert (orientation[0] == true);
                    assert (intinfo2[0] == 4);
                    assert (intinfo2[1] == 4);
                    break;
                case 1: // 2 domains, homogeneous BC, left patch is extended to the right
                    assert (orientation[1] == true);
                    assert (intinfo2[1] == 4);
                    
                    if ( (temp_mu.p() == 1)|| (temp_nu.p() == 1) )
                    {
                        assert (centerpatchnumber == 1);
                        assert (geometry_type == 0);
                        assert (orientation[0] == (temp_mu.p() == 1));
                        
                        if (intinfo2[0] == 4)
                        {
                            assert (temp_mu.p() == 1);
                            assert (temp_nu.p() == 1);
                        }
                        else
                        {
                            if (intinfo2[0] == 7) // RM
                            {
                                assert (temp_nu.p() == 0);
                                assert (temp_mu.p() == 1);
                                assert (     ( (temp_nu.e()[0] == 0) 
                                              && (temp_nu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0) ))
                                             )
                                            || 
                                             ( (temp_nu.e()[0] == 1)
                                              && ( temp_nu.k()[0] > (1<<(temp_nu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                             )
                                        );
                            }
                            else
                            {
                                if (intinfo2[0] == 5) // MR
                                {
                                    assert (temp_nu.p() == 1);
                                    assert (temp_mu.p() == 0);
                                    //cout << temp_nu << endl;
                                    //cout << temp_mu << endl;
                                    //cout << nu_and_mu_intersect << endl;
                                    //cout << temp_b2 << endl;
                                    
                                    //cout << qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0) ) << endl;
                                    //cout << ( temp_mu.k()[0] > (1<<(temp_mu.j()[0]) ) -1 - qtbasis.get_numofbw() ) << endl;
                                    assert (     ( (temp_mu.e()[0] == 0) 
                                                  && (temp_mu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0) ))
                                                 )
                                                || 
                                                 ( (temp_mu.e()[0] == 1)
                                                  && ( temp_mu.k()[0] > (1<<(temp_mu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                                 )
                                            );
                                }
                                else
                                {
                                    abort();
                                }
                            }
                        }
                    }
                    else 
                    {
                        // patch == 0
                        assert (centerpatchnumber == 0);
                        assert (orientation[0] == true);
                        assert (intinfo2[1] == 4); // MM
                        //qtbasis.get_j0(0,0);
                        //(((qtbasis.get_bc()[0][0])?0:2) + ((qtbasis.get_bc()[0][1])?0:1));
                        //qtbasis.get_bases_infact()[0][0];
                        //qtbasis.get_bases_infact()[0][(((qtbasis.get_bc()[0][0])?0:2) + ((qtbasis.get_bc()[0][1])?0:1))].DeltaRmax(qtbasis.get_j0(0,0));
                        // nui_basisnum in x direction: 
                        //temp_i1 = (((qtbasis.get_bc()[temp_nu.p()][0])?0:2) + ((qtbasis.get_bc()[temp_nu.p()][1])?0:1));
                        //mui_basisnum = (((qtbasis.get_bc()[temp_mu.p()][2*i])?0:2) + ((qtbasis.get_bc()[temp_mu.p()][2*i+1])?0:1));
                        if (false)
                        //if ((n==273) && (m==732))
                        {
                            cout << "temp_mu = " << temp_mu << endl;
                            cout << ( (temp_mu.e()[0] == 0) 
                                         && (temp_mu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0)) )
                                       ) << endl;
                            cout << "(1<<(temp_mu.j()[0]) ) = " << (1<<(temp_mu.j()[0]) ) << endl;
                            cout << "qtbasis.get_numofbw() = " << qtbasis.get_numofbw() << endl;
                            cout << ( (temp_mu.e()[0] == 1)
                                         && ( temp_mu.k()[0] > (1<<(temp_mu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                       ) << endl;
                            cout << "temp_nu = " << temp_nu << endl;
                            cout << ( (temp_nu.e()[0] == 0) 
                                         && (temp_nu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0)) )
                                       ) << endl;
                            cout << ( (temp_nu.e()[0] == 1)
                                         && ( temp_nu.k()[0] > (1<<(temp_nu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                       ) << endl;
                        }
                        
                                
                        if  (  (
                                   ( (temp_mu.e()[0] == 0) 
                                     && (temp_mu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0)) )
                                   )
                               ||  ( (temp_mu.e()[0] == 1)
                                     && ( temp_mu.k()[0] > (1<<(temp_mu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                   )
                               )
                            && (
                                   ( (temp_nu.e()[0] == 0) 
                                     && (temp_nu.k()[0] == qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0)) )
                                   )
                               ||  ( (temp_nu.e()[0] == 1)
                                     && ( temp_nu.k()[0] > (1<<(temp_nu.j()[0]) ) -1 - qtbasis.get_numofbw() ) 
                                   )
                               )
                            )
                        {
                            // mui, nui are both extended
                            assert (geometry_type == 1);
                            assert (centerpatchnumber == 0);
                            assert (intinfo2[0] == 8); // RR
                        }
                        else
                        {
                            assert (geometry_type == 0);
                            assert (centerpatchnumber == 0);
                            assert (intinfo2[0] == 4); // MM
                        }
                    }
                    break;
                case 7: // 2 domains, homogeneous BC, patch 0 is extended northwards to patch 1 (right in the y direction for DIM=2)
                    
                    assert (orientation[0] == true);
                    assert (intinfo2[0] == 4);
                    if ( (temp_mu.p() == 1)|| (temp_nu.p() == 1) )
                    {
                        assert (centerpatchnumber == 1);
                        assert (geometry_type == 0);
                        assert (orientation[1] == (temp_mu.p() == 1));
                        
                        if (intinfo2[1] == 4)
                        {
                            assert (temp_mu.p() == 1);
                            assert (temp_nu.p() == 1);
                        }
                        else
                        {
                            if (intinfo2[1] == 7) // RM
                            {
                                assert (temp_nu.p() == 0);
                                assert (temp_mu.p() == 1);
                                assert (     ( (temp_nu.e()[1] == 0) 
                                              && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaRmax(qtbasis.get_j0(0,1)))
                                             )
                                            || 
                                             ( (temp_nu.e()[1] == 1)
                                              && ( temp_nu.k()[1] > (1<<(temp_nu.j()[1]) ) -1 - qtbasis.get_numofbw() ) 
                                             )
                                        );
                            }
                            else
                            {
                                if (intinfo2[1] == 5) // MR
                                {
                                    assert (temp_nu.p() == 1);
                                    assert (temp_mu.p() == 0);
                                    assert (     ( (temp_mu.e()[1] == 0) 
                                                  && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaRmax(qtbasis.get_j0(0,1) ))
                                                 )
                                                || 
                                                 ( (temp_mu.e()[1] == 1)
                                                  && ( temp_mu.k()[1] > (1<<(temp_mu.j()[1]) ) -1 - qtbasis.get_numofbw() ) 
                                                 )
                                            );
                                }
                                else
                                {
                                    assert(false);
                                }
                            }
                        }
                    }
                    else
                    {
                        // patch == 0
                        if  (  (
                                   ( (temp_mu.e()[1] == 0) 
                                     && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaRmax(qtbasis.get_j0(0,1) ) )
                                   )
                               ||  ( (temp_mu.e()[1] == 1)
                                     && (temp_mu.k()[1] > (1<<(temp_mu.j()[1]) ) -1 - qtbasis.get_numofbw() )
                                   )
                               )
                            && (
                                   ( (temp_nu.e()[1] == 0) 
                                     && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaRmax(qtbasis.get_j0(0,1) ) )
                                   )
                               ||  ( (temp_nu.e()[1] == 1)
                                     && (temp_nu.k()[1] > (1<<(temp_nu.j()[1]) ) -1 - qtbasis.get_numofbw() )
                                   )
                               )
                            )
                        {
                            // mui, nui are both extended
                            assert (geometry_type == 2);
                            assert (centerpatchnumber == 0);
                            assert (orientation[1] == true);
                            assert (intinfo2[1] == 8); // RR
                        }
                        else
                        {
                            assert (geometry_type == 0);
                            assert (centerpatchnumber == 0);
                            assert (orientation[1] == true);
                            assert (intinfo2[1] == 4); // MM
                        }
                    }
                    break;
                case 8: // 2 domains, homogeneous BC, right patch is extended to the left
                    assert (orientation[1] == true);
                    assert (intinfo2[1] == 4);
                    
                    if ( (temp_mu.p() == 0)|| (temp_nu.p() == 0) )
                    {
                        assert (centerpatchnumber == 0);
                        assert (geometry_type == 0);
                        assert (orientation[0] == (temp_mu.p() == 0));
                        
                        if (intinfo2[0] == 4)
                        {
                            assert (temp_mu.p() == 0);
                            assert (temp_nu.p() == 0);
                        }
                        else
                        {
                            if (intinfo2[0] == 1) // LM
                            {
                                assert (temp_nu.p() == 1);
                                assert (temp_mu.p() == 0);
                                assert (     ( (temp_nu.e()[0] == 0) 
                                              && (temp_nu.k()[0] == qtbasis.bases() [1] [0] -> DeltaLmin())
                                             )
                                            || 
                                             ( (temp_nu.e()[0] == 1)
                                              && ( temp_nu.k()[0] < qtbasis.get_numofbw() ) 
                                             )
                                        );
                            }
                            else
                            {
                                if (intinfo2[0] == 3) // ML
                                {
                                    assert (temp_nu.p() == 0);
                                    assert (temp_mu.p() == 1);
                                    //cout << temp_nu << endl;
                                    //cout << temp_mu << endl;
                                    //cout << nu_and_mu_intersect << endl;
                                    //cout << temp_b2 << endl;
                                    
                                    //cout << qtbasis.bases() [0] [0] -> DeltaRmax(qtbasis.get_j0(0,0) ) << endl;
                                    //cout << ( temp_mu.k()[0] > (1<<(temp_mu.j()[0]) ) -1 - qtbasis.get_numofbw() ) << endl;
                                    assert (     ( (temp_mu.e()[0] == 0) 
                                                  && (temp_mu.k()[0] == qtbasis.bases() [1] [0] -> DeltaLmin())
                                                 )
                                                || 
                                                 ( (temp_mu.e()[0] == 1)
                                                  && ( temp_mu.k()[0] < qtbasis.get_numofbw() ) 
                                                 )
                                            );
                                }
                                else
                                {
                                    abort();
                                }
                            }
                        }
                    }
                    else 
                    {
                        // patch == 1
                        if  (  (
                                   ( (temp_mu.e()[0] == 0) 
                                     && (temp_mu.k()[0] == qtbasis.bases() [1] [0] -> DeltaLmin() )
                                   )
                               ||  ( (temp_mu.e()[0] == 1)
                                     && ( temp_mu.k()[0] < qtbasis.get_numofbw() ) 
                                   )
                               )
                            && (
                                   ( (temp_nu.e()[0] == 0) 
                                     && (temp_nu.k()[0] == qtbasis.bases() [1] [0] -> DeltaLmin() )
                                   )
                               ||  ( (temp_nu.e()[0] == 1)
                                     && ( temp_nu.k()[0] < qtbasis.get_numofbw() ) 
                                   )
                               )
                            )
                        {
                            // mui, nui are both extended
                            assert (geometry_type == 1);
                            assert (orientation[0] == false);
                            assert (centerpatchnumber == 0);
                            assert (intinfo2[0] == 0); // LL
                        }
                        else
                        {
                            assert (geometry_type == 0);
                            assert (orientation[0] == true);
                            assert (centerpatchnumber == 1);
                            assert (intinfo2[0] == 4); // MM
                        }
                    }
                    break;
                case 2: // 2 domains, homogeneous BC, patch 1 is extended southwards to patch 0 (left in the y direction for DIM=2)
                    assert (orientation[0] == true);
                    assert (intinfo2[0] == 4);
                    if ( (temp_mu.p() == 0)|| (temp_nu.p() == 0) )
                    {
                        assert (centerpatchnumber == 0);
                        assert (geometry_type == 0);
                        assert (orientation[1] == (temp_mu.p() == 0));
                        
                        if (intinfo2[1] == 4)
                        {
                            assert (temp_mu.p() == 0);
                            assert (temp_nu.p() == 0);
                        }
                        else
                        {
                            if (intinfo2[1] == 1) // LM
                            {
                                assert (temp_nu.p() == 1);
                                assert (temp_mu.p() == 0);
                                assert (     ( (temp_nu.e()[1] == 0) 
                                              && (temp_nu.k()[1] == qtbasis.bases() [1] [1] -> DeltaLmin())
                                             )
                                            || 
                                             ( (temp_nu.e()[1] == 1)
                                              && ( temp_nu.k()[1] < qtbasis.get_numofbw() ) 
                                             )
                                        );
                            }
                            else
                            {
                                if (intinfo2[1] == 3) // ML
                                {
                                    assert (temp_nu.p() == 0);
                                    assert (temp_mu.p() == 1);
                                    assert (     ( (temp_mu.e()[1] == 0) 
                                                  && (temp_mu.k()[1] == qtbasis.bases() [1] [1] -> DeltaLmin() )
                                                 )
                                                || 
                                                 ( (temp_mu.e()[1] == 1)
                                                  && ( temp_mu.k()[1] < qtbasis.get_numofbw() ) 
                                                 )
                                            );
                                }
                                else
                                {
                                    assert(false);
                                }
                            }
                        }
                    }
                    else
                    {
                        // patch == 1
                        if  (  (
                                   ( (temp_mu.e()[1] == 0) 
                                     && (temp_mu.k()[1] == qtbasis.bases() [1] [1] -> DeltaLmin() )
                                   )
                               ||  ( (temp_mu.e()[1] == 1)
                                     && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                   )
                               )
                            && (
                                   ( (temp_nu.e()[1] == 0) 
                                     && (temp_nu.k()[1] == qtbasis.bases() [1] [1] -> DeltaLmin() )
                                   )
                               ||  ( (temp_nu.e()[1] == 1)
                                     && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                   )
                               )
                            )
                        {
                            // mui, nui are both extended
                            assert (geometry_type == 2);
                            assert (centerpatchnumber == 0);
                            assert (orientation[1] == false);
                            assert (intinfo2[1] == 0); // LL
                        }
                        else
                        {
                            assert (geometry_type == 0);
                            assert (centerpatchnumber == 1);
                            assert (orientation[1] == true);
                            assert (intinfo2[1] == 4); // MM
                        }
                    }
                    break;
                case 3: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 1
                    switch (geometry_type)
                    {
                        case 1:
                            assert (centerpatchnumber == 1);
                            assert (intinfo2[1] == 4); // MM
                            switch (intinfo2[0])
                            {
                                case 0: // LL
                                    assert (temp_nu.p() == 2);
                                    assert (temp_mu.p() == 2);
                                    assert ( ( (temp_mu.e()[0] == 0) 
                                              && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_mu.e()[0] == 1)
                                              && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert ( ( (temp_nu.e()[0] == 0) 
                                              && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_nu.e()[0] == 1)
                                              && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert (orientation[0] == false);
                                    assert (orientation[1] == true);
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        case 2:
                            assert (centerpatchnumber == 1);
                            assert (intinfo2[0] == 4);
                            switch (intinfo2[1])
                            {
                                case 0: // LL
                                    assert (temp_nu.p() == 0);
                                    assert (temp_mu.p() == 0);
                                    assert ( ( (temp_mu.e()[1] == 0) 
                                              && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_mu.e()[1] == 1)
                                              && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert ( ( (temp_nu.e()[1] == 0) 
                                              && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_nu.e()[1] == 1)
                                              && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert (orientation[0] == true);
                                    assert (orientation[1] == false);
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        case 0:
                            switch (centerpatchnumber)
                            {
                                case 0:
                                    assert (temp_nu.p() == 0);
                                    assert (temp_mu.p() == 0);
                                    assert (intinfo2[0] == 4);
                                    assert (intinfo2[1] == 4);
                                    assert (orientation[0] == true);
                                    assert (orientation[1] == true);
                                    temp_b = ( ( (temp_mu.e()[1] == 0) 
                                              && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_mu.e()[1] == 1)
                                              && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    temp_b = temp_b &&
                                            (
                                            ( ( (temp_nu.e()[1] == 0) 
                                              && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_nu.e()[1] == 1)
                                              && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                             )
                                           )
                                           );
                                    assert (!temp_b); // otherwise LL in y direction
                                    break;
                                case 1:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert ( ( (temp_nu.e()[1] == 0) 
                                                    && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                    )
                                                    || 
                                                    ( (temp_nu.e()[1] == 1)
                                                    && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                                    )
                                                    );
                                            if (temp_mu.p() == 1)
                                            {
                                                assert(intinfo[0] == 4);
                                                assert(intinfo[1] == 1); // LM
                                                assert(orientation[0] == true);
                                                assert(orientation[1] == true);
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 2)
                                                {
                                                    assert(intinfo[0] == 3); // ML
                                                    assert(intinfo[1] == 1); // LM
                                                    assert(orientation[0] == false);
                                                    assert(orientation[1] == true);
                                                    assert ( ( (temp_mu.e()[0] == 0) 
                                                            && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                                            )
                                                            || 
                                                            ( (temp_mu.e()[0] == 1)
                                                            && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                                            )
                                                            );
                                                }
                                                else
                                                {
                                                    abort();
                                                }
                                            }
                                            break;
                                        case 1:
                                            if (temp_mu.p() == 0)
                                            {
                                                assert(intinfo2[0] == 4);
                                                assert(intinfo2[1] == 3); // ML
                                                assert(orientation[0] == true);
                                                assert(orientation[1] == false);
                                                assert ( ( (temp_mu.e()[1] == 0) 
                                                        && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                        )
                                                        || 
                                                        ( (temp_mu.e()[1] == 1)
                                                        && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                                        )
                                                        );
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 1)
                                                {
                                                    assert(intinfo2[0] == 4); // MM
                                                    assert(intinfo2[1] == 4); // MM
                                                    assert(orientation[0] == true);
                                                    assert(orientation[1] == true);
                                                }
                                                else
                                                {
                                                    assert (temp_mu.p() == 2);
                                                    assert(intinfo2[0] == 3); // ML
                                                    assert(intinfo2[1] == 4); // MM
                                                    assert(orientation[0] == false);
                                                    assert(orientation[1] == true);
                                                    assert ( ( (temp_mu.e()[0] == 0) 
                                                            && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                                            )
                                                            || 
                                                            ( (temp_mu.e()[0] == 1)
                                                            && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                                            )
                                                            );
                                                }
                                            }
                                            break;
                                        case 2:
                                            assert ( ( (temp_nu.e()[0] == 0) 
                                                && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                                )
                                                || 
                                                ( (temp_nu.e()[0] == 1)
                                                && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                                )
                                                );
                                            if (temp_mu.p() == 0)
                                            {
                                                assert(intinfo2[0] == 1); // LM
                                                assert(intinfo2[1] == 3); // ML
                                                assert(orientation[0] == true);
                                                assert(orientation[1] == false);
                                                assert ( ( (temp_mu.e()[1] == 0) 
                                                        && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                        )
                                                        || 
                                                        ( (temp_mu.e()[1] == 1)
                                                        && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                                        )
                                                        );
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 1)
                                                {
                                                    assert(intinfo2[0] == 1); // MM
                                                    assert(intinfo2[1] == 4); // MM
                                                    assert(orientation[0] == true);
                                                    assert(orientation[1] == true);
                                                }
                                                else
                                                {
                                                    abort();
                                                }
                                            }
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 2:
                                    assert (temp_nu.p() == 2);
                                    assert (temp_mu.p() == 2);
                                    assert (intinfo2[0] == 4);
                                    assert (intinfo2[1] == 4);
                                    assert (orientation[0] == true);
                                    assert (orientation[1] == true);
                                    temp_b = ( ( (temp_mu.e()[0] == 0) 
                                              && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_mu.e()[0] == 1)
                                              && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    temp_b = temp_b && ( ( (temp_nu.e()[0] == 0) 
                                              && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_nu.e()[0] == 1)
                                              && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert (!temp_b); // otherwise LL in x direction
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        default:
                            abort();
                            break;
                    }
                    break;
                case 4: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 0
                    switch (geometry_type)
                    {
                        case 1:
                            assert (centerpatchnumber == 0);
                            assert (intinfo2[1] == 4); // MM
                            assert (intinfo2[0] == 0); // LL
                            
                            assert (temp_nu.p() == 2);
                            assert (temp_mu.p() == 2);
                            assert ( ( (temp_mu.e()[0] == 0) 
                                      && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                     )
                                    || 
                                     ( (temp_mu.e()[0] == 1)
                                      && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                     )
                                   );
                            assert ( ( (temp_nu.e()[0] == 0) 
                                      && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                     )
                                    || 
                                     ( (temp_nu.e()[0] == 1)
                                      && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                     )
                                   );
                            assert (orientation[0] == false);
                            assert (orientation[1] == true);
                            break;
                        case 2:
                            assert (centerpatchnumber == 1);
                            assert (intinfo2[0] == 4);
                            assert (intinfo2[1] == 0); // LL
                            assert (temp_nu.p() == 0);
                            assert (temp_mu.p() == 0);
                            assert ( ( (temp_mu.e()[1] == 0) 
                                      && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                     )
                                    || 
                                     ( (temp_mu.e()[1] == 1)
                                      && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                     )
                                   );
                            assert ( ( (temp_nu.e()[1] == 0) 
                                      && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                     )
                                    || 
                                     ( (temp_nu.e()[1] == 1)
                                      && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                     )
                                   );
                            assert (orientation[0] == true);
                            assert (orientation[1] == false);
                            break;
                        case 0:
                            switch (centerpatchnumber)
                            {
                                case 0:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            if (temp_mu.p() == 0)
                                            {
                                                assert(intinfo2[0] == 4);
                                                assert(intinfo2[1] == 4);
                                                assert(orientation[0] == true);
                                                assert(orientation[1] == true);
                                                temp_b = ( ( (temp_nu.e()[1] == 0) 
                                                        && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                        )
                                                        || 
                                                        ( (temp_nu.e()[1] == 1)
                                                        && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                                        )
                                                        );
                                                temp_b = temp_b && ( ( (temp_mu.e()[1] == 0) 
                                                        && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                        )
                                                        || 
                                                        ( (temp_mu.e()[1] == 1)
                                                        && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                                        )
                                                        );
                                                assert (!temp_b);
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 2)
                                                {
                                                    assert(intinfo2[0] == 3); // ML
                                                    assert(intinfo2[1] == 4); // MM
                                                    assert(orientation[0] == false);
                                                    assert(orientation[1] == true);
                                                    assert ( ( (temp_mu.e()[0] == 0) 
                                                            && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                                            )
                                                            || 
                                                            ( (temp_mu.e()[0] == 1)
                                                            && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                                            )
                                                            );
                                                }
                                                else
                                                {
                                                    abort();
                                                }
                                            }
                                            break;
                                        case 1:
                                            abort();
                                            break;
                                        case 2:
                                            assert ( ( (temp_nu.e()[0] == 0) 
                                                && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                                )
                                                || 
                                                ( (temp_nu.e()[0] == 1)
                                                && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                                )
                                                );
                                            assert (temp_mu.p() == 0);
                                            
                                            assert(intinfo2[0] == 1); // LM
                                            assert(intinfo2[1] == 4); // MM
                                            assert(orientation[0] == true);
                                            assert(orientation[1] == true);
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 1:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 1);
                                            assert(intinfo2[0] == 4);
                                            assert(intinfo2[1] == 1); // LM
                                            assert(orientation[0] == true);
                                            assert(orientation[1] == true);
                                            assert( ( (temp_nu.e()[1] == 0) 
                                                    && (temp_nu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                    )
                                                    || 
                                                    ( (temp_nu.e()[1] == 1)
                                                    && (temp_nu.k()[1] < qtbasis.get_numofbw()) 
                                                    )
                                                    );
                                            break;
                                        case 1:
                                            if (temp_mu.p() == 0)
                                            {
                                                assert(intinfo2[0] == 4);
                                                assert(intinfo2[1] == 3); // ML
                                                assert(orientation[0] == true);
                                                assert(orientation[1] == false);
                                                assert( ( (temp_mu.e()[1] == 0) 
                                                        && (temp_mu.k()[1] == qtbasis.bases() [0] [1] -> DeltaLmin() )
                                                        )
                                                        || 
                                                        ( (temp_mu.e()[1] == 1)
                                                        && (temp_mu.k()[1] < qtbasis.get_numofbw()) 
                                                        )
                                                        );
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 1)
                                                {
                                                    assert(intinfo2[0] == 4);
                                                    assert(intinfo2[1] == 4);
                                                    assert(orientation[0] == true);
                                                    assert(orientation[1] == true);
                                                }
                                                else
                                                {
                                                    abort();
                                                }
                                            }
                                            break;
                                        case 2:
                                            abort();
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 2:
                                    assert (temp_nu.p() == 2);
                                    assert (temp_mu.p() == 2);
                                    assert (intinfo2[0] == 4);
                                    assert (intinfo2[1] == 4);
                                    assert (orientation[0] == true);
                                    assert (orientation[1] == true);
                                    temp_b = ( ( (temp_mu.e()[0] == 0) 
                                              && (temp_mu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_mu.e()[0] == 1)
                                              && (temp_mu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    temp_b = temp_b && ( ( (temp_nu.e()[0] == 0) 
                                              && (temp_nu.k()[0] == qtbasis.bases() [2] [0] -> DeltaLmin() )
                                             )
                                            || 
                                             ( (temp_nu.e()[0] == 1)
                                              && (temp_nu.k()[0] < qtbasis.get_numofbw()) 
                                             )
                                           );
                                    assert (!temp_b); // otherwise LL in x direction
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        default:
                            abort();
                            break;
                    }
                    break;
                    
                    
                    
                    
                    

                case 6: // Cube, decomposed into 4 subcubes, homogeneous BC
                    // Patch 0 is extended to the east to patch 1 and north to patch 3.
                    // Patch 1 is extended to the north to patch 2
                    // Patch 2 is extended to the west to patch 3
                    
                    for (unsigned int i = 0; i<DIM; ++i)
                    {
                        // mu is extended left?
                        mu_lbw[i] = ( ( (temp_mu.e()[i] == 0) 
                                          && (temp_mu.k()[i] == qtbasis.bases() [temp_mu.p()] [i] -> DeltaLmin() )
                                         )
                                        || 
                                         ( (temp_mu.e()[i] == 1)
                                          && (temp_mu.k()[i] < qtbasis.get_numofbw()) 
                                         )
                                       );
                        // nu is extended left?
                        nu_lbw[i] = ( ( (temp_nu.e()[i] == 0) 
                                          && (temp_nu.k()[i] == qtbasis.bases() [temp_nu.p()] [i] -> DeltaLmin() )
                                         )
                                        || 
                                         ( (temp_nu.e()[i] == 1)
                                          && (temp_nu.k()[i] < qtbasis.get_numofbw()) 
                                         )
                                       );
                        // mu is extended right?
                        mu_rbw[i] = ( ( (temp_mu.e()[i] == 0) 
                                          && (temp_mu.k()[i] == qtbasis.bases() [temp_mu.p()] [i] -> DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) ))
                                         )
                                        || 
                                         ( (temp_mu.e()[i] == 1)
                                          && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw())
                                         )
                                       );
                        // nu is extended right?
                        nu_rbw[i] = ( ( (temp_nu.e()[i] == 0) 
                                          && (temp_nu.k()[i] == qtbasis.bases() [temp_nu.p()] [i] -> DeltaRmax( qtbasis.get_j0(temp_nu.p(),i) ))
                                         )
                                        || 
                                         ( (temp_nu.e()[i] == 1)
                                          && (temp_nu.k()[i] > (1<<temp_nu.j()[i])- 1 - qtbasis.get_numofbw())
                                         )
                                       ); 
                        assert (! (nu_lbw[i] && nu_rbw[i]));
                        assert (! (mu_lbw[i] && mu_rbw[i]));
                    }
                    
                    switch (geometry_type)
                    {
                        case 0: // decomposition 6; geometry type 0
                            switch (centerpatchnumber)
                            {
                                case 0:
                                    assert (!nu_rbw[0] || !mu_rbw[0]);
                                    assert (!nu_rbw[1] || !mu_rbw[1]);
                                    assert (temp_nu.p() == 0);
                                    assert (temp_mu.p() == 0);
                                    assert (intinfo2[0] == 4);
                                    assert (intinfo2[1] == 4);
                                    assert (orientation[0] == true);
                                    assert (orientation[1] == true);
                                    break;
                                case 1:
                                    if (temp_nu.p() == 1)
                                    {
                                        if (temp_mu.p() == 1)
                                        {
                                            assert (!nu_rbw[1] || !mu_rbw[1]);
                                            assert (intinfo2[0] == 4);
                                            assert (intinfo2[1] == 4);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                        }
                                        else
                                        {
                                            assert (temp_mu.p() == 0);
                                            assert (mu_rbw[0] == true);
                                            assert (!nu_rbw[1] || !mu_rbw[1]);
                                            assert (intinfo2[0] == 5);
                                            assert (intinfo2[1] == 4);
                                            assert (orientation[0] == false);
                                            assert (orientation[1] == true);
                                        }
                                    }
                                    else
                                    {
                                        assert ( (temp_nu.p() == 0) && (temp_mu.p() == 1) );
                                        assert (nu_rbw[0] == true);
                                        assert (!nu_rbw[1] || !mu_rbw[1]);
                                        assert (intinfo2[0] == 7);
                                        assert (intinfo2[1] == 4);
                                        assert (orientation[0] == true);
                                        assert (orientation[1] == true);
                                    }
                                    break;
                                case 2:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 2);
                                            assert(nu_rbw[0] == true);
                                            assert(nu_rbw[1] == true);
                                            assert(mu_lbw[0] == false);
                                            assert (intinfo2[0] == 7);
                                            assert (intinfo2[1] == 7);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                            break;
                                        case 1:
                                            assert (nu_rbw[1] == true);
                                            assert (temp_mu.p() == 2);
                                            assert (intinfo2[0] == 4);
                                            assert (intinfo2[1] == 7);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                            break;
                                        case 2:
                                            if (temp_mu.p() == 0)
                                            {
                                                assert(mu_rbw[0] == true);
                                                assert(mu_rbw[1] == true);
                                                assert(nu_lbw[0] == false);
                                                assert (intinfo2[0] == 5);
                                                assert (intinfo2[1] == 5);
                                                assert (orientation[0] == false);
                                                assert (orientation[1] == false);    
                                            }
                                            else
                                            {
                                                if (temp_mu.p() == 1)
                                                {
                                                    assert (mu_rbw[1] == true);
                                                    assert (intinfo2[0] == 4);
                                                    assert (intinfo2[1] == 5);
                                                    assert (orientation[0] == true);
                                                    assert (orientation[1] == false);
                                                }
                                                else
                                                {
                                                    assert (temp_mu.p() == 2);
                                                    assert (!nu_lbw[0] || !mu_lbw[0]);
                                                    assert (intinfo2[0] == 4);
                                                    assert (intinfo2[1] == 4);
                                                    assert (orientation[0] == true);
                                                    assert (orientation[1] == true);
                                                }
                                            }
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 3:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert(nu_rbw[1] == true);
                                            if (temp_mu.p() == 2)
                                            {
                                                assert (nu_rbw[0] == false);
                                                assert (intinfo2[0] == 3);
                                                assert (intinfo2[1] == 7);
                                                assert (orientation[0] == false);
                                                assert (orientation[1] == true);
                                            }
                                            else
                                            {                                        
                                                assert (temp_mu.p() == 3);
                                                assert (intinfo2[0] == 4);
                                                assert (intinfo2[1] == 7);
                                                assert (orientation[0] == true);
                                                assert (orientation[1] == true);
                                            }
                                            break;
                                        case 2:
                                            assert (nu_lbw[0] == true);
                                            if (temp_mu.p() == 0)
                                            {
                                                assert (mu_rbw[0] == false);
                                                assert (mu_rbw[1] == true);
                                                assert (intinfo2[0] == 1);
                                                assert (intinfo2[1] == 5);
                                                assert (orientation[0] == true);
                                                assert (orientation[1] == false);
                                            }
                                            else
                                            {
                                                assert (temp_mu.p() == 3);
                                                assert (intinfo2[0] == 1);
                                                assert (intinfo2[1] == 4);
                                                assert (orientation[0] == true);
                                                assert (orientation[1] == true);
                                            }
                                            break;
                                        case 3:
                                            switch (temp_mu.p())
                                            {
                                                case 0:
                                                    break;
                                                    assert(mu_rbw[1] == true);
                                                    assert (intinfo2[0] == 4);
                                                    assert (intinfo2[1] == 5);
                                                    assert (orientation[0] == true);
                                                    assert (orientation[1] == false);    
                                                case 2:
                                                    assert (mu_lbw[0] == true);
                                                    assert (intinfo2[0] == 3);
                                                    assert (intinfo2[1] == 4);
                                                    assert (orientation[0] == false);
                                                    assert (orientation[1] == true);
                                                    break;
                                                case 3:
                                                    assert (intinfo2[0] == 4);
                                                    assert (intinfo2[1] == 4);
                                                    assert (orientation[0] == true);
                                                    assert (orientation[1] == true);
                                                    break;
                                                default:
                                                    abort();
                                                    break;
                                            }
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        case 1: // decomposition 6; geometry type 1
                            switch (centerpatchnumber)
                            {
                                case 0:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 0);
                                            assert (nu_rbw[0] == true);
                                            assert (mu_rbw[0] == true);
                                            assert (!nu_rbw[1] || !mu_rbw[1]);
                                            assert (intinfo2[0] == 8);
                                            assert (intinfo2[1] == 4);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 3:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 2);
                                            assert (nu_rbw[0] == true);
                                            assert (nu_rbw[1] == true);
                                            assert (mu_lbw[0] == true);
                                            assert (intinfo2[0] == 6);
                                            assert (intinfo2[1] == 7);
                                            assert (orientation[0] == false);
                                            assert (orientation[1] == true);
                                            break;
                                        case 2:
                                            
                                            /*                                           
                                                if (mu_lbw[0] != false)
                                            {
                                                cout << "temp_nu = " << temp_nu << endl;
                                                cout << "temp_mu = " << temp_mu << endl;
                                                cout << "nu_rbw = " << nu_rbw << endl;
                                                cout << "nu_lbw = " << nu_lbw << endl;
                                                cout << "mu_rbw = " << mu_rbw << endl;
                                                cout << "mu_lbw = " << mu_lbw << endl;
                                                for (unsigned int i = 0; i < DIM; ++i)
                                                {
                                                    
                                                    cout << "temp_mu.e()[i] = " << temp_mu.e()[i] << endl;
                                                    cout << "temp_mu.k()[i] = " << temp_mu.k()[i] << endl;
                                                    cout << "qtbasis.bases() [temp_mu.p()] [i] -> DeltaLmin() = " << qtbasis.bases() [temp_mu.p()] [i] -> DeltaLmin() << endl;
                                                    cout << ( (temp_mu.e()[i] == 0) 
                                                            && (temp_mu.k()[i] == qtbasis.bases() [temp_mu.p()] [i] -> DeltaLmin() )
                                                            ) << endl;
                                                    cout << ( (temp_mu.e()[i] == 1)
                                                            && (temp_mu.k()[i] < qtbasis.get_numofbw()) 
                                                            ) << endl;
                                                }
                                                cout << "orientation = " << orientation << endl;
                                                cout << "intinfo2 = " << intinfo2 << endl;
                                            }
                                                  */ 
                                            
                                            
                                            if (temp_mu.p() == 0)
                                            {
                                                assert (nu_lbw[0] == true);
                                                assert (mu_rbw[0] == true);
                                                assert (mu_rbw[1] == true);
                                                assert (intinfo2[0] == 2);
                                                assert (intinfo2[1] == 5);
                                                assert (orientation[0] == true);
                                                assert (orientation[1] == false);
                                            }
                                            else
                                            {
                                                assert (temp_mu.p() == 2);
                                                assert (nu_lbw[0] == true);
                                                assert (mu_lbw[0] == true);
                                                assert (intinfo2[0] == 0);
                                                assert (intinfo2[1] == 4);
                                                assert (orientation[0] == false);
                                                assert (orientation[1] == true);
                                            }
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        case 2: // decomposition 6; geometry type 2
                            switch (centerpatchnumber)
                            {
                                case 0:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 0);
                                            assert (!nu_rbw[0] || !mu_rbw[0]);
                                            assert (nu_rbw[1] == true);
                                            assert (mu_rbw[1] == true);
                                            assert (intinfo2[0] == 4);
                                            assert (intinfo2[1] == 8);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                case 1:
                                    switch (temp_nu.p())
                                    {
                                        case 0:
                                            assert (temp_mu.p() == 1);
                                            assert (nu_rbw[0] == true);
                                            assert (nu_rbw[1] == true);
                                            assert (mu_rbw[1] == true);
                                            assert (intinfo2[0] == 7);
                                            assert (intinfo2[1] == 8);
                                            assert (orientation[0] == true);
                                            assert (orientation[1] == true);
                                            break;
                                        case 1:
                                            if (temp_mu.p() == 0)
                                            {
                                                assert (nu_rbw[1] == true);
                                                assert (mu_rbw[0] == true);
                                                assert (mu_rbw[1] == true);
                                                assert (intinfo2[0] == 5);
                                                assert (intinfo2[1] == 8);
                                                assert (orientation[0] == false);
                                                assert (orientation[1] == true);
                                            }
                                            else
                                            {
                                                assert (temp_mu.p() == 1);
                                                assert (nu_rbw[1] == true);
                                                assert (mu_rbw[1] == true);
                                                assert (intinfo2[0] == 4);
                                                assert (intinfo2[1] == 8);
                                                assert (orientation[0] == true);
                                                assert (orientation[1] == true);
                                            }
                                            break;
                                        default:
                                            abort();
                                            break;
                                    }
                                    break;
                                default:
                                    abort();
                                    break;
                            }
                            break;
                        case 7: // decomposition 6; geometry type 7
                            assert (temp_nu.p() == 0);
                            assert (temp_mu.p() == 0);
                            assert (nu_rbw[0] == true);
                            assert (nu_rbw[1] == true);
                            assert (mu_rbw[0] == true);
                            assert (mu_rbw[1] == true);
                            assert (intinfo2[0] == 8);
                            assert (intinfo2[1] == 8);
                            assert (orientation[0] == true);
                            assert (orientation[1] == true);
                    } // end of decomposition 6; switch (geometry type)
                    break;
                default:
                    cout << "main:: error! test not implemented for this decomposition!" << endl;
                    abort();
                    break;
            }
        } // end of loop over m
    } // end of loop over n
    cout << "done" << endl;
#else
    cout << "skipping test of routines from qtbasis"  << endl;
#endif
    
    

    
    cout << "main:: Initialize f, a and q coeffs, cached qt problem ..." << endl;
    
    char rhs_filename[250];
//    const char* path = "/import/shared/friedrich/source/precomputed/";
    const char* path = "output/";
    sprintf(rhs_filename, "%sf_prec_unso_rhs%d_dec%d_off%d_%d_%d",path,_RHS_NUMBER,_DECOMP,_MAXOFFSET,_PRIMAL_ORDER_D,_DUAL_ORDER_DT);

    char logfile[250];
    sprintf(logfile,"%slog_rhs%d_dec%d_off%d_%d_%d_hmax%d.m",path,_RHS_NUMBER,_DECOMP,_MAXOFFSET,_PRIMAL_ORDER_D,_DUAL_ORDER_DT,_ONEDIMHAARCOUNT);
    
    
    
    /*
     * if DIM == 2: plot coeffs of the solution with options "true" (primary) and "scaled"
     */
#if (_DIMENSION == 2)    
    char plotfile[250];
    sprintf(plotfile,"%splot_rhs%d_dec%d_off%d_%d_%d_hmax%d.m",path,_RHS_NUMBER,_DECOMP,_MAXOFFSET,_PRIMAL_ORDER_D,_DUAL_ORDER_DT,_ONEDIMHAARCOUNT);
#endif
    
      
//    char ucoeffsfile[250];
//    sprintf(ucoeffsfile,"%su_coeffs_rhs%d_dec%d_off%d_%d_%d_hmax%d.m",path,_RHS_NUMBER,_DECOMP,_MAXOFFSET,_PRIMAL_ORDER_D,_DUAL_ORDER_DT,_ONEDIMHAARCOUNT);
    
            
            

    Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> > agencoeffs, qgencoeffs;
    agencoeffs.resize(num_of_patches);
    qgencoeffs.resize(num_of_patches);
    for (unsigned int p=0; p< num_of_patches; ++p)
    {
        agencoeffs[p] = FixedMatrix<double, _ONEDIMHAARCOUNT> (1.0);
        // qgencoeffs[p] = FixedMatrix<double, _ONEDIMHAARCOUNT> (0.0); // -> poisson!
        qgencoeffs[p] = FixedMatrix<double, _ONEDIMHAARCOUNT> (1.0); // -> helmholtz!
// CLEANUP        
        for (unsigned int i=0; i< _ONEDIMHAARCOUNT; ++i)
        {
            for (unsigned int j=0; j< _ONEDIMHAARCOUNT; ++j)
            {
                assert (agencoeffs[p].get_entry(i,j) == 1.0);
                assert (qgencoeffs[p].get_entry(i,j) == 1.0);
                //agencoeffs[p].set_entry(i,j,1.0);
                //qgencoeffs[p].set_entry(i,j,1.0);
            }
        }
        /*
         */
    }
    cout << "done initializing a,q" << endl;

#if _COMPUTE_RHS
    TestRHS<DIM, _RHS_NUMBER> f;
    
    // uncoment for constant RHS ...
    //ConstantFunction<DIM> constant_rhs(Vector<double>(1, "2.0"));
    //PoissonBVP<DIM> poisson(&constant_rhs);

    // uncomment in 2 space dimensions
    //TestRHS<2> rhs;
    //PoissonBVP<dim> poisson(&rhs);

    
    //InfiniteVector<double, int> acoeffs,qcoeffs;
    
    //TestRHS<0> f; // zero
    //TestRHS<4> f; // sin * sin
    //TestRHS<5> f; // 3.0
    //TestRHS<6> f; // 3.0 on (0,1)^2
    

    
    //CachedQTProblem<Basis,_ONEDIMHAARCOUNT> qtproblem();
    //CachedQTProblem<Basis> qtproblem(&qtbasis,acoeffs,qcoeffs,&f,1.0,1.0);
    
    cout << "done initializing f" << endl;
    CachedQTProblem<Basis,_ONEDIMHAARCOUNT> cqtproblem (&qtbasis,
            agencoeffs,
            qgencoeffs,
            &f,
            rhs_filename,
            5.0,
            10.0);
#else
    Function<DIM>* f = 0;
    CachedQTProblem<Basis,_ONEDIMHAARCOUNT> cqtproblem (&qtbasis,
            agencoeffs,
            qgencoeffs,
            f,
            rhs_filename,
            5.0,
            10.0);
    cout << "rhs is not initialized as a function. Rather the coeffs are loaded from disc. Filename = " << rhs_filename << endl;
#endif
    
// CLEANUP

    for (unsigned int i=0; i < qtbasis.degrees_of_freedom(); ++i)
    {
        cout << "i = " << i << " = " << *qtbasis.get_wavelet(i) << "; F_precond[" << i << "] = " << cqtproblem.f(i) << endl;
    }
    
    //abort();
    
    cout << "done initializing cached qt problem" << endl;
    
#if ( (_COMPARE_QT_a_AND_T_a || _CALLS_TO_CQTPROBLEM_f) || ( _COMPARE_QT_RHS_AND_T_RHS ) )
    cout << "for comparison with tbasis: " << endl << "initialize poisson problems and cached t problems on each patch" << endl;
    
    //TBasis basisBC(bc)
    typedef TensorBasis<Basis1d,DIM> TBasis;
    typedef TBasis::Index TIndex;
    typedef TensorEquation<Basis1d,DIM,TBasis> TEquation;
    typedef CachedTProblem<TEquation> CTEquation;
    /*
    FixedArray1D<TBasis*,num_of_patches> TBasisArray;
    for (unsigned int p=0; p<num_of_patches; ++p)
    {
        *TBasisArray[p] = TBasis(bc_bool[p]);
        TBasisArray[p]->set_jmax(multi_degree(qtbasis.j0()[0])+maxleveloffset);
    }
     */
    
    ConstantFunction<DIM> constant_rhs(Vector<double>(1, "2.0"));
    PoissonBVP<DIM> poisson(&constant_rhs);
    
    cout << "done initializing poisson problem with constant rhs" << endl;
    FixedArray1D<TEquation*,num_of_patches> TEqArray;
    FixedArray1D<CTEquation*,num_of_patches> CTEqArray;
    for (unsigned int p=0; p<num_of_patches; ++p)
    {
        TEqArray[p] = new TEquation(&poisson, bc_bool[p]);
        TEqArray[p]->set_jmax(multi_degree(qtbasis.j0()[0])+maxleveloffset);
        CTEqArray[p] = new CTEquation(TEqArray[p],5.5,35.0);
    }
    
    cout << "done initializing cached t equations on each patch" << endl;
#else
    cout << "no comparison with TBasis routines" << endl;
#endif    
    //TensorEquation<Basis1d,dim,TBasis> eq(&poisson, bc);
    //eq.set_jmax(multi_degree(eq.basis().j0())+offset); //calls setup_full_collection
  
    //CachedTProblem<TensorEquation<Basis1d,dim,TBasis> > ctproblem(&eq,5.5,35.0);

    MultiIndex<int,DIM> level_j;
    double temp_d1, temp_d2, temp_d3, temp_d4, maxdif;
    unsigned int level_p, maxdif_nunum, maxdif_munum;
    
    typedef CachedQTProblem<Basis,_ONEDIMHAARCOUNT>::ColumnCache ColumnCache;
    typedef CachedQTProblem<Basis,_ONEDIMHAARCOUNT>::Column Column;
    typedef CachedQTProblem<Basis,_ONEDIMHAARCOUNT>::Block Block;
    
    ColumnCache* cachePointer;
    Column* cachePointerGen;
    



#if _COMPARE_QT_a_AND_T_a
    cout << "comparing qtbasis::a and tbasis::a" <<endl;
    //cout << "(levelnum, nunum) =";
    /*
    cout << "qtbasis.get_numberoflevels() = " << qtbasis.get_numberoflevels() << endl;
    for (unsigned int levelnum = 0; levelnum < qtbasis.get_numberoflevels(); ++levelnum)
    {
        cout << "levelnum == " << levelnum << endl;
        Index first_wav = qtbasis.first_wavelet(levelnum);
        //Index temp_ind = qtbasis.first_wavelet(levelnum);
        cout << "first_wav = " << first_wav << endl;
        qtbasis.get_level(levelnum,level_j,level_p);
        Index last_wav = qtbasis.last_wavelet(levelnum);
        cout << "last_wav = " << last_wav << endl;
    }
    */
    for (unsigned int levelnum = 0; levelnum < qtbasis.get_numberoflevels(); ++levelnum)
    {
        cout << "\nlevelnum == " << levelnum << endl;
        Index first_wav = qtbasis.first_wavelet(levelnum);
        //Index temp_ind = qtbasis.first_wavelet(levelnum);
        qtbasis.get_level(levelnum,level_j,level_p);
        
        Index last_wav = qtbasis.last_wavelet(levelnum);

        for (unsigned int nunum = first_wav.number(); nunum <= last_wav.number(); ++nunum)
        {
            temp_nu = qtbasis.get_wavelet(nunum);
            //cout << "nunum = " << nunum << "; temp_nu = " << temp_nu << endl;
            cout << " (" << levelnum << "," << nunum << ")";
            TIndex temp_ti1(temp_nu.j(),temp_nu.e(), temp_nu.k(), &(CTEqArray[first_wav.p()]->basis()));
            for (unsigned int munum = first_wav.number(); munum <= last_wav.number(); ++munum)
            {
                temp_mu = qtbasis.get_wavelet(munum);
                //cout << "nunum = " << nunum << "; temp_nu = " << temp_nu << "; munum = " << munum << "; temp_mu = " << temp_mu << endl;
                TIndex temp_ti2(temp_mu.j(),temp_mu.e(), temp_mu.k(), & (CTEqArray[first_wav.p()]->basis()));
                
                temp_d1 = cqtproblem.a(munum,nunum);
                temp_d2 = CTEqArray[first_wav.p()]->a(temp_ti2,temp_ti1);
                temp_d3 = cqtproblem.a(munum,nunum);
                temp_d4 = CTEqArray[first_wav.p()]->a(temp_ti2,temp_ti1);
                assert (temp_d1 == temp_d3);
                assert (temp_d2 == temp_d4);
                // if nu and mu are both extended in one direction, temp_d2 is too small!
                // check LMR info and multiply temp_d2 *2 in each extension direction
                temp_b = qtbasis.get_LMR_info(nunum, level_j, level_p, intinfo);
                // include information about mu into LMR info
                // -> store update in intinfo2
                for (unsigned int i=0; i<DIM; ++i)
                {
                    // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                    switch (intinfo[i])
                    {
                        case 0:
                            // LL: check whether mu is a left boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                                  ||
                                  ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw() ) ) ) )
                            {
                                // mu is not extended left
                                intinfo[i] = 4; // MM
                            }

                            break;
                        case 2:
                            // LR: check whether mu is a right boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                            {
                                intinfo[i] = 1; // LM
                            }
                            break;
                        case 3:
                            // ML : check whether mu is a left boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                                  ||
                                  ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw()) ) ) )
                            {
                                // mu is not extended left. There is no intersection of mu and nu!
                                if (temp_mu.e()[i] == 0)
                                {
                                    assert(temp_mu.k()[i] > kmaxgen[i]);
                                }
                                else
                                {
                                    assert(temp_mu.k()[i] < kmaxwav[i]);
                                }
                            }
                            break;
                        case 5:
                            // MR: check whether mu is a right boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                            {
                                // mu is not extended right. There is no intersection of mu and nu!
                                if (temp_mu.e()[i] == 0)
                                {
                                    assert(temp_mu.k()[i] < kmingen[i]);
                                }
                                else
                                {
                                    assert(temp_mu.k()[i] < kminwav[i]);
                                }
                            }
                            break;
                        case 6:
                            // RL : check whether mu is a left boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                                  ||
                                  ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] < qtbasis.get_numofbw()) ) ) )
                            {
                                // mu is not extended left. There is no intersection of mu and nu!
                                intinfo[i] = 7; // RM
                            }
                            break;
                        case 8:
                            // RR: check whether mu is a right boundary gen/wav
                            if (!( ( (temp_mu.e()[i] == 0) && (temp_mu.k()[i] == qtbasis.get_bases(temp_mu.p(),i)->DeltaRmax( qtbasis.get_j0(temp_mu.p(),i) )) ) 
                              ||
                              ( (temp_mu.e()[i] == 1) && (temp_mu.k()[i] > (1<<temp_mu.j()[i])- 1 - qtbasis.get_numofbw()) ) ))
                            {
                                // nu is not extended right
                                intinfo[i] = 4; // MM
                            }
                            break;
                        case 1:
                        case 4:
                        case 7:
                            break;
                        default:
                            abort();
                            break;
                    }
                }
                for (unsigned int i = 0; i<DIM; ++i)
                {
                    if ( ((intinfo[i]%2 == 0) && (intinfo[i] !=4)) )
                    // if (intinfo[i] != 4)
                    {
                        temp_d2 *= 2;
                    }
                }
                if (maxdif < abs(temp_d1-temp_d2))
                {
                    maxdif = abs(temp_d1-temp_d2);
                    maxdif_nunum = nunum;
                    maxdif_munum = munum;
                }
                if (abs (temp_d1-temp_d2) > 1e-6)
                {
                    cout << "error!" << endl;
                    cout << "levelnum = " << levelnum << endl;
                    cout << "level_j = " << level_j << endl;
                    cout << "level_p = " << level_p << endl;
                    cout << "nunum = " << nunum << endl;
                    cout << "temp_nu = " << temp_nu << endl;
                    cout << "temp_ti1 = " << temp_ti1 << endl;
                    cout << "munum = " << munum << endl;
                    cout << "temp_mu = " << temp_mu << endl;
                    cout << "temp_ti2 = " << temp_ti2 << endl;
                    cout << "temp_d1 = " << temp_d1 << endl;
                    cout << "temp_d2 = " << temp_d2 << endl;
                    cout << "difference = " << abs (temp_d1-temp_d2) << endl;
                    temp_d1 = cqtproblem.a(munum,nunum);
                    temp_d2 = CTEqArray[first_wav.p()]->a(temp_ti2,temp_ti1);
                    //assert (temp_d1 == temp_d2);
                }
            }
        }
    }
    cout << "maximal difference of cqtproblem.a and CTEqArray[p]->a occurred at" 
            << endl << "numum = " << maxdif_nunum 
            << "; munum = " << maxdif_munum 
            << "; difference = " << maxdif << endl;
    assert(maxdif < 1e-10);
    
#else
    cout << "skipping comparison of qtbasis::a and tbasis::a" << endl;
#endif
  
    
#if _SYMMETRY_OF_CACHED_A
    cout << "testing symmetry of cached a ..."<<endl;
    maxdif = 0;
//    cout << "n =";
    for (unsigned int n=0; n < qtbasis.degrees_of_freedom(); ++n)
    {
        //cout << "n= " << n << endl;
//        cout << " " << n;
        for (unsigned int m=0; m<qtbasis.degrees_of_freedom(); ++m)
        {
            if (m%500 == 0)
            {
                cout << "n = " << n << "; m= " << m << endl;
            }
            temp_d1 = cqtproblem.a(m,n);
            temp_d2 = cqtproblem.a(n,m);
            if (maxdif < abs(temp_d1-temp_d2))
            {
                maxdif = abs(temp_d1-temp_d2);
                maxdif_nunum = n;
                maxdif_munum = m;
                cout << "n = " << maxdif_nunum 
                        << "; m = " << maxdif_munum 
                        << "; difference = " << maxdif << endl;
                temp_d3 = cqtproblem.a(m,n);
                temp_d4 = cqtproblem.a(n,m);
                assert (temp_d1 == temp_d3);
                assert (temp_d2 == temp_d4);
                temp_d3 = cqtproblem.a(m,n);
                temp_d4 = cqtproblem.a(n,m);
            }
            //assert (temp_d1 == temp_d2);
        }
    }
    cout << "symmetry test for a()" << endl
            << "maximal difference occoured at" 
            << endl << "n = " << maxdif_nunum 
            << "; m = " << maxdif_munum 
            << "; difference = " << maxdif << endl;
    assert(maxdif < 1e-10);
  
#else
  cout << "skipping symmetry test for cached a" << endl;
#endif
  

#if _CALLS_TO_CQTPROBLEM_f
    cout << "test calls to f(lambdanum)" << endl;
    cout << "ensure that c qt problem and c t problem yield the same as long as lambda is not extended" << endl;
    cout << "comparison only works if f is independent of the corner! (because in tbasis the corner is assumed to be 0)" << endl;
    
    
    Point <DIM> temp_poi;
    for (unsigned int i = 0; i>DIM; ++i)
        temp_poi[i] = cqtproblem.basis()->get_corners(0)[i];
    temp_d1 = cqtproblem.f_->value(temp_poi);
    //cout << "x = " << temp_poi << "; temp_d1 = " << temp_d1 << endl;
        
    for (unsigned int p=1; p< cqtproblem.basis()->get_nop(); ++p)
    {
        for (unsigned int i = 0; i>DIM; ++i)
        {
            temp_poi[i] = cqtproblem.basis()->get_corners(p)[i];
        }
        temp_d2 = cqtproblem.f_->value(temp_poi);
        //cout << "x = " << temp_poi << "; temp_d2 = " << temp_d2 << endl;
        assert (cqtproblem.f_->value(temp_poi) == temp_d1);
    }
    maxdif = 0;
    maxdif_nunum = 0;
    for (unsigned int n = 64; n < cqtproblem.basis()->degrees_of_freedom(); ++n)
    {
        temp_d1 = cqtproblem.basis()->integrate(cqtproblem.f_, n);
        temp_nu = cqtproblem.basis()->get_wavelet(n);
        cout << "n = " << n << "; temp_nu = " << temp_nu << endl;
        temp_b = qtbasis.get_LMR_info(n, temp_nu.j(), temp_nu.p(), intinfo);
        assert(temp_b);
        exit = (intinfo[0] != 4);
        for (unsigned int i=1; i<DIM; ++i)
        {
            exit = exit || (intinfo[i] != 4);
        }
           
        if (exit)
        {
            // wavelet is extended. in general: ctproblem.f != cqtproblem.f ->continue!
            //continue; 
            // for constant RHS: 2 * ctproblem.f == cqtproblem.f (reflection!)
            temp_d1 /= 2;
        }
        TIndex temp_ti(temp_nu.j(),temp_nu.e(), temp_nu.k(), &(CTEqArray[temp_nu.p()]->basis()));
        temp_d2 = CTEqArray[temp_nu.p()]->basis().integrate(cqtproblem.f_, temp_ti);
        
        if (maxdif < abs(temp_d1-temp_d2))
        {
            maxdif = abs(temp_d1-temp_d2);
            maxdif_nunum = n;
            cout << "n= " << n <<"; temp_nu = " << temp_nu << "; temp_ti = " << temp_ti << endl;
            cout << "temp_d1 = " << temp_d1 << "; temp_d2 = " << temp_d2 << "; difference = " << abs(temp_d1-temp_d2) << endl;
            temp_d1 = cqtproblem.basis()->integrate(cqtproblem.f_, n);
            temp_d2 = CTEqArray[temp_nu.p()]->basis().integrate(cqtproblem.f_, temp_ti);
        }
    }
    cout << "maximal difference of cqtproblem.f and CTEqArray[p]->f occurred at" 
        << endl << "n = " << maxdif_nunum 
        << "; difference = " << maxdif << endl;
    assert(maxdif < 1e-10);
    
    for (unsigned int n=0; n<cqtproblem.basis()->degrees_of_freedom(); n++)
    {
        temp_d1 = cqtproblem.basis()->integrate(cqtproblem.f_, n);
    }
    
#else
    cout << "skip tests to f(lambdanum)" << endl;
#endif
  
 
#if _TEST_NORM_A_AND_NORM_A_INV 
    /*
    MultiIndex<int,DIM> temp_j, temp_e, temp_k;
    temp_j[0] = temp_j[1] = 3;
    temp_e[0] = temp_e[1] = 1;
    temp_k[0] = temp_e[1] = 6;
    Index temp_ind1(temp_j,temp_e,temp_k,3);
    temp_e[0] = temp_e[1] = 0;
    temp_k[0] = temp_e[1] = 3;
    Index temp_ind2(temp_j,temp_e,temp_k,0);
    cqtproblem.a(temp_ind2.number(), temp_ind1.number());
    cqtproblem.a(temp_ind2.number(), temp_ind1.number());
     */
    /*
    typedef std::map<int, entries> Block;
    typedef std::map<int, Block> Column;
    typedef std::map<int, Column> ColumnCache;
    typedef FixedArray1D<ColumnCache,16> typeIcache;
    typedef FixedArray1D<ColumnCache,12> typeIIcache;
    typedef FixedArray1D<ColumnCache,12> typeIIIcache;
    typedef FixedArray1D<Column,16> typeIcachedGens; // there is only one level with generators, so we can omit one "map"
    typedef FixedArray1D<Column,12> typeIIcachedGens;
    typedef FixedArray1D<Column,12> typeIIIcachedGens;
    typeIcache typeIcache_;
    typeIIcache typeIIcache_;
    typeIIIcache typeIIIcache_;
    typeIcachedGens typeIcachedGens_;
    typeIIcachedGens typeIIcachedGens_;
    typeIIIcachedGens typeIIIcachedGens_;
      */  
 
    cout << "testing norm_A, norm_Ainv, normtest" << endl;
    temp_d1 = cqtproblem.norm_A();
    temp_d2 = cqtproblem.norm_Ainv();
    cout << "norm_A = " << temp_d1 << "; norm_A_inv = " << temp_d2 << "; Kondition = " << (temp_d1*temp_d2) << endl;
    for (unsigned int offset = 0; offset <= maxleveloffset; ++offset)
    {
        cqtproblem.normtest(offset);
    }
#else
    cout << "skipping testing of norm_A"<<endl;
#endif

        
#if _CHECK_CACHED_A_BY_HAND
    cout << "call a() for 2 distinct wavelets and check the cached values by hand" << endl;
    int nunum(cqtproblem.basis()->degrees_of_freedom()-1), munum(cqtproblem.basis()->degrees_of_freedom()-1);
    //int nunum(1), munum(1);
    temp_nu = cqtproblem.basis()->get_wavelet(nunum);
    temp_mu = cqtproblem.basis()->get_wavelet(munum);
    int nu_p(temp_nu.p()), mu_p(temp_mu.p());
    MultiIndex<int,DIM> mu_j(temp_mu.j());
    MultiIndex<int,DIM> mu_e(temp_mu.e());
    MultiIndex<int,DIM> mu_k(temp_mu.k());
    MultiIndex<int,DIM> nu_j(temp_nu.j());
    MultiIndex<int,DIM> nu_e(temp_nu.e());
    MultiIndex<int,DIM> nu_k(temp_nu.k());
    
    cout << "Indices of interest:" << endl;
    
    cout << "nu = " << *cqtproblem.basis()->get_wavelet(nunum) << "; mu = " << *cqtproblem.basis()->get_wavelet(munum) << endl;
    
    cout << "look at the cache!" << endl;
   
    temp_b = cqtproblem.basis()->get_LMR_info(nunum,mu_j,mu_p, intinfo);
    for (unsigned int i=0; i<DIM; ++i)
    {
        nui_basisnum = (((cqtproblem.basis()->get_bc()[nu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[nu_p][2*i+1])?0:1));
        mui_basisnum = (((cqtproblem.basis()->get_bc()[mu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[mu_p][2*i+1])?0:1));
        bool mu_min_type_i = (mu_j[i]==cqtproblem.basis()->j0()[mu_p][i]) ? true:false;
        switch (intinfo[i])
        {
            case 0:
            case 4:
            case 8: // LL MM RR -> type 1
                cachePointer = & cqtproblem.typeIcache_[4*nui_basisnum + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIcachedGens_[4*nui_basisnum + mui_basisnum];
                break;
            case 1:
            case 7: 
            case 2:
            case 6: // LM RM LR RL -> type 2
                cachePointer = &cqtproblem.typeIIcache_[4*(nui_basisnum-1) + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIcachedGens_[4*(nui_basisnum-1) + mui_basisnum];
                break;
            case 3:
            case 5: // ML MR -> type 3
                cachePointer = &cqtproblem.typeIIIcache_[3*nui_basisnum + (mui_basisnum-1)];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIIcachedGens_[3*nui_basisnum + (mui_basisnum-1)];
                break;
            default:
                abort();
        }
        
        
        Basis::IntervalBasis::Index nui(nu_j[i],nu_e[i],nu_k[i],cqtproblem.basis()->get_bases_infact()[nui_basisnum]);
        unsigned int nui_num = nui.number();
        ColumnCache::iterator col_lb(cachePointer->lower_bound(nui_num));
        ColumnCache::iterator col_it(col_lb);
        if (col_lb == cachePointer->end() ||
            cachePointer->key_comp()(nui_num, col_lb->first))
        {
            cout << "the nui-th column has not been requested so far" << endl;
            //typedef typename ColumnCache::value_type value_type;
            //col_it = cachePointer->insert(col_lb, value_type(nui_num, Column()));
        }
        else
        {
            cout << "col_it points to the column of \nu_i" << endl;
            Column& col(col_it->second);
            int mui_levelnum (mu_j[i] - cqtproblem.basis()->j0()[mu_p][i]);
            Column::iterator lb(col.lower_bound(mui_levelnum));
            Column::iterator it(lb);

            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                cout << "no entries have ever been computed for this column and this level (wavCache)" << endl;
            }
            else
            {
                cout << "column and levelblock already exist (wavCache)" << endl;
                if (mu_e[i] == 0)
                {
                    Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    if (lb2 == cachePointerGen->end() ||
                        cachePointer->key_comp()(nui_num, lb2->first))
                    {
                        cout << "the nui-th column has not been requested so far (genCache). this should not happen!" << endl;
                        abort();
                    }
                    Block& block(lb2->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // 
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
                else
                {
                    Block& block(it->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // level exists, but in row 'mui_k' no entry is available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
            }
        }
    }
    cout << "done looking at the cache" << endl;
    cout << "calling a with:" << endl;
    cout << "nu = " << *cqtproblem.basis()->get_wavelet(nunum) << "; mu = " << *cqtproblem.basis()->get_wavelet(munum) << endl;
    double temp_d (cqtproblem.a(munum, nunum));
    cout << "yields " << temp_d << endl;
    cout << "look at the cache!" << endl;
    temp_b = cqtproblem.basis()->get_LMR_info(nunum,mu_j,mu_p, intinfo);
    for (unsigned int i=0; i<DIM; ++i)
    {
        nui_basisnum = (((cqtproblem.basis()->get_bc()[nu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[nu_p][2*i+1])?0:1));
        mui_basisnum = (((cqtproblem.basis()->get_bc()[mu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[mu_p][2*i+1])?0:1));
        bool mu_min_type_i = (mu_j[i]==cqtproblem.basis()->j0()[mu_p][i]) ? true:false;
        switch (intinfo[i])
        {
            case 0:
            case 4:
            case 8: // LL MM RR -> type 1
                cachePointer = & cqtproblem.typeIcache_[4*nui_basisnum + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIcachedGens_[4*nui_basisnum + mui_basisnum];
                break;
            case 1:
            case 7: 
            case 2:
            case 6: // LM RM LR RL -> type 2
                cachePointer = &cqtproblem.typeIIcache_[4*(nui_basisnum-1) + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIcachedGens_[4*(nui_basisnum-1) + mui_basisnum];
                break;
            case 3:
            case 5: // ML MR -> type 3
                cachePointer = &cqtproblem.typeIIIcache_[3*nui_basisnum + (mui_basisnum-1)];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIIcachedGens_[3*nui_basisnum + (mui_basisnum-1)];
                break;
            default:
                abort();
        }
        
        
        Basis::IntervalBasis::Index nui(nu_j[i],nu_e[i],nu_k[i],cqtproblem.basis()->get_bases_infact()[nui_basisnum]);
        unsigned int nui_num = nui.number();
        ColumnCache::iterator col_lb(cachePointer->lower_bound(nui_num));
        ColumnCache::iterator col_it(col_lb);
        if (col_lb == cachePointer->end() ||
            cachePointer->key_comp()(nui_num, col_lb->first))
        {
            cout << "the nui-th column has not been requested so far" << endl;
            //typedef typename ColumnCache::value_type value_type;
            //col_it = cachePointer->insert(col_lb, value_type(nui_num, Column()));
        }
        else
        {
            cout << "col_it points to the column of nu_i" << endl;
            Column& col(col_it->second);
            int mui_levelnum (mu_j[i] - cqtproblem.basis()->j0()[mu_p][i]);
            Column::iterator lb(col.lower_bound(mui_levelnum));
            Column::iterator it(lb);

            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                cout << "no entries have ever been computed for this column and this level (wavCache)" << endl;
            }
            else
            {
                cout << "column and levelblock already exist (wavCache)" << endl;
                if (mu_e[i] == 0)
                {
                    Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    if (lb2 == cachePointerGen->end() ||
                        cachePointer->key_comp()(nui_num, lb2->first))
                    {
                        cout << "the nui-th column has not been requested so far (genCache). this should not happen!" << endl;
                        abort();
                    }
                    Block& block(lb2->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // 
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
                else
                {
                    Block& block(it->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // level exists, but in row 'mui_k' no entry is available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
            }
        }
    }
    cout << "done looking at the cache" << endl;
    cout << "Calling a again!" << endl;
    cqtproblem.a(munum, nunum);
    cout << "look at the cache!" << endl;
    temp_b = cqtproblem.basis()->get_LMR_info(nunum,mu_j,mu_p, intinfo);
    for (unsigned int i=0; i<DIM; ++i)
    {
        nui_basisnum = (((cqtproblem.basis()->get_bc()[nu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[nu_p][2*i+1])?0:1));
        mui_basisnum = (((cqtproblem.basis()->get_bc()[mu_p][2*i])?0:2) + ((cqtproblem.basis()->get_bc()[mu_p][2*i+1])?0:1));
        bool mu_min_type_i = (mu_j[i]==cqtproblem.basis()->j0()[mu_p][i]) ? true:false;
        switch (intinfo[i])
        {
            case 0:
            case 4:
            case 8: // LL MM RR -> type 1
                cachePointer = & cqtproblem.typeIcache_[4*nui_basisnum + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIcachedGens_[4*nui_basisnum + mui_basisnum];
                break;
            case 1:
            case 7: 
            case 2:
            case 6: // LM RM LR RL -> type 2
                cachePointer = &cqtproblem.typeIIcache_[4*(nui_basisnum-1) + mui_basisnum];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIcachedGens_[4*(nui_basisnum-1) + mui_basisnum];
                break;
            case 3:
            case 5: // ML MR -> type 3
                cachePointer = &cqtproblem.typeIIIcache_[3*nui_basisnum + (mui_basisnum-1)];
                if (mu_min_type_i)
                    cachePointerGen = & cqtproblem.typeIIIcachedGens_[3*nui_basisnum + (mui_basisnum-1)];
                break;
            default:
                abort();
        }
        
        
        Basis::IntervalBasis::Index nui(nu_j[i],nu_e[i],nu_k[i],cqtproblem.basis()->get_bases_infact()[nui_basisnum]);
        unsigned int nui_num = nui.number();
        ColumnCache::iterator col_lb(cachePointer->lower_bound(nui_num));
        ColumnCache::iterator col_it(col_lb);
        if (col_lb == cachePointer->end() ||
            cachePointer->key_comp()(nui_num, col_lb->first))
        {
            cout << "the nui-th column has not been requested so far" << endl;
            //typedef typename ColumnCache::value_type value_type;
            //col_it = cachePointer->insert(col_lb, value_type(nui_num, Column()));
        }
        else
        {
            cout << "col_it points to the column of \nu_i" << endl;
            Column& col(col_it->second);
            int mui_levelnum (mu_j[i] - cqtproblem.basis()->j0()[mu_p][i]);
            Column::iterator lb(col.lower_bound(mui_levelnum));
            Column::iterator it(lb);

            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                cout << "no entries have ever been computed for this column and this level (wavCache)" << endl;
            }
            else
            {
                cout << "column and levelblock already exist (wavCache)" << endl;
                if (mu_e[i] == 0)
                {
                    Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    if (lb2 == cachePointerGen->end() ||
                        cachePointer->key_comp()(nui_num, lb2->first))
                    {
                        cout << "the nui-th column has not been requested so far (genCache). this should not happen!" << endl;
                        abort();
                    }
                    Block& block(lb2->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // 
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
                else
                {
                    Block& block(it->second);
                    Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    Block::iterator block_it(block_lb);
                    // level exists, but in row 'mui_k' no entry is available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        cout << "intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero" << endl;
                    }
                    else 
                    {
                        cout << "'mui_k' is available" << endl;
                        cout << "entries.first = " << (block_it->second).first << endl;
                        cout << "entries.second = " << (block_it->second).second << endl;
                    }
                }
            }
        }
    }
    cout << "done looking at the cache" << endl;
#else
    cout << "do not call a() for 2 distinct wavelets and check the cached values by hand" << endl;
#endif
  
#if _TEST_ADDBALL
    cout << "testing add_ball"<<endl;
        
    /*
    typedef std::pair<FixedArray1D<double,_ONEDIMHAARCOUNT>, FixedArray1D<double,_ONEDIMHAARCOUNT> > entries;
    typedef map<int, entries> Block;
    Block temp_B;
    Block* testmap = & temp_B;
    FixedArray1D<double,_ONEDIMHAARCOUNT> gram, der;
    typedef typename Block::value_type value_type_block;
    
    cout << "?? =  " << testmap->size() << endl;
    cout << "?? =  " << testmap->begin()->first << endl;
    cout << "?? =  " << testmap->rbegin()->first << endl;
    testmap->insert(testmap->end(), value_type_block(13, make_pair(gram,der)));
    cout << "?? =  " << testmap->size() << endl;
    cout << "?? =  " << testmap->begin()->first << endl;
    cout << "?? =  " << testmap->rbegin()->first << endl;
*/
    Vector<double> temp_w;
    temp_w.resize(cqtproblem.basis()->degrees_of_freedom());
// cleanup
    for (unsigned int m=0; m< temp_w.size(); ++m)
    {
        assert (temp_w[m] == 0);
    }
    Index add_ball_nu, add_ball_mu;
    unsigned int maxdif_n(0), maxdif_m(0);
    double add_ball_maxdif(1e-14);
    
    unsigned int temp_i;
    double factor = 1.7;
    bool precond(true);
    
//    int m=257, n=256;
//    add_ball_nu = cqtproblem.basis()->get_wavelet(n);
//    add_ball_mu = cqtproblem.basis()->get_wavelet(m);
//    cout << "N = " << n << "; add_ball_nu = " << add_ball_nu << "; M = " << m << "; add_ball_mu = " << add_ball_mu << endl;
//    double temp_d = cqtproblem.a(m,n);
//    cqtproblem.add_ball(n, temp_w, 1, factor, 99, tensor_simple, precond);
//    cout << "  cqtproblem.a(m,n) = " << cqtproblem.a(m,n) << "; cqtproblem.D(m) = " << cqtproblem.D(m) << "; cqtproblem.D(n) = " << cqtproblem.D(n) << endl;
//    cout << "  cqtproblem.a(m,m) = " << cqtproblem.a(m,m) << "; cqtproblem.a(n,n) = " << cqtproblem.a(n,n)  << endl;
//    cout << "  factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) = " << factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) << endl;
//    cout << "  temp_w[m] = " << temp_w[m] << endl;
//    if (precond)
//    {
//        cout << "  difference = " << abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) -temp_w[m]) << endl;
//        add_ball_maxdif = abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) -temp_w[m]);
//    }
//    else
//    {
//        cout << "  difference = " << abs(factor*cqtproblem.a(m,n) -temp_w[m]) << endl;
//        add_ball_maxdif = abs(factor*cqtproblem.a(m,n) -temp_w[m]);
//    }
//    if (precond)
//    {
//        assert (abs(temp_w[m] - factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) )<1e-14);
//    }
//    else
//    {
//        assert (abs(temp_w[m] - factor*cqtproblem.a(m,n))<1e-14);
//    }
        
    
    add_ball_maxdif = 1e-15; //
    for (unsigned int levelrange = 0; levelrange <= cqtproblem.basis()->get_jmax() - multi_degree(cqtproblem.basis()->first_wavelet(0).j()); ++ levelrange)
    {
        for (unsigned int n = 0; n < cqtproblem.basis()->degrees_of_freedom(); ++n)
        {
            cout << "N = " << n << "; levelrange = " << levelrange << endl;
            add_ball_nu = cqtproblem.basis()->get_wavelet(n);
            temp_w = 0;
            cqtproblem.add_ball(n, temp_w, levelrange, factor, 99, tensor_simple, precond);
            
            /*
            for (unsigned int m=0; m< temp_w.size(); ++m)
            {
                if (temp_w[m] != 0)
                    cout << "  w["<<m<<"] = " << temp_w[m] << endl;
            }
            */
            
            for (unsigned int m = 0; m < cqtproblem.basis()->degrees_of_freedom(); ++m)
            {
                add_ball_mu = cqtproblem.basis()->get_wavelet(m);
//                if (m == M_)
//                {
//                    cout << "closer look at m = " << m << "; add_ball_mu = " << add_ball_mu << endl;
//                }
                
                //cout << "N = " << n << "; M = " << m << endl;
                temp_i = 0;
                for (unsigned int i=0; i<DIM; ++i)
                {
                    temp_i += abs (add_ball_nu.j()[i]-add_ball_mu.j()[i]);
                }
                if (temp_i <= levelrange)
                {
                    //if 
                    if (precond ? (add_ball_maxdif < abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) -temp_w[m])) 
                                : (add_ball_maxdif < abs(factor*cqtproblem.a(m,n) -temp_w[m])) )
                    {
                        cout << "N = " << n << "; add_ball_nu = " << add_ball_nu << "; M = " << m << "; add_ball_mu = " << add_ball_mu << endl;
                        cout << "  cqtproblem.a(m,n) = " << cqtproblem.a(m,n) << "; cqtproblem.D(m) = " << cqtproblem.D(m) << "; cqtproblem.D(n) = " << cqtproblem.D(n) << endl;
                        cout << "  cqtproblem.a(m,m) = " << cqtproblem.a(m,m) << "; cqtproblem.a(n,n) = " << cqtproblem.a(n,n)  << endl;
                        cout << "  factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) = " << factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) << endl;
                        cout << "  temp_w[m] = " << temp_w[m] << endl;
                        if (precond)
                        {
                            cout << "  difference = " << abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) -temp_w[m]) << endl;
                            cout << "  quotient f*a(" << m << ", " << n << ")/D(" << m<< ")/D(" << n << ") / temp_w[" << m << "] = " << abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n)/temp_w[m]) << endl;
                            add_ball_maxdif = abs(factor*cqtproblem.a(m,n)/cqtproblem.D(m)/cqtproblem.D(n) -temp_w[m]);
                        }
                        else
                        {
                            cout << "  difference = " << abs(factor*cqtproblem.a(m,n) -temp_w[m]) << endl;
                            cout << "  quotient f*a(" << m << ", " << n << ") / temp_w[" << m << "] = " << abs(factor*cqtproblem.a(m,n)/temp_w[m]) << endl;
                            add_ball_maxdif = abs(factor*cqtproblem.a(m,n) -temp_w[m]);
                        }
                        maxdif_n = n;
                        maxdif_m = m;
                        //assert (false);
                        //abort();
                    }
                }
                else
                {
                    assert (temp_w[m] == 0);
                }
            }
            
        }
    }
//                    cout << "error!" << endl;
//                    cout << "levelnum = " << levelnum << endl;
//                    cout << "level_j = " << level_j << endl;
//                    cout << "level_p = " << level_p << endl;
//                    cout << "nunum = " << nunum << endl;
//                    cout << "temp_nu = " << temp_nu << endl;
//                    cout << "temp_ti1 = " << temp_ti1 << endl;
//                    cout << "munum = " << munum << endl;
//                    cout << "temp_mu = " << temp_mu << endl;
//                    cout << "temp_ti2 = " << temp_ti2 << endl;
//                    cout << "temp_d1 = " << temp_d1 << endl;
//                    cout << "temp_d2 = " << temp_d2 << endl;
//                    cout << "difference = " << abs (temp_d1-temp_d2) << endl;
//                    temp_d1 = cqtproblem.a(munum,nunum);
//                    temp_d2 = CTEqArray[first_wav.p()]->a(temp_ti2,temp_ti1);
    
    cout << "maximal difference of cqtproblem.add_ball and cqtproblem.a occurred at" 
            << endl << "n = " << maxdif_n 
            << "; m = " << maxdif_m 
            << "; difference = " << add_ball_maxdif << endl;
    assert(add_ball_maxdif < 1e-10);
#else
    cout << "skipping test of add_ball"<<endl;
#endif
  
#if _TEST_CDD2
    cout << "Test setup and solution with CDD2_SOLVE" << endl;
    
    tstart = clock();
    
    assert (cqtproblem.basis()->first_generator() == cqtproblem.basis()->first_wavelet(0));
    
//    set<Index> Lambda_galerkin;
//    for (Index lambda ( cqtproblem.basis()->first_generator()), itend(cqtproblem.basis()->get_wavelet(cqtproblem.basis()->degrees_of_freedom()));; ++lambda) 
//    {
//        Lambda_galerkin.insert(lambda);
//        if (lambda == itend) 
//            break;
//    }
    set<int> Lambda_galerkin;
    //for (unsigned int n=0; n<cqtproblem.basis()->degrees_of_freedom();++n)
    for (unsigned int n=0; n<15;++n)
    {
        Lambda_galerkin.insert(n);
    }
    SparseMatrix<double> A_galerkin;
    //A_apply.resize(Lambda_apply.size(), Lambda_apply.size());
    cout << "degrees_of_freedom = " << cqtproblem.basis()->degrees_of_freedom() << endl;
    cout << "main:: setting up stiffness matrix ..." << endl;
    setup_stiffness_matrix(cqtproblem, Lambda_galerkin, A_galerkin);
    cout << "  done" << endl;
    
    //cout << "A_galerkin = " << endl << A_galerkin << endl;

    cout << "main:: setting up RHS" << endl;
    InfiniteVector<double, Index> F_eta;
    cqtproblem.RHS(1e-6, F_eta);
    cout << "l2norm(f) = " << l2_norm(F_eta) << endl;
    cout << "main:: calling norm_Ainv" << endl;
    //const double nAinv = cqtproblem.norm_Ainv();
    //cout << "main:: calling norm_A" << endl;
    //const double nA = ctproblem.norm_A();
    //cout << "main:: calling l2_norm" << endl;
    //const double l2normF = l2_norm(F_eta);
    //cout << "l2norm(f_eta) = " << l2normF << " normA = "<< nA << " normAinv = " << nAinv << " kond = "<<(nA*nAinv)<<endl;
    double nu = cqtproblem.norm_Ainv() * l2_norm(F_eta);
    cout << "nu = " << nu << endl;

    InfiniteVector<double, int> u_epsilon;
    cout << "calling CDD2_SOLVE_TENSOR" << endl;
    //unsigned int maxlevel = ctproblem.basis().get_jmax();

    //nu = min(nu,2.0);
    //nu = 3000;
    
    InfiniteVector<double,int> ff, v, Av;
    cqtproblem.RHS(1e-6, ff);
    APPLY_TENSOR(cqtproblem, v, 1e-4, Av, 99, tensor_simple, true);
    cout << "Number of degrees of freedom " << Av.size() << endl;
    cout << "current residual error ||f-Av||=" << l2_norm(ff - Av) << endl;
    v += 0.5 * (ff - Av);
    cout << "v.size() = " << v.size() << endl;
    APPLY_TENSOR(cqtproblem, v, 1e-4, Av, 99, tensor_simple, true);
    cout << "Number of degrees of freedom " << Av.size() << endl;
    cout << "current residual error ||f-Av||=" << l2_norm(ff - Av) << endl;
    v += 0.5 * (ff - Av);
    cout << "v.size() = " << v.size() << endl;
    
    CDD2_SOLVE(cqtproblem, nu, 1e-4, u_epsilon,99);
// TODO    
//    cout << "  ... done. Saving output in file 'u_adaptCDD2.m'" << endl;
//    cout << "main:: u_epsilon.size() = " << u_epsilon.size() << endl;
//    u_epsilon.scale(&cqtproblem, -1);
//    cout << "main:: u_epsilon.size() = " << u_epsilon.size() << endl;
//    SampledMapping<DIM> s(evaluate(cqtproblem.basis(), u_epsilon, true, d+DIM+maxleveloffset+1));
//    std::ofstream u_stream("u_adaptCDD2.m");
//    s.matlab_output(u_stream);
//    u_stream.close();
//    //  cout << "  ... done" << endl;
//    cout << "  ... done. Coarsening with tol = 1e-6. Saving output in file 'u_adaptCDD2c.m'" << endl;
//
//    InfiniteVector<double, Index> u_epsilon_coarse;
//    u_epsilon.COARSE(1.0e-6, u_epsilon_coarse);
//    cout << "main:: u_epsilon_coarse.size() = " << u_epsilon_coarse.size() << endl;    
//    SampledMapping<DIM> s2(evaluate(eq.basis(), u_epsilon_coarse, true, d+DIM+maxleveloffset+1));
//    std::ofstream u_stream2("u_adaptCDD2c.m");
//    s2.matlab_output(u_stream2);
//    u_stream2.close();
//    cout << "  ... done" << endl;
//
//    if (DIM == 2)
//    {
//        std::ofstream plotstream;
//        plotstream.open("coefficient_plotCDD2.m");
//        cout << "Writing active coefficients to coefficient_plotCDD2.m" << endl;
//        plot_indices(cqtproblem.basis(), u_epsilon, maxleveloffset, plotstream, "jet", true, true);
//        plotstream.close();
//    }

    //cout << "u_epsilon = " << u_epsilon << endl;

    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
    cout << "Skipping test of CDD2" << endl;
#endif
    
#if _COMPARE_QT_RHS_AND_T_RHS
    cout << "comparing qtbasis::RHS and tbasis::RHS" <<endl;
    
    TestRHS<5> f2; // 3.0
//    int res(8);
//    Grid<2> grid3(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<res, 1<<res);
//    SampledMapping<2> abbildung(grid3, f2);
//    ofstream datei("Constant_Function.m");
//    abbildung.matlab_output(datei);
//    datei.close();
    
    FixedArray1D<bool,2*DIM> bc;
    bc[0] = false;
    bc[1] = false;
    bc[2] = false;
    bc[3] = false;
    TBasis tbasis(bc);
    //Basis1d basis1d(false,false);
    //FixedArray1D<Basis1d*,2> bases1d;
    //bases1d[0] = bases1d[1] = &basis1d;
    //TBasis tbasis(bases1d);
    
    const int maxrange = 0;
    MultiIndex<int,2> jmax;
    jmax[0]=tbasis.bases()[0]->j0()+maxrange;
    jmax[1]=tbasis.bases()[1]->j0();
    tbasis.set_jmax(jmax);
    InfiniteVector<double,TIndex> coeffs3;
    tbasis.expand(&f2, false, jmax, coeffs3);
    cout << coeffs3 << endl;
    SampledMapping<DIM> plot_patch2 (evaluate (tbasis, coeffs3, false,6)); // correct!
    plot_patch2.matlab_output(cout); // <-- spam, but helps
    //plot_patch2 = evaluate (tbasis, coeffs3, true,3); // incorrect!
    //plot_patch2.matlab_output(cout);
    std::ofstream plotstream_RHStest;
    plotstream_RHStest.open("plot_F_t2.m");
    plot_patch2.matlab_output(plotstream_RHStest);
    plotstream_RHStest.close();
    
    InfiniteVector<double, int> F_qt;
    cqtproblem.RHS(1e-6, F_qt);
    F_qt.scale(&cqtproblem,1);
    
    Array1D<InfiniteVector<double,int> > F_t;
    F_t.resize(cqtproblem.basis()->get_nop());
    for (unsigned int p=0; p< F_t.size(); ++p)
    {
        CTEqArray[p]->RHS(1e-6, F_t[p]); 
        // F_t[p].scale(CTEqArray[p],1); // <- only implemented for "Index"
        for (InfiniteVector<double,int>::const_iterator it(F_t[p].begin()), itend(F_t[p].end()); it!=itend; ++it)
        {
//            cout << "it.index() = " << it.index();
//            cout << "; CTEqArray[p]->D(CTEqArray[p]->basis().get_wavelet(it.index()) = " << CTEqArray[p]->D(CTEqArray[p]->basis().get_wavelet(it.index())) << endl;
            temp_d1 = *it * CTEqArray[p]->D(CTEqArray[p]->basis().get_wavelet(it.index()));
            F_t[p].set_coefficient(it.index(),temp_d1);
        }
    }
    
    
    cout << "plotting F_qt to disk" << endl;
    Array1D<SampledMapping<DIM> > plot_sample (qtbasis.sampled_output(F_qt,
                false,
                5));
    
    plotstream_RHStest.open("plot_F_qt.m");
    matlab_output(plot_sample, plotstream_RHStest);
    plotstream_RHStest.close();
    
    SampledMapping<DIM> plot_patch (evaluate (CTEqArray[0]->basis(), F_t[0], false,5));
    plotstream_RHStest.open("plot_F_t.m");
    plot_patch.matlab_output(plotstream_RHStest);
    plotstream_RHStest.close();
    
#else
    cout << "skipping comparison of qtbasis::RHS and tbasis::RHS" << endl;
#endif
    
    
#if _TEST_CDD1
    /*
    int kmingeni, kmaxgeni, kminwavi, kmaxwavi;
    qtbasis.get_onedim_intersections(3,
                        3,
                        0,
                        3,
                        0,
                        true,
                        3,
                        3,
                        kmingeni,
                        kmaxgeni,
                        kminwavi,
                        kmaxwavi,
                        gen_intersection_i,
                        wav_intersection_i);
    cout << "intinfoi = 3" << endl;
    cout << "kmingeni = " << kmingeni << "; "
            << "kmaxgeni = " << kmaxgeni << "; "
            << "kminwavi = " << kminwavi << "; "
            << "kmaxwavi = " << kmaxwavi << "; "
            << "gen_intersection_i = " << gen_intersection_i << "; "
            << "wav_intersection_i = " << wav_intersection_i << endl;
    
    qtbasis.get_onedim_intersections(5,
                        3,
                        0,
                        3,
                        0,
                        true,
                        3,
                        3,
                        kmingeni,
                        kmaxgeni,
                        kminwavi,
                        kmaxwavi,
                        gen_intersection_i,
                        wav_intersection_i);
    cout << "intinfoi = 5" << endl;
    cout << "kmingeni = " << kmingeni << "; "
            << "kmaxgeni = " << kmaxgeni << "; "
            << "kminwavi = " << kminwavi << "; "
            << "kmaxwavi = " << kmaxwavi << "; "
            << "gen_intersection_i = " << gen_intersection_i << "; "
            << "wav_intersection_i = " << wav_intersection_i << endl;
    
   */  
    cout << "Testing adaptive wavelet-Galerkin solution with CDD1_SOLVE ..." << endl;
    
    std::ofstream logstream;
    logstream.open(logfile);
           
    tstart = clock();
    cout << "main:: setting up RHS" << endl;
    InfiniteVector<double, int> F_cdd1;
    cqtproblem.RHS(1e-6, F_cdd1);
    cout << "l2norm(f) (preconditioned) = " << l2_norm(F_cdd1) << endl;
    
    logstream << "%[Degrees_of_freedom l2_Residual l2_Diag_times_residual]" << endl;
    logstream << "erg = [0 " << l2_norm(F_cdd1) << " ";
    
    F_cdd1.scale(&cqtproblem, 1);
    cout << "l2norm(f) (unpreconditioned) = " << l2_norm(F_cdd1) << endl;
    logstream << l2_norm(F_cdd1) << "];" << endl;

    
    cout << "degrees_of_freedom = " << cqtproblem.basis()->degrees_of_freedom() << endl;
    InfiniteVector<double, int> u_cdd1;
    cout << "calling CDD1_SOLVE" << endl;

    // CDD1_SOLVE(cqtproblem, 1e-1, u_cdd1,99,tensor_simple);
    CDD1_SOLVE_LOGGED(logstream, cqtproblem, 1e-6, u_cdd1,99,tensor_simple);
    cout << "  ... done. u_epsilon.size() = " << u_cdd1.size() << endl;

#if (_DIMENSION == 2)
    std::ofstream plotstream;
    
//    const int NN = 67;
//    cout << "n = " << NN << "; " << *cqtproblem.basis()->get_wavelet(NN) << endl;
//    Array1D<SampledMapping<DIM> > ergebnis (qtbasis.sampled_output(NN, 2.3,
//                true,
//                5));
//    cout << "plot " << NN << " to disk" << endl;
//    plotstream.open("plot_dummy_t.m");
//    matlab_output(ergebnis, plotstream);
//    plotstream.close();
//    ergebnis = qtbasis.sampled_output(NN, 2.3,
//                false,
//                5);
//    cout << "plot " << NN << " to disk" << endl;
//    plotstream.open("plot_dummy_f.m");
//    matlab_output(ergebnis, plotstream);
//    plotstream.close();
    cout << "plotting u_cdd1 to disk" << endl;
    Array1D<SampledMapping<DIM> > ergebnis2;
//    ergebnis2 = qtbasis.sampled_output(u_cdd1,
//                true,
//                6);
//    plotstream.open("plot_u_cdd1_t.m");
//    matlab_output(ergebnis2, plotstream);
//    plotstream.close();
//    cout << "done" << endl;
    
//    ergebnis2 = qtbasis.sampled_output(u_cdd1,
//                false,
//                6);
//    plotstream.open("plot_u_cdd1_f.m");
//    matlab_output(ergebnis2, plotstream);
//    plotstream.close();
//    cout << "done" << endl;
    u_cdd1.scale(&cqtproblem, -1);
    ergebnis2 = qtbasis.sampled_output(u_cdd1,
                true,
                6);
    
    
    plotstream.open(plotfile);
    matlab_output(ergebnis2, plotstream);
    plotstream.close();
    cout << "done" << endl;
    
//    ergebnis2 = qtbasis.sampled_output(u_cdd1,
//                false,
//                6);
//    plotstream.open("plot_u_cdd1_fs.m");
//    matlab_output(ergebnis2, plotstream);
//    plotstream.close();
//    cout << "done." << endl;
    
    
//
//    if (dim == 2)
//    {
//        std::ofstream plotstream2;
//        plotstream2.open("coefficient_plotCDD1.m");
//        cout << "Writing active coefficients to coefficient_plotCDD1.m" << endl;
//        plot_indices(cqtproblem.basis(), u_epsilon3, offset, plotstream2, "jet", true, true);
//        plotstream2.close();
//    }
//
//    cout << "saving computed solution to u_adaptCDD1.m" << endl;
//    SampledMapping<dim> s3(evaluate(eq.basis(), u_epsilon3, true, d+dim+offset+1));
//    std::ofstream u_stream3("u_adaptCDD1.m");
//    s3.matlab_output(u_stream3);
//    u_stream3.close();
//
//    cout << "  ... done. Coarsening u_epsilon with tol = 1e-6. " << endl;
//
//    InfiniteVector<double, int> u_cdd1_coarse;
//    u_cdd1.COARSE(1.0e-6, u_cdd1_coarse);
//    cout << "main:: u_cdd1_coarse.size() = " << u_cdd1_coarse.size() << endl;
//    cout << "Saving output in file 'u_adaptCDD1c.m'" << endl;
//    SampledMapping<dim> s4(evaluate(eq.basis(), u_cdd1_coarse, true, d+dim+offset+1));
//    std::ofstream u_stream4("u_adaptCDD1c.m");
//    s4.matlab_output(u_stream4);
//    u_stream4.close();
//    cout << "  ... done" << endl;
#endif
    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    
    
//    std::ofstream ucdd1stream;
//    ucdd1stream.open(ucoeffsfile);
//    ucdd1stream.close();
    

//    const int temp_iii(qtbasis.degrees_of_freedom());
//    FixedVector<double, temp_iii > u_cdd1_vector;
//    transform_IVTofixed(u_cdd1, u_cdd1_vector);
//    fixed_vector_writeToFile(ucoeffsfile,u_cdd1_vector);
//    cout << "u_cdd1 (coeffs) has been written to " << ucoeffsfile << endl;
    
#if (_DIMENSION == 2)
    cout << "Plots have been written to the plotfile:" << endl << plotfile << endl;
#endif
    cout << "time needed in total = " << time << " seconds" << endl;
    logstream << "time = " << time << ";" << endl;
    logstream.close();
    cout << "Output has been written to the logfile:" << endl << logfile << endl;
    
    cout << "Plot the results with MATLAB" << endl;
    cout << "A plot script is contained in test_cached_qtproblem" << endl;
#if 0 
% Matlab code to plot output of test_cached_qtproblem
% First half for dec3, code for L-shaped domain below
clear
figure

run('log_rhs1_dec3_off4_2_2_hmax1'); % selected
erg(1,1) = 1;
loglog(erg(:,1), erg(:,2)/erg(1,2),'k-');
clear('erg');
hold
run('log_rhs1_dec13_off1_2_2_hmax1'); % selected
erg(1,1) = 1;
loglog(erg(:,1), erg(:,2)/erg(1,2),'k--');
clear('erg');
% add triangle
% if d==2 use this code
line ([5000 50000],[0.1 0.01],'color','k')
line ([5000 50000],[0.1 0.1],'color','k')
line ([50000 50000],[0.1 0.01],'color','k')
set(gca,'Xtick',[1 10 100 1000 10000]);
set(gca,'Ytick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]);
set(gca,'FontSize',16)
% else, if d==3
line ([5000 50000],[0.1 0.001],'color','k')
line ([5000 50000],[0.1 0.1],'color','k')
line ([50000 50000],[0.1 0.001],'color','k')

% Plot the solution
% L-shaped domain:
run('plot_rhs1_dec3_off4_2_2_hmax1');
figure;
hold;
 [M N] =size(x);
 for i=1:N
     surf(x{i},y{i},z{i});
     % colorscale gray;
     shading interp;
 end
#endif // Matlab code

           
#else
    cout << "Skipping test of CDD2" << endl;
#endif
    
#if _TEST_ADAPTIVE_STUFF
    
    tstart = clock();
    
    
    set<int> Lambda;
    for (unsigned int n=0; n< cqtproblem.basis()->degrees_of_freedom(); ++n)
    {
        Lambda.insert(n);
    }

    SparseMatrix<double> A_apply;
    //A_apply.resize(Lambda_apply.size(), Lambda_apply.size());
    cout << "  setting up stiffness matrix ..." << endl;
    setup_stiffness_matrix(cqtproblem, Lambda, A_apply);
    cout << "  done" << endl;
    Vector<double> x(cqtproblem.basis()->degrees_of_freedom()), 
          res_apply1(cqtproblem.basis()->degrees_of_freedom()), 
          res_apply2(cqtproblem.basis()->degrees_of_freedom());
    Vector<double> apply_diff;
    double temp_apply1(0), temp_apply2(0);
    ctmax = 0;
    eqmax = 0;
    for (unsigned int n=0; n<cqtproblem.basis()->degrees_of_freedom(); ++n)
    //for (set<Index>::iterator it1(Lambda_apply.begin()), itend(Lambda_apply.end());it1!=itend;it1++)
    {
        x[n] = 1;
        A_apply.apply(x,res_apply1);
        cqtproblem.apply(Lambda,x,res_apply2);
        apply_diff = res_apply1;
        apply_diff.subtract(res_apply2);
        temp_apply1 = linfty_norm(apply_diff);
        temp_apply2 = l2_norm(apply_diff);
        if (temp_apply1 > ctmax)
        {
            ctmax = temp_apply1;
            it1max = n;
        }
        if (temp_apply2 > eqmax)
        {
            eqmax = temp_apply2;
            it2max = n;
        }
        if (abs(temp_apply1)+abs(temp_apply2)>1e-15)
        {
            cout << "n = " << n << "; " << *cqtproblem.basis()->get_wavelet(n) << ";  Apply(x) - Ax. linfty error=" << temp_apply1 << "; l2 error = " << temp_apply2 << endl;
        }
        x[n] = 0;
    }
    if (abs(ctmax)+abs(eqmax)>1e-16)
    {
        cout << "largest errors in apply:" <<endl;
        cout << "linfty = "<< ctmax << " at index ("<< it1max << ") = " << *cqtproblem.basis()->get_wavelet(it1max) << endl;
        cout << "l2     = "<< eqmax << " at index ("<< it2max << ") = " << *cqtproblem.basis()->get_wavelet(it2max) << endl;
    }
    else
    {
        cout << " all errors below tolerance of 1e-16" << endl;
    }
    cout << " done"<<endl;
    
#else
    cout << "skip testing of _TEST_ADAPTIVE_STUFF" << endl;
#endif
    
    
    
    
    
    
    
    cout << "cleaning up the memory" << endl;
    for (unsigned int p=0; p<num_of_patches;++p)
    {
        //cout << TEqArray[p] -> basis().degrees_of_freedom() << endl;
        //cout << CTEqArray[p] -> basis().degrees_of_freedom() << endl;
#if ( (_COMPARE_QT_a_AND_T_a || _CALLS_TO_CQTPROBLEM_f) || ( _COMPARE_QT_RHS_AND_T_RHS ) )
        delete TEqArray[p];
        delete CTEqArray[p];
#endif
        //cout << TEqArray[p] -> basis().degrees_of_freedom() << endl;
        //cout << CTEqArray[p] -> basis().degrees_of_freedom() << endl;
    }
    
#if 0          

  
#if 0
  cout << "testing D ..."<<endl;;
#else
  cout << "skipping testing of D"<<endl;
#endif
  
#if 0
  cout << "testing calls to the cache ..."<<endl;
#else
  cout << "skipping testing of cache calls"<<endl;
#endif

#if 0 // I forgot to implement&call set_jmax for the underlying equation!
#endif

#if 0 // use only if apply test is active
  cout << "testing apply for specific indices" << endl;
#endif

#if 0
  cout << "testing add_ball"<<endl;
#else
  cout << "skipping test of add_ball" << endl;
#endif

#if 0
  // computation of normA with a Lanczos iteration with tolerance tol produces only results up to tolerance tol !!
#endif
  
# if 0
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD2_SOLVE ..." << endl;
#else
  cout << "skipping CDD2 test"<<endl;
#endif
  
# if 0
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD1_SOLVE ..." << endl;
#else
  cout << "skipping CDD1 test"<<endl;
#endif
  
#endif


         
    return 0;
}

#if 0 
  return 0;
}
#endif





#if 0
    typedef PBasis<d,dT> Basis1d;
    typedef Basis1d::Index Index1d;
    typedef CubeBasis<Basis1d,DIM> Basis;
    typedef Basis::Index Index;
    typedef Basis::Support Support;
    TestRHS<5> f2; // 3.0
//    int res(8);
//    Grid<2> grid3(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<res, 1<<res);
//    SampledMapping<2> abbildung(grid3, f2);
//    ofstream datei("Constant_Function.m");
//    abbildung.matlab_output(datei);
//    datei.close();
    
    FixedArray1D<bool,2*DIM> bc;
    bc[0] = false;
    bc[1] = false;
    bc[2] = false;
    bc[3] = false;
    Basis cbasis(bc);
    //Basis1d basis1d(false,false);
    //FixedArray1D<Basis1d*,2> bases1d;
    //bases1d[0] = bases1d[1] = &basis1d;
    //TBasis tbasis(bases1d);
    
    const int maxrange = 0;
    //MultiIndex<int,2> jmax;
    //jmax[0]=tbasis.bases()[0]->j0()+maxrange;
    //jmax[1]=tbasis.bases()[1]->j0();
    
    int jmax = cbasis.j0()+maxrange;
    cbasis.set_jmax(jmax);
    InfiniteVector<double,Index> coeffs;
    cbasis.expand(&f2, false, jmax, coeffs);
    cout << coeffs << endl;
    
    SampledMapping<DIM> plot_patch (evaluate (cbasis, coeffs, false,6));
    plot_patch.matlab_output(cout);
    
//    std::ofstream plotstream_RHStest;
//    plotstream_RHStest.open("plot_F_t2.m");
//    plot_patch2.matlab_output(plotstream_RHStest);
//    plotstream_RHStest.close();
//    
    
#endif  
    
    
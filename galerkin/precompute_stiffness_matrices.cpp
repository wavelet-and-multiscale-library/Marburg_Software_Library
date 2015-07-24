/* 
 * This file is for precomputing stiffness matrices and storing them on the hard disk
 * 
 * Motivation: You are interested in the matrix A_sigma = (\int sigma grad psi_lambda grad psi_mu)_{lambda,mu} where
 * - psi_lambda, psi_mu are 2d wavelets as given by TBasis
 * - sigma = sum_eta=0^N c_eta h_eta, c_eta\in\R, h_eta = Haar wavelet number eta
 * You need A_sigma for several sigma.
 * 
 * Strategy: Compute all Haar matrices (\int h_eta grad psi_lambda grad psi_mu)_{lambda,mu}, eta=0,...,N
 * Then compose (update) A_sigma with the help of these matrices.
 * 
 * Example 1: Consider the solution of a parabolic PDE
 *      u' + A(t)u = f(t,u)
 * where A depends on the time
 * 
 * Example 2: Solve the inverse problem related to
 *      u' + A(\sigma, t)u = f(t,u).
 *      That is: \sigma is a parameter with an unknown decomposition into Haar wavelets.
 *      The decomposition is approximated by an iteration scheme. 
 *      Consequently lots of parabolic PDEs have to be solved and 
 *      consequently lots and lots of stiffness matrices A(\sigma_n, t_j) have to be computed.
 * 
  * 
 * Strategy 1: Compute and store all Haar matrices
 * PRO: Computing the matrices for each Haar wavelet means less computing when they are actually used, e.g., by a parabolic solver
 * CON: This approach needs A LOT of disk space!
 * Strategy 2: Compute only matrices that correspond to one space dimension:
 *      \int h_eta grad psi_lambda grad psi_mu
 *      is an integral over \Omega, e.g. the cube in N dimensions. 
 *      Since all wavelets are tensors of one dimensional functions it is possible 
 *      to only consider one dimensional integrals. The full matrices can be assembled 
 *      from the one dimensional values whenever needed (here it would be possible to use 
 *      a cache for the full N-dimensional entries)
 * PRO: Highly reduced need for disk space when compared to Strategy 1.
 * CON: Values of the final matrix now need to be composed from the one dimensional matrices == more computation effort. 
 * If the N-dimensional matrix is stored in memory, then the total amount of memory consumption is even higher than that of Strategy 1.
 * 
 * Strategy 3: Don not use any precomputation. 
 * Compute the values that are needed for the one dimensional integrals and store them in a cache. 
 * Do not use a cache for the N-dimensional matrices (as this costs too much memory).
 * PRO: Least amount of used memory. No requirements on disk space.
 * CON: More computation effort than other strategies. However, no IO!
 * Strategy 3 is realized in cached_QTproblem
 * 
 * 
 * 
 * "spatial wavelets":
 * 1d wavelet basis: Primbs
 * 2d wavelet basis: TBasis
 * 
 * - computed matrices:
 * gramian: \int psi_lambda psi_mu  // psi_lambda, psi_mu are 1 or 2d wavelets
 * laplacian: \int grad psi_lambda grad psi_mu
 * haar_gen_gramian: \int g_eta psi_lambda psi_mu // g_eta is a Haar generator on some level with number eta;
 * haar_wav_gramian: \int h_eta psi_lambda psi_mu // h_eta is the Haar wavelet with number eta
 * haar_gen_laplacian: \int g_eta grad psi_lambda grad psi_mu
 * haar_wav_laplacian: \int h_eta grad psi_lambda grad psi_mu
 * 
 * - as helper matrices: 
 * ij-subblocks of the above matrices for 1d spatial wavelets. The subblocks correspond to all wavelets psi_lambda, psi_mu with |lambda|=i, |mu|=j
 * 
 * - for the project with bremen: compute a specific RHS corresponding to the function 1/2
 * 
 * - transition matrices (<h_eta,psi_mu>)_(eta,mu)
 * 
 * 
 * If desired output is preconditioned with the diagonal.
 * 
 * Usage:
 * Program is controlled by makros.
 * ...OFFSET_START X, ...OFFSET_END Y: 
 *      compute Y-X+1 matrices. They correspond to all wavelets from level j0 up to jo+X, ..., j0 up to j0+Y
 *      Y < X : Skip this computation (if not otherwise specified)
 * 
 * ...HAAR_LEVEL_START, ...HAAR_LEVEL_END: 
 *      consider Haar wavelets or generators between this values
 *      ...End < ...START: Skip this computation
 *      These makros are only relevant if _HAAR_WAV_MATRIX_FIRST_ETA is negative.
 * 
 * ..._FIRST_ETA, ..._LAST_ETA:
 *      Computation of the 2d matrices for all Haar wavelets may take a while.
 *      Parallel computation of the matrices with different instances of this program is advised.
 *      This makro makes it possible to balance the workload.
 *      Instead of computing the matrices for all Haar wavelets on a specified level in a serial way it is possible to restrict the computation to eta between some bounds.
 *      This is in partivcular useful for larger levels where simply more Haar wavelets exist than on lower levels
 *      For packages with similar workloads see below
 * 
 * computation of matrices for 1d primbs wavelets:
 * specify the wanted makros and run this program.
 * 
 * computation of matrices for 2d TBasis wavelets:
 * 2d matrices are computed as tensor products of subblocks of matrices that correspond to 1d prombs wavelets.
 * Therefore you need to compute the correct 1d matrices first. E.g. for boundary conditions tttf (tt in x, tf in y direction):
 * 1. Compute tt and tf 1d matrices.
 * 2. Compute subblocks of the 1d matrices
 *      Note: In order to compute a 2d Haar_wav matrix you will need 1d Haar_wav and 1d Haar_gen matrices!
 * 3. Compute tttf 2d Matrix. 
 *      Note: The minimal level for a 2d matrix is given by ||j0||_1 = j0[x] + j0[y]
 *      All 2d wavelets psi_lambda, psi_mu with ||lambda-j0||_1,||mu-j0||_1 <= OFFSET are considered
 * It is not possible to compute 2d Haar_gen matrices at the moment. They are simply not needed!
 * 
 * Example parameters
 * 
 *      Spatial offset:
 *      1d: 
 *      d=dt=2,3 : 7-10 (maxlevel 10,13)
 *      2d:
 *      d=dt=2 : 4 (maxlevel 8)
 *      d=dt=3 : 3 (maxlevel 9)
 *      Haar levels
 *      1d start =0, end =5 => 64 Haar wavelets
 *      2d start =0, end =3,4,5 => 256,1024,4096 Haar wavelets.
 *      
 *      2d, d=dt=2, bc=ffff, jmax=4,5,6,7,8 => spatial dof = 81, 225, 577, 1409, 3329
 *      2d, d=dt=3, bc=ffff, jmax=6,7,8,9   => spatial dof = 324,900,2308,5636
 * 
 *      2d: d=dt=2, bc=ffff, jmax=8 => each haar_wav_gramian matrix uses 10-30M bytes!!
 *      2d: d=dt=3, bc=ffff, jmax=8 => each haar_wav_gramian matrix uses 6-43M bytes!!
 *      2d: d=dt=3, bc=ffff, jmax=9 => each haar_wav_gramian matrix uses 30-120M bytes!!
 * 
 *      2d: d=dt=2, spatial maxoffset 3, haarlevel 4 seems balanced?!
 *      2d: d=dt=3, spatial maxoffset 2, haarlevel 4 seems balanced?!
 * 
 * Caution:
 * - Numbering: wavelet numbering begins at 0. 0th wavelet is indeed a generator.
 * Numbering of generators in 1d: level 0: No 1. level 1: No 2,3. level 2: 4,5,6,7 ...
 * Numbering of generators in 2d: analogous!
 * 
 * - Depending on the parameters the program computes a LOT matrices.
 * Especially in 2d this may fill up the harddrive!
 * Example: 1d Haar wavelets on level 0 up to 5: 64 wavelets
 * Blocks up to an offset of 7: ij = {00,01,...,07,10,...,17,...,77} : 64
 * For all 2d haar_wav_laplacian matrix with spatial offset 7, BC tttf (not even up to 7!):
 * different BC in each dim: *2
 * 1d _gen and _wav needed : *2
 * 1d _gramian and _laplacian needed in each dim: *2
 * 64*64*2*2*2 = 16384 subblock matrices!!
 * 
 * I suggest to delete the subblock matrices as soon as they are not longer needed
 * 
 */

// more output for the cached problem (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controlled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0

#define _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY 0 // more output for compose_matrix

// switch between isotropic and anisotropic Wavelets (in cdd1.h)
#define _WAVELETTL_USE_TBASIS 1

#define _DIMENSION 2
#define _D 2
#define _DT 2
#define _PRECONDITIONED 0 // 0,1 compute preconditioned versions of the matrices?

#if _DIMENSION == 1
#define _BOUNDARY_CONDITIONS 3 // 1D :: 0 = tt, 1 = tf, 2 = ft, 3 = ff
#else
#define _BOUNDARY_CONDITIONS 1 // 2D :: 0 = tttt, 1=ffff
#endif

#define _NO_JPMAX_grt_70_WARNING 1 // turn off warning in APPLY_TENSOR

#define _LAPLACIAN_MATRIX_OFFSET_START 0 // minimal offset for computation of stiffnessmatrix (laplacian)
#define _LAPLACIAN_MATRIX_OFFSET_END -2 // maximal offset for computation of stiffnessmatrix (laplacian)

#define _GRAMIAN_MATRIX_OFFSET_START 0// same for \int \psi_\lambda\psi_\mu
#define _GRAMIAN_MATRIX_OFFSET_END -2

#define _RHS_OFFSET_START 3// same for \int 1/2 \psi_\lambda
#define _RHS_OFFSET_END 2

// \int haar_gen_eta psi_lambda psi_mu
#define _HAAR_GEN_MATRIX_SPATIAL_OFFSET_START 0// same for \int \h_\eta \psi_\lambda \psi_\mu
#define _HAAR_GEN_MATRIX_SPATIAL_OFFSET_END 2
#define _HAAR_GEN_MATRIX_HAAR_LEVEL_START 0 // range of Haar generators considered
#define _HAAR_GEN_MATRIX_HAAR_LEVEL_END 2

// \int haar_wav_eta psi_lambda psi_mu
// _HAAR_WAV_MATRIX_FIRST_ETA < 0 means: compute values for all haar wavelets (from HAAR_LEVEL_START to _END)!
// This will take some time (in 2d)! Therefore it is possible to specify FIRST/LAST_ETA
#define _HAAR_WAV_MATRIX_SPATIAL_OFFSET_START 0// family of matrices is computed for jmin,...,jmin+offset. same for \int \h_\eta \psi_\lambda \psi_\mu
#define _HAAR_WAV_MATRIX_SPATIAL_OFFSET_END 2 //2D: ddt22=>4, ddt33=>3, sonst platzt die festplatte
#define _HAAR_WAV_MATRIX_FIRST_ETA -1 // any nonnegative value disables the level functionality // 1D,j=5:: 0,...,2^(j+1)-1;  2D,j=5 :: 0,...,2^(2j+2)-1
#define _HAAR_WAV_MATRIX_LAST_ETA -7 // packets with similar workload: 0_650, 651_1350, 1351_2100, 2101_3100, 3101_4095
#define _HAAR_WAV_MATRIX_HAAR_LEVEL_START 0 // range of Haar wavelets considered. only active if _HAAR_WAV_MATRIX_FIRST_ETA < 0
#define _HAAR_WAV_MATRIX_HAAR_LEVEL_END 2 // 2D: 5 = 4096 wavelets => 20 sec per matrix <=> 1 day


// \int haar_gen_eta grad psi_lambda grad psi_mu
#define _HAAR_GEN_LAPLACE_SPATIAL_OFFSET_START 0
#define _HAAR_GEN_LAPLACE_SPATIAL_OFFSET_END 2
#define _HAAR_GEN_LAPLACE_HAAR_LEVEL_START 0 
#define _HAAR_GEN_LAPLACE_HAAR_LEVEL_END 2

// \int haar_wav_eta grad psi_lambda grad psi_mu
#define _HAAR_WAV_LAPLACE_SPATIAL_OFFSET_START 0
#define _HAAR_WAV_LAPLACE_SPATIAL_OFFSET_END 2
#define _HAAR_WAV_LAPLACE_FIRST_ETA -513 
#define _HAAR_WAV_LAPLACE_LAST_ETA -1023 // packets with similar workload: 0_650, 651_1350, 1351_2100, 2101_3100, 3101_4095
#define _HAAR_WAV_LAPLACE_HAAR_LEVEL_START 0
#define _HAAR_WAV_LAPLACE_HAAR_LEVEL_END -2

#define _BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_START 0 // precompute basistransitionmatrices (<h_eta,psi_mu>)_(eta,mu)
#define _BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END -4
#define _BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_START 0 
#define _BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_END -5

// Blocks = submatrices of matrices corresponding to 1 dimensional primbs wavelets
// computation: load the matrix corresponding to all wavelets up to jmax, then cut out the blocks.
// Blocks are used in a second step to compose matrices corresponding to anisotropic wavelets in two dimensions

//#define _LAPLACIAN_MATRIX_BLOCK_START 0
#define _LAPLACIAN_MATRIX_BLOCK_MAXOFFSET -2 // compute blocks of the 1D stiffness matrix: (<grad psi_lambda, grad psi_mu>)_(lambda,mu) with |mu|=j0+const1, |lambda|=j0+const2 and range of const1,2 given by MAXOFFSET

//#define _GRAMIAN_MATRIX_BLOCK_START 0
#define _GRAMIAN_MATRIX_BLOCK_MAXOFFSET -2 // same for <psi_lambda,psi_mu>


// \int haar_gen_eta psi_lambda psi_mu 
//#define _HAAR_GEN_MATRIX_BLOCK_SPATIAL_START 0
#define _HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET -2 // for each fixed Haar generator compute the blocks <h_eta psi_lambda,psi_mu> for lambda and mu from on level j0 up to j0+MAXOFFSET.
#define _HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_START 0 // The levels considered for the Haar generators are defined by HAARLEVEL_START and HAARLEVEL_END
#define _HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_END -2

// \int haar_wav_eta psi_lambda psi_mu 
//#define _HAAR_WAV_MATRIX_BLOCK_SPATIAL_START 0
#define _HAAR_WAV_MATRIX_BLOCK_SPATIAL_MAXOFFSET -2 // same for 1d Haar wavelets
#define _HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_START 0
#define _HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_END -2

// \int haar_gen_eta grad psi_lambda grad_psi mu
#define _HAAR_GEN_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET -2 // for each fixed Haar generator compute the blocks <h_eta psi_lambda,psi_mu> for lambda and mu from on level j0 up to j0+MAXOFFSET.
#define _HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_START 0 // The levels considered for the Haar generators are defined by HAARLEVEL_START and HAARLEVEL_END
#define _HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_END -2

// \int haar_wav_eta grad psi_lambda grad_psi mu
#define _HAAR_WAV_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET -2 // same for 1d Haar wavelets
#define _HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_START 0
#define _HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_END -2


#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
//#include <interval/ds_basis.h>
#include <interval/p_basis.h>

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tbasis_indexplot.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
 * -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega
 *
 * a = 0, q =h_eta (Haar wavelet)
 *
 * h_eta (x)=
 * 1 2^(-j)k <= x <= 2^(-j)(k+1/2)
 * -1 2^(-j)(k+1/2) < x <= 2^(-j)(k+1)
 *
 * where j= level = floor (log_2(eta)), k = eta - 2^j and eta >= 1
 * h_0 == 1 denotes the generator. (propably not used
 *
 * For dimensions > 1 we use isotropic tensor Haar wavelets. The resultung matrices
 * are tensors of the 1D matrices and need additionally generators on highre levels.
 * -> HaarGeneratorBVP
 */
class HaarWaveletBVP
: public EllipticBVP<1>
{
public:
    /*!
      constructor with given right-hand side
    */
    HaarWaveletBVP(const Function<1>* f = 0, const unsigned int eta = 1)
    : EllipticBVP<1> (f, f, f), j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta > 0);
    }
    /*!
      diffusion coefficient a
     */
    const double a(const Point<1>& x) const { return 0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<1>& x) const
    {
        return (
        (1.0/(1<<j) *k <= x[0])?
            ( (1.0/(1<<j) *(k+0.5) < x[0]) ?
                ( ((1.0/(1<<j) *(k+1)) >= x[0])? -1.0 : 0)
                                           : 1.0
            )
                               : 0.0
               )* (double) twotothejhalf(j);
    }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return false; }
    
//protected:
    const int j;
    const int k;
};


class HaarGeneratorBVP
: public EllipticBVP<1>
{
public:
    /*!
      constructor with given right-hand side
    */
    HaarGeneratorBVP(const Function<1>* f = 0, const unsigned int eta = 1)
    : EllipticBVP<1> (f, f, f), j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta>0);
    }
    /*!
      diffusion coefficient a
     */
    const double a(const Point<1>& x) const { return 0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<1>& x) const
    {
        return (double) twotothejhalf(j)*( ( (1.0/(1<<j) *k <= x[0]) &&  ( ( (1.0/(1<<j) *(k+1)) > x[0]) || ((k+1) == (1<<j) ) )  )?1.0:0.0);
    }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return false; }

//protected:
    const int j;
    const int k;
};




class HaarWaveletPoissonBVP
: public EllipticBVP<1>
{
public:
    /*!
      constructor with given right-hand side
    */
    HaarWaveletPoissonBVP(const Function<1>* f = 0, const unsigned int eta = 1)
    : EllipticBVP<1> (f, f, f), j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta > 0);
    }
    

    /*!
      diffusion coefficient a
    */
    const double a(const Point<1>& x) const
    {
        return (
        (1.0/(1<<j) *k <= x[0])?
            ( (1.0/(1<<j) *(k+0.5) < x[0]) ?
                ( ((1.0/(1<<j) *(k+1)) >= x[0])? -1.0 : 0)
                                           : 1.0
            )
                               : 0.0
               )* (double) twotothejhalf(j);
    }
    
    /*!
      reaction coefficient q
     */
    const double q(const Point<1>& x) const { return 0; }
    

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return false; }
    
//protected:
    const int j;
    const int k;
};


class HaarGeneratorPoissonBVP
: public EllipticBVP<1>
{
public:
    /*!
      constructor with given right-hand side
    */
    HaarGeneratorPoissonBVP(const Function<1>* f = 0, const unsigned int eta = 1)
    : EllipticBVP<1> (f, f, f), j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta>0);
    }
    /*!
      diffusion coefficient a
     */
    const double a(const Point<1>& x) const
    {                   
        return (double) twotothejhalf(j)*( ( (1.0/(1<<j) *k <= x[0]) &&  ( ( (1.0/(1<<j) *(k+1)) > x[0]) || ((k+1) == (1<<j) ) )  )?1.0:0.0);
    }

    /*!
      reaction coefficient q
    */
    const double q(const Point<1>& x) const { return 0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return false; }

//protected:
    const int j;
    const int k;
};

/* gives the Haar wavelet in 1D*/
class HaarWavelet1D
: public Function<1>
{
public:
    HaarWavelet1D(const unsigned int eta = 1)
    : j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta > 0);
    }
    virtual ~HaarWavelet1D() {};
    double value(const Point<1>& x, const unsigned int component = 0) const
    {
        return (
        (1.0/(1<<j) *k <= x[0])?
            ( (1.0/(1<<j) *(k+0.5) < x[0]) ?
                ( ((1.0/(1<<j) *(k+1)) >= x[0])? -1.0 : 0)
                                           : 1.0
            )
                               : 0.0
               )* (double) twotothejhalf(j);
    }
    void vector_value(const Point<1>& p, Vector<double>& values) const
    {
        values[0] = value(p);
    }
    //protected:
    const int j,k;
};

/* gives the Haar wavelet in 2D*/
class HaarWavelet2D
: public Function<2>
{
public:
    HaarWavelet2D(const int eta = 1)
    //: j((int) floor(log2(eta))), k(  eta - (1<<(int) floor(log2(eta))) )
    {
        assert(eta > 0);
        int n(eta),l=0;
        while ((n>>=1)>0) ++l; // log_2 of eta
        j = l/2; // 2d level j of eta
        int type = eta / (1<<(j<<1)); // 2d type: 1 = gen x wav, 2 = wav x gen, 3 = wav x wav
        int k = eta - type*(1<<(j<<1)); // number in current type
        xnum = (k / (1<<j));
        ynum = (k % (1<<j)); // number of gen or wav in x and y direction in the current level j. range is {0,...,2^j-1}. caution: numbering in HaarGenBVP begins at 1 (as does it in HaarWavBVP)
        //xnum+=(1<<j);ynum+=(1<<j); // number (first gen or wav on first level has number 1)
        switch (type)
        {
            case 1:
                xtype=0; //gen
                ytype=1; //wav
                break;
            case 2:
                xtype=1;
                ytype=0;
                break;
            case 3:
                xtype=1;
                ytype=1;
                break;
            default:
                abort();
                break;
        }
    }
    virtual ~HaarWavelet2D() {};
    double value(const Point<2>& x, const unsigned int component = 0) const
    {
        return (xtype == 1 ?
                 (
        (1.0/(1<<j) *xnum <= x[0])?
            ( (1.0/(1<<j) *(xnum+0.5) < x[0]) ?
                ( ((1.0/(1<<j) *(xnum+1)) >= x[0])? -1.0 : 0)
                                           : 1.0
            )
                               : 0.0
               )
                           : (((1.0/(1<<j) *xnum <= x[0]) && ((1.0/(1<<j) *(xnum+1)) > x[0]))?1.0:0.0)
                ) *
               (ytype == 1 ?
                 (
        (1.0/(1<<j) *ynum <= x[1])?
            ( (1.0/(1<<j) *(ynum+0.5) < x[1]) ?
                ( ((1.0/(1<<j) *(ynum+1)) >= x[1])? -1.0 : 0)
                                           : 1.0
            )
                               : 0.0
               )
                           : (((1.0/(1<<j) *ynum <= x[1]) && ((1.0/(1<<j) *(ynum+1)) > x[1]))?1.0:0.0)
                ) * (1<<j);
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const
    {
        values[0] = value(p);
    }
    //protected:
    int j,xnum,ynum,xtype,ytype;
};

/* A Method to compose the matrix with entries
 *   c_{lambda mu} = a_{lambda_1 mu_1}*b_{lambda_2 mu_2}
 * where lambda and mu are 2d anisotropic tensor wavelet indices (tbasis)
 * and B,C are given matrices corresponding to 1D wavelet indices.
 * We assume that the rectangular matrices B_ij, C_kl corresponding of all pairs
 * of wavelets with |lambda_1|=i and |mu_1|=j (and C analog) exist as precomputed
 * files.
 * Entries up to ||lambda-j0||_1,||mu-j0||_1 <= offset are considered.
 *
 * example input:
 * a_namestart   = /a_path/laplacian_primbs_d_dt_3_3_bc_
 * a_bcstring     = tt
 * a_nameend     = npcd_block_ //blocknumbers are added as they are needed
 * b_namestart   = /a_path/laplacian_primbs_d_dt_3_3_bc_
 * b_bcstring     = tf
 * b_nameend     = npcd_block_
 * out_namestart = /out_path/laplacian_primbs_d_dt_3_3_bc_
 * out_nameend   = jmax_3_npcd
 * offset        = 0
 * numberofgens  = 8 //primbs tt ddT=3
 *
 * gives output:
 * out_namestart_a_bcstringb_bcstring_out_nameend
 * = /outpath/laplacian_primbs_d_dt_3_3_bc_tttf_jmax_6_npcd
 */
void compose_matrix(SparseMatrix<double> & A,
        const char* a_namestart, const char* a_bcstring, const char* a_nameend,
        const char* b_namestart, const char* b_bcstring, const char* b_nameend,
        //const char* out_namestart, const char* out_nameend,
        const unsigned int offset,
        const unsigned int number_of_gens)//, const unsigned int size)
{
    cout << "compose_matrix :: A.column_dimension() = " << A.column_dimension() << " A.row_dimension() = " << A.row_dimension() << endl;
    //SparseMatrix<double> A(size,size);
     /*!
      write access to a subblock;
      if the "reflect" flag is set, rows and columns of the block are reflected before writing
    */
    //template <class MATRIX>
    //void set_block(const size_type firstrow, const size_type firstcolumn,
//		   const MATRIX& M,
		   //const bool reflect = false);
    unsigned int startrow(0), startcolumn; // row and columnindex where the current block will be inserted
    //unsigned int startrow, startcolumn; // indices where a_(i,j)*BBlock_{bi,bj} is inserted
    // iterate over the level blocks, eg all lambda and mu indices have the same norm
    for (unsigned int rowoffset=0; rowoffset <= offset; ++rowoffset)
    {
        for (unsigned int ai = 0; ai <= rowoffset; ++ai)
        {
            //iterate over the sublevels (g0w7 w0w7 w1w6 w2w5 ... w7g0 w7w0 are the sublevels of the level 7)
            // Observe ai+bi == current_rowoffset, aj+bj == current_columnoffset
            startcolumn = 0;
            for (unsigned int columnoffset=0; columnoffset <= offset; ++columnoffset) // maybe offset-rowoffset is enough
            {
                // bi == rowoffset - ai;
                //startcolumn = startcolumn;
                for (unsigned int aj = 0; aj <= columnoffset; ++aj)
                {

#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "adding block (ro,co) = (" << rowoffset << ", " << columnoffset << "), (ai,aj) = (" << ai << ", " << aj << ") at position (lro,lco) = (" << startrow << ", " << startcolumn << ")" << endl;
#endif
                    // bj == columnoffset - aj
                    // insert the current sublevel block.
                    // load matrices with all 1d integrals
                    SparseMatrix<double> ABlock, BBlock; // it is inefficient load so many matrices again and again
                    char matrix_filename[250];
                    sprintf(matrix_filename, "%s%s%s%d_%d",a_namestart,a_bcstring,a_nameend,ai,aj);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "loading file: " << matrix_filename << endl;
#endif
                    ABlock.matlab_input(matrix_filename);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "ABlock.column_dimension() = " << ABlock.column_dimension() << " ABlock.row_dimension() = " << ABlock.row_dimension() << endl;
#endif
                    //cout << "ABlock = " << endl;
                    //ABlock.print(cout, 4, 4);

                    sprintf(matrix_filename, "%s%s%s%d_%d",b_namestart,b_bcstring,b_nameend, rowoffset-ai, columnoffset-aj);// bi,bj);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "loading file: " << matrix_filename << endl;
#endif
                    BBlock.matlab_input(matrix_filename);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "BBlock.column_dimension() = " << BBlock.column_dimension() << " BBlock.row_dimension() = " << BBlock.row_dimension() << endl;
#endif

                    // special care needs to be taken to handle the fact that generators and wavelets on the first level are stored
                    // in a single 1d-file
                    if (rowoffset == 0)
                    {
                        if (columnoffset == 0)
                        {
                            // first subblock on first level!
                            // A1xB1 A1xB2 A2xB1 A2xB2
                            // A1xB3 A1xB4 A2xB3 A2xB4
                            // A3xB1 A3xB2 A4xB1 A4xB2
                            // A3xB3 A3xB4 A4xB3 A4xB4

                            //void get_block(const size_type firstrow, const size_type firstcolumn,
                            //               const size_type rows, const size_type columns,
                            //               MATRIX& M) const;
                            SparseMatrix<double> A1,A2,A3,A4,B1,B2,B3,B4;
                            ABlock.get_block(0             ,0             ,number_of_gens                          ,number_of_gens                          ,A1);
                            ABlock.get_block(0             ,number_of_gens,number_of_gens                          ,ABlock.column_dimension()-number_of_gens,A2);
                            ABlock.get_block(number_of_gens,0             ,ABlock.column_dimension()-number_of_gens,number_of_gens                          ,A3);
                            ABlock.get_block(number_of_gens,number_of_gens,ABlock.column_dimension()-number_of_gens,ABlock.column_dimension()-number_of_gens,A4);

                            BBlock.get_block(0             ,0             ,number_of_gens                          ,number_of_gens                          ,B1);
                            BBlock.get_block(0             ,number_of_gens,number_of_gens                          ,BBlock.column_dimension()-number_of_gens,B2);
                            BBlock.get_block(number_of_gens,0             ,BBlock.column_dimension()-number_of_gens,number_of_gens                          ,B3);
                            BBlock.get_block(number_of_gens,number_of_gens,BBlock.column_dimension()-number_of_gens,BBlock.column_dimension()-number_of_gens,B4);

// CLEANUP
                            /*
                            int a1rd = A1.row_dimension();
                            int a2rd = A2.row_dimension();
                            int a3rd = A3.row_dimension();
                            int a4rd = A4.row_dimension();
                            int a1cd = A1.column_dimension();
                            int a2cd = A2.column_dimension();
                            int a3cd = A3.column_dimension();
                            int a4cd = A4.column_dimension();
                            int b1rd = B1.row_dimension();
                            int b2rd = B2.row_dimension();
                            int b3rd = B3.row_dimension();
                            int b4rd = B4.row_dimension();
                            int b1cd = B1.column_dimension();
                            int b2cd = B2.column_dimension();
                            int b3cd = B3.column_dimension();
                            int b4cd = B4.column_dimension();
                            cout << "a1rd = " << a1rd << " a1cd = " << a1cd << endl
                                 << "a2rd = " << a2rd << " a2cd = " << a2cd << endl
                                 << "a3rd = " << a3rd << " a3cd = " << a3cd << endl
                                 << "a4rd = " << a4rd << " a4cd = " << a4cd << endl
                                 << "b1rd = " << b1rd << " b1cd = " << b1cd << endl
                                 << "b2rd = " << b2rd << " b2cd = " << b2cd << endl
                                 << "b3rd = " << b3rd << " b3cd = " << b3cd << endl
                                 << "b4rd = " << b4rd << " b4cd = " << b4cd << endl;
                            */

                            assert ((A1.row_dimension() == A2.row_dimension())
                                    && (A3.row_dimension() == A4.row_dimension())
                                    && (A1.column_dimension() == A3.column_dimension())
                                    && (A2.column_dimension() == A4.column_dimension())
                                    && (B1.row_dimension() == B2.row_dimension())
                                    && (B3.row_dimension() == B4.row_dimension())
                                    && (B1.column_dimension() == B3.column_dimension())
                                    && (B2.column_dimension() == B4.column_dimension()));

                            for (unsigned int i=0; i < A1.row_dimension(); ++i)
                            {
                                for (unsigned int j=0; j < A1.column_dimension(); ++j)
                                {
                                    A.set_block(i*B1.row_dimension()                                      ,                                            j*B1.column_dimension(),B1, false, A1.get_entry(i,j));
                                    A.set_block(i*B2.row_dimension()                                      ,A1.column_dimension()*B1.column_dimension()+j*B2.column_dimension(),B2, false, A1.get_entry(i,j));
                                    A.set_block(A1.row_dimension()*B1.row_dimension()+i*B3.row_dimension(),                                            j*B3.column_dimension(),B3, false, A1.get_entry(i,j));
                                    A.set_block(A1.row_dimension()*B2.row_dimension()+i*B4.row_dimension(),A1.column_dimension()*B3.column_dimension()+j*B4.column_dimension(),B4, false, A1.get_entry(i,j));
                                }
                                for (unsigned int j=0; j < A2.column_dimension(); ++j)
                                {
                                    A.set_block(i*B1.row_dimension()                                      ,j*B1.column_dimension() + A1.column_dimension()*BBlock.column_dimension()                                            ,B1, false, A2.get_entry(i,j));
                                    A.set_block(i*B2.row_dimension()                                      ,j*B2.column_dimension() + A1.column_dimension()*BBlock.column_dimension()+A2.column_dimension()*B1.column_dimension(),B2, false, A2.get_entry(i,j));
                                    A.set_block(A2.row_dimension()*B1.row_dimension()+i*B3.row_dimension(),j*B3.column_dimension() + A1.column_dimension()*BBlock.column_dimension()                                            ,B3, false, A2.get_entry(i,j));
                                    A.set_block(A2.row_dimension()*B2.row_dimension()+i*B4.row_dimension(),j*B4.column_dimension() + A1.column_dimension()*BBlock.column_dimension()+A2.column_dimension()*B3.column_dimension(),B4, false, A2.get_entry(i,j));
                                    
                                }
                            }
                            for (unsigned int i=0; i< A3.row_dimension(); ++i)
                            {
                                for (unsigned int j=0; j < A3.column_dimension(); ++j)
                                {
                                    A.set_block(i*B1.row_dimension() + A1.row_dimension()*BBlock.row_dimension()                                      ,j*B1.column_dimension()                                              ,B1, false, A3.get_entry(i,j));
                                    A.set_block(i*B2.row_dimension() + A1.row_dimension()*BBlock.row_dimension()                                      ,j*B2.column_dimension() + A3.column_dimension()*B1.column_dimension(),B2, false, A3.get_entry(i,j));
                                    A.set_block(i*B3.row_dimension() + A1.row_dimension()*BBlock.row_dimension()+A3.row_dimension()*B1.row_dimension(),j*B3.column_dimension()                                              ,B3, false, A3.get_entry(i,j));
                                    A.set_block(i*B4.row_dimension() + A1.row_dimension()*BBlock.row_dimension()+A3.row_dimension()*B2.row_dimension(),j*B4.column_dimension() + A3.column_dimension()*B3.column_dimension(),B4, false, A3.get_entry(i,j));
                                }
                                for (unsigned int j=0; j < A4.column_dimension(); ++j)
                                {
                                    A.set_block(i*B1.row_dimension() + A2.row_dimension()*BBlock.row_dimension()                                      ,j*B1.column_dimension() + A3.column_dimension()*BBlock.column_dimension()                                            ,B1, false, A4.get_entry(i,j));
                                    A.set_block(i*B2.row_dimension() + A2.row_dimension()*BBlock.row_dimension()                                      ,j*B2.column_dimension() + A3.column_dimension()*BBlock.column_dimension()+A4.column_dimension()*B1.column_dimension(),B2, false, A4.get_entry(i,j));
                                    A.set_block(i*B3.row_dimension() + A2.row_dimension()*BBlock.row_dimension()+A4.row_dimension()*B1.row_dimension(),j*B3.column_dimension() + A3.column_dimension()*BBlock.column_dimension()                                            ,B3, false, A4.get_entry(i,j));
                                    A.set_block(i*B4.row_dimension() + A2.row_dimension()*BBlock.row_dimension()+A4.row_dimension()*B2.row_dimension(),j*B4.column_dimension() + A3.column_dimension()*BBlock.column_dimension()+A4.column_dimension()*B3.column_dimension(),B4, false, A4.get_entry(i,j));
                                }
                            }
                        } // end of 00-block
                        else // 0j-block, rowoffset = 0, columnoffset > 0
                        {
                            // A1xB1 A2xB1 | C1xD1 | ... | E1xF1 E1xF2
                            // A1xB2 A2xB2 | C1xD2 | ... | E1xF3 E1xF4
                            // A3xB1 A4xB1 | C2xD1 | ... | E2xF1 E2xF2
                            // A3xB2 A4xB2 | C2xD2 | ... | E2xF3 E2xF4
                            // (A1 A2)xB1 | same | ... | same
                            // (A1 A2)xB2 | same | ... | same
                            // (A3 A4)xB1 | same | ... | same
                            // (A3 A4)xB2 | same | ... | same
                            // => j=0 : merge first two column-operations. effekt: cases aj=0,..,aj=columnoffset-1 can be handled the same way
                            if (aj < columnoffset)
                            {
                                // A00 = A^gg A^gw = A1
                                //       A^wg A^ww   A2
                                // able to process as rows
                                SparseMatrix<double> A1,A2,B1,B2;
                                ABlock.get_block(0             ,0,number_of_gens                       ,ABlock.column_dimension(),A1);
                                ABlock.get_block(number_of_gens,0,ABlock.row_dimension()-number_of_gens,ABlock.column_dimension(),A2);
                                BBlock.get_block(0             ,0,number_of_gens                       ,BBlock.column_dimension(),B1);
                                BBlock.get_block(number_of_gens,0,BBlock.row_dimension()-number_of_gens,BBlock.column_dimension(),B2);
                                
                                assert ((A1.column_dimension() == A2.column_dimension())
                                        && (B1.column_dimension() == B2.column_dimension())
                                        );
                                assert (startrow == 0);
                                for (unsigned int j=0; j < A1.column_dimension(); ++j)
                                {
                                    for (unsigned int i=0; i < A1.row_dimension(); ++i)
                                    {
                                        A.set_block(i*B1.row_dimension()                                                                                    ,startcolumn + j*B1.column_dimension(),B1, false, A1.get_entry(i,j));
                                        A.set_block(i*B2.row_dimension() + A1.row_dimension()*B1.row_dimension()                                            ,startcolumn + j*B2.column_dimension(),B2, false, A1.get_entry(i,j));
                                    }
                                    for (unsigned int i=0; i < A2.row_dimension(); ++i)
                                    {
                                        A.set_block(i*B1.row_dimension() + A1.row_dimension()*BBlock.row_dimension()                                        ,startcolumn + j*B1.column_dimension(),B1, false, A2.get_entry(i,j));
                                        A.set_block(i*B2.row_dimension() + A1.row_dimension()*BBlock.row_dimension() + A2.row_dimension()*B1.row_dimension(),startcolumn + j*B2.column_dimension(),B2, false, A2.get_entry(i,j));
                                    }
                                }
                            }
                            else // aj == columnoffset
                            {
                                // need to split BBlock in 4 parts and ABlock in 2
                                // inserting
                                // A1xB1 A1xB2
                                // A1xB3 A1xB4
                                // A2xB1 A2xB2
                                // A2xB3 A2xB4
                                SparseMatrix<double> A1,A2,B1,B2,B3,B4;
                                ABlock.get_block(0,0,number_of_gens,ABlock.column_dimension(),A1);
                                ABlock.get_block(number_of_gens,0,ABlock.row_dimension()-number_of_gens,ABlock.column_dimension(),A2);
                                BBlock.get_block(0             ,0             ,number_of_gens                          ,number_of_gens                          ,B1);
                                BBlock.get_block(0             ,number_of_gens,number_of_gens                          ,BBlock.column_dimension()-number_of_gens,B2);
                                BBlock.get_block(number_of_gens,0             ,BBlock.column_dimension()-number_of_gens,number_of_gens                          ,B3);
                                BBlock.get_block(number_of_gens,number_of_gens,BBlock.column_dimension()-number_of_gens,BBlock.column_dimension()-number_of_gens,B4);
                                assert (
                                    (A1.column_dimension() == A2.column_dimension())
                                    && (B1.row_dimension() == B2.row_dimension())
                                    && (B3.row_dimension() == B4.row_dimension())
                                    && (B1.column_dimension() == B3.column_dimension())
                                    && (B2.column_dimension() == B4.column_dimension()));
                                assert (startrow == 0);

                                for (unsigned int j=0; j < A1.column_dimension(); ++j)
                                {
                                    for (unsigned int i=0; i < A1.row_dimension(); ++i)
                                    {
                                        A.set_block(i*B1.row_dimension()                                                                                    ,startcolumn + j*B1.column_dimension()                                              ,B1, false, A1.get_entry(i,j));
                                        A.set_block(i*B2.row_dimension()                                                                                    ,startcolumn + j*B2.column_dimension() + A1.column_dimension()*B1.column_dimension(),B2, false, A1.get_entry(i,j));
                                        A.set_block(i*B3.row_dimension() + A1.row_dimension()*B1.row_dimension()                                            ,startcolumn + j*B3.column_dimension()                                              ,B3, false, A1.get_entry(i,j));
                                        A.set_block(i*B4.row_dimension() + A1.row_dimension()*B2.row_dimension()                                            ,startcolumn + j*B4.column_dimension() + A1.column_dimension()*B3.column_dimension(),B4, false, A1.get_entry(i,j));
                                    }
                                    for (unsigned int i=0; i < A2.row_dimension(); ++i)
                                    {
                                        A.set_block(i*B1.row_dimension() + A1.row_dimension()*BBlock.row_dimension()                                        ,startcolumn + j*B1.column_dimension()                                              ,B1, false, A2.get_entry(i,j));
                                        A.set_block(i*B2.row_dimension() + A1.row_dimension()*BBlock.row_dimension()                                        ,startcolumn + j*B2.column_dimension() + A2.column_dimension()*B1.column_dimension(),B2, false, A2.get_entry(i,j));
                                        A.set_block(i*B3.row_dimension() + A1.row_dimension()*BBlock.row_dimension() + A2.row_dimension()*B1.row_dimension(),startcolumn + j*B3.column_dimension()                                              ,B3, false, A2.get_entry(i,j));
                                        A.set_block(i*B4.row_dimension() + A1.row_dimension()*BBlock.row_dimension() + A2.row_dimension()*B2.row_dimension(),startcolumn + j*B4.column_dimension() + A2.column_dimension()*B3.column_dimension(),B4, false, A2.get_entry(i,j));
                                    }
                                }
                            }
                        } // end of 0j-block
                    } // end of if rowoffset == 0
                    else
                    { // rowoffset > 0
                        if (columnoffset == 0)
                        {   // i0-block
                            // for ai != rowoffset A can be decomposed into two parts
                            if (ai < rowoffset)
                            {
                                // A00 = A^gg A^gw = A1 A2
                                //       A^wg A^ww
                                // able to process as columns
                                // inserting
                                // A1xB1 A1xB2 A2xB1 A2xB2
                                SparseMatrix<double> A1,A2,B1,B2;
                                ABlock.get_block(0,0             ,ABlock.row_dimension(),number_of_gens,A1);
                                ABlock.get_block(0,number_of_gens,ABlock.row_dimension(),ABlock.column_dimension()-number_of_gens,A2);
                                BBlock.get_block(0,0             ,BBlock.row_dimension(),number_of_gens,B1);
                                BBlock.get_block(0,number_of_gens,BBlock.row_dimension(),BBlock.column_dimension()-number_of_gens,B2);

                                assert ((A1.row_dimension() == A2.row_dimension())
                                        && (B1.row_dimension() == B2.row_dimension()));
                                assert (startcolumn == 0);

                                for (unsigned int i=0; i < A1.row_dimension(); ++i)
                                {
                                    for (unsigned int j=0; j < A1.column_dimension(); ++j)
                                    {
                                        A.set_block(startrow + i*B1.row_dimension(),j*B1.column_dimension()                                                                                                ,B1, false, A1.get_entry(i,j));
                                        A.set_block(startrow + i*B2.row_dimension(),j*B2.column_dimension() + A1.column_dimension()*B1.column_dimension()                                                  ,B2, false, A1.get_entry(i,j));
                                    }
                                    for (unsigned int j=0; j < A2.column_dimension(); ++j)
                                    {
                                        A.set_block(startrow + i*B1.row_dimension(),j*B1.column_dimension() + A1.column_dimension()*BBlock.column_dimension()                                              ,B1, false, A2.get_entry(i,j));
                                        A.set_block(startrow + i*B2.row_dimension(),j*B2.column_dimension() + A1.column_dimension()*BBlock.column_dimension() + A2.column_dimension()*B1.column_dimension(),B2, false, A2.get_entry(i,j));
                                    }
                                }
                            }
                            else // i0-block, ai = rowoffset
                            { 
                                // need to split BBlock in 4 parts and ABlock in 2
                                // inserting
                                // A1xB1 A1xB2 A2xB1 A2xB2
                                // A1xB3 A1xB4 A2xB3 A2xB4

                                SparseMatrix<double> A1,A2,B1,B2,B3,B4;
                                ABlock.get_block(0,0,ABlock.row_dimension(),number_of_gens,A1);
                                ABlock.get_block(0,number_of_gens,ABlock.row_dimension(),ABlock.column_dimension()-number_of_gens,A2);
                                BBlock.get_block(0             ,0             ,number_of_gens                          ,number_of_gens                          ,B1);
                                BBlock.get_block(0             ,number_of_gens,number_of_gens                          ,BBlock.column_dimension()-number_of_gens,B2);
                                BBlock.get_block(number_of_gens,0             ,BBlock.column_dimension()-number_of_gens,number_of_gens                          ,B3);
                                BBlock.get_block(number_of_gens,number_of_gens,BBlock.column_dimension()-number_of_gens,BBlock.column_dimension()-number_of_gens,B4);

                                assert ((A1.row_dimension() == A2.row_dimension())
                                        && (B1.row_dimension() == B2.row_dimension())
                                        && (B3.row_dimension() == B4.row_dimension())
                                        && (B1.column_dimension() == B3.column_dimension())
                                        && (B2.column_dimension() == B4.column_dimension()));
                                assert (startcolumn == 0);
                                for (unsigned int i=0; i < A1.row_dimension(); ++i)
                                {
                                    for (unsigned int j=0; j < A1.column_dimension(); ++j)
                                    {
                                        A.set_block(startrow + i*B1.row_dimension()                                        , j*B1.column_dimension()                                                                                             ,B1, false, A1.get_entry(i,j));
                                        A.set_block(startrow + i*B2.row_dimension()                                        , j*B2.column_dimension() + A1.column_dimension()*B1.column_dimension()                                               ,B2, false, A1.get_entry(i,j));
                                        A.set_block(startrow + i*B3.row_dimension() + A1.row_dimension()*B1.row_dimension(), j*B3.column_dimension()                                                                                             ,B3, false, A1.get_entry(i,j));
                                        A.set_block(startrow + i*B4.row_dimension() + A1.row_dimension()*B2.row_dimension(), j*B4.column_dimension() + A1.column_dimension()*B3.column_dimension()                                               ,B4, false, A1.get_entry(i,j));
                                    }
                                    for (unsigned int j=0; j < A2.column_dimension(); ++j)
                                    {
                                        A.set_block(startrow + i*B1.row_dimension()                                        , j*B1.column_dimension() + A1.column_dimension()*BBlock.column_dimension()                                              ,B1, false, A2.get_entry(i,j));
                                        A.set_block(startrow + i*B2.row_dimension()                                        , j*B2.column_dimension() + A1.column_dimension()*BBlock.column_dimension() + A2.column_dimension()*B1.column_dimension(),B2, false, A2.get_entry(i,j));
                                        A.set_block(startrow + i*B3.row_dimension() + A2.row_dimension()*B1.row_dimension(), j*B3.column_dimension() + A1.column_dimension()*BBlock.column_dimension()                                              ,B3, false, A2.get_entry(i,j));
                                        A.set_block(startrow + i*B4.row_dimension() + A2.row_dimension()*B2.row_dimension(), j*B4.column_dimension() + A1.column_dimension()*BBlock.column_dimension() + A2.column_dimension()*B3.column_dimension(),B4, false, A2.get_entry(i,j));
                                    }
                                }
                            }

                        } // end of columnoffset == 0
                        else
                        { // columnoffset > 0, ij-block
                            if (ai == rowoffset)
                            {
                                if (aj == columnoffset)
                                {
                                    // last row & last column
                                    // splitting BBlock in 4 parts necessary
                                    // inserting
                                    // AxB1 AxB2
                                    // AxB3 AxB4
                                    SparseMatrix<double> B1,B2,B3,B4;
                                    BBlock.get_block(0             ,0             ,number_of_gens                       ,number_of_gens                          ,B1);
                                    BBlock.get_block(0             ,number_of_gens,number_of_gens                       ,BBlock.column_dimension()-number_of_gens,B2);
                                    BBlock.get_block(number_of_gens,0             ,BBlock.row_dimension()-number_of_gens,number_of_gens                          ,B3);
                                    BBlock.get_block(number_of_gens,number_of_gens,BBlock.row_dimension()-number_of_gens,BBlock.column_dimension()-number_of_gens,B4);

                                    assert ((B1.row_dimension() == B2.row_dimension())
                                         && (B3.row_dimension() == B4.row_dimension())
                                         && (B1.column_dimension() == B3.column_dimension())
                                         && (B2.column_dimension() == B4.column_dimension()));

                                    for (unsigned int i=0; i < ABlock.row_dimension(); ++i)
                                    {
                                        for (unsigned int j=0; j < ABlock.column_dimension(); ++j)
                                        {
                                            A.set_block(startrow + i*B1.row_dimension()                                            , startcolumn + j*B1.column_dimension()                                                   ,B1, false, ABlock.get_entry(i,j));
                                            A.set_block(startrow + i*B2.row_dimension()                                            , startcolumn + j*B2.column_dimension() + ABlock.column_dimension()*B1.column_dimension() ,B2, false, ABlock.get_entry(i,j));
                                            A.set_block(startrow + i*B3.row_dimension() + ABlock.row_dimension()*B1.row_dimension(), startcolumn + j*B3.column_dimension()                                                   ,B3, false, ABlock.get_entry(i,j));
                                            A.set_block(startrow + i*B4.row_dimension() + ABlock.row_dimension()*B2.row_dimension(), startcolumn + j*B4.column_dimension() + ABlock.column_dimension()*B3.column_dimension() ,B4, false, ABlock.get_entry(i,j));
                                        }
                                    }
                                }
                                else
                                {
                                    // last row (without the last box)
                                    // vertical splitting of BBlock necessary
                                    // inserting
                                    // AxB1
                                    // AxB2
                                    SparseMatrix<double> B1,B2;
                                    BBlock.get_block(0,0,             number_of_gens,BBlock.column_dimension(),B1);
                                    BBlock.get_block(number_of_gens,0,BBlock.row_dimension()-number_of_gens,BBlock.column_dimension(),B2);

                                    assert(B1.column_dimension() == B2.column_dimension());
                                    for (unsigned int i=0; i < ABlock.row_dimension(); ++i)
                                    {
                                        for (unsigned int j=0; j < ABlock.column_dimension(); ++j)
                                        {
                                            A.set_block(startrow + i*B1.row_dimension()                                            , startcolumn + j*B1.column_dimension(),B1, false, ABlock.get_entry(i,j));
                                            A.set_block(startrow + i*B2.row_dimension() + ABlock.row_dimension()*B1.row_dimension(), startcolumn + j*B2.column_dimension(),B2, false, ABlock.get_entry(i,j));
                                        }
                                    }
                                }
                            }
                            else
                            {
                                if (aj == columnoffset)
                                {
                                    // last column, 
                                    // BBlock needs to be splitted horizontally
                                    // inserting
                                    // AxB1 AxB2
                                    SparseMatrix<double> B1,B2;
                                    BBlock.get_block(0,0,BBlock.row_dimension(), number_of_gens,B1);
                                    BBlock.get_block(0,number_of_gens,BBlock.row_dimension(),BBlock.column_dimension()-number_of_gens,B2);
                                    assert(B1.row_dimension() == B2.row_dimension());
                                    for (unsigned int i=0; i < ABlock.row_dimension(); ++i)
                                    {
                                        for (unsigned int j=0; j < ABlock.column_dimension(); ++j)
                                        {
                                            A.set_block(startrow + i*B1.row_dimension() , startcolumn + j*B1.column_dimension()                                                  ,B1, false, ABlock.get_entry(i,j));
                                            A.set_block(startrow + i*B2.row_dimension() , startcolumn + j*B2.column_dimension() + ABlock.column_dimension()*B1.column_dimension(),B2, false, ABlock.get_entry(i,j));
                                        }
                                    }
                                }
                                else
                                {
                                    // a normal box, no splitting necessary
                                    for (unsigned int i=0; i < ABlock.row_dimension(); ++i)
                                    {
                                        for (unsigned int j=0; j < ABlock.column_dimension(); ++j)
                                        {
                                            A.set_block(startrow + i*BBlock.row_dimension() , startcolumn + j*BBlock.column_dimension(),BBlock, false, ABlock.get_entry(i,j));
                                        }
                                    }
                                }
                            }
                        }
                    } // end of rowoffset > 0

#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                    cout << "added block (ro,co) = (" << rowoffset << ", " << columnoffset << "), (ai,aj) = (" << ai << ", " << aj << ") at position (lro,lco) = (" << startrow << ", " << startcolumn << ")" << endl;
#endif

                    startcolumn+=ABlock.column_dimension()*BBlock.column_dimension();
                    if ((aj == columnoffset)&&(columnoffset == offset))
                    {
                        // we have added the last block on the current ai line and will increase ai
                        startrow += ABlock.row_dimension()*BBlock.row_dimension();
                    }
                } // end of for(aj ...)
                // iteration über aj beendet
                //startcolumn
            }
            //startrow=startrow; // hiernach wird ai erhoeht!
            //startrow+=
        }
    }
};

//using MathTL::SimpleSturmBVP;
//using MathTL::CG;

#if 0
int main()
{
    const int d  = _D;
    const int dT = _DT;

    const int dim = _DIMENSION;
    //const int offset = _OFFSET; // Radius for the index set Lambda


#if _PRECONDITIONED
    const char* precondition_string = "pcd";
    const char* precondition_path = "preconditioned";
#else
    const char* precondition_string = "npcd";
    const char* precondition_path = "unpreconditioned";
#endif
    //typedef DSBasis<d,dT> Basis1d;
    typedef PBasis<d,dT> Basis1d;
    const char* basis_string="primbs";

    typedef TensorBasis<Basis1d,dim> Basis;
    typedef Basis::Index Index;
    typedef Index::level_type index_lt;

    PoissonBVP<dim> poisson(0);
    const char* problem_string="laplacian";
    
    IdentityBVP<dim> identity(0);
    const char* identity_string="gramian";

    FixedArray1D<bool,(2*dim)> bc;
#if _DIMENSION == 1
#if _BOUNDARY_CONDITIONS == 0
    const char* bc_string ={"tt"};
    bc[0] = bc[1] = true;
#elif _BOUNDARY_CONDITIONS == 1
    const char* bc_string ={"tf"};
    bc[0] = true; bc[1] = false;
#elif _BOUNDARY_CONDITIONS == 2
    const char* bc_string ={"ft"};
    bc[0] = false; bc[1] = true;
#elif _BOUNDARY_CONDITIONS == 3
    const char* bc_string ={"ff"};
    bc[0] = bc[1] = false;
#endif
#else
#if _BOUNDARY_CONDITIONS == 0
    bc[0] = bc[1] = bc[2] = bc[3] = true;
    const char* bc_string = {"tttt"};
    const char* a_bcstring ={"tt"};
    const char* b_bcstring ={"tt"};
#elif _BOUNDARY_CONDITIONS == 1
    bc[0] = bc[1] = bc[2] = bc[3] = false;
    const char* bc_string = {"ffff"};
    const char* a_bcstring ={"ff"};
    const char* b_bcstring ={"ff"};
#endif
#endif

    const char* gramianStorageFolder = {"/import/shared/friedrich/source/precomputed/gramian"};
    const char* laplacianStorageFolder = {"/import/shared/friedrich/source/precomputed/laplacian"};
    const char* haar_wav_gramianStorageFolder = {"/import/shared/friedrich/source/precomputed/haar_wav_gramian"};
    const char* haar_gen_gramianStorageFolder = {"/import/shared/friedrich/source/precomputed/haar_gen_gramian"};
    //const char* haar_wav_laplacianStorageFolder = {"/home/friedrich/source/precomputed_matrices/haar_wav_laplacian"};
    //const char* haar_gen_laplacianStorageFolder = {"/home/friedrich/source/precomputed_matrices/haar_gen_laplacian"};
    const char* haar_wav_laplacianStorageFolder = {"/import/shared/friedrich/source/precomputed/haar_wav_laplacian"};
    const char* haar_gen_laplacianStorageFolder = {"/import/shared/friedrich/source/precomputed/haar_gen_laplacian"};
    const char* transitionMatrixStorageFolder = {"/import/shared/friedrich/source/precomputed/transition_matrices"};
    const char* matrixBlocksStorageFolder = {"/import/shared/friedrich/source/precomputed/matrix_blocks"};
    const char* rhsstorageFolder = {"/import/shared/friedrich/source/precomputed/functions"};
    //const char* rhsstorageFolder = {"/home/friedrich/source/precomputed/functions"};

    clock_t tstart, tend;
    double time;

    cout << "main:: setting up equation" << endl;
    TensorEquation<Basis1d,dim,Basis> te_pois(&poisson, bc, false);
    TensorEquation<Basis1d,dim,Basis> te_gram(&identity, bc, false);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_pois(&te_pois);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_gram(&te_gram);
/*
    eq.basis_.set_jmax(10);
    for (unsigned int i=0; i < 4000; ++i)
        cout << " i = " << i << " index = " << *(eq.basis().get_wavelet(i)) << endl;
    //cout << " 100  = " << *(eq.basis().get_wavelet(100)) << endl
    //     << " 1044 = " << *(eq.basis().get_wavelet(1044)) << endl
    //     << " 4000 = " << *(eq.basis().get_wavelet(4000)) << endl;
    abort();
  */
    for (int offset=_LAPLACIAN_MATRIX_OFFSET_START;offset <= _LAPLACIAN_MATRIX_OFFSET_END; ++offset)
    {
        const int jmax = multi_degree(ctp_pois.basis().j0()) + offset;
        //eq.basis_.set_jmax(jmax); // nur in 1D nötig
        te_pois.basis_.set_jmax(jmax);

        char matrix_filename[250];
        //laplacian_primbs_d_dt_3_3_bc_tttf_jmax_3
        
        sprintf(matrix_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",laplacianStorageFolder,precondition_path,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);

        //cout << "main:: setting up problem" << endl;
        //CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem(&eq,5.5,35.0);

        cout << "main:: compute stiffness matrix int grad psi_lambda grad psi_mu" << endl;

        tstart = clock();

#if _DIMENSION == 1
        /* setup and fill index set */
        
        set<Index> Lambda_IndexSet;
        for (Index lambda = ctp_pois.basis().first_generator();; ++lambda)
        {
            Lambda_IndexSet.insert(lambda);
            if (lambda == ctp_pois.basis().last_wavelet(jmax)) break;
        }
        //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());

        cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;

        /* setup stiffness matrix A */
        SparseMatrix<double> A;        /* stiffness matrix*/
        setup_stiffness_matrix(ctp_pois, Lambda_IndexSet, A, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */
#else // dim == 2
        cout << "compose laplacian matrix out of 1D matrices" << endl;
        unsigned int dof = te_pois.basis().last_wavelet(jmax).number()+1;
        cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << dof << endl;
        char a_namestart [250], b_namestart[250], a_nameend[250], b_nameend[250];
        SparseMatrix<double> A(dof,dof);
        //Blockfilename = laplacian_primbs_d_dt_3_3_bc_tt_npcd_block_2_10
        //                gramian_primbs_d_dt_3_3_bc_tf_npcd_block_2_10
        sprintf(a_namestart,"%s/%s/%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,problem_string,basis_string,d,dT);
        sprintf(b_namestart,"%s/%s/%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
        sprintf(a_nameend,"_%s_block_",precondition_string);
        sprintf(b_nameend,"_%s_block_",precondition_string);
        compose_matrix(A,a_namestart,a_bcstring,a_nameend,
                          b_namestart,b_bcstring,b_nameend,
                          offset,
                          te_pois.basis().bases()[0]->Deltasize(te_pois.basis().j0()[0]));

        SparseMatrix<double> A2(dof,dof);
        //Blockfilename = gramian_primbs_d_dt_3_3_bc_tt_npcd_block_2_10
        //                laplacian_primbs_d_dt_3_3_bc_tf_npcd_block_2_10
        sprintf(a_namestart,"%s/%s/%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
        sprintf(b_namestart,"%s/%s/%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,problem_string,basis_string,d,dT);
        sprintf(a_nameend,"_%s_block_",precondition_string);
        sprintf(b_nameend,"_%s_block_",precondition_string);
        compose_matrix(A2,a_namestart,a_bcstring,a_nameend,
                          b_namestart,b_bcstring,b_nameend,
                          offset,
                          te_pois.basis().bases()[0]->Deltasize(te_pois.basis().j0()[0]));

        //cout << "precompute_stuff:: laplacian makes trouble when adding matrices!" << endl;
        //cout << "A = " << A << endl;
        //cout << "A2 = " << A2 << endl;
        
        A.add(1.0,A2);

#endif
        cout << "write to file: " << matrix_filename << endl;
        A.matlab_output(matrix_filename, "A", 1);
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;

    }
    

    for (int offset=_GRAMIAN_MATRIX_OFFSET_START;offset <= _GRAMIAN_MATRIX_OFFSET_END; ++offset)
    {
        const int jmax = multi_degree(ctp_gram.basis().j0()) + offset;
        te_gram.basis_.set_jmax(jmax);
        char matrix_filename[250];
        //gramian_primbs_d_dt_3_3_bc_tttf_jmax_3
        sprintf(matrix_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",gramianStorageFolder,precondition_path,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);
        cout << "main:: compute stiffness matrix int psi_lambda psi_mu" << endl;

        tstart = clock();
#if _DIMENSION == 1
        /* setup and fill index set */
        set<Index> Lambda_IndexSet;
        for (Index lambda = ctp_gram.basis().first_generator();; ++lambda)
        {
            Lambda_IndexSet.insert(lambda);
            if (lambda == ctp_gram.basis().last_wavelet(jmax)) break;
        }
        //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());

        cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;

        /* setup stiffness matrix A */
        SparseMatrix<double> A;        /* stiffness matrix*/
        setup_stiffness_matrix(ctp_gram, Lambda_IndexSet, A, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */

#else // dim == 2
        cout << "compose gramian matrix out of 1D matrices" << endl;
        unsigned int dof = te_gram.basis().last_wavelet(jmax).number()+1;
        cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << dof << endl;
        char a_namestart [250], a_nameend[250];
        SparseMatrix<double> A(dof,dof);
        //Blockfilename = gramian_primbs_d_dt_3_3_bc_tt_npcd_block_2_10
        //                gramian_primbs_d_dt_3_3_bc_tf_npcd_block_2_10
        sprintf(a_namestart,"%s/%s/%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
        sprintf(a_nameend,"_%s_block_",precondition_string);
        
        // caution! a_name=b_name, but only in this case!
        compose_matrix(A,a_namestart,a_bcstring,a_nameend,
                         a_namestart,a_bcstring,a_nameend,
                         offset,
                         te_gram.basis().bases()[0]->Deltasize(te_gram.basis().j0()[0]));
#endif
        cout << "write to file: " << matrix_filename << endl;
        A.matlab_output(matrix_filename, "A", 1);
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;
    }
    
    //_RHS_OFFSET_START
    //for (int offset=_RHS_OFFSET_START;offset <= _RHS_OFFSET_END; ++offset)
    if (_RHS_OFFSET_START <= _RHS_OFFSET_END)
    {
        cout << "main:: compute vector int 1/2 psi_lambda. Basis functions are unpreconditioned." << endl;
        if (_RHS_OFFSET_START < _RHS_OFFSET_END)
        {
            cout << "observe that only generators are active in the expansion of a constant function. Thus the RHS is computed only once." << endl;
        }
        //const int jmax = multi_degree(ctproblem.basis().j0()) + offset;
        const int jmax = multi_degree(ctp_gram.basis().j0());
        const char* function_string = "onehalf";
        te_gram.basis_.set_jmax(jmax);

        char vector_filename[250];
        //onehalf_primbs_d_dt_3_3_bc_tttf_jmax_3
        //observe: only generator coefficients are active in the expansion of 1/2! (at least for primbs 3 3)

        
        //sprintf(vector_filename, "%s/%s_%s_d_dt_%d_%d_bc_%s_jmax_%d.iv",storageFolder,function_string,basis_string,d,dT,bc_string,jmax);
        sprintf(vector_filename, "%s/%s_%s_d_dt_%d_%d_bc_%s.iv",rhsstorageFolder,function_string,basis_string,d,dT,bc_string);

        tstart = clock();

        ConstantFunction<dim> constant_fkt(Vector<double>(1, "0.5"));
        InfiniteVector<double,Index> coeffs;
        //eq.basis().expand(&constant_fkt, false, multi_degree(eq.basis().j0())+offset, coeffs);
        te_gram.basis().expand(&constant_fkt, false, multi_degree(te_gram.basis().j0()), coeffs);

        //cout << "number of coefficients before coarsening: " << coeffs.size() << endl;
        cout << "coarsening ..." << endl;
        coeffs.compress(1e-13);
        cout << "number of coefficients after coarsening: " << coeffs.size() << endl;
        cout << "write to file: " << vector_filename << endl;
        cout << "coeffs = " << coeffs << endl;
        cout << "" << endl;
        
        cout << " THERE IS NO WRITE TO FILE METHOD AT THE MOMENT" << endl;
        abort();
        //coeffs.writeToFile(vector_filename);
        
        
        //coeffs.readFromFile(vector_filename);
        //cout << "coeffs = " << coeffs << endl;
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;
    }




    for (int haar_level=_HAAR_GEN_MATRIX_HAAR_LEVEL_START;haar_level <= _HAAR_GEN_MATRIX_HAAR_LEVEL_END; ++haar_level)
    {
#if _DIMENSION == 1
        // Generators, numbering and count is for the 1D case!
        for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
        {

            HaarGeneratorBVP haarbvp(0,eta);
            cout << "HaarGeneratorBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
            //Point<1> poi(9.0/16.0);
            //cout << "haarbvp.q(9/16) = " << haarbvp.q(poi) << endl;
            /*
            for (int k=0; k<= 4*(1<< haarbvp.j);++k)
            {
                cout << ((double)k/(4*(1<< haarbvp.j))) << " ";
            }
            cout << endl;
            for (int k=0; k<= 4*(1<< haarbvp.j);++k)
            {
                Point<1> poi((double)k/(4*(1<< haarbvp.j)));
                cout << haarbvp.q(poi) << " ";
            }
            cout << endl;
            continue;
            */
            // const char* identity_string="haar"; // =gramian. see matrix_filename below
            TensorEquation<Basis1d,dim,Basis> te_haargram(&haarbvp, bc, false);
            CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_haargram(&te_haargram);
            cout << "main:: compute stiffness matrices int g_eta psi_lambda psi_mu" << endl;
            for (int spatial_offset=_HAAR_GEN_MATRIX_SPATIAL_OFFSET_START;spatial_offset <= _HAAR_GEN_MATRIX_SPATIAL_OFFSET_END; ++spatial_offset)
            {
                const int jmax = multi_degree(ctp_haargram.basis().j0()) + spatial_offset;
                //te_gram.basis_.set_jmax(jmax);
                te_haargram.basis_.set_jmax(jmax);

                char matrix_filename[250];
                //1D haar_gen_2_gramian_primbs_d_dt_3_3_bc_tt_jmax_3_npcd
                //2D haar_gen_2_gramian_primbs_d_dt_3_3_bc_tttt_jmax_3_npcd
                
                sprintf(matrix_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_gen_gramianStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                tstart = clock();

                /* setup and fill index set */
                set<Index> Lambda_IndexSet;
                for (Index lambda = ctp_haargram.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == ctp_haargram.basis().last_wavelet(jmax)) break;
                }
                //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());

                cout << "eta = " << eta << " dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;

                /* setup stiffness matrix A */
                SparseMatrix<double> A;        /* stiffness matrix*/

                setup_stiffness_matrix(ctp_haargram, Lambda_IndexSet, A, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */

                cout << "write to file: " << matrix_filename << endl;
                A.matlab_output(matrix_filename, "A", 1);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
            }
        }
#else
        cout << "This computation doesn't make sense and hence isn't done" << endl;
        // matrices are used nowhere
#endif
    }

#if _HAAR_WAV_MATRIX_FIRST_ETA < 0
    for (int haar_level=_HAAR_WAV_MATRIX_HAAR_LEVEL_START;haar_level <= _HAAR_WAV_MATRIX_HAAR_LEVEL_END; ++haar_level)
    {
        cout << "main:: begin computation of haar_gramian on haar_level = " << haar_level << endl;
        clock_t tstart2 = clock();
        clock_t tend2;
        int firsteta = (_DIMENSION == 1 ? (1<<haar_level) : (1<<(haar_level<<1)));
        int lasteta = (_DIMENSION == 1 ? (1<<(haar_level+1)) : (1<<((haar_level<<1) + 2 ))) -1;
        if (haar_level == 0)
            firsteta = 0;
#else
        int firsteta = _HAAR_WAV_MATRIX_FIRST_ETA;
        int lasteta = _HAAR_WAV_MATRIX_LAST_ETA;
#endif
        // Wavelets, numbering and count is for the 1D case!
        //cout << " hl= " << haar_level << " hl<<1 = " << (haar_level<<1) << " (hl<<1) +2 = " << ((haar_level<<1) + 2 ) << " " << (1<<((haar_level<<1) + 2 )) << endl;
        for (int eta = firsteta; eta <=lasteta ;++eta)
        {
#if _DIMENSION == 1
            HaarWaveletBVP haarbvp(0,eta);
            cout << "HaarWaveletBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
            // Point<1> poi(9.0/16.0);
            // cout << "haarbvp.q(9/16) = " << haarbvp.q(poi) << endl;
            /*
            for (int k=0; k<= 4*(1<< haarbvp.j);++k)
            {
                cout << ((double)k/(4*(1<< haarbvp.j))) << " ";
            }
            cout << endl;
            for (int k=0; k<= 4*(1<< haarbvp.j);++k)
            {
                Point<1> poi((double)k/(4*(1<< haarbvp.j)));
                cout << haarbvp.q(poi) << " ";
            }
            cout << endl;
            continue;
            */
            // const char* identity_string="haar"; // =gramian. see matrix_filename below
            TensorEquation<Basis1d,dim,Basis> te_haarwavgram(&haarbvp, bc, false);
            CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_haargram(&te_haarwavgram);
            cout << "main:: compute stiffness matrix int h_eta psi_lambda psi_mu" << endl;
#else
            
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
            cout << "compose haar_wav stiffness matrices int h_eta psi_lambda psi_mu out of 1D matrices" << endl;
#endif
            char a_namestart [250], a_nameend[250], b_namestart [250], b_nameend[250];
            //Blockfilename = haar_wav_17_gramian_primbs_d_dt_3_3_bc_ff_nprc_block_0_2
            //                haar_gen_23_gramian_primbs_d_dt_3_3_bc_ff_nprc_block_0_2
            if (eta == 0)  // gen x gen on lvl 0
            {
                sprintf(a_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
                sprintf(b_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "generator!" << endl;
#endif
            }
            else
            {
                int n(eta),l=0;
                while ((n>>=1)>0) ++l; // log_2 of eta
                int j = l/2; // 2d level j of eta
                int type = eta / (1<<(j<<1)); // 2d type: 1 = gen x wav, 2 = wav x gen, 3 = wav x wav
                int k = eta - type*(1<<(j<<1)); // number in current type
                int xnum (k / (1<<j)), ynum (k % (1<<j)); // number of gen or wav in x and y direction in the current level j. range is {0,...,2^j-1}. caution: numbering in HaarGenBVP begins at 1 (as does it in HaarWavBVP)
                xnum+=(1<<j);ynum+=(1<<j); // number (first gen or wav on first level has number 1)
                switch (type)
                {
                    case 1:
                        sprintf(a_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        break;
                    case 2:
                        sprintf(a_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        break;
                    case 3:
                        sprintf(a_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        break;
                }
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "2d eta = " << eta << " results in j = " << j << " type = " << type << " (xnum, ynum) = (" << xnum << ", " << ynum << ")" << endl;
#endif
            }
            sprintf(a_nameend,"_%s_block_",precondition_string);
            sprintf(b_nameend,"_%s_block_",precondition_string);
#endif
            for (int spatial_offset=_HAAR_WAV_MATRIX_SPATIAL_OFFSET_START;spatial_offset <= _HAAR_WAV_MATRIX_SPATIAL_OFFSET_END; ++spatial_offset)
            {
                //cout << "main:: setting up problem" << endl;
                //CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem(&eq,5.5,35.0);
                tstart = clock();
#if _DIMENSION == 1
                const int jmax = multi_degree(ctp_haargram.basis().j0()) + spatial_offset;
                te_haarwavgram.basis_.set_jmax(jmax);
                //haarwa.basis_.set_jmax(jmax);
                /* setup and fill index set */
                set<Index> Lambda_IndexSet;
                for (Index lambda = ctp_haargram.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == ctp_haargram.basis().last_wavelet(jmax)) break;
                }
                //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());

                cout << "eta = " << eta << " dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;

                /* setup stiffness matrix A */
                SparseMatrix<double> A;        /* stiffness matrix*/

                setup_stiffness_matrix(ctp_haargram, Lambda_IndexSet, A, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */
#else
                const int jmax = multi_degree(ctp_pois.basis().j0()) + spatial_offset;
                te_pois.basis_.set_jmax(jmax);
                unsigned int dof = te_pois.basis().last_wavelet(jmax).number()+1;
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << dof << endl;
#endif
                
                SparseMatrix<double> A(dof,dof);
                //Blockfilename = haar_gen_48_gramian_primbs_d_dt_3_3_bc_tt_npcd_block_5_10
                //
                // compute type of haar-wavelet number eta
                compose_matrix(A,a_namestart,a_bcstring,a_nameend,
                                 b_namestart,b_bcstring,b_nameend,
                                 spatial_offset,
                                 te_pois.basis().bases()[0]->Deltasize(te_pois.basis().j0()[0]));
#endif
                char matrix_filename[250];
                //haar_wav_2_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                sprintf(matrix_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_wav_gramianStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);

                cout << "write to file: " << matrix_filename << endl;
                A.matlab_output(matrix_filename, "A", 1);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
            }
        }
#if _HAAR_WAV_MATRIX_FIRST_ETA < 0
        tend2 = clock();
        time = (double)(tend2-tstart2)/CLOCKS_PER_SEC;
        cout << "  computation for haar_level = " << haar_level << " done. Time needed: " << time << " seconds" << endl;
    }
#else
#endif
    
    // -----------------------------------------
    // -----------------------------------------
    // \int haar_gen_eta grad psi_lambda grad psi_mu
    // analogous to the previous blocks
    for (int haar_level=_HAAR_GEN_LAPLACE_HAAR_LEVEL_START;haar_level <= _HAAR_GEN_LAPLACE_HAAR_LEVEL_END; ++haar_level)
    {
#if _DIMENSION == 1
        // Generators, numbering and count is for the 1D case!
        for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
        {
            HaarGeneratorPoissonBVP haarbvp(0,eta);
            cout << "HaarGeneratorPoissonBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
            // const char* identity_string="haar"; // =gramian. see matrix_filename below
            TensorEquation<Basis1d,dim,Basis> te_haargenlap(&haarbvp, bc, false);
            CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_haargenlap(&te_haargenlap);
            cout << "main:: compute stiffness matrices int g_eta grad psi_lambda grad psi_mu" << endl;
            for (int spatial_offset=_HAAR_GEN_LAPLACE_SPATIAL_OFFSET_START;spatial_offset <= _HAAR_GEN_LAPLACE_SPATIAL_OFFSET_END; ++spatial_offset)
            {
                const int jmax = multi_degree(ctp_haargenlap.basis().j0()) + spatial_offset;
                //te_pois.basis_.set_jmax(jmax);
                te_haargenlap.basis_.set_jmax(jmax);

                char matrix_filename[250];
                //1D haar_gen_2_laplacian_primbs_d_dt_3_3_bc_tt_jmax_3_npcd
                //2D haar_gen_2_laplacian_primbs_d_dt_3_3_bc_tttt_jmax_3_npcd
                                
                sprintf(matrix_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_gen_laplacianStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                tstart = clock();

                /* setup and fill index set */
                set<Index> Lambda_IndexSet;
                for (Index lambda = te_haargenlap.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == te_haargenlap.basis().last_wavelet(jmax)) break;
                }
                //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());
                cout << "eta = " << eta << " dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;
                /* setup stiffness matrix A */
                SparseMatrix<double> A;        /* stiffness matrix*/
                setup_stiffness_matrix(ctp_haargenlap, Lambda_IndexSet, A, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */
                cout << "write to file: " << matrix_filename << endl;
                A.matlab_output(matrix_filename, "A", 1);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
            }
        }
#else
        cout << "This computation doesn't make sence and hence isn't done" << endl;
        // matrices are used nowhere
#endif
    }
    
    // -----------------------------------------
    // -----------------------------------------
    // same with Haar wavelets

    
#if _HAAR_WAV_LAPLACE_FIRST_ETA < 0
    for (int haar_level=_HAAR_WAV_LAPLACE_HAAR_LEVEL_START;haar_level <= _HAAR_WAV_LAPLACE_HAAR_LEVEL_END; ++haar_level)
    {
        cout << "main:: begin computation of haar_laplacian on haar_level = " << haar_level << endl;
        clock_t tstart2 = clock();
        clock_t tend2;
        int firsteta = (_DIMENSION == 1 ? (1<<haar_level) : (1<<(haar_level<<1)));
        int lasteta = (_DIMENSION == 1 ? (1<<(haar_level+1)) : (1<<((haar_level<<1) + 2 ))) -1;
        if (haar_level == 0)
            firsteta = 0;
#else
        int firsteta = _HAAR_WAV_LAPLACE_FIRST_ETA;
        int lasteta = _HAAR_WAV_LAPLACE_LAST_ETA;
#endif
        // Wavelets, numbering and count is for the 1D case!
        //cout << " hl= " << haar_level << " hl<<1 = " << (haar_level<<1) << " (hl<<1) +2 = " << ((haar_level<<1) + 2 ) << " " << (1<<((haar_level<<1) + 2 )) << endl;
        for (int eta = firsteta; eta <=lasteta ;++eta)
        {
#if _DIMENSION == 1
            HaarWaveletPoissonBVP haarbvp(0,eta);
            cout << "HaarWaveletPoissonBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
            // const char* identity_string="haar"; // =gramian. see matrix_filename below
            TensorEquation<Basis1d,dim,Basis> te_haarwavlap(&haarbvp, bc, false);
            CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctp_haarwavlap(&te_haarwavlap);
            cout << "main:: compute stiffness matrix int h_eta psi_lambda psi_mu" << endl;
#else
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
            cout << "compose haar_wav stiffness matrices int h_eta grad psi_lambda grad psi_mu out of 1D matrices" << endl;
#endif
            char a1_namestart [250], a2_namestart [250], a_nameend[250], b1_namestart [250], b2_namestart [250], b_nameend[250];
            //Blockfilename = haar_wav_17_gramian_primbs_d_dt_3_3_bc_ff_nprc_block_0_2
            //                haar_gen_23_gramian_primbs_d_dt_3_3_bc_ff_nprc_block_0_2
            if (eta == 0)  // gen x gen on lvl 0
            {
                sprintf(a1_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,problem_string,basis_string,d,dT);
                sprintf(b1_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
                sprintf(a2_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT);
                sprintf(b2_namestart,"%s/%s/haar_gen_1_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,problem_string,basis_string,d,dT);
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "generator!" << endl;
#endif
            }
            else
            {
                int n(eta),l=0;
                while ((n>>=1)>0) ++l; // log_2 of eta
                int j = l/2; // 2d level j of eta
                int type = eta / (1<<(j<<1)); // 2d type: 1 = gen x wav, 2 = wav x gen, 3 = wav x wav
                int k = eta - type*(1<<(j<<1)); // number in current type
                int xnum (k / (1<<j)), ynum (k % (1<<j)); // number of gen or wav in x and y direction in the current level j. range is {0,...,2^j-1}. caution: numbering in HaarGenBVP begins at 1 (as does it in HaarWavBVP)
                xnum+=(1<<j);ynum+=(1<<j); // number (first gen or wav on first level has number 1)
                //cout << "eta = " << eta << ":: j=" << j << "; type=" << type << "; k=" << k << "; xnum=" << xnum << "; ynum=" << ynum << endl;
                switch (type)
                {
                    case 1:
                        sprintf(a1_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,problem_string,basis_string,d,dT);
                        sprintf(b1_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        sprintf(a2_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b2_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,problem_string,basis_string,d,dT);
                        break;
                    case 2:
                        sprintf(a1_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,problem_string,basis_string,d,dT);
                        sprintf(b1_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        sprintf(a2_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b2_namestart,"%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,problem_string,basis_string,d,dT);
                        break;
                    case 3:
                        sprintf(a1_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,problem_string,basis_string,d,dT);
                        sprintf(b1_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,identity_string,basis_string,d,dT);
                        sprintf(a2_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,xnum,identity_string,basis_string,d,dT);
                        sprintf(b2_namestart,"%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_",matrixBlocksStorageFolder,precondition_path,ynum,problem_string,basis_string,d,dT);
                        break;
                }
                cout << "2d eta = " << eta << " results in j = " << j << " type = " << type << " (xnum, ynum) = (" << xnum << ", " << ynum << ")" << endl;
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "2d eta = " << eta << " results in j = " << j << " type = " << type << " (xnum, ynum) = (" << xnum << ", " << ynum << ")" << endl;
#endif
            }
            sprintf(a_nameend,"_%s_block_",precondition_string);
            sprintf(b_nameend,"_%s_block_",precondition_string);
#endif
            for (int spatial_offset=_HAAR_WAV_LAPLACE_SPATIAL_OFFSET_START;spatial_offset <= _HAAR_WAV_LAPLACE_SPATIAL_OFFSET_END; ++spatial_offset)
            {
                
                 //cout << "main:: setting up problem" << endl;
                //CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem(&eq,5.5,35.0);
                tstart = clock();
#if _DIMENSION == 1
                const int jmax = multi_degree(ctp_haarwavlap.basis().j0()) + spatial_offset;
                te_haarwavlap.basis_.set_jmax(jmax);
                //haarwa.basis_.set_jmax(jmax);
                /* setup and fill index set */
                set<Index> Lambda_IndexSet;
                for (Index lambda = ctp_haarwavlap.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == ctp_haarwavlap.basis().last_wavelet(jmax)) break;
                }
                //assert(eq.basis().degrees_of_freedom() == Lambda_IndexSet.size());
                cout << "eta = " << eta << " dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << Lambda_IndexSet.size() << endl;
                /* setup stiffness matrix A */
                SparseMatrix<double> A1;        /* stiffness matrix*/
                setup_stiffness_matrix(ctp_haarwavlap, Lambda_IndexSet, A1, _PRECONDITIONED);   /* in galerkin_utils (this will take some time: jmax 17 takes 1h) */
#else
                const int jmax = multi_degree(ctp_pois.basis().j0()) + spatial_offset;
                te_pois.basis_.set_jmax(jmax);
                unsigned int dof = te_pois.basis().last_wavelet(jmax).number()+1;
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 1
                cout << "dim = " << dim << " jmax = " << jmax << " degrees_of_freedom = " << dof << endl;
#endif
                SparseMatrix<double> A1(dof,dof),A2(dof,dof);
                //Blockfilename = haar_gen_48_gramian_primbs_d_dt_3_3_bc_tt_npcd_block_5_10
                //
                // compute type of haar-wavelet number eta
                compose_matrix(A1,a1_namestart,a_bcstring,a_nameend,
                                 b1_namestart,b_bcstring,b_nameend,
                                 spatial_offset,
                                 te_pois.basis().bases()[0]->Deltasize(te_pois.basis().j0()[0]));
                compose_matrix(A2,a2_namestart,a_bcstring,a_nameend,
                                 b2_namestart,b_bcstring,b_nameend,
                                 spatial_offset,
                                 te_pois.basis().bases()[0]->Deltasize(te_pois.basis().j0()[0]));
                A1.add(1.0,A2);
#endif
                char matrix_filename[250];
                //haar_wav_2_laplacian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                //sprintf(matrix_filename, "/home/friedrich/source/precomputed_matrices/haar_wav_laplacian/unpreconditioned/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",eta,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                sprintf(matrix_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_wav_laplacianStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
               
                cout << "write to file: " << matrix_filename << endl;
                A1.matlab_output(matrix_filename, "A", 1);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
            }
        }
#if _HAAR_WAV_LAPLACE_FIRST_ETA < 0
        tend2 = clock();
        time = (double)(tend2-tstart2)/CLOCKS_PER_SEC;
        cout << "  computation for haar_level = " << haar_level << " done. Time needed: " << time << " seconds" << endl;
    }
#else
#endif
    // -----------------------------------------
    // -----------------------------------------
    
    
    
    
    // basis transition matrices

    if ( (_BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_END >= _BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_START)
       &&(_BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END >= _BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_START) )
    {

        assert (_PRECONDITIONED == 0);


        int numberofhaarwavelets, numberofspatialwavelets;


        if (_DIMENSION == 1)
        {
            numberofhaarwavelets = 1<<(_BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_END+1);
        }
        else
        {
            numberofhaarwavelets = (1<<((_BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_END<<1) + 2 ));
        }
        numberofspatialwavelets = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+_BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END).number()+1;
        /* setup the biggest matrix A */
        SparseMatrix<double> A( numberofhaarwavelets ,  numberofspatialwavelets  );
        cout << "A.row_dimension() = " << A.row_dimension() << " A.column_dimension() = " << A.column_dimension() << endl;
        //A.set_entry(A.row_dimension(), A.column_dimension(), 11 );
        cout << "main:: compute basis transition matrix int h_eta psi_mu" << endl;
#if _PRECOMPUTE_STUFF_COMPOSE_VERBOSITY > 0
        cout << " numberofhaarwavelets = " << numberofhaarwavelets << " numberofspatialwavelets = " << numberofspatialwavelets << endl;
#endif
        tstart = clock();
        for (int eta = 0; eta < numberofhaarwavelets;++eta)
        {
            // compose_matrix assumes anisotropic tensor structure for both bases, here this is only fulfilled for one basis...

            InfiniteVector<double,Index> row;
            if (eta == 0)
            {
                //cout << "eta = " << eta << " doing generator row" << endl;
                ConstantFunction<dim> constant_fkt(Vector<double>(1, "1"));
                te_pois.basis().expand(&constant_fkt, false, multi_degree(te_pois.basis().j0())+_BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END, row);
                //cout << "coarsening ..." << endl;
                row.compress(1e-13);
            }
            else
            {
#if _DIMENSION == 1
                //cout << "eta = " << eta << " doing wavelet row" << endl;
                HaarWavelet1D hwav(eta);
#else
                HaarWavelet2D hwav(eta);
#endif
                te_pois.basis().expand(&hwav, false, multi_degree(te_pois.basis().j0())+_BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END, row);
                //cout << "coarsening ..." << endl;
                row.compress(1e-13);
            }
            // set the etath row to coeff
            // brute force
            for (InfiniteVector<double,Index>::const_iterator it(row.begin()),itend(row.end());it!=itend;++it)
            {
                //cout << "setting entry a("<<eta<<", " << it.index().number() << ") = " << (*it) << endl;
                if (abs((*it)) > 1e-13)
                {
                    A.set_entry(eta,it.index().number(),(*it));
                }
            }
        }
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;

        cout << "store pruned versions of the matrix. maybe only the biggest version is used in the end, but who knows?" << endl;
        tstart = clock();
        for (int haar_level=_BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_END;haar_level >= _BASIS_TRANSITION_MATRIX_HAARMAXLEVEL_START; --haar_level)
        {
            for (int spatial_offset=_BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_END;spatial_offset >= _BASIS_TRANSITION_MATRIX_SPATIAL_OFFSET_START; --spatial_offset)
            {
                char matrix_filename[250];
                //transition_haar_jmax_5_primbs_d_dt_3_3_bc_tf_jmax_3
                int haar_jmax = haar_level;
                int spatialjmax = spatial_offset+multi_degree(te_pois.basis().j0());
                sprintf(matrix_filename, "%s/%s/transition_haar_jmax_%d_%s_d_dt_%d_%d_bc_%s_jmax_%d",transitionMatrixStorageFolder,precondition_path,haar_jmax,basis_string,d,dT,bc_string,spatialjmax);

                cout << "write to file: " << matrix_filename << endl;
#if _DIMENSION == 1
                int number_of_last_row = (1<<(haar_level+1))-1; // row and column numbering starts at 0
#else
                int number_of_last_row = (1<<((haar_level<<1) + 2 )) -1;
#endif
                int levelnorm = multi_degree(te_pois.basis().j0());
                int number_of_last_column = te_pois.basis().last_wavelet(levelnorm+spatial_offset).number();
                A.matlab_output(matrix_filename, "A", 1,number_of_last_row,number_of_last_column);
            }
        }
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;

    }

    // Laplacian blocks
    if (0 <= _LAPLACIAN_MATRIX_BLOCK_MAXOFFSET)
    {
#if _DIMENSION == 1
        SparseMatrix<double> A;
        unsigned int jmax = multi_degree(te_pois.basis().j0())+_LAPLACIAN_MATRIX_BLOCK_MAXOFFSET;
        te_pois.basis_.set_jmax(jmax);
        char matrix_filename[250];
        //laplacian_primbs_d_dt_3_3_bc_tt_jmax_3_npcd
        tstart = clock();
        cout << "compute blocks of the laplacian stiffness matrix" << endl;
        sprintf(matrix_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",laplacianStorageFolder,precondition_path,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
        cout << "loading file: " << matrix_filename << endl;
        A.matlab_input(matrix_filename);
//CLEANUP
        //cout << "A = " << endl;
        //A.print(cout,6,4);

        for (unsigned int i=0; i <= _LAPLACIAN_MATRIX_BLOCK_MAXOFFSET; ++i)
        {
            for (unsigned int j=0; j <= _LAPLACIAN_MATRIX_BLOCK_MAXOFFSET; ++j)
            {
                // store <psi'_lambda,psi'_mu> where |lambda| == j0+i, |mu| == j0+j
                char block_filename[250];
                //laplacian_primbs_d_dt_3_3_bc_tt_npcd_block_2_10
                sprintf(block_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,problem_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                cout << "writing block (" << i << ", " << j << ") to file: " << block_filename << endl;
                int rowstart = (i==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                int rowend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                int columnstart = (j==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                int columnend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
            }
        }
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
        cout << "branch not implemented" << endl;
#endif
    }

    // gramian blocks
    if (0 <= _GRAMIAN_MATRIX_BLOCK_MAXOFFSET)
    {
#if _DIMENSION == 1
        SparseMatrix<double> A;
        unsigned int jmax = multi_degree(te_gram.basis().j0())+_GRAMIAN_MATRIX_BLOCK_MAXOFFSET;
        te_gram.basis_.set_jmax(jmax);
        char matrix_filename[250];
        tstart = clock();
        cout << "compute blocks of the gramian stiffness matrix" << endl;
        //gramian_primbs_d_dt_3_3_bc_tt_jmax_3_npcd
        sprintf(matrix_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",gramianStorageFolder,precondition_path,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);
        cout << "loading file: " << matrix_filename << endl;
        A.matlab_input(matrix_filename);
        for (unsigned int i=0; i <= _GRAMIAN_MATRIX_BLOCK_MAXOFFSET; ++i)
        {
            for (unsigned int j=0; j <= _GRAMIAN_MATRIX_BLOCK_MAXOFFSET; ++j)
            {
                // store <psi_lambda,psi_mu> where |lambda| == j0+i, |mu| == j0+j
                char block_filename[250];
                //gramian_primbs_d_dt_3_3_bc_tt_npcd_block_2_10
                sprintf(block_filename, "%s/%s/%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,identity_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                cout << "writing block (" << i << ", " << j << ") to file: " << block_filename << endl;
                int rowstart = (i==0)?0:te_gram.basis().first_wavelet(multi_degree(te_gram.basis().j0())+i).number();
                int rowend = te_gram.basis().last_wavelet(multi_degree(te_gram.basis().j0())+i).number();
                int columnstart = (j==0)?0:te_gram.basis().first_wavelet(multi_degree(te_gram.basis().j0())+j).number();
                int columnend = te_gram.basis().last_wavelet(multi_degree(te_gram.basis().j0())+j).number();
                A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
            }
        }
        tend = clock();
        time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        cout << "  ... done, time needed: " << time << " seconds" << endl;
#else
        cout << "branch not implemented" << endl;
#endif
    }

    // blocks of (<g_eta psi_lambda , psi_mu>)_(lambda,mu) where g_eta is the Haar generator number eta
    if ((_HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_END >= _HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_START)
        && (_HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET >= 0) )
    {
        cout << "Compute blocks of the matrices <g_eta psi_lambda, psi_mu>" << endl;
#if _DIMENSION == 1
        for (int haar_level=_HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_START;haar_level <= _HAAR_GEN_MATRIX_BLOCK_HAARLEVEL_END; ++haar_level)
        {
            // Generators
            // caution: first generator on level 0 gets number 1 (level 1 consists of number 2 and 3, ...)
            // but observe that in the 1d basis the first element ist the first generator and that this time it gets the number 0 (because the first wavelet gets number 1)
            for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
            {

                SparseMatrix<double> A;
                unsigned int jmax = multi_degree(te_pois.basis().j0())+_HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET;
                te_pois.basis_.set_jmax(jmax);
                tstart = clock();
                cout << "setting up the whole matrix" << endl;

                // loading matrix from file
                char matrix_filename[250];
                //haar_gen_2_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                sprintf(matrix_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_gen_gramianStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                cout << "loading file: " << matrix_filename << endl;
                A.matlab_input(matrix_filename);

                /*
                unsigned int jmax = multi_degree(ctproblem.basis().j0())+_HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET;
                SparseMatrix<double> A;
                // setup and fill index set
                set<Index> Lambda_IndexSet;
                for (Index lambda = ctproblem.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == ctproblem.basis().last_wavelet(jmax)) break;
                }
                cout << "setting up the whole matrix" << endl;
                tstart = clock();
                HaarGeneratorBVP haarbvp(0,eta);
                cout << "HaarGeneratorBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
                TensorEquation<Basis1d,dim,Basis> haarqe(&haarbvp, bc, false);
                CachedTProblem<TensorEquation<Basis1d,dim,Basis> > cthaargramian(&haarqe);
                setup_stiffness_matrix(cthaargramian, Lambda_IndexSet, A, _PRECONDITIONED);   // in galerkin_utils (this will take some time: jmax 17 takes 1h)
                */
                //cout << "write to file: " << matrix_filename << endl;
                //A.matlab_output(matrix_filename, "A", 1,(1<<(haar_level+1))-1,eq.basis().last_wavelet(multi_degree(eq.basis().j0())+spatil_offset).number());
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
                cout << "storing matrix blocks in separate files" << endl;
                for (unsigned int i = 0; i<=_HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET; ++i)
                {
                    for (unsigned int j = 0; j<=_HAAR_GEN_MATRIX_BLOCK_SPATIAL_MAXOFFSET; ++j)
                    {
                        char block_filename[250];
                        //haar_gen_2_gramian_primbs_d_dt_3_3_bc_tt_npcd_block_2_7
                        sprintf(block_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                        cout << "writing " << block_filename << endl;
                        int rowstart = (i==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int rowend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int columnstart = (j==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        int columnend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
                    }
                }
            }
        }
#else
        cout << "branch not implemented" << endl;
#endif
    }



    // blocks of (<h_eta psi_lambda , psi_mu>)_(lambda,mu) where h_eta is the Haar wavelet number eta
    if ((_HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_END >= _HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_START)
        && (_HAAR_WAV_MATRIX_BLOCK_SPATIAL_MAXOFFSET >= 0) )
    {
        cout << "Compute blocks of the matrices <h_eta psi_lambda, psi_mu>" << endl;
#if _DIMENSION == 1
        for (int haar_level=_HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_START;haar_level <= _HAAR_WAV_MATRIX_BLOCK_HAARLEVEL_END; ++haar_level)
        {
            // Wavelets

            // caution: first generator on level 0 gets number 1 (level 1 consists of number 2 and 3, ...)
            // but observe that in the 1d basis the first element ist the first generator and that this time it gets the number 0 (because the first wavelet gets number 1)
            for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
            {
                SparseMatrix<double> A;
                unsigned int jmax = multi_degree(te_pois.basis().j0())+_HAAR_WAV_MATRIX_BLOCK_SPATIAL_MAXOFFSET;
                te_pois.basis_.set_jmax(jmax);
                tstart = clock();
                cout << "setting up the whole matrix" << endl;

                // loading matrix from file
                char matrix_filename[250];
                //haar_wav_2_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                sprintf(matrix_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_wav_gramianStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                cout << "loading file: " << matrix_filename << endl;
                A.matlab_input(matrix_filename);
                //

                /*
                // compute matrix
                // setup and fill index set
                set<Index> Lambda_IndexSet;
                for (Index lambda = ctproblem.basis().first_generator();; ++lambda)
                {
                    Lambda_IndexSet.insert(lambda);
                    if (lambda == ctproblem.basis().last_wavelet(jmax)) break;
                }
                                
                HaarWaveletBVP haarbvp(0,eta);
                cout << "HaarWaveletBVP:: eta = " << eta << " j = " << haarbvp.j << " k = " << haarbvp.k << endl;
                TensorEquation<Basis1d,dim,Basis> haarqe(&haarbvp, bc, false);
                CachedTProblem<TensorEquation<Basis1d,dim,Basis> > cthaargramian(&haarqe);
                setup_stiffness_matrix(cthaargramian, Lambda_IndexSet, A, _PRECONDITIONED);   // in galerkin_utils (this will take some time: jmax 17 takes 1h)

                */
                //cout << "write to file: " << matrix_filename << endl;
                //A.matlab_output(matrix_filename, "A", 1,(1<<(haar_level+1))-1,eq.basis().last_wavelet(multi_degree(eq.basis().j0())+spatil_offset).number());
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
                cout << "storing matrix blocks in separate files" << endl;
                for (unsigned int i = 0; i<=_HAAR_WAV_MATRIX_BLOCK_SPATIAL_MAXOFFSET; ++i)
                {
                    for (unsigned int j = 0; j<=_HAAR_WAV_MATRIX_BLOCK_SPATIAL_MAXOFFSET; ++j)
                    {
                        char block_filename[250];
                        //haar_wav_2_gramian_primbs_d_dt_3_3_bc_tt_npcd_block_2_7
                        sprintf(block_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,eta,identity_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                        cout << "writing " << block_filename << endl;
                        int rowstart = (i==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int rowend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int columnstart = (j==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        int columnend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
                    }
                }
            }
        }
#else
        cout << "branch not implemented" << endl;
#endif
    }
    
    /**************************/
    
    // blocks of (<g_eta grad psi_lambda , grad psi_mu>)_(lambda,mu) where g_eta is the Haar generator number eta
    if ((_HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_END >= _HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_START)
        && (_HAAR_GEN_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET >= 0) )
    {
#if _DIMENSION == 1
        cout << "Compute blocks of the matrices <g_eta grad psi_lambda, grad psi_mu> (one dimensional primbs wavelets)" << endl;
        for (int haar_level=_HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_START;haar_level <= _HAAR_GEN_LAPLACIAN_BLOCK_HAARLEVEL_END; ++haar_level)
        {
            // Generators
            // caution: first generator on level 0 gets number 1 (level 1 consists of number 2 and 3, ...)
            // but observe that in the 1d basis the first element ist the first generator and that this time it gets the number 0 (because the first wavelet gets number 1)
            for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
            {
                SparseMatrix<double> A;
                unsigned int jmax = multi_degree(te_pois.basis().j0())+_HAAR_GEN_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET;
                te_pois.basis_.set_jmax(jmax);
                tstart = clock();
                cout << "setting up the whole matrix" << endl;
                // loading matrix from file
                char matrix_filename[250];
                //haar_gen_2_laplacian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                sprintf(matrix_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_gen_laplacianStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                cout << "loading file: " << matrix_filename << endl;
                A.matlab_input(matrix_filename);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
                cout << "storing matrix blocks in separate files" << endl;
                for (unsigned int i = 0; i<=_HAAR_GEN_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET; ++i)
                {
                    for (unsigned int j = 0; j<=_HAAR_GEN_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET; ++j)
                    {
                        char block_filename[250];
                        //haar_gen_2_laplacian_primbs_d_dt_3_3_bc_tt_npcd_block_2_7
                        sprintf(block_filename, "%s/%s/haar_gen_%d_%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                        cout << "writing " << block_filename << endl;
                        int rowstart = (i==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int rowend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int columnstart = (j==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        int columnend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
                    }
                }
            }
        }
#else
        cout << "branch not implemented" << endl;
#endif
    }

    // blocks of (<h_eta grad psi_lambda , grad psi_mu>)_(lambda,mu) where h_eta is the Haar wavelet number eta, psi_lambda, psi_mu are one dimensional primbs wavelets
    if ((_HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_END >= _HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_START)
        && (_HAAR_WAV_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET >= 0) )
    {
        cout << "Compute blocks of the matrices <h_eta grad psi_lambda, grad psi_mu>" << endl;
#if _DIMENSION == 1
        for (int haar_level=_HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_START;haar_level <= _HAAR_WAV_LAPLACIAN_BLOCK_HAARLEVEL_END; ++haar_level)
        {
            // Wavelets
            // caution: first generator on level 0 gets number 1 (level 1 consists of number 2 and 3, ...)
            // but observe that in the 1d basis the first element ist the first generator and that this time it gets the number 0 (because the first wavelet gets number 1)
            for (int eta = (1<<haar_level); eta < (1<<(haar_level+1));++eta)
            {
                SparseMatrix<double> A;
                unsigned int jmax = multi_degree(te_pois.basis().j0())+_HAAR_WAV_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET;
                te_pois.basis_.set_jmax(jmax);
                tstart = clock();
                cout << "setting up the whole matrix" << endl;

                // loading matrix from file
                char matrix_filename[250];
                //haar_wav_2_laplacin_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                sprintf(matrix_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_jmax_%d_%s",haar_wav_laplacianStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,jmax,precondition_string);
                cout << "loading file: " << matrix_filename << endl;
                A.matlab_input(matrix_filename);
                tend = clock();
                time = (double)(tend-tstart)/CLOCKS_PER_SEC;
                cout << "  ... done, time needed: " << time << " seconds" << endl;
                cout << "storing matrix blocks in separate files" << endl;
                for (unsigned int i = 0; i<=_HAAR_WAV_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET; ++i)
                {
                    for (unsigned int j = 0; j<=_HAAR_WAV_LAPLACIAN_BLOCK_SPATIAL_MAXOFFSET; ++j)
                    {
                        char block_filename[250];
                        //haar_wav_2_laplacian_primbs_d_dt_3_3_bc_tt_npcd_block_2_7
                        sprintf(block_filename, "%s/%s/haar_wav_%d_%s_%s_d_dt_%d_%d_bc_%s_%s_block_%d_%d",matrixBlocksStorageFolder,precondition_path,eta,problem_string,basis_string,d,dT,bc_string,precondition_string,i,j);
                        cout << "writing " << block_filename << endl;
                        int rowstart = (i==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int rowend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+i).number();
                        int columnstart = (j==0)?0:te_pois.basis().first_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        int columnend = te_pois.basis().last_wavelet(multi_degree(te_pois.basis().j0())+j).number();
                        A.matlab_output(block_filename,"A",1,rowstart, rowend,columnstart,columnend);
                    }
                }
            }
        }
#else
        cout << "branch not implemented" << endl;
#endif
    }
    
    /*******************************/
    return 0;
}
#endif

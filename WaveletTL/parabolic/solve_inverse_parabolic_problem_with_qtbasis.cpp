// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Ulrich Friedrich                                                   |
// +--------------------------------------------------------------------+

/*
 * Main file for the treatment of the inverse problem in my project
 * - uses the qtbasis
  */

// code related flags:
// define verbosity and tbasis:
// more output for the cached problemm (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0
#define _COMPUTE_UPDATE_W_VERBOSITY 0 // output in parabolic_tools::compute_update_w
#define _NO_JPMAX_grt_70_WARNING 1 // turn off warning in APPLY_TENSOR
#define _SKIP_MINIMIZATION_PROBLEM 1 // replace the costly minimization in APPLY by always taking jtilde = jmax
// switch between isotropic and anisotropic Wavelets (in cdd1.h)
#define _WAVELETTL_USE_TBASIS 1
// switch between some code. Use 0 every time!
#define _EXPANSIONTYPE_F_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_f; 1 liefert Unterschied falls f konstant ist und man TIME_CONSTANT =0,1 wählt
#define _EXPANSIONTYPE_FT_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_ft


/* implemented Decompositions:
 * 0 :: 1 patch, 2d or 3d
 * 1,2,7,8 :: 2 patches, homogeneous BC
 * 1: extension to the right
 * 2: extension to the south
 * 7: extension to the north
 * 8: extension to the left
 * 3,4 :: L-shaped domain, homogeneous BC
 * (!) 5:: L-shaped. corner domain is extended in 2 directions. inhomogeneous BC at reentrant corner. L-shaped supports are not implemented in compute_sum_over_patches !!
 * 6: big cube decomposed into 4 subcubes. Homogeneous BC
 * 9: Donut; 8 domains
 * 10: 3 domains in a line. middle is extended left and right
 * 11: Slit domain (0,2)^2 \ {1}x(1,2)  (4 sub-domains)
 * 12: Big square in 3d (8 sub-squares)
 * 13:  (0,2)^3 \ (1,2)^3 = Fichera corner domain (7 sub-squares)
 */
#define _DECOMPOSITION 0
#define _NOP 1
#define _MAXOFFSET 0
//#define _RHS_NUMBER 5
#define _DIMENSION 2
#define _PRIMAL_ORDER_D 3
#define _DUAL_ORDER_DT 3
// #define _COMP_COMPUTE_RHSUTE_RHS 1 
//#define _ONEDIMHAARCOUNT 1


// define mathematical setting
//#define _MAX_INVERSE_ITERATIONS 30000 // now given as an argument
//#define _SPATIAL_JMAX 3 // maximal 13 für 1D, 9 im 2D Fall; Gramian Matrizen nicht weiter vorberechnet
#define _HAAR_JMAX 4

//#define _D 2// primal polynomial degree; only used for the initialization of the pbasis used to define Index
//#define _DT 2// dual polynomial degree; only used for the initialization of the pbasis used to define Index

#define _FORWARD_PROBLEM_NO 1 // 13,14 problem number for uexact (time forward case). Determines initial value.


// we use NOTS many time steps for the parabolic equation (time discretization has this number+1 many points) (equidistant)
// must fit with the definition of wtrue
#define _NUMBER_OF_TIME_STEPS 8

// Initialization of the coupling Matrix W
// The (unknown) true coupling Matrix W is specified as a parameter from the command line. Its number determines many output filenames as well as utrue (if it is computed at run time)
// If the computation is not resumed, i.e., starts at iteration 0, an initial guess for the matrix W has to be specified. This is usual 0, but does not have to
#define _SECOND_SETUP_COUPLING_MATRIX_W 10 // W=empty 10 = 1 patch, 30 = 3 patches;

// TODO If RESUME is specified: do not initialize W. Only use the number for the right filename
/* if wanted, one could initialize W used in the iteration with something different than 0.
 * For instance one could resume a previous iteration (there is a little bit of
 * code to be implemented for this, but the big parts are already there)
 * If one wanted to do this, the old logfile had to be resumed. Also the following two makros have to be set accordingly
 */
//#define _INVERSE_ITERATION_COUNT 0 // you need to change this number if you want to resume an iteration!

// enable/disable some debugging routines
#define _TEST_PARABOLIC_TOOLS 0


        
/*
 * u0 = initial value. It determines together with W u_true and v_true, the solutions of the time forward and backward parabolic problem.
 * (With 1/2 and u_true as the righthand sides, respectively)
 * 
 * If _COMPUTE_u0 == 1 is specified u0 is computed (and stored on disk). 
 * u0 is determined by u0_function given by _FORWARD_PROBLEM_NO)
 */
#define _COMPUTE_u0 1

/*
 * If _COMPUTE_SIMULATED_DATA == 1 is specified (and _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 is NOT)
 * u_true, v_true are computed as solutions of the parabolic problem. Afterwards they are stored in files.
 * If the flag is not set (== 0) u_true, v_true are loaded from files.
 */
#define _COMPUTE_SIMULATED_DATA 1 // compute u_true, v_true (true as in "corresponding to the looked for W")

/*
 * If _COMPUTE_SIMULATED_DATA == 1 is specified and _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 == 1 is specified as well,
 * then we assume that u_exact_function given by _FORWARD_PROBLEM_NO describes the solution of the parabolic problem at all time steps.
 * Coefficients are computed directly using the function
 */
#define _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 0 // 1 means uexact describes the solution for all time points. 0 = use uexact only for u0. this determines whether uexact can be used to store u_problem

/*
 * compute the coefficients of f = 1/2 and store them to disk
 * OR
 * load them from disk
 */
#define _COMPUTE_ONEHALF_COEFFS 1 // compute u0, u_true, v_true (true == corresponding to the looked for W

// deactivate the inverse solver for debugging (everything is initialized, but no iteration computed)
#define _SOLVE_INVERSE_PROBLEM 0

/* activate the following switches to
 * - fill the logfile with interesting numbers
 * - print the l2 error of the current forward/backward iterates ||u_i-utrue|| and ||v_i-vtrue|| to cout
 * utrue,vtrue are computed if _COMPUTE_SIMULATED_DATA is set.
 * If they are _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 determines whether they are
 * given by _FORWARD_PROBLEM_NO, _BACKWARD_PROBLEM_NO or if the
 * solution coefficients computed by the the parabolic solver are used.
 */

/* _COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL = 1 :: include errors and norms ue2,u0,u2 in logfile*/
#define _COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL 1


/* if _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL is set, the program will compare
 * the current parameter estimation with either w_true, or with a given
 * approximation of the limit of the iteration (w_limit).
 * => include we1,we2,w0,w1,w2 in logfile
 */
#define _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL 1 // 1 == l2 difference of current coefficients for W and w_true

#define _SAVE_SPATIAL_SOLUTIONS_TO_FILE 1 // save parabolic.forward_solution_ & .backward_solution_ to disk (after the last iteration)
#define _SAVE_UPDATE_U_H_TO_FILE 1 // save output of w_update ( u,v and unpostprocessed w_update=u*h) to plotfile. This is only done if ((max_inverse_terations > 0) && (iteration_count == max_inverse_iterations))
// W is stored on HDD (.iv). This allows to resume the computation

// TODO: add a nice description of the usage of tolerances:
#define _TOLERANCE 1e-5
#define _INCREMENT_TOLERANCE 1e-6 // tolerance for the call of increment (parabolic time step)
/* If the program is called with max_inverse_iterations=0 the relative error of
 * the l1-Norm of the update (of W) is used to determine whether the iteration
 * should stop. The tolerance for this step is given here:*/
#define _INVERSE_ITERATION_STOPPING_TOLERANCE 1e-5 // tolerance for inverse iteration stopping criterion

#define _DIFFUSION_COEFFICIENT 1

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
#include <interval/p_evaluate.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <geometry/sampled_mapping.h>

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tbasis_indexplot.h>

#include <general_domain/qtbasis.h>
#include <general_domain/qtbasis_index.h>

#include <parabolic/example_parabolic_problems.cpp>
#include <parabolic/aff_lin_par_eq.h>

#include <parabolic/lin_par_eq_tensor.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>

#include <parabolic/parabolic_tools.h>
#include <io/infinite_vector_io.h>
#include <io/fixed_matrix_io.h>

#include <algebra/fixed_vector.h>
#include <algebra/fixed_matrix.h>


using namespace std;
using namespace MathTL;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

//int main()
int main(int argc, char *argv[])
{
    cout << "This program realizes the iterated soft shrinkage for the linearized PDE in " << _DIMENSION << " space dimensions" << endl;
#if 0
    // ./solve_PDE 30000 3 5 0 10000
    // ./solve_PDE 500 3 5 0 10000
    if (argc != 9)
    {
        cout << "Command line arguments are: " << endl
             << "[1]: number_of_inverse_iteration" << endl << " = inverse iteration stops after this amount of iterations. Specifying 0 activates a stop condition based upon the relative error. In this case the iteration stops if the ell2 norms of the update of the coupling matrix W relative to the norm of W is less or equal tolerance = " << _INVERSE_ITERATION_STOPPING_TOLERANCE << "." << endl
             << "[2]: spatial wavelet maxlevel spatial_jmax" << endl << " = maximum level of wavelet considered for the forward solver. Needs to be greater than " << (_PRIMAL_ORDER_D*_DIMENSION) << "." << endl
             << "[3]: initialization number of W" << endl << " = the initialization choice that is or was used to compute ydata." << endl
             << "[4]: number of best available approximation to W" << endl << " = output of the program includes the norm difference of the current approximation of W and either W_true (if you specify 0) or W_N (value after N iterations, if you specify N). The idea is to run the program twice and specify N=max_iterations for the second run, to get convergence rates W_k->W_limit(alpha)." << endl
             << "[5]: shrinkage parameter alpha as integer factor times 1e-7" << endl
             << "[6]: compute noise (1) or load from file (0)" << endl << " = if you want the exact same noise and test for different alpha specify 0." << endl
             << "[7]: size of noise delta as integer factor times 1e-5" << endl << " = specify the magnitude of the noise as N*1e-5. Do specify a value even if compute_noise = 0, because otherwise the program cannot load the right noise file." << endl
             << "[8]: number of first iteration. 0 for new computation"
             << endl << "Hardcoded parameters:" << endl
                << "Constant diffusion coeficient of magnitude " << _DIFFUSION_COEFFICIENT << "." << endl
             << "Number of spatial dimensions considered is " << _DIMENSION << "." << endl
             << "Domain decomposition number is " << _DECOMPOSITION << "." << endl;
        cout << "Wavelet basis considered: Primbs" << endl;
        cout << "Primal and dual order of the wavelets _PRIMAL_ORDER_D = " << _PRIMAL_ORDER_D << ", _DUAL_ORDER_DT = " << _DUAL_ORDER_DT << "." << endl;
        cout << "Haar wavelet maximal level used for discretization of W is " << _HAAR_JMAX << "." << endl
             << "Relative error stopping tolerance for the ell1 stopping criterion (active only if max_inverse_iterations==0) is " << _INVERSE_ITERATION_STOPPING_TOLERANCE << "." << endl;
        if (_PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 == 0)
        {
            if (_COMPUTE_SIMULATED_DATA == 1)
                cout << "ydata is computed with W specified above and u0 given by _PROBLEM_NUMBER " << _FORWARD_PROBLEM_NO << " (see test_problems.cpp)." << endl;
            else
                cout << "ydata is loaded from disk. The files are specified by W given above and u0 given by _PROBLEM_NUMBER " << _FORWARD_PROBLEM_NO << " (see test_problems.cpp)." << endl;
        }
        else
        {
            cout << "ydata is determined by uexact (time forward case) given by _FORWARD_PROBLEM_NUMBER = " << _FORWARD_PROBLEM_NO << ","
                 << "see test_problems.cpp." << endl;
            if (_COMPUTE_SIMULATED_DATA == 0)
                cout << "ydata is loaded from disk." << endl;
        }
        if (_SECOND_SETUP_COUPLING_MATRIX_W == 0)
            cout << "W for the inverse iteration is initialized with 0." << endl;
        else
            cout << "W for the inverse iteration is initialized with setting " << _SECOND_SETUP_COUPLING_MATRIX_W << "." << endl;        
        cout << "Number of timesteps used for discretization is " << _NUMBER_OF_TIME_STEPS << "." << endl;
        cout << endl << "sample call:" << endl << "./solve_PDE 500 6 6 0 1000 1 100 0" << endl;
        return 0;
    }
    const int max_inverse_iterations = atoi(argv[1]) + atoi(argv[8]);
    const int spatial_jmax = atoi(argv[2]);
    const int first_setup_coupling_matrix_w = atoi(argv[3]);
    const int w_limit_number = atoi(argv[4]);
    const int alpha_times = atoi(argv[5]);
    const int compute_noise = atoi(argv[6]);
    const int delta_times = atoi(argv[7]);
    const int start_iteration = atoi(argv[8]);
    
#else
    //                                     it jmax winit wbestapprox alpha
    // ./solve_PDE_dim1_ddT2_tsteps10_hmax5 20 8 5 0 100
    // ./solve_PDE_dim1_ddT3_tsteps10_hmax5 20 8 5 0 100
    // ./solve_PDE_dim1_ddT2_tsteps10_hmax5 5 10 5 0 100 zu teuer
    // ./solve_PDE_dim1_ddT2_tsteps10_hmax6 10 8 5 0 100
    const int max_inverse_iterations = 1; // needs to be greater than start_iteration !
    const int spatial_jmax = 3*_DIMENSION;
    const int first_setup_coupling_matrix_w = 11;
    const int w_limit_number = 0;
    const int alpha_times = 100;
    const int compute_noise = 1;
    const int delta_times = 500000;
    const int start_iteration = 0;
#endif
    const double alpha = alpha_times * 1e-7;
    const double delta = delta_times * 1e-5;
    cout << "called with" << endl
         << "max_inverse_iterations = " << max_inverse_iterations << endl
         << "spatial wavelet maxlevel jmax = " << spatial_jmax << endl
         << "initialization number of W = " << first_setup_coupling_matrix_w << endl
         << "iteration number corresponding to reference reconstruction Wlimit = " << w_limit_number << endl
         << "shrinkage_parameter alpha = " << alpha << endl
         << "compute noise (1) or load from file (0) = " << compute_noise << endl
         << "noise level delta = " << delta << endl
         << "start_iteration = " << start_iteration << endl;

    //abort();
    clock_t tstart, tend;
    double time;

    const double diff_coeff(_DIFFUSION_COEFFICIENT);
    const int dim = _DIMENSION;
    const int d  = _PRIMAL_ORDER_D;
    const int dT = _DUAL_ORDER_DT;
    const int haar_jmax = _HAAR_JMAX;
    const int precise_evaluate_granularity = (1<< (_HAAR_JMAX+2));
    const int onedimhaarcount = (1<<_HAAR_JMAX);
    const int decomposition = _DECOMPOSITION;
    unsigned int num_of_patches;
    //const int number_of_timesteps (1<<_TIME_DISCRETIZATION_GRANULARITY);
    const unsigned int number_of_timesteps = _NUMBER_OF_TIME_STEPS;
    

// TODO: add a nice description of the usage of tolerances:
    const double tolerance_APPLY = _TOLERANCE;
    const double increment_tolerance = _INCREMENT_TOLERANCE; // tolerance for the call of increment (parabolic time step)
    //const double alpha = _SHRINKAGE_PARAMETER; // this line could be moved somewhere else, parameter could be changed at each iteration ...
    const double stopping_tolerance = _INVERSE_ITERATION_STOPPING_TOLERANCE;


    typedef PBasis<d,dT> Basis1d;

    //typedef TensorBasis<Basis1d,dim> TBasis;
    typedef QTBasis<Basis1d,dim> QTBasis;
    typedef QTBasis::Index Index;
    const unsigned int forward_problem_number=_FORWARD_PROBLEM_NO;
    //const unsigned int backward_problem_number=_BACKWARD_PROBLEM_NO;

    const int resolution = d+dim+((dim == 1)? (spatial_jmax-3):(spatial_jmax-6)) +1; // just some number that guarantees that the sampling grid of the matlab output (in sampledmapping) is fine enough
    
    std::ofstream resultstream; // ostream used for matlab output by SampledMapping
    std::ofstream logstream, plotstream;


    // construct the basis
    Array1D<Point<dim,int> > corners;
    Array1D<FixedArray1D<int,2*dim> > neighbours;
    Array1D<FixedArray1D<bool,2*dim> > bc_bool;
    Array1D<FixedArray1D<int,2*dim> > bc_int;


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
            assert (dim == 2);
            break;
        case 4: // L-shaped domain, homogeneous BC, patch 0 is extended to the south to patch 1. Patch 2 is extended to the left to patch 0
            num_of_patches = 3;
            assert (dim == 2);
            break;
        case 5: // L-shaped domain, homogeneous BC everywhere except around the reentrant corner
            // patch 0 is extended to right to patch 1 and to the north to patch 2
            num_of_patches = 3;
            assert (dim == 2);
            break;
        case 6: // Big square, decomposed into 4 sub-squares
            // patch 0 is extended to right to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            num_of_patches = 4;
            assert (dim == 2);
            break;
        case 9: // Donut! Decomposed into 8 sub-squares
            // 6 5 4
            // 7 x 3
            // 0 1 2
            // 
            // patches 1,3,5,7 are extended to their neighbours
            // homogeneous BC wherever possible
            num_of_patches = 8;
            assert (dim == 2);
            break;
        case 10: // Line
            // 0 1 2
            // 
            // patches 1 is extended left and right
            // homogeneous BC wherever possible
            num_of_patches = 3;
            assert (dim == 2);
            break;
        case 11: // Slit domain (0,2)^2 \ {1}x(1,2)  (4 sub-domains)
            // patch 0 is extended to right to patch 1
            // patch 1 is extended nowhere
            // patch 2 is extended to the south to patch 1. But is is not extended to the patch left of it (number 3)
            // patch 3 is extended to the south to patch 0. But is is not extended to the patch right of it (number 2)
            num_of_patches = 4;
            assert (dim == 2);
            break;
        case 12: // Big square in 3d (8 sub-squares)
            // patches 0,1,2,3 and 4,5,6,7 are extended like in decomposition 6 to form bases for (0,2)^2x (0,1) and (0,2)^2x(1,2), respectively.
            // patches 0,1,2,3 are extended up in z direction.
            // E.g. :
            // patch 0 is extended to the east to patch 1 and to the north to patch 3
            // patch 1 is extended to the north to patch 2
            // patch 2 is extended to the left to patch 3
            num_of_patches = 8;
            assert (dim == 3);
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
    corners.resize(num_of_patches);
    neighbours.resize(num_of_patches);
    bc_bool.resize(num_of_patches);
    bc_int.resize(num_of_patches);

    for (unsigned int p = 0; p < num_of_patches; ++p)
    {
        for (unsigned int i = 0; i < dim; ++i)
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
//            bc_bool[0][0] = false;
            bc_bool[0][1] = false;
//            bc_bool[0][2] = false;
            bc_bool[0][3] = false;
//            bc_int[0][1] = 0;
//            bc_int[0][3] = 0;
            break;
        case 1: // 2 domains, homogeneous BC, left patch is extended to the right
            assert (num_of_patches == 2);
            if (dim == 1)
            {
                corners[1][0] = 1;
            }
            else if (dim == 2)
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
            if (dim == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (dim == 2)
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
            if (dim == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (dim == 2)
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
            if (dim == 1)
            {
                corners[0][0] = 19;
                corners[1][0] = 20;
            }
            else if (dim == 2)
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 2);
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
            assert (dim == 3);
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
            assert (dim == 3);
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
    
    QTBasis qtbasis(corners, neighbours,bc_bool);
    qtbasis.set_jmax(spatial_jmax);
    
#if 0
//CLEANUP
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
#elif _BOUNDARY_CONDITIONS == 1
    bc[0] = bc[1] = bc[2] = bc[3] = false;
    const char* bc_string = {"ffff"};
#endif
#endif

#endif
    
    
    
    
    
    
    
    
    // Set up the PDE and discretization

    

    

    // set up the filenames

#if _TEST_PARABOLIC_TOOLS == 1
    //const char* temp_filename = "/home/friedrich/source/precomputed/delete_on_sight";
    const char* temp_filename = "/import/shared/friedrich/source/precomputed/delete_on_sight";
#endif
    // file containing the coefficients of the constant driving term ( R*1/2 with R=1)
    //const char* path = "/home/friedrich/source/precomputed";
    const char* path = "/import/shared/friedrich/source/precomputed/"; 
    //const char* path = "/import/shared/friedrich/source/X"; // for testing the parallel shell skript
    
    // auf Marc: (computescript3.sh)
    //char path[100];
    //sprintf(path, "/scratch/friedrich/d%d_dt%d_w%d_a%d_d%d",d,dT,first_setup_coupling_matrix_w,alpha_times,delta_times);
    
    
    
    
    char onehalf_filename[250]; // filename of IV that contains the coeffs of onehalf_function
    char onehalf_precond_gramian_filename[250]; // Vector<double> of onehalf, divided by sqrt of the diagonal of gramian 
    FixedArray1D< FixedArray1D<char*, number_of_timesteps+1>,2> onehalf_mode_filename; // divided by sqrt of T and alphaI-T
    
    sprintf(onehalf_filename, "%sfunctions/onehalf_primbs_d_dt_%d_%d_dec%d.iv",path,d,dT,decomposition);
    sprintf(onehalf_precond_gramian_filename, "%sfunctions/onehalf_primbs_d_dt_%d_%d_dec%d_gramian.iv",path,d,dT,decomposition);
    
    char temp_filename[250];
    for (int mode = 0; mode < 2 ; ++mode)
    {
        for (int t=0; t<=number_of_timesteps; ++t)
        {
            sprintf(temp_filename, "%sfunctions/onehalf_primbs_d_dt_%d_%d_dec%d_mode_%d_t_%d_%d.iv",path,d,dT,decomposition,mode,number_of_timesteps,t);
            onehalf_mode_filename[mode][t] = temp_filename;
        }
    }

    // the matrices \int h_eta psi_lambda psi_mu
//    char haar_gramian_filename_start[250];
//    sprintf(haar_gramian_filename_start, "%sstiffness_matrices/unpreconditioned/haar_wav_",path);
//    char haar_gramian_filename_end[250];
//    sprintf(haar_gramian_filename_end, "_gramian_primbs_d_dt_%d_%d_bc_%s_jmax_%d_npcd",d,dT,bc_string,spatial_jmax);
//    char laplacian_filename[250];
//    sprintf(laplacian_filename, "%sstiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_%d_%d_bc_%s_jmax_%d_npcd",path,d,dT,bc_string,spatial_jmax);
    char u0_filename[250];
    // initial value file. Contains an InfiniteVector of wavelet coefficients w.r.t. the primbs basis
    // u0_problem_13_wnum_0_primbs_d_dt_3_3_bc_tt_jmax_7.iv
    sprintf(u0_filename, "%sinitial_values/u0_problem_%d_wnum_%d_primbs_d_dt_%d_%d_dec%d_jmax_%d.iv",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition,spatial_jmax);
    // files that hold solutions of the parabolic forward and backward problem corresponding to the looked for parameters. Matlab format.
    // u_problem_13_wnum_0_primbs_d_dt_3_3_bc_tt_jmax_7_tstep_10_4.m
    char true_forward_solution_filename_start[250];
    sprintf(true_forward_solution_filename_start,"%sydata/u_p%d_w%d_p%d%d_dec_%d_j%d_t%d_",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition,spatial_jmax, number_of_timesteps);
    const char * true_forward_solution_filename_matlab_end = ".m";
    const char * true_forward_solution_filename_coeffs_end = ".iv";
    //char true_backward_solution_filename_start[250];
    // we do not use _BACKWARD_PROBLEM_NO since the output is determined by ydata and W (_FORWARD_PROBLEM_NO and wnum)
    //sprintf(true_backward_solution_filename_start,"%sydata/v_p%d_w%d_p%d%d_bc%s_j%d_t%d_",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,bc_string,spatial_jmax, number_of_timesteps);
    //const char * true_backward_solution_filename_matlab_end = ".m";
    //const char * true_backward_solution_filename_coeffs_end = ".iv";

    // files that hold current solutions of the parabolic forward and backward problem. Matlab format.
    // u_p13_w0_p33_bctt_j7_a6800_d0_t10_4.m
    char forward_solution_filename_start[250];
    sprintf(forward_solution_filename_start,"%soutput/u_p%d_w%d_p%d%d_dec%d_j%d_a%d_d%d_t%d_",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition,spatial_jmax, alpha_times, delta_times, number_of_timesteps);
    const char * forward_solution_filename_end = ".m";
    char backward_solution_filename_start[250];
    sprintf(backward_solution_filename_start,"%soutput/v_p%d_w%d_p%d%d_dec%d_j%d_a%d_d%d_t%d_",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition,spatial_jmax, alpha_times, delta_times, number_of_timesteps);
    const char * backward_solution_filename_end = ".m";
    // file that holds a solution of the inverse problem (coeffs of W wrt Haar basis)
    // the Haar basis does have free boundary conditions, but the properties of the spatial basis do influence the output!
    // depending on the iteration number the loaded file will contain
    // wlimit = the reconstruction to be considered "final" that is used to compute reconstruction errors
    // woutput = where the computed reconstruction is stored (important for resuming or for computation of a better wlimit)
    // wresume = the construction used for initializing the inverse iteration (0 in the first iteration)
    //char wlimit_coefficients_filename_start[250];
    //sprintf(wlimit_coefficients_filename_start,"%sydata/w_p%d_w%d_p%d%d_bc%s_j%d_h%d_a%d_d%d_i",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,bc_string, spatial_jmax, haar_jmax, alpha_times, delta_times);
    //const char * wlimit_coefficients_filename_matlab_end = ".m";
    //const char * wlimit_coefficients_filename_end = ".iv"; // this InfiniteVector contains coefficients indexed with integers
    char w_coefficients_filename_start[250];
    sprintf(w_coefficients_filename_start,"%soutput/w_p%d_w%d_p%d%d_dec%d_j%d_h%d_a%d_d%d_i",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition, spatial_jmax, haar_jmax, alpha_times, delta_times);
    //const char * w_coefficients_filename_matlab_end = ".m";
    const char * w_coefficients_filename_end = ".iv"; // this InfiniteVector contains coefficients indexed with integers
    
    // old: logfile_problem_13_wnum_0_primbs_d_dt_3_3_bc_tt_jmax_9_haar_jmax_5_tstep_10.m
    // log_p_13_wnum_0_p_d_dt_3_3_bc_tt_jmax_9_hmax_5_alp_10000_t_10.m // old names got too long for matlab!
    // log_p13_w0_p33_bctt_j9_h5_a10000_d5000_t10.m
    char logfile[250];
    sprintf(logfile,"%soutput/log_p%d_w%d_p%d%d_dec%d_j%d_h%d_a%d_d%d_t%d.m",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition, spatial_jmax, haar_jmax, alpha_times, delta_times, number_of_timesteps);
    
    // for debugging: file to store haar_generator coefficients of forward solution and backward solution (u,h) as well as their product w_update
    // currently only one iteration is stored in this file
    char plotfile[250];
    sprintf(plotfile,"%soutput/plot_p%d_w%d_p%d%d_dec%d_j%d_h%d_a%d_d%d_t%d.m",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition, spatial_jmax, haar_jmax, alpha_times, delta_times, number_of_timesteps);

    // files that contains noise for each timestep
    // the noise does not depend on the problem at hand, but we do not want to use different noises while changing alpha
    char noise_filename_start[250];
    sprintf(noise_filename_start,"%snoise/noise_p%d_w%d_p%d%d_dec%d_j%d_d%d_t%d_",path,_FORWARD_PROBLEM_NO,first_setup_coupling_matrix_w,d,dT,decomposition,spatial_jmax, delta_times,number_of_timesteps);

    tstart = clock();
    
    cout << "main:: Initialize parabolic model." << endl;

    // Compute & store the parabolic solution of the problem determined by the
    // initial value u0 and the coupling matrix W:
    // u' = D*Delta u - Wu + 1/2

#if 0    
    // Initialize cached gramian matrix of the problem
    // needed to transform coefficient vectors (dual into primal)
    char matrix_filename[250];
    // haar_wav_0_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
    sprintf(matrix_filename, "%s%d%s",haar_gramian_filename_start, 0, haar_gramian_filename_end);
    cout << "main:: initialize CompressedProblemFromMatrix(gram)" << endl;
    IdentityBVP<dim> identity_bvp(NULL); //function reference is neverr used ...
    TensorEquation<Basis1d,dim,Basis> gram_problem (&identity_bvp, bc, false); 
    gram_problem.set_jmax(spatial_jmax,false); //false == do not call compute_rhs()
    // pbasis is initialized. pbasis has some internal maxlevel.
    // quite some memory is consumed depending on it. standard: about 21 mb memory!!!!
    // memory and time consumption is doubled each level step. jmax=14 => 43 mb
    // see: p_basis.h::L.390:

    SparseMatrix<double> gram_matrix(1);
    cout << "main:: initialize gram_matrix with file = " << matrix_filename << endl;
    gram_matrix.matlab_input(matrix_filename);
// TODO: add a nice description for the parameters
    CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> > ctgramian(&gram_problem, &gram_matrix, 5.5, 35.0, false); // false == matrix not preconditioned
#endif
        
    // set up the coupling matrix W
    //Array1D<InfiniteVector<double, int> > w_true; // for each timestep: active Haar wavelet coefficients of W
    //w_true.resize(number_of_timesteps+1);
    FixedArray1D<Array1D<FixedMatrix<double, onedimhaarcount> >, number_of_timesteps +1> w_true;
    //switch (_FIRST_SETUP_COUPLING_MATRIX_W)

    //initialize_coupling_matrix<_HAAR_JMAX>(w_true, first_setup_coupling_matrix_w, d);
    initialize_coupling_matrix(w_true, first_setup_coupling_matrix_w);
    
    // set up time discretization
    // so far we use a fixed time_discretization
    //Array1D<double> time_discretization;
    //time_discretization.resize(number_of_timesteps+1);
    FixedArray1D<double, number_of_timesteps +1 > time_discretization;
    for (unsigned int i = 0; i <= number_of_timesteps; ++i)
    {
        time_discretization[i] = i / (double)(number_of_timesteps);
    }
    
    
    // Initialize initial value u0. Simulated Data is either computed on run time or loaded from files
#if _COMPUTE_u0 || _COMPUTE_SIMULATED_DATA
#if _DIMENSION == 1
    Exact_Sol1D<forward_problem_number> uexact_function;
#else
    Exact_Sol2D<forward_problem_number> uexact_function;
#endif
#endif
    
#if _COMPUTE_ONEHALF_COEFFS
    ConstantFunction<dim,double> onehalf_function(Vector<double>(1,"0.5"));
#endif
    
    cout << "calling constructor AffLinParEq_qtbasis" << endl;
    tstart = clock();
    

CachedQTProblem<QTBasis, onedimhaarcount,1>(&qtbasis, 
                      5.5, 40.0); 
cout << Basis1d::primal_polynomial_degree() << endl;
AffLinParEq_qtbasis<number_of_timesteps, QTBasis, precise_evaluate_granularity, onedimhaarcount>(&qtbasis,
        diff_coeff, 
                             onehalf_filename,
                             5.5,
                              40
                              );
//Dummy<QTBasis, onedimhaarcount> temp_dummy(&qtbasis, onehalf_filename, 5.5,40);


    AffLinParEq_qtbasis<number_of_timesteps, QTBasis, precise_evaluate_granularity, onedimhaarcount> parabolic_problem
                              (&qtbasis,
#if _COMPUTE_u0
                              &uexact_function,
#else
                              NULL,
#endif
                              u0_filename,
                              diff_coeff,
                              w_true,
                              time_discretization,
#if _COMPUTE_ONEHALF_COEFFS
                              &onehalf_function,
#else
                              NULL, //const Function<DIM>* f, // from CQTProblem: rhs as a function
#endif
                              onehalf_filename, // filename of IV that contains the coeffs of onehalf_function
                              onehalf_precond_gramian_filename, 
                              onehalf_mode_filename,
                              5.5, // from CQTProblem: estimate for \|dI+W\|
                              40.0, // from CQTProblem:  estimate for \|(dI+W)^{-1}\|
                              tolerance_APPLY
                              );
    
//        return 0;
//}
//#if 0

    
//    InfiniteVector<double,int> u0T;


    // We have to compute the coefficients of u0 w.r.t. the primal basis
    // The following code is equivalent to the (not implemented!!) call
    // ctgramian.basis().expand(&uexact, TRUE, spatial_jmax, u0T);
// CLEANUP
    // cout << "main:: compute u0T" << endl;
    //ctgramian.basis().expand(&uexact, false, spatial_jmax, u0T);
    //cout << "main :: ctgramian.basis().expand(&uexact, false, spatial_jmax, u0T):: u0T = " << endl << u0T << endl;
    //u0T.clear();

    
    //cout << "main:: u0T after CDD1_SOLVE :: " << endl << u0T << endl;
//#else
//    // load data
//    cout << "main:: load u0T from file " << u0_filename << endl;
//    //u0T.readFromFile(u0_filename);
//    std::ifstream u0_ifstream(u0_filename, std::ifstream::binary);
//    u0_ifstream.open(u0_filename);
//    readIVFromFile(u0T.basis(), u0T, u0_ifstream);
//    u0_ifstream.close();
//    cout << "main:: done." << endl;
//#endif

    

    

//    AffLinParEq_Precomputed<TensorEquation<Basis1d,dim,QTBasis>, number_of_timesteps > parabolic
//                            (&gram_problem, // only the basis of this problem is used, no calls to f are made
//                            u0T,
//                            diff_coeff,
//                            w_true,
//                            time_discretization,
//                            f_filename,
//                            //true_forward_solution_filename_start,
//                            //true_forward_solution_filename_coeffs_end,
//                            haar_gramian_filename_start,
//                            haar_gramian_filename_end,
//                            laplacian_filename,
//                            haar_jmax,
//                            spatial_jmax,
//                            tolerance
//                            );
    
    {
// TODO: Viel zu aufwendig!
        
        // äquivalent zu parabolic_problem.backward_solution_[0].clear();
        
        parabolic_problem.set_time_direction(false);  // time still counts from 0 to 1. The change of direction is handled by aff_lin_par_eq
        InfiniteVector<double,int> temp_iv;
        temp_iv.clear();
        parabolic_problem.set_current_timestep(0); //0 means the last timestep!
        parabolic_problem.set_solution(temp_iv,tolerance_APPLY);
        parabolic_problem.set_time_direction(true);
    }

    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "main:: initialization of parabolic_problem completed. (true_forward_solution_preprocessed is not initialized jet).\nTime needed sofar = " << time << " seconds" << endl;

    cout << "main:: initialize ROW method" << endl;
    //ROWMethod<InfiniteVector<double,int> > method(WMethod<InfiniteVector<double,int> >::ROS2);
    //ROWMethod<InfiniteVector<double,int>, AbstractCachedIVP<InfiniteVector<double,int> > > method(WMethod<InfiniteVector<double,int>, AbstractCachedIVP<InfiniteVector<double,int> > >::ROS2);
    ROWMethod<InfiniteVector<double,int>, 
              AffLinParEq_qtbasis<number_of_timesteps, QTBasis, precise_evaluate_granularity, onedimhaarcount>
            > method(WMethod<InfiniteVector<double,int>, 
                             AffLinParEq_qtbasis<number_of_timesteps, QTBasis, precise_evaluate_granularity, onedimhaarcount>
                            >::ROS2);
    
    
    method.set_preprocessor(&parabolic_problem);
    cout << "  done." << endl;
    
    FixedArray1D<InfiniteVector<double,int>, number_of_timesteps+1 > utrue_coeffs; // solution of the time forward case using the unknown true coupling matrix Wtrue

#if _COMPUTE_SIMULATED_DATA == 1
    // simulate utrue
    cout << "main:: begin computation of simulated data." << endl;

    tstart = clock();
    cout << "main:: time forward case ..." << endl;
    
    
    if (_PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 == 1)
    {
        // uexact is actually meaningful, there is no need to compute any data
        InfiniteVector<double,int> temp_iv;
        for (unsigned int i=0; i<=number_of_timesteps;++i)
        {
            uexact_function.set_time(i/(double)number_of_timesteps);
            parabolic_problem.basis()->expand(&uexact_function, false, temp_iv);
            parabolic_problem.set_another_rhs_(&temp_iv);
            parabolic_problem.set_problem_mode(2);
            //ctgramian.set_f(&uexact);
            CDD1_SOLVE(parabolic_problem, tolerance_APPLY, temp_iv, spatial_jmax, tensor_simple);
            temp_iv.scale(&parabolic_problem,-1);
            temp_iv.compress(1e-14);
            cout << "main:: compute ydata_matlab_output forward case." << endl;
            cout << "       at t=" << i/(double)number_of_timesteps << " ell_2 error = " << l2_norm(parabolic_problem.forward_solution_[i] - temp_iv) << endl;
            parabolic_problem.forward_solution_[i]=temp_iv;
        }
    }
    else
    {
// TODO: comment on use of this tolerance 
//        dummy(parabolic_problem);
//        InfiniteVector<double,int> temp_iv, result, error_estimate;
//        unsigned int i=0;
//        method.increment(&parabolic_problem, (double)(i)/(double)number_of_timesteps, temp_iv, 1.0/number_of_timesteps, result, error_estimate, increment_tolerance);
        solve_parabolic_problem(parabolic_problem, method, true, increment_tolerance, tolerance_APPLY);
    }
    utrue_coeffs = parabolic_problem.forward_solution_;
    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "done. time needed = " << time << "seconds" << endl;
    
    cout << "write coefficients and matlab output of utrue HDD" << endl;
    // used for checking correctness and postprocessing, l2-error estimates of current forward solution.
    // possibly shift this code to parabolic_tools to increase readability
    store_uexact(parabolic_problem, true_forward_solution_filename_start, true_forward_solution_filename_coeffs_end, true_forward_solution_filename_matlab_end, resolution);

    /* // Same for backward problem. Currently not used!
#if _DIMENSION == 1
    Exact_Sol1D<backward_problem_number> vexact;
#else
    Exact_Sol2D<backward_problem_number> vexact;
#endif

    { 
        cout << "determine vtrue and write coeffs and matlab output to HDD" << endl;
        InfiniteVector<double,Index> temp_iv;
        char filename[250];
        for (unsigned int i=0; i<=number_of_timesteps;++i)
        {
            if (_PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 == 1) // if uexact is meaningful, we should use it! => overwrite parabolic.forward_solution_ for determining ydata!
            {
                vexact.set_time(i/(double)number_of_timesteps);
                //ctgramian.basis().expand(&vexact, false, spatial_jmax, temp_iv);
                ctgramian.set_f(&vexact);
                CDD1_SOLVE(ctgramian, tolerance, temp_iv, spatial_jmax, tensor_simple);
                temp_iv.scale(&ctgramian,-1);
                temp_iv.compress(1e-14);
                cout << "main:: compute ydata_matlab_output backward case." << endl;
                cout << "       at t=" << i/(double)number_of_timesteps << " ell_2 error = " << l2_norm(parabolic.backward_solution_[i] - temp_iv) << endl;
                parabolic.backward_solution_[i]=temp_iv;
            }

            ostringstream output_filename2;
            output_filename2 << true_backward_solution_filename_start << i << true_backward_solution_filename_matlab_end;
            resultstream.open(output_filename2.str().c_str());
            SampledMapping<dim> vi_plot(evaluate(ctgramian.basis(), parabolic.backward_solution_[i], true, resolution));
            vi_plot.matlab_output(resultstream);
            resultstream.close();
            sprintf(filename,"%s%d%s",true_backward_solution_filename_start,i,true_backward_solution_filename_coeffs_end);
            cout << "main:: write vtrue_coeffs[" << i << "] to file = " << endl << filename << endl;
            parabolic.backward_solution_[i].writeToFile(filename);
        }
        //vtrue_coeffs = parabolic.backward_solution_;
    }
*/


    
#if _TEST_PARABOLIC_TOOLS == 1

#if _DIMENSION == 1
    FixedVector<double, (1<<(_HAAR_JMAX+1))> data_haar_gen_coeffs;
#else
    FixedMatrix<double, (1<<(_HAAR_JMAX+1)), (1<<(_HAAR_JMAX+1))> data_haar_gen_coeffs;
#endif
    // TEST some of the parabolic_tools
    // data_haar_gen_coeffs should contain the Haar generator coefficients of parabolic.forward_solution_[number_of_timesteps]

    cout << "######################################" << endl;
    cout << "######################################" << endl;
    cout << "begin tests of parabolic_tools" << endl;
    cout << "######################################" << endl;
    cout << "######################################" << endl;

    //cout << "output of precise_evaluate(parabolic.forward_solution_[" << number_of_timesteps << "]) == " << endl;
#if _DIMENSION == 1

    Exact_Sol1D<forward_problem_number> temp_fkt;
#else
    Exact_Sol2D<31> temp_fkt;
#endif
    data_haar_gen_coeffs.scale(0.0);
    //data_haar_gen_coeffs.set_entry(4,2,1.0); // wavelet #57
    //data_haar_gen_coeffs.set_entry(5,2,-1.0);
    //data_haar_gen_coeffs.set_entry(4,3,-1.0);
    //data_haar_gen_coeffs.set_entry(5,3,1.0);
    
    FixedMatrix<double, (1<<(_HAAR_JMAX+1)), (1<<(_HAAR_JMAX+1))> temp_haar_wav_coeffs;
    InfiniteVector<double, int> temp_haar_wav_coeffs_iv;
    
    cout << "test precise_evaluate (*fkt)" << endl;
    precise_evaluate(&temp_fkt, d, data_haar_gen_coeffs);
    cout << "data_haar_gen_coeffs" << endl << data_haar_gen_coeffs << endl;
    cout << "plot_haar_gen_coeffs(data_haar_gen_coeffs, " << haar_jmax+1 << ") =" << endl;
    plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
    
    cout << "test haar_wavelet_transform" << endl;
    haar_wavelet_transform(data_haar_gen_coeffs,temp_haar_wav_coeffs);
    cout << "temp_haar_wav_coeffs" << endl << temp_haar_wav_coeffs << endl;
    
    cout << "test transform_fixedToIV" << endl;
    transform_fixedToIV (temp_haar_wav_coeffs, temp_haar_wav_coeffs_iv);
    cout << "temp_haar_wav_coeffs_iv" << endl << temp_haar_wav_coeffs_iv << endl;
    
    cout << "test transform_IVTofixed" << endl;
    transform_IVTofixed (temp_haar_wav_coeffs_iv, temp_haar_wav_coeffs);
    cout << "temp_haar_wav_coeffs" << endl << temp_haar_wav_coeffs << endl;
    
    cout << "test inverse_haar_wavelet_transform" << endl;
    inverse_haar_wavelet_transform(temp_haar_wav_coeffs,data_haar_gen_coeffs);
    cout << "data_haar_gen_coeffs" << endl << data_haar_gen_coeffs << endl;
    cout << "plot_haar_gen_coeffs(data_haar_gen_coeffs, " << haar_jmax+1 << ") =" << endl;
    plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
    
    //cout << "data_haar_gen_coeffs = " << endl << data_haar_gen_coeffs << endl;
    
    // CLEANUP ydata does no longer exists as a Haar wavelet/generator coefficient vector
    //cout << "testing transform_fixedToIV(data_haar_gen_coeffs,ydata["<<number_of_timesteps << "])" << endl;
    //cout << "ydata[" << number_of_timesteps << "] = " << ydata[number_of_timesteps] << endl;
    //cout << "testing transform_IVTofixed(ydata["<<number_of_timesteps<<"], data_haar_gen_coeffs)"<<endl;
    //transform_IVTofixed(ydata[number_of_timesteps], data_haar_gen_coeffs);
    //cout << "data_haar_gen_coeffs = " << endl << data_haar_gen_coeffs << endl;

    cout << "compare the matlab output of parabolic.forward_solution_[" << number_of_timesteps << "], located in .../ydata/u_problem_?_wnum_?_primbs_d_dt_3_3_bc_?_jmax_?_tstep_?_?.m"
            " with the following:" << endl;
    cout << "plot_haar_gen_coeffs(data_haar_gen_coeffs, " << haar_jmax+1 << ") =" << endl;
    plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);

    cout << "testing Haar wavelet transform" << endl;

#if _DIMENSION == 1
    FixedVector<double, (1<<(haar_jmax+1))> data_haar_wav_coeffs;
#else
    FixedMatrix<double, (1<<(haar_jmax+1)), (1<<(haar_jmax+1))> data_haar_wav_coeffs;
#endif
    cout << "data_haar_gen_coeffs = " << endl << data_haar_gen_coeffs << endl;
    cout << "haar_wavelet_transform(data_haar_gen_coeffs, data_haar_wav_coeffs) => " << endl;
    haar_wavelet_transform(data_haar_gen_coeffs, data_haar_wav_coeffs);
    cout << "data_haar_wav_coeffs = " << endl << data_haar_wav_coeffs << endl;
    cout << "inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs) => " << endl;
    inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs);
    cout << "data_haar_gen_coeffs = " << endl << data_haar_gen_coeffs << endl;


    cout << "Testing fixed_vector_READ/WRITE routines" << endl;
    cout << "writing data_haar_gen_coeffs to filename " << endl << temp_filename << endl;
#if _DIMENSION == 1
    fixed_vector_writeToFile(temp_filename, data_haar_gen_coeffs);
    fixed_vector_readFromFile(temp_filename, data_haar_wav_coeffs);
#else
    fixed_matrix_writeToFile(temp_filename, data_haar_gen_coeffs);
    fixed_matrix_readFromFile(temp_filename, data_haar_wav_coeffs);
#endif
    data_haar_gen_coeffs.scale(-1);
    data_haar_wav_coeffs.add(data_haar_gen_coeffs);
    
    cout << "difference between original and stored vector = " << endl;
    cout << "diff = " << endl << data_haar_wav_coeffs << endl;

    cout << "testing compute_update_w with parabolic.forward_solution_, parabolic.backward_solution_ and ydata" << endl;
    Array1D<InfiniteVector<double, int> > temp_w_update;
    temp_w_update.resize(number_of_timesteps+1);

    cout << "p.fs,p.bs,ydata" << endl;
    compute_update_w<Basis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(&(parabolic.problem_->basis_), parabolic.forward_solution_, parabolic.backward_solution_, temp_w_update);
    for (unsigned int i = 0; i<= number_of_timesteps; ++i)
    {
        cout << "w_update["<<i<<"] = " << endl;
        transform_IVTofixed(temp_w_update[i], data_haar_wav_coeffs);
        inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs);
        plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
    }
    Array1D<InfiniteVector<double,Index> > empty;
    empty.resize(number_of_timesteps+1);
    cout << "0,p.bs,ydata" << endl;
    compute_update_w<Basis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(&(parabolic.problem_->basis_), empty, parabolic.backward_solution_, temp_w_update);
    for (unsigned int i = 0; i<= number_of_timesteps; ++i)
    {
        cout << "w_update[" << i << "] = " << endl;
        //cout << endl;
        //cout << "temp_w_update["<<i<<"] = " << endl << temp_w_update[i] << endl;
        transform_IVTofixed(temp_w_update[i], data_haar_wav_coeffs);
        //cout << "IVTofixed => " << endl << data_haar_wav_coeffs << endl;
        inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs);
        //cout << "inverse_haar_wavelet_transform => " << endl << data_haar_gen_coeffs << endl;
        //cout << "plot data_haar_gen_coeffs" << endl;
        plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
    }
    cout << "p.fs,0,ydata" << endl;
    compute_update_w<Basis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(&(parabolic.problem_->basis_), parabolic.forward_solution_, empty, temp_w_update);
    for (unsigned int i = 0; i<= number_of_timesteps; ++i)
    {
        cout << "main::w_update["<<i<<"] = " << endl;
        transform_IVTofixed(temp_w_update[i], data_haar_wav_coeffs);
        inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs);
        plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
        //cout << "main:: temp_w_update[" << i << "] = " << endl << temp_w_update[i] << endl;
        //cout << "main:: IVTofixed => data_haar_wav_coeffs = " << endl << data_haar_wav_coeffs << endl;
        //cout << "main:: invHWT => data_haar_gen_coeffs = " << endl << data_haar_gen_coeffs << endl;
    }
    //Array1D<InfiniteVector<double,int> > empty2;
    //empty2.resize(number_of_timesteps+1);
    cout << "p.fs,p.bs,0" << endl;
    compute_update_w<Basis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(&(parabolic.problem_->basis_), parabolic.forward_solution_, parabolic.backward_solution_, temp_w_update);
    for (unsigned int i = 0; i<= number_of_timesteps; ++i)
    {
        cout << "w_update["<<i<<"] = " << endl;
        transform_IVTofixed(temp_w_update[i], data_haar_wav_coeffs);
        inverse_haar_wavelet_transform(data_haar_wav_coeffs, data_haar_gen_coeffs);
        plot_haar_gen_coeffs(data_haar_gen_coeffs, haar_jmax+1);
    }

    cout << "test create_noise" << endl;


    Array1D<InfiniteVector<double, Index> > temp_noise_coeffs;
    temp_noise_coeffs.resize(number_of_timesteps+1);
    create_noise(1,1e-4,&(parabolic.problem_->basis_), temp_noise_coeffs);
    for (unsigned int i=0; i< noise_coeffs.size(); ++i)
    {
        cout << "noise_coeffs[i] = " << endl << noise_coeffs[i] << endl;
    }

    Array1D<double> temp_noise_norms_sqr;
    temp_noise_norms_sqr.resize(temp_noise_coeffs.size());
    double temp_double;

    for (unsigned int i=0; i< temp_noise_coeffs.size(); ++i)
    {
        temp_double+= l2_norm_sqr(temp_noise_coeffs[i]);
    }
    temp_double = sqrt(temp_double);
    cout << "norm of computed noise is = " << temp_double << " whereas delta = " << 1e-4 << endl;




    cout << "######################################" << endl;
    cout << "######################################" << endl;
    cout << "end tests of parabolic_tools" << endl;
    cout << "######################################" << endl;
    cout << "######################################" << endl;

#endif

        cout << "main:: computation of simulated data completed." << endl;

#else
    // _COMPUTE_SIMULATED_DATA == 0
    cout << "main:: load simulated data (uexact_coeffs) from file " << true_forward_solution_filename_start << "i" << true_forward_solution_filename_coeffs_end << endl;
    {
        char filename[250];
        //cout << "main:: load u_true, v_true from files " << endl;
        //cout << true_forward_solution_filename_start << "i" << true_forward_solution_filename_coeffs_end << endl;
        //cout << true_backward_solution_filename_start << "i" << true_backward_solution_filename_coeffs_end << endl;

        for (unsigned int i=0; i <= number_of_timesteps; ++i)
        {
            sprintf(filename,"%s%d%s",true_forward_solution_filename_start,i,true_forward_solution_filename_coeffs_end);
            utrue_coeffs[i].readFromFile(filename);
            //sprintf(filename,"%s%d%s",true_backward_solution_filename_start,i,true_backward_solution_filename_coeffs_end);
            //vtrue_coeffs[i].readFromFile(filename);
        }
    }
    cout << "main:: data has been loaded from disk" << endl;
#endif

   
// compute&store noise OR load noise from file
    FixedArray1D<InfiniteVector<double, int>, number_of_timesteps+1 > noise_coeffs;
    //noise_coeffs.resize(number_of_timesteps+1);
    {
        char filename[250];
        if (compute_noise == 1)
        {
            // compute & store noise
            create_noise(1.0,delta,parabolic_problem.basis(), noise_coeffs);
            for (unsigned int i=0; i<= number_of_timesteps; ++i)
            {
                sprintf(filename,"%s%d%s",noise_filename_start,i,".iv");
                cout << "main:: write noise_coeffs[" << i << "] to file = " << endl << filename << endl;
                writeIVToFile(noise_coeffs[i],filename);
            }
        }
        else
        {
            for (unsigned int i=0; i<= number_of_timesteps; ++i)
            {
                sprintf(filename,"%s%d%s",noise_filename_start,i,".iv");
                cout << "main:: load noise_coeffs[" << i << "] from file = " << endl << filename << endl;
                readIVFromFile(noise_coeffs[i],filename);
            }
        }
    }

    //  add noise where it is needed
    cout << "main:: adding noise" << endl;
    FixedArray1D<InfiniteVector<double,int>, number_of_timesteps+1 > uobserved_coeffs; // utrue+noise=observed data
    for (unsigned int i=0; i<= number_of_timesteps;++i)
    {
        uobserved_coeffs[i]=utrue_coeffs[i];
        uobserved_coeffs[i].add(noise_coeffs[i]);
    }
    parabolic_problem.u0_.add(noise_coeffs[0]);
    parabolic_problem.set_true_forward_solution_preprocessed(uobserved_coeffs, tolerance_APPLY);
    // side effect: parabolic_forward_solution_ was used as a temporary storage! == is currently not meaningful

    // compute Norm of utrue_coeffs. This is done to get an idea of the relative magnitude of the noise. Result is written to standard output and the logfile. Value does not influence further computation
    double utrue_coeffs_norm(0), uobserved_coeffs_norm(0);
    {
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            utrue_coeffs_norm+=l2_norm_sqr(utrue_coeffs[i]);
            uobserved_coeffs_norm+=l2_norm_sqr(uobserved_coeffs[i]);
        }
        utrue_coeffs_norm = sqrt(utrue_coeffs_norm);
        uobserved_coeffs_norm = sqrt(uobserved_coeffs_norm);
        cout << "Magnitude of noise delta = " << delta << endl;
        cout << "l2 norm of ydata (without noise) is = " << utrue_coeffs_norm << "  => relative noise level = " << (delta / utrue_coeffs_norm) << endl;
        cout << "l2 norm of ydata (with noise) is = " << uobserved_coeffs_norm << "  => relative noise level = " << (delta / uobserved_coeffs_norm) << endl;
    }

#if _SOLVE_INVERSE_PROBLEM == 0
    cout << "main :: skipping inverse iteration. program ends." << endl;
    abort();
#endif

    cout << "main:: begin inverse iteration" << endl;
    //Array1D<InfiniteVector<double, int> > w, w_limit;
    //w.resize(number_of_timesteps+1);
    //w_limit.resize(number_of_timesteps+1);
    FixedArray1D<Array1D<FixedMatrix<double, onedimhaarcount> >, number_of_timesteps +1> w, w_limit;
    //FixedArray1D<InfiniteVector<double, int>, number_of_timesteps+1 > w_update;
    FixedArray1D<Array1D<FixedMatrix<double, onedimhaarcount> >, number_of_timesteps +1> w_next;
    
    // wlimit either contains the (unknown) true value wtrue, or it contains the limit of the iterative process (denoted wlimit). We have access to it by a previous computation
    // wlimit is used to compute the error of approximaton ||w_current - w_limit||_x
    if (w_limit_number == 0)
    {
        cout << "main:: w_limit_number == 0. Use wtrue as value of wlimit." << endl;
        for (unsigned int i=0; i<=number_of_timesteps; ++i)
        {
            w_limit[i] = w_true[i];
        }
    }
    else
    {
        cout << "main:: load wlimit from disk" << endl;
        char filename[250];
        for (unsigned int i=0; i<=number_of_timesteps; ++i)
        {
            w_limit[i].resize(qtbasis.get_nop());
            for (int patch = 0; patch < qtbasis.get_nop(); ++patch)
            {
                sprintf(filename,"%s%d_tstep%d_%d_patch_%d%s",w_coefficients_filename_start,w_limit_number,number_of_timesteps,i,patch,w_coefficients_filename_end);
                cout << "loading wlimit[" << i << "][" << patch << "] from file" << endl << filename << endl;
                fixed_matrix_readFromFile(w_limit[i][patch],filename);
                //w_limit[i].readFromFile(filename);
            }
        }
    }

     
    

    bool done = false; // determines whether the main iteration should halt

    // local variables for the log file
    //double relative_error(0); // store relative error ||w_update||_2/||w||_2 ; unused at the moment.
    int iteration_count = start_iteration;
    
    clock_t tstart2, tend2,tend3, tend4;
    //double temp_d;
    //unsigned int temp_i;
    
    tstart2=clock();

    if (start_iteration == 0)
    {
        cout << "main:: initialize coupling matrix W for the iteration" << endl;
        initialize_coupling_matrix(w, _SECOND_SETUP_COUPLING_MATRIX_W);
        //w_next = w;
        logstream.open(logfile);
        initialize_logstream(logstream, compute_noise, first_setup_coupling_matrix_w, max_inverse_iterations, stopping_tolerance, spatial_jmax, w_limit_number, alpha, delta, utrue_coeffs_norm, uobserved_coeffs_norm, number_of_timesteps);
        log_reconstruction_errors(logstream, w, w_limit);
    }
    else
    {
        char filename[250];
        cout << "main:: load coupling matrix W from disk" << endl;
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            w[i].resize(qtbasis.get_nop());
            for (int patch = 0; patch < qtbasis.get_nop(); ++patch)
            {
                sprintf(filename,"%s%d_t%d_%d_patch_%d%s",w_coefficients_filename_start,start_iteration,number_of_timesteps,i,patch,w_coefficients_filename_end);
                cout << "loading W[" << i << "][" << patch << "] from file" << endl << filename << endl;
                fixed_matrix_readFromFile(w[i][patch],filename);
            }
        }
        logstream.open(logfile, ios_base::app);
        logstream << "% RESUMING at Iteration " << start_iteration << " for additional " << (max_inverse_iterations - start_iteration) << " iterations\n";
    }
    logstream.close();
    
    cout << "main:: begin main computation" << endl;
    while (!done)
    {
        tstart = clock();
        //w = w_next;
        parabolic_problem.set_W(w);
        tend = clock();
        solve_parabolic_problem(parabolic_problem, method, true, increment_tolerance, tolerance_APPLY);
        tend2 = clock();
        //cout << "main:: solve_parabolic_problem:: time backward case" << endl;
        solve_parabolic_problem(parabolic_problem, method, false, increment_tolerance, tolerance_APPLY);
        tend3 = clock();
        
        // shrinkage_iteration computes w_update = -h*u
#if _SAVE_UPDATE_U_H_TO_FILE == 1
        // plot stuff in the last iteration. for efficience reasons this is done in the compute_update routine
        if ((max_inverse_iterations > 0) && (iteration_count == (max_inverse_iterations-1) ))
        {
            // this is the last iteration!
            //cout << "main:: last iteration. plotting current u,h and w_update (=u*h) to file =" << endl << plotfile << endl;
            plotstream.open(plotfile);
            shrinkage_iteration(parabolic_problem, w, alpha, w_next, &plotstream);
            //shrinkage_iteration<QTBasis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(parabolic_problem, w, alpha, w_next, &plotstream);
            plotstream.close();
        }
        else
        {
            shrinkage_iteration(parabolic_problem, w, alpha, w_next, NULL);
            //shrinkage_iteration<QTBasis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(parabolic_problem, w, alpha, w_next, NULL);
        }
#else
        // we are not interested in plots of u,h,w_update
        // do not store w_update, u,h in their Haar wavelet expansion:
        shrinkage_iteration<QTBasis, (1<< (_HAAR_JMAX+1)), _DIMENSION>(parabolic_problem, w, alpha, w_next, NULL);
#endif
        
        w = w_next;
        tend4 = clock();
        
        // check stopping criterion.
        // done = check_stopping_criterion(wtrue,w,w_update,w_new)
        // do not use the unshrinked w_update! You may run into loops

        /*
        // stopping criterion is only active if max_inverse_iterations is not specified
        w_update = w;
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            w_update[i].subtract(w_next[i]);
        }
        if ((max_inverse_iterations == 0) && (iteration_count > 0))
        {            
            //done = relative_change_criterion_ell1(w, w_update, stopping_tolerance, relative_error);
            if (iteration_count == 1)
            {
                //logstream << "r = " << relative_error << "; % relative error as computed by relative_change_criterion_ell1 \n";
            }
            else
            {
                //logstream << "r = [r; " << relative_error << "];\n";
            }
        }
        */
        ++iteration_count;
        cout << "Iteration " << iteration_count << " completed. time needed = " << (double)(tend4-tstart)/CLOCKS_PER_SEC << " seconds" << endl;
        logstream.open(logfile, ios_base::app);
        //total_time = (double)(tend4-tstart)/CLOCKS_PER_SEC;
        //assemble_time = (double)(tend-tstart)/CLOCKS_PER_SEC;
        //forward_time = (double)(tend2-tend)/CLOCKS_PER_SEC;
        //backward_time = (double)(tend3-tend2)/CLOCKS_PER_SEC;
        //update_time = (double)(tend4-tend3)/CLOCKS_PER_SEC;

        // Note: we use w here (since u,h are computed with it) and NOT w_next
        // w_next in the last iteration is treated in finalize_logstream
        /*
        update_logstream<Basis, _NUMBER_OF_TIME_STEPS>(logstream, iteration_count, 
                         (double)(tend4-tstart)/CLOCKS_PER_SEC, 
                         (double)(tend-tstart)/CLOCKS_PER_SEC, 
                         (double)(tend2-tend)/CLOCKS_PER_SEC, 
                         (double)(tend3-tend2)/CLOCKS_PER_SEC, 
                         (double)(tend4-tend3)/CLOCKS_PER_SEC, 
                         w, w_limit, parabolic.forward_solution_, parabolic.backward_solution_,utrue_coeffs);
        */
        logstream << "% Iteration " << iteration_count << "\n";
        // write computation time to logfile
        // logstream << "t = [t; " << total_time << " " << assemble_time << " " << forward_time << " " << backward_time << " " << update_time << " " << shrinkage_time << "];\n";
        logstream << "t = [t; " << (double)(tend4-tstart)/CLOCKS_PER_SEC
                << " " << (double)(tend-tstart)/CLOCKS_PER_SEC
                << " " << (double)(tend2-tend)/CLOCKS_PER_SEC
                << " " << (double)(tend3-tend2)/CLOCKS_PER_SEC
                << " " << (double)(tend4-tend3)/CLOCKS_PER_SEC 
                << "];\n";
/*
        log_times<Basis, _NUMBER_OF_TIME_STEPS>(logstream, iteration_count, 
                         (double)(tend4-tstart)/CLOCKS_PER_SEC, 
                         (double)(tend-tstart)/CLOCKS_PER_SEC, 
                         (double)(tend2-tend)/CLOCKS_PER_SEC, 
                         (double)(tend3-tend2)/CLOCKS_PER_SEC, 
                         (double)(tend4-tend3)/CLOCKS_PER_SEC, 
 * */
#if _COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL == 1
        log_solution_errors(logstream, parabolic_problem.forward_solution_, uobserved_coeffs);
#endif
#if _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL == 1
        log_reconstruction_errors(logstream,w, w_limit);
#endif
        if ((max_inverse_iterations > 0) && (max_inverse_iterations == iteration_count))
        {
            done = true;
            //logstream << "% Iteration maximum reached. No error criterion was checked\n";
            logstream << "% Iteration complete" << endl;
            
        }
        logstream.close(); // writes logstream to HDD
    } // end of while(!done)
    
    tend = clock();
    time = (double)(tend-tstart2)/CLOCKS_PER_SEC; // total
    
    logstream.open(logfile, ios_base::app);
    //finalize_logstream<_NUMBER_OF_TIME_STEPS>(logstream, time,max_inverse_iterations, w_next, w_limit);
    if (start_iteration == 0)
    {
        logstream << "time_total = " << time << ";\n";
    }
    else
    {
        logstream << "time_total = time_total + " << time << ";\n";
    }
    logstream.close();
    
    plotstream.open(plotfile, ios_base::app);
    plot_solutions(plotstream, w, w_true, w_limit);
    plotstream.close();
    // plot stuff ...

    // store solutions on hard drive
#if _SAVE_SPATIAL_SOLUTIONS_TO_FILE == 1
// CLEANUP
    cout << "main :: save the forward and backward solutions (for every t_step i) to files" << endl
         << forward_solution_filename_start << "i" << forward_solution_filename_end << endl
         << backward_solution_filename_start << "i" << backward_solution_filename_end << endl;
    Array1D<SampledMapping<dim> > sol_i_sampling;
    for (unsigned int i=0; i<= number_of_timesteps;++i)
    {
        ostringstream output_filename;
        output_filename << forward_solution_filename_start << i << forward_solution_filename_end;
        resultstream.open(output_filename.str().c_str());
        
//         ergebnis2 = qtbasis.sampled_output(u_cdd1,
//                true,
//                6);
//         
        sol_i_sampling = qtbasis.sampled_output(parabolic_problem.forward_solution_[i],
                true,
                resolution);
        matlab_output(sol_i_sampling, resultstream);
        resultstream.close();

        ostringstream output_filename2;
        output_filename2 << backward_solution_filename_start << i << backward_solution_filename_end;
        resultstream.open(output_filename2.str().c_str());
        sol_i_sampling = qtbasis.sampled_output(parabolic_problem.backward_solution_[i],
                true,
                resolution);
        matlab_output(sol_i_sampling, resultstream);
        resultstream.close();
    }
#endif // finished storing solutions to HDD
        
    /* store W on HDD
     * - to make computation of error ||W*-W_n|| possible
     * - to make resume of inverse iteration possible
     * one could store matlab output in every iteration and produce a movie
     */
    char filename[250];
    for (unsigned int i=0; i<= number_of_timesteps; ++i)
    {
        //sprintf(filename,"%s%d_t%d_%d%s",w_coefficients_filename_start,max_inverse_iterations,number_of_timesteps,i,w_coefficients_filename_end);
        //cout << "writing solution to file" << endl << filename << endl;
        //w_next[i].writeToFile(filename);
        for (int patch = 0; patch < qtbasis.get_nop(); ++patch)
        {
            sprintf(filename,"%s%d_t%d_%d_patch_%d%s",w_coefficients_filename_start,max_inverse_iterations,number_of_timesteps,i,patch,w_coefficients_filename_end);
            cout << "writing W[" << i << "][" << patch << "] to file" << endl << filename << endl;
            fixed_matrix_writeToFile(w[i][patch],filename);
        }
    }
    if (w_limit_number <  max_inverse_iterations)
    {
        cout << "main:: computed approximation to W is superior!" << endl;
    }
    if (iteration_count < max_inverse_iterations)
    {
        cout << "main:: Caution! Observe that the inverse iteration stopped after " << iteration_count << " iterations!" << endl;
    }
    
    cout << "Output has been written to the logfile:" << endl << logfile << endl;
    cout << "Plots have been written to the plotfile:" << endl << plotfile << endl;
    cout << "time needed in total = " << time << " seconds" << endl;
} // end of main

//#endif

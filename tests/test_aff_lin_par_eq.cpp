// define verbosity and tbasis:
// more output for the cached problemm (in normA())
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 0
// normA uses setup_stiffness_matrix. here the verbosity of the call is controled:
#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
// for verbose output of CDD1
#define _WAVELETTL_CDD1_VERBOSITY 0
#define _NO_JPMAX_grt_70_WARNING 0 // turn off warning in APPLY_TENSOR
// switch between isotropic and anisotropic Wavelets (in cdd1.h)
#define _WAVELETTL_USE_TBASIS 1

// define mathematical setting
#define _DIMENSION 2 // only used for the initialization of the pbasis used to define Index
#define _D 3 // only used for the initialization of the pbasis used to define Index
#define _DT 3 // only used for the initialization of the pbasis used to define Index

#define _PROBLEM_NO 14
#define _BOUNDARY_CONDITIONS 1 // 1D :: 0 = tt, 1 = tf, 2 = ft, 3 = ff; 2D :: 0 = tttt, 1=ffff
#define _PAR_JMAX 3
#define _SPATIAL_JMAX 6 // maximal 9 im 3D Fall, da Gramian Matrizen nicht weiter vorberechnet
#define _TIME_DISCRETIZATION_GRANULARITY 3 // we use 1<<_T_D_G many time steps for the parabolic equation (time discretization has this number+1 many points)
//#define _TIME_FORWARD_CASE 1 // 1=time forward case (prob, 0= time backward case
//#define _OFFSET 0

// define behaiviour of this test file
#define _TEST_INITIALIZATION 0
//#define _COMPARE_WITH_LINPAREQ 1
#define _TEST_LIN_PAR_EQ_TENSOR 0
#define _TEST_AFF_LIN_PAR_EQ 1

#define _EXPANSIONTYPE_F_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_f; 1 liefert Unterschied falls f konstant ist und man TIME_CONSTANT =0,1 wählt
#define _EXPANSIONTYPE_FT_ 0 // 1 = Thorsten, 0 = Ulli in lin_par_eq_tensor::evaluate_ft

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

#include <parabolic/example_parabolic_problems.cpp>
#include <parabolic/aff_lin_par_eq.h>

#include <parabolic/lin_par_eq_tensor.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>

#include <algebra/fixed_vector.h>
#include <algebra/fixed_matrix.h>

#include <parabolic/parabolic_tools.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

int main()
{
    cout << "Testing::  affine linear parabolic equation" << endl;

    clock_t tstart, tend;
    double time;

    const int d  = _D;
    const int dT = _DT;
    const int dim = _DIMENSION;
    //const int offset = _OFFSET; // Radius for the index set Lambda
    typedef PBasis<d,dT> Basis1d;

    typedef TensorBasis<Basis1d,dim> Basis;
    typedef Basis::Index Index;
    //typedef Index::level_type index_lt;
    const unsigned int N=_PROBLEM_NO;
#if _DIMENSION == 1
    Exact_Sol1D<N> uexact;
#else
    Exact_Sol2D<N> uexact;
#endif

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


    const double diff_coeff(1);
    const int number_of_timesteps (1<<_TIME_DISCRETIZATION_GRANULARITY);
    Array1D<InfiniteVector<double, int> > w;
    w.resize(number_of_timesteps+1);
    Array1D<double> time_discretization;
    time_discretization.resize(number_of_timesteps+1);
    for (unsigned int i = 0; i <= number_of_timesteps; ++i)
    {
        w[i].clear(); // an empty w vector should result in the heat equation (for d=1)
        //w[i].set_coefficient(3,3.24);
        //cout << "main:: w[i] = " << w[i] << endl;
        time_discretization[i] = i / (double)(number_of_timesteps);
    }
    /*
    w[0].set_coefficient(3,120);
    w[0].set_coefficient(5,120);
    w[1].set_coefficient(5,120);
    w[3].set_coefficient(7,120);
    */
    InfiniteVector<double, int> empty_iv;
    cout << "compare with an empty InfiniteVector<double, int>: empty = " << empty_iv << endl;
    empty_iv.clear();
    cout << "compare with an empty InfiniteVector<double, int>: empty = " << empty_iv << endl;

    //const char* f_filename = "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_ffff.iv";
    //const char* haar_gramian_filename_start = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_";
    //const char* haar_gramian_filename_end = "_gramian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd";
    //const char* laplacian_filename = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd";

    char f_filename[250];
    //sprintf(f_filename, "/home/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_%s.iv",bc_string);
    sprintf(f_filename, "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_%s.iv",bc_string);
    //const char* f_filename = "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_ff.iv";

    //const char* haar_gramian_filename_start = "/home/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_";
    const char* haar_gramian_filename_start = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_";
    char haar_gramian_filename_end[250];
    //const char* haar_gramian_filename_end = "_gramian_primbs_d_dt_3_3_bc_ff_jmax_13_npcd";
    sprintf(haar_gramian_filename_end, "_gramian_primbs_d_dt_3_3_bc_%s_jmax_%d_npcd",bc_string,_SPATIAL_JMAX);

    char laplacian_filename[250];
    //sprintf(laplacian_filename, "/home/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_%s_jmax_%d_npcd",bc_string,_SPATIAL_JMAX);
    sprintf(laplacian_filename, "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_%s_jmax_%d_npcd",bc_string,_SPATIAL_JMAX);
    //const char* laplacian_filename = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_ff_jmax_13_npcd";

    const int par_jmax = _PAR_JMAX;
    const int spatial_jmax = _SPATIAL_JMAX;
    const double tolerance = 1e-6;
    cout << "main:: compute u0T" << endl;
    InfiniteVector<double,Index> u0T;
    // gleichbedeutend mit :: compute ctproblem.basis().expand(&uexact, TRUE, jmax, u0T);
    char matrix_filename[250];
    //haar_wav_0_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
    sprintf(matrix_filename, "%s%d%s",haar_gramian_filename_start, 0, haar_gramian_filename_end);
    //cout << "assemble_W:: this is stored in memory: raw_cache_["<< it.index() << "] = " << (raw_cache_[it.index()]) << endl;
    cout << "main:: initialize CompressedProblemFromMatrix(gram)" << endl;


    //ConstantFunction<dim> zero_function(Vector<double>(1, "0.0"));
    IdentityBVP<dim> identity_bvp(&uexact);
    TensorEquation<Basis1d,dim,Basis> gram_problem (&identity_bvp, bc, false); // 21 mb speicher!!!!, jmax=14 => 43 mb

    // teste Größe durch kopieren:
    //Array1D<std::pair<Index,double> > fcoeffs2(gram_problem.fcoeffs); // unwesentlicher Speicherbedarf

    //Basis temp_basis(bc); // 43 mb, dh +22 mb speicher! bzw 85mb für level 14

    // was macht die Basis so groß?
    //Array1D<Index> full_collection2(temp_basis.full_collection); // kaum etwas

    //Basis temp_basis2(temp_basis.bases()); // werden nur pointer kopiert. kein zusätzlicher Speicheraufwand
    //Basis1d temp_bas1d(bc[0],bc[1]); // 64 mb speicher, +21 mb. Das scheint der normale Aufwand für eine 1D Basis zu sein! p_basis.h::Z.390:jmax=14=> 128mb



    gram_problem.set_jmax(spatial_jmax,true);



// CLEANUP another problem?
    //cout << "gram_problem.basis_.degrees_of_freedom() = " << gram_problem.basis_.degrees_of_freedom() << endl;

    SparseMatrix<double> gram_matrix(1);
    cout << "main:: initialize gram_matrix with file = " << matrix_filename << endl;
    gram_matrix.matlab_input(matrix_filename);
// CLEANUP
    //cout << gram_problem.a(gram_problem.basis_.get_wavelet(12),gram_problem.basis_.get_wavelet(12)) << endl;
    CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> > ctgramian(&gram_problem, &gram_matrix, 5.5, 35.0, false);
    //ctgramian.entries_cache.matlab_input(matrix_filename);
    cout << "main:: use CDD1 to compute u0T" << endl;


    /*
// CLEANUP another problem?
    cout << "ctgramian.basis().degrees_of_freedom() = " << ctgramian.basis().degrees_of_freedom() << endl;
    cout << "ctgramian.problem->basis().degrees_of_freedom() = " << ctgramian.problem->basis().degrees_of_freedom() << endl;

    cout << "ctgramian.basis().first_generator() = " << ctgramian.basis().first_generator() << endl;
    cout << "ctgramian.basis().last_generator() = " << ctgramian.basis().last_generator() << endl;
    cout << "ctgramian.basis().first_wavelet(3) = " << ctgramian.basis().first_wavelet(3) << endl;
    cout << "ctgramian.basis().last_wavelet(3) = " << ctgramian.basis().last_wavelet(3) << endl;

    cout << "ctgramian.basis().bases()[0]->Deltasize(3) = " << ctgramian.basis().bases()[0]->Deltasize(3) << endl;
    cout << "ctgramian.basis().bases()[0]->Nablasize(3) = " << ctgramian.basis().bases()[0]->Nablasize(3) << endl;
    //ctgramian.s
 */

    uexact.set_time(0);
    ctgramian.set_f(&uexact);
// CLEANUP
    //cout << "main :: u0T = " << endl << u0T << endl;
    //cout << ctgramian.entries_cache->get_entry(12,12) << endl;;
    CDD1_SOLVE(ctgramian, tolerance, u0T, spatial_jmax, tensor_simple);
    //cout << "main:: u0T after CDD1_SOLVE :: " << endl << u0T << endl;
    u0T.scale(&ctgramian,-1);
    u0T.compress(1e-14);

    const char* temp_filename = "/import/shared/friedrich/source/precomputed/delete_on_sight";
    cout << "main:: writing u0T to matlab_output. filename = " << endl << temp_filename << ".m" << endl;

    const int resolution = d+dim+((dim == 1)? (_SPATIAL_JMAX-3):(_SPATIAL_JMAX-6)) +1;
    std::ofstream resultstream; // ostream used for matlab output by SampledMapping

    /*
    ostringstream output_filename;
    output_filename << temp_filename << "_false.m";
    resultstream.open(output_filename.str().c_str());
    SampledMapping<dim> u0T_plot(evaluate(ctgramian.basis(), u0T, false, resolution));
    u0T_plot.matlab_output(resultstream);
    resultstream.close();

    ostringstream output_filename2;
    output_filename2 << temp_filename << "_true.m";
    resultstream.open(output_filename2.str().c_str());
    SampledMapping<dim> u0T_plot2(evaluate(ctgramian.basis(), u0T, true, resolution));
    u0T_plot2.matlab_output(resultstream);
    resultstream.close();
     */


/* von weiter unten kopiert:
    InfiniteVector<double,Index> u0T2;
    uexact.set_time(0);
    ctgramian2.set_f(&uexact);
    CDD1_SOLVE(ctgramian2, tolerance, u0T2, spatial_jmax, tensor_simple);
    u0T2.scale(&ctgramian2,-1);
    u0T2.compress(1e-14);
*/
    /*
// CLEANUP
    PoissonBVP<dim> elliptic(&uexact);
    TensorEquation<Basis1d,dim,Basis> eq(&elliptic, bc, false);
    eq.set_jmax(spatial_jmax, false);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem(&eq,5.5,35.0);
    cout << "ctproblem.basis().degrees_of_freedom() = " << ctproblem.basis().degrees_of_freedom() << endl;
*/

    cout << "main:: w = " << w << endl;
#if _DIMENSION == 1
    Exact_Sol1D<0> will_abort_if_called;
#else
    Exact_Sol2D<0> will_abort_if_called;
#endif

    AbortBVP<dim> messy_bvp(&will_abort_if_called);
    TensorEquation<Basis1d,dim,Basis> unused_problem (&messy_bvp, bc, false); // 21 mb speicher!!!!, jmax=14 => 43 mb
    unused_problem.set_jmax(spatial_jmax,false);

    cout << "calling constructor AffLinParEq_Precomputed" << endl;

    tstart = clock();
    cout << "main :: insert dummy or true_forward_solution_filename into constructor!" << endl;
    abort();
    AffLinParEq_Precomputed<TensorEquation<Basis1d,dim,Basis>, number_of_timesteps > parabolic1
                            (&unused_problem,
                            u0T,
                            diff_coeff,
                            w,
                            time_discretization,
                            f_filename,
            //f_filename,
            //f_filename,
                            haar_gramian_filename_start,
                            haar_gramian_filename_end,
                            laplacian_filename,
                            par_jmax,
                            spatial_jmax,
                            tolerance
                            );
    /*
    cout << "dateinamen::" << endl << "f_filename = " << endl << f_filename << endl
            << "haar_gramian_filename_start =" << endl << haar_gramian_filename_start << endl
            << "haar_gramian_filename_end =" << endl << haar_gramian_filename_end << endl
            << "laplacian_filename =" << endl << laplacian_filename << endl;
    */

    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "main:: initialization of parabolic completed. time needed = " << time << " seconds" << endl;
#if _TEST_INITIALIZATION == 1 // test initialisation of parabolic
    // test stored variables:
    // protected variables have been set to public! cached_tproblem.h, aff_lin_par_eq.h
    cout << "Test initialization" << endl;
    tstart = clock();
    cout << "parabolic->identity_->entries_cache" << endl << (parabolic1.identity_->entries_cache) << endl;
    cout << "u0T = " << u0T << endl;
    cout << "forward_solution_ " << endl << parabolic1.forward_solution_ << endl;
    cout << "forward_solution_preprocessed_ " << endl << parabolic1.forward_solution_preprocessed_ << endl;

    InfiniteVector<double, Index> temp_u0;
    unused_problem.basis().expand(&uexact, false, spatial_jmax, temp_u0);

    temp_u0.compress(1e-13);
    cout << "expand(&uexact, false) = " << endl << temp_u0 << endl;

    cout << "backward_solution_ " << endl << parabolic1.backward_solution_ << endl;
    cout << "f_" << endl << parabolic1.f_ << endl;
    cout << "d_" << endl << parabolic1.d_ << endl;
    cout << "time_discretization_" << endl << parabolic1.time_discretization_ << endl;

    cout << "time_direction_" << endl << parabolic1.time_direction_ << endl;
    cout << "current_timestep_" << endl << parabolic1.current_timestep_ << endl;
    cout << "raw_cache_" << endl << parabolic1.raw_cache_ << endl;
    cout << "scaled_laplacian_matrix_" << endl << parabolic1.scaled_laplacian_matrix_ << endl;
    cout << "haar_gramian_filename_start_" << endl << parabolic1.haar_gramian_filename_start_ << endl;
    cout << "haar_gramian_filename_end_" << endl << parabolic1.haar_gramian_filename_end_ << endl;
    cout << "spatial_jmax_" << endl << parabolic1.spatial_jmax_ << endl;
    cout << "assembled_problems_" << endl;
    for (unsigned int i = 0; i<= number_of_timesteps; ++i)
    {
        cout << "   i = " << i << " *parabolic.assembled_problems_[i].entries_cache " << endl << *parabolic1.assembled_problems_[i].entries_cache << endl;
    }

    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "test initialization finished. time needed = " << time << "seconds" << endl;
#endif

//#if _COMPARE_WITH_LINPAREQ == 1

#if _TEST_LIN_PAR_EQ_TENSOR == 1
    //cout << "Compare aff_lin_par_eq and lin_par_eq" << endl;

    tstart = clock();

    ConstantFunction<dim> one_half_function(Vector<double>(1, "0.5"));

    IdentityBVP<dim> identity(&one_half_function); // rechte Seite muss gesetzt werden, um dualen coeff vektor berechnen zu koennen
    double eta = 1;

    PoissonBVP<dim> elliptic(&one_half_function); // -u''
    //PertubedBVP<dim>elliptic(&one_half_function,1); // -u''+u

    TensorEquation<Basis1d,dim,Basis> eq(&elliptic, bc, false);
    TensorEquation<Basis1d,dim,Basis> gram(&identity, bc, false);
    //jmax = multi_degree(eq.basis_.j0())+offset;
    eq.set_jmax(spatial_jmax, false);
    gram.set_jmax(spatial_jmax, false);

    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctproblem2(&eq,5.5,35.0);
    CachedTProblem<TensorEquation<Basis1d,dim,Basis> > ctgramian2(&gram,5.5,35.0);

    InfiniteVector<double,Index> u0T2;
    uexact.set_time(0);
    ctgramian2.set_f(&uexact);
    CDD1_SOLVE(ctgramian2, tolerance, u0T2, spatial_jmax, tensor_simple);
    cout << "main:: u0T2 after CDD1_SOLVE :: " << endl << u0T2 << endl;
    u0T2.scale(&ctgramian2,-1);
    u0T2.compress(1e-14);

#if 1
    // ctgramian2 and ctgramian seem to be quite different!!!!
    // compare them:

    cout << "u0T = " << endl << u0T << endl;
    cout << "u0T2 = " << endl << u0T2 << endl;

    set<Index> Lambda1, Lambda2;
    for (Index lambda ( ctgramian.basis().first_generator() ), itend(ctgramian.basis().last_wavelet(spatial_jmax));; ++lambda)
    {
        Lambda1.insert(lambda);
        if (lambda == itend) break;
    }
    for (Index lambda ( ctgramian2.basis().first_generator() ), itend(ctgramian2.basis().last_wavelet(spatial_jmax));; ++lambda)
    {
        Lambda2.insert(lambda);
        if (lambda == itend) break;
    }
    // "cout << "Lambda1 = " << endl << Lambda1 << endl << "Lambda2 = " << endl << Lambda2 << endl;"
    cout << "Lambda1" << endl;
    for (set<Index>::iterator it(Lambda1.begin()), itend(Lambda1.end()); it != itend; ++it)
    {
        cout << *it << " ";
    }
    cout << endl << "Lambda2" << endl;
    for (set<Index>::iterator it(Lambda1.begin()), itend(Lambda1.end()); it != itend; ++it)
    {
        cout << *it << " ";
    }
    cout << endl;
    cout << "compare ctgramian and ctgramian2:: check bilinearform" << endl;
    for (set<Index>::iterator it1(Lambda1.begin()), itend(Lambda1.end()); it1 != itend; ++it1)
    {
        for (set<Index>::iterator it2(Lambda1.begin()); it2 != itend; ++it2)
        {
            if (abs(ctgramian.a(*it1,*it2) - ctgramian2.a(*it1,*it2)) > 1e-14)
            {
                cout << "difference detected: diff = " << abs(ctgramian.a(*it1,*it2) - ctgramian2.a(*it1,*it2)) << " ctgramian.a(" << *it1 << ", " << *it2 << ") = " << ctgramian.a(*it1,*it2) << " ctgramian2.a(" << *it1 << ", " << *it2 << ") = " << ctgramian2.a(*it1,*it2) << endl;
            }
        }
    }
    cout << "done." << endl;
    cout << "compare ctgramian and ctgramian2:: check f" << endl;
    for (set<Index>::iterator it1(Lambda1.begin()), itend(Lambda1.end()); it1 != itend; ++it1)
    {
        if (abs(ctgramian.f(*it1) - ctgramian2.f(*it1)) > 1e-14)
        {
            cout << "difference detected: diff = " << abs(ctgramian.f(*it1) - ctgramian2.f(*it1)) << " ctgramian.f(" << *it1 << ") = " << ctgramian.f(*it1) << " ctgramian2.f(" << *it1 << ") = " << ctgramian2.f(*it1) << endl;
        }
    }
    cout << "done" << endl;
     cout << "compare ctgramian and ctgramian2:: check D" << endl;
    for (set<Index>::iterator it1(Lambda1.begin()), itend(Lambda1.end()); it1 != itend; ++it1)
    {
        if (abs(ctgramian.D(*it1) - ctgramian2.D(*it1)) > 1e-14)
        {
            cout << "difference detected: diff = " << abs(ctgramian.D(*it1) - ctgramian2.D(*it1)) << endl;
            cout << "   ctgramian.D(" << *it1 << ") = " << ctgramian.D(*it1) << endl;
            cout << "   ctgramian2.D(" << *it1 << ") = " << ctgramian2.D(*it1) << endl;
            //cout << " compare with abs(ctgramian.a("<< *it1 << ", " << *it1 << ") - ctgramian2.a(" << *it1 << ", " << *it1 << ")) = " << abs(ctgramian.a(*it1,*it1) - ctgramian2.a(*it1,*it1));
            cout << "   ctgramian.a(" << *it1 << ", " << *it1 << ") = " << ctgramian.a(*it1,*it1) << endl;
            cout << "   ctgramian2.a(" << *it1 << ", " << *it1 << ") = " << ctgramian2.a(*it1,*it1) << endl;
        }
    }
    cout << "done" << endl;
    cout << "compare ctgramian and ctgramian2:: RHS" << endl;
    InfiniteVector<double,Index> temp_rhs1, temp_rhs2;
    ctgramian.RHS(1e-5,temp_rhs1);
    ctgramian2.RHS(1e-5,temp_rhs2);
    temp_rhs1.subtract(temp_rhs2);
    if(l2_norm(temp_rhs1) > 1e-14)
    {
        cout << "difference detected: l2_norm(diff) = " << l2_norm(temp_rhs1);
        ctgramian.RHS(1e-5,temp_rhs1);
        cout << " ctgramian.RHS = " << temp_rhs1 << endl << "ctgramian2.RHS = " << temp_rhs2 << endl;
    }
    cout << "done" << endl;
    cout << "compare ctgramian and ctgramian2:: F_norm" << endl;
    if (abs(ctgramian.F_norm() - ctgramian2.F_norm()) > 1e-14)
    {
        cout << "difference detected: diff = " << abs(ctgramian.F_norm() - ctgramian2.F_norm()) << " ctgramian.F_norm = " << ctgramian.F_norm() << " ctgramian2.F_norm = " << ctgramian2.F_norm() << endl;
    }
    cout << "done" << endl;

#endif

#if 0 // really compare something...
    LinearParabolicEquationTensor<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > parabolic2(&ctproblem2, &ctgramian2, u0T2, &one_half_function, spatial_jmax);
#else // just run some tests with the old code
#if _DIMENSION == 1
    TestRHS1D<_PROBLEM_NO> driving_f;
#else
    TestRHS2D<_PROBLEM_NO> driving_f;
#endif
    LinearParabolicEquationTensor<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > parabolic2(&ctproblem2, &ctgramian2, u0T2, &driving_f, spatial_jmax);
#endif


// ###############################
    cout << "solving the parabolic problem with lin_par_eq_tensor" << endl;
    ROWMethod<InfiniteVector<double,Index> > method2(WMethod<InfiniteVector<double,Index> >::ROS2);
    method2.set_preprocessor(&parabolic2);
    OneStepScheme<InfiniteVector<double,Index> >* scheme2 = &method2;
    IVPSolution<InfiniteVector<double,Index> > results2;
    InfiniteVector<double,Index> temp2, result2, error_estimate2, temp_f2;

    //for (int expo = 6; expo <= 6; expo++)
    {
        int expo = _TIME_DISCRETIZATION_GRANULARITY;
        temp2 = parabolic2.u0;
        const int resolution = d+dim+((dim == 1)? (spatial_jmax-3):(spatial_jmax-6)) +1;
        std::ofstream resultstream;

        //SampledMapping<dim> u0_plot(evaluate(ctproblem2.basis(), temp2, true, resolution));
        SampledMapping<dim> u0_plot(evaluate(ctgramian.basis(), temp2, true, resolution));
        resultstream.open("u0.m");
        u0_plot.matlab_output(resultstream);
        resultstream.close();

        cout << "main:: calling parabolic2.eval_f with argument temp2 = " << endl << temp2 << endl;
        parabolic2.evaluate_f(0,temp2,1e-6,temp_f2);
        SampledMapping<dim> f0_plot(evaluate(ctproblem2.basis(), temp_f2, true, resolution));
        resultstream.open("f0.m");
        f0_plot.matlab_output(resultstream);
        resultstream.close();

        int N = 1<<expo;
        double h = 1.0/N;
        cout << "  h = " << h << ":" << endl;

        results2.t.push_back(0);
        results2.u.push_back(temp2);

        for (int i = 1; i <= N; i++)
        //for (int i = 1; i <= 0; i++)
        {
            // cout << "---------------- before increment() -----------------------" << endl;
            cout << "main:: increment with i = " << i << endl;
            method2.increment(&parabolic2, (i-1)*h, temp2, h, result2, error_estimate2, 1e-5);
            //abort();
            //scheme->increment(&parabolic2, (i-1)*h, temp, h, result, error_estimate, 1e-5);
            results2.t.push_back(i*h);
            results2.u.push_back(result2);

            ostringstream output_filename, output_filename2;
            output_filename << "u" << i << ".m";
            output_filename2 << "f" << i << ".m";

            resultstream.open(output_filename.str().c_str());
            SampledMapping<dim> ui_plot(evaluate(ctproblem2.basis(), result2, true, resolution));
            ui_plot.matlab_output(resultstream);
            resultstream.close();

            //eval_ft or parts thereof

            parabolic2.evaluate_f(i*h,result2,1e-6,temp_f2);

            resultstream.open(output_filename2.str().c_str());
            SampledMapping<dim> fi_plot(evaluate(ctproblem2.basis(), temp_f2, true, resolution));
            fi_plot.matlab_output(resultstream);
            resultstream.close();

//       cout << "---------------- after increment() -----------------------" << endl;

            InfiniteVector<double,Index> uexact_coeffs;
            uexact.set_time(i*h);
            ctproblem2.basis().expand(&uexact, false, spatial_jmax, uexact_coeffs);

            ctgramian2.set_f(&uexact);
            CDD1_SOLVE(ctgramian2, tolerance, uexact_coeffs, spatial_jmax, tensor_simple);
            uexact_coeffs.scale(&ctgramian2,-1);
            uexact_coeffs.compress(1e-14);

            cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result2 - uexact_coeffs) << endl;
            temp2 = result2;
        }
    }
    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "finished run with lin_par_eq_tensor. time needed = " << time << "seconds" << endl;

#endif




// ##########################
#if _TEST_AFF_LIN_PAR_EQ == 1
    tstart = clock();
    cout << "solving the parabolic problem with aff_lin_par_eq" << endl;
    ROWMethod<InfiniteVector<double,Index> > method1(WMethod<InfiniteVector<double,Index> >::ROS2);
    method1.set_preprocessor(&parabolic1);
    OneStepScheme<InfiniteVector<double,Index> >* scheme1 = &method1;
    IVPSolution<InfiniteVector<double,Index> > results1;
#if _PROBLEM_NO == 13
    cout << "heat equation. no need to change the already assembled W" << endl;
#else
#if _PROBLEM_NO == 14
    for (unsigned int i=0; i<= number_of_timesteps; ++i)
    {
        w[i].set_coefficient(0,1); // this adds the Gram matrix, resulting in A= -delta+Id
    }
    cout << "calling assemble_W" << endl;
    parabolic1.assemble_W(w);
#else
    //cout << "Problem number should be 13 or 14, but it is actually " << _PROBLEM_NO << ". aborting!";
    //abort();
#endif
#endif

    /*
    // the helper of aff_lin_par_eq makes problems!
    double alpha = 2;
    AffLinParEq_Precomputed_StageEquationHelper<CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> > > helper1(alpha, &parabolic1.assembled_problems_[0], parabolic1.identity_, u0T);
    LinParEqTenROWStageEquationHelper<CachedTProblem<TensorEquation<Basis1d,dim,Basis> > > helper2(alpha, parabolic2.elliptic, parabolic2.identity, u0T);

    InfiniteVector<double, Index> temp_iv;
    CDD1_SOLVE(helper2, tolerance, temp_iv, spatial_jmax);
    cout << "main:: helper2=>temp_iv = " << endl << temp_iv << endl;
    CDD1_SOLVE(helper1, tolerance, temp_iv, spatial_jmax);
    cout << "main:: helper1=>temp_iv = " << endl << temp_iv << endl;
*/

    cout << "time forward case. Problem No = " << _PROBLEM_NO << endl;
    InfiniteVector<double,Index> temp1, result1, error_estimate1, temp_f1;
    //for (int expo = 6; expo <= 6; expo++)
    {
        int expo = _TIME_DISCRETIZATION_GRANULARITY;
        temp1 = parabolic1.u0;
        const int resolution = d+dim+((dim == 1)? (spatial_jmax-3):(spatial_jmax-6)) +1;
        std::ofstream resultstream;


        SampledMapping<dim> u0_plot(evaluate(ctgramian.basis(), temp1, true, resolution));
        resultstream.open("u0.m");
        u0_plot.matlab_output(resultstream);
        resultstream.close();
        parabolic1.set_current_timestep(0);
        parabolic1.set_solution(temp1,tolerance);

        //cout << "main:: calling parabolic1.eval_f with argument temp1 = " << endl << temp1 << endl;
        parabolic1.evaluate_f(0,temp1,1e-6,temp_f1);
        //cout << "main:: result of eval_f = " << endl << temp_f1 << endl;
        //abort();
        SampledMapping<dim> f0_plot(evaluate(ctgramian.basis(), temp_f1, true, resolution));
        resultstream.open("f0.m");
        f0_plot.matlab_output(resultstream);
        resultstream.close();

        /*
        set <Index> temp_set;
        temp1.support(temp_set);
        cout << "test_aff_lin_par_eq: support of temp1" << endl;
        for ( set< CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> >::WaveletBasis::Index>::const_iterator it(temp_set.begin());
                it != temp_set.end(); ++it)
            cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
        temp_set.clear();

        temp_f1.support(temp_set);
        cout << "test_aff_lin_par_eq: support of temp_f1" << endl;
        for ( set< CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> >::WaveletBasis::Index>::const_iterator it(temp_set.begin());
                it != temp_set.end(); ++it)
            cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
        temp_set.clear();
        abort();
*/
        int N = 1<<expo;
        double h = 1.0/N;
        cout << "  h = " << h << ":" << endl;

        results1.t.push_back(0);
        results1.u.push_back(temp1);

        for (int i = 1; i <= N; i++)
        //for (int i = 1; i <= 0; i++)
        {
            // cout << "---------------- before increment() -----------------------" << endl;
            cout << "main:: increment with i = " << i << endl;
/*
            cout << "debugging ..." << endl;
            cout << "eval_f for some arguments " << endl;
            parabolic1.evaluate_f(0,temp1,1e-6,temp_f1);
            parabolic1.evaluate_f(h,temp1,1e-6,temp_f1);
            parabolic1.evaluate_f(0.34*h,temp1,1e-6,temp_f1);

            //parabolic1.solve_ROW_stage_equation(t,infvec,alph,y,tol,res)
            InfiniteVector<double,Index> temp_res;
            cout << "die rechte seite des helper problems macht das problem!!!" << endl;
            cout << "geht nicht: temp f1 = " << endl << temp_f1 << endl;
            cout << "geht wohl: u0T = " << endl << u0T << endl;
            AffLinParEq_Precomputed_StageEquationHelper<CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> > >  helperXYC(0.1, &parabolic1.assembled_problems_[parabolic1.current_timestep_], parabolic1.identity_, u0T);
            CDD1_SOLVE(helperXYC, tolerance, temp_res, spatial_jmax); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
            cout << "#############################" << endl;
            cout << "#############################" << endl;
            cout << "#############################" << endl;
            cout << "#############################" << endl;
            cout << "#############################" << endl << endl << endl << endl;

            SparseMatrix<double> A_Lambda;

            set<Index> Problem_Set;
            cout << "inserting " << *ctgramian.basis().get_wavelet(1) << endl;
            Problem_Set.insert(ctgramian.basis().get_wavelet(1));
            cout << "inserting " << *ctgramian.basis().get_wavelet(2) << endl;
            Problem_Set.insert(ctgramian.basis().get_wavelet(2));
            cout << "inserting " << *ctgramian.basis().get_wavelet(7) << endl;
            Problem_Set.insert(ctgramian.basis().get_wavelet(7));
            cout << "inserting " << *ctgramian.basis().get_wavelet(8) << endl;
            Problem_Set.insert(ctgramian.basis().get_wavelet(8));
            AffLinParEq_Precomputed_StageEquationHelper<CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> > >  helperXYD(0.1, &parabolic1.assembled_problems_[parabolic1.current_timestep_], parabolic1.identity_, temp_f1);
            setup_stiffness_matrix(helperXYD, Problem_Set, A_Lambda);
            cout << "main:: A_problem = " << endl << A_Lambda << endl;

            InfiniteVector<double,Index> F;
            helperXYD.RHS(1e-3, F);
            //set< CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> >::WaveletBasis::Index> temp_set;
            //F.support(temp_set);
            //cout << "test_aff_lin_par_eq: support of P.RHS("<<1e-3<<")" << endl;
            //for ( set< CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> >::WaveletBasis::Index>::const_iterator it(temp_set.begin());
            //        it != temp_set.end(); ++it)
            //    cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
 */
/*
 * set <Index> temp_set;
            helperXYC.y.support(temp_set);
            cout << "test_aff_lin_par_eq: support of y" << endl;
            for ( set< CompressedProblemFromMatrix<TensorEquation<Basis1d,dim,Basis> >::WaveletBasis::Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            */


/*
            CDD1_SOLVE(helperXYD, tolerance, temp_res, spatial_jmax); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
            temp_res.scale(&helperXYC, -1); // Dx -> x
            parabolic1.solve_ROW_stage_equation((i-1)*h,parabolic1.forward_solution_[0],0.1,temp_f1,tolerance,temp_res);
*/

            //cout << "time forward: pre increment: i = " << i << " temp1 = " << endl << temp1 << endl << "result1 = " << endl << result1 << endl;
            method1.increment(&parabolic1, (i-1)*h, temp1, h, result1, error_estimate1, 1e-5);
            //cout << "time forward: post increment: i = " << i << " temp1 = " << endl << temp1 << endl << "result1 = " << endl << result1 << endl;
            //scheme->increment(&parabolic2, (i-1)*h, temp, h, result, error_estimate, 1e-5);
            results1.t.push_back(i*h);
            results1.u.push_back(result1);

            ostringstream output_filename, output_filename2;
            //output_filename << "alpe_u" << i << ".m";
            //output_filename2 << "alpe_f" << i << ".m";

            output_filename << "u" << i << ".m";
            output_filename2 << "f" << i << ".m";

            resultstream.open(output_filename.str().c_str());
            SampledMapping<dim> ui_plot(evaluate(ctgramian.basis(), result1, true, resolution));
            ui_plot.matlab_output(resultstream);
            resultstream.close();

            parabolic1.evaluate_f(i*h,result1,1e-6,temp_f1);

            resultstream.open(output_filename2.str().c_str());
            SampledMapping<dim> fi_plot(evaluate(ctgramian.basis(), temp_f1, true, resolution));
            fi_plot.matlab_output(resultstream);
            resultstream.close();

//       cout << "---------------- after increment() -----------------------" << endl;

            parabolic1.set_current_timestep(i);
            parabolic1.set_solution(result1,tolerance);
            InfiniteVector<double,Index> uexact_coeffs;
            uexact.set_time(i*h);

            ctgramian.basis().expand(&uexact, false, spatial_jmax, uexact_coeffs);

            ctgramian.set_f(&uexact);
            CDD1_SOLVE(ctgramian, tolerance, uexact_coeffs, spatial_jmax, tensor_simple);
            uexact_coeffs.scale(&ctgramian,-1);
            uexact_coeffs.compress(1e-14);

            cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result1 - uexact_coeffs) << endl;
            temp1 = result1;
        }
    }
    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "finished aff_lin_par_eq time forward run. time needed = " << time << "seconds" << endl;

    cout << "begin aff_lin_par_eq time backwards run." << endl;
    // time still counts from 0 to 1. The change of direction is handled by aff_lin_par_eq
    IVPSolution<InfiniteVector<double,Index> > results_back;
#if _DIMENSION == 1
    Exact_Sol1D<N+2> uexact_b;
#else
    Exact_Sol2D<N+2> uexact_b;
#endif

    // time backward case
    parabolic1.set_time_direction(false);
    {
        int expo = _TIME_DISCRETIZATION_GRANULARITY;
        const int resolution = d+dim+((dim == 1)? (spatial_jmax-3):(spatial_jmax-6)) +1;

        std::ofstream resultstream;

        temp1.clear();
        SampledMapping<dim> v0_plot(evaluate(ctgramian.basis(), temp1, true, resolution));
        resultstream.open("v0.m");
        v0_plot.matlab_output(resultstream);
        resultstream.close();
        parabolic1.set_current_timestep(0); //0 means the last timestep!
        parabolic1.set_solution(temp1,tolerance);

        cout << "main:: calling parabolic1.eval_f with argument temp1 = " << endl << temp1 << endl;
        parabolic1.evaluate_f(0,temp1,1e-6,temp_f1);
        SampledMapping<dim> g0_plot(evaluate(ctgramian.basis(), temp_f1, true, resolution));
        resultstream.open("g0.m");
        g0_plot.matlab_output(resultstream);
        resultstream.close();

        int N = 1<<expo;
        double h = 1.0/N;
        cout << "  h = " << h << ":" << endl;

        results_back.t.push_back(0);
        results_back.u.push_back(temp1);

        for (int i = 1; i <= N; i++)
        //for (int i = 1; i <= 0; i++)
        {
            // cout << "---------------- before increment() -----------------------" << endl;
            cout << "main:: increment with i = " << i << endl;
            method1.increment(&parabolic1, (i-1)*h, temp1, h, result1, error_estimate1, 1e-5);
            //scheme->increment(&parabolic2, (i-1)*h, temp, h, result, error_estimate, 1e-5);
            results_back.t.push_back(i*h);
            results_back.u.push_back(result1);

            ostringstream output_filename, output_filename2;
            //output_filename << "alpe_v" << i << ".m";
            //output_filename2 << "alpe_g" << i << ".m";

            output_filename << "v" << i << ".m";
            output_filename2 << "g" << i << ".m";

            resultstream.open(output_filename.str().c_str());
            SampledMapping<dim> vi_plot(evaluate(ctgramian.basis(), result1, true, resolution));
            vi_plot.matlab_output(resultstream);
            resultstream.close();

            parabolic1.evaluate_f(i*h,result1,1e-6,temp_f1);

            resultstream.open(output_filename2.str().c_str());
            SampledMapping<dim> gi_plot(evaluate(ctgramian.basis(), temp_f1, true, resolution));
            gi_plot.matlab_output(resultstream);
            resultstream.close();

//       cout << "---------------- after increment() -----------------------" << endl;

            parabolic1.set_current_timestep(i);
            parabolic1.set_solution(result1,tolerance);
            InfiniteVector<double,Index> uexact_b_coeffs;
            uexact_b.set_time(i*h);

            ctgramian.basis().expand(&uexact, false, spatial_jmax, uexact_b_coeffs);

            ctgramian.set_f(&uexact_b);
            CDD1_SOLVE(ctgramian, tolerance, uexact_b_coeffs, spatial_jmax, tensor_simple);
            uexact_b_coeffs.scale(&ctgramian,-1);
            uexact_b_coeffs.compress(1e-14);

            cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result1 - uexact_b_coeffs) << endl;
            temp1 = result1;
        }
    }
    clock_t tend2 = clock();
    time = (double)(tend2-tend)/CLOCKS_PER_SEC;
    cout << "finished aff_lin_par_eq time backward run. time needed = " << time << "seconds" << endl;
    time = (double)(tend2-tstart)/CLOCKS_PER_SEC;
    cout << "aff_lin_par_eq time needed in total = " << time << "seconds" << endl;
#endif
// ##########################

/*
    char transition_matrix_filename[250];
    sprintf(transition_matrix_filename, "/import/shared/friedrich/source/precomputed/transition_matrices/unpreconditioned/transition_haar_jmax_%d_primbs_d_dt_3_3_bc_%s_jmax_%d",par_jmax,bc_string,spatial_jmax);
    cout << "main:: loading base transition from file " << endl << transition_matrix_filename << endl;
    SparseMatrix<double> transition_matrix(1);
    transition_matrix.matlab_input(transition_matrix_filename);
    Vector<double> u_coeffs, v_coeffs;
    u_coeffs.resize(ctgramian.basis().degrees_of_freedom(),true);
    v_coeffs.resize(ctgramian.basis().degrees_of_freedom(),true);
    for (unsigned int i=0; i<= number_of_timesteps; ++i)
    {
        for (InfiniteVector<double,Index>::const_iterator it(parabolic1.forward_solution_[i].begin()), itend(parabolic1.forward_solution_[i].end());it!=itend;++it)
        {
            //it.index().number();
            u_coeffs[it.index().number()]=*it;
        }
        for (InfiniteVector<double,Index>::const_iterator it(parabolic1.backward_solution_[i].begin()), itend(parabolic1.backward_solution_[i].end());it!=itend;++it)
        {
            //it.index().number();
            v_coeffs[it.index().number()]=*it;
        }
    }
    Vector<double> u_haar_gen_coeffs, v_haar_gen_coeffs;
    u_haar_gen_coeffs.resize((1<<(par_jmax+1)),true);
    v_haar_gen_coeffs.resize((1<<(par_jmax+1)),true);
 */

#if _DIMENSION == 1
    FixedVector<double, (1<<(par_jmax+1))> u_haar_gen_coeffs, v_haar_gen_coeffs;
    Array1D<FixedVector<double, (1<<(par_jmax+1))> > product_over_time;

    //cout << "u_haar_gen_coeffs[0] = " << u_haar_gen_coeffs[0] << endl;
    //cout << "u_haar_gen_coeffs[1] = " << u_haar_gen_coeffs[1] << endl;
    //cout << "u_haar_gen_coeffs[1] = " << u_haar_gen_coeffs[1] << endl;
    //cout << "u_haar_gen_coeffs[1] = " << u_haar_gen_coeffs[1] << endl;
#else
    FixedMatrix<double, (1<<(par_jmax+1)), (1<<(par_jmax+1))> u_haar_gen_coeffs, v_haar_gen_coeffs;
    Array1D<FixedMatrix<double, (1<<(par_jmax+1)), (1<<(par_jmax+1))> > product_over_time;
#endif
    product_over_time.resize(number_of_timesteps+1);

    // at each timestep: multiply generator coefficients => coefficients of the product
    // dont forget the weights at the Haar generators (they do evaluate to 0 or 1/L_2norm = 2^((level*dim)/2)
    for (unsigned int i=0; i<= number_of_timesteps; ++i)
    {
        /*
        dummy_function1(&ctgramian.basis());
        dummy_function2(parabolic1.forward_solution_[i]);
        dummy_function12(&ctgramian.basis(), parabolic1.forward_solution_[i]);
        dummy_function3(u_haar_gen_coeffs);

        dummy_function4<Basis, 64, 2>(&ctgramian.basis(), parabolic1.forward_solution_[i], u_haar_gen_coeffs);
*/
        //cout << "parabolic1.forward_solution_[" << i << "] = " << endl << parabolic1.forward_solution_[i] << endl;
        precise_evaluate<Basis, (1<< (_PAR_JMAX+1)), _DIMENSION>(&ctgramian.basis(),parabolic1.forward_solution_[i], u_haar_gen_coeffs);
        //cout << "parabolic1.backward_solution_[" << (number_of_timesteps - i) << "] = " << endl << parabolic1.backward_solution_[number_of_timesteps - i ] << endl;
        precise_evaluate<Basis, (1<< (_PAR_JMAX+1)), _DIMENSION>(&ctgramian.basis(),parabolic1.backward_solution_[number_of_timesteps - i ], v_haar_gen_coeffs);
#if _DIMENSION == 1
        for (unsigned int j = 0; j<(1<<(par_jmax+1)) ;++j)
        {
            product_over_time[i][j] = u_haar_gen_coeffs[j] * v_haar_gen_coeffs[j] * twotothejhalf(par_jmax+1);
        }
#else
        for (unsigned int j = 0; j<(1<<(par_jmax+1)) ;++j)
        {
            for (unsigned int k = 0; k<(1<<(par_jmax+1)) ;++k)
            {
                product_over_time[i].set_entry(j,k, u_haar_gen_coeffs.get_entry(j,k) * v_haar_gen_coeffs.get_entry(j,k) * (1<<(par_jmax+1)) );
            }
        }
#endif
        cout << "main:: plotting u_haar_gen_coeffs" << endl;
        plot_haar_gen_coeffs(u_haar_gen_coeffs,par_jmax+1);
        cout << "main:: plotting v_haar_gen_coeffs" << endl;
        plot_haar_gen_coeffs(v_haar_gen_coeffs,par_jmax+1);
        cout << "main:: plotting product_over_time[" << i << "]" << endl;
        plot_haar_gen_coeffs(product_over_time[i],par_jmax+1);
        // inverse Haar wavelet transform
        // compute sums and differences for the (1<< (par_jmax+1) ) many generators
#if _DIMENSION == 1
        FixedVector<double, (1<<(par_jmax+1))> temp_fv;
        const double sqrt2(sqrt(2));
        for (int j=par_jmax; j>=0; --j)
        {
            for (unsigned int k=0; k< (1<< j); ++k)
            {
                temp_fv[k]        = (product_over_time[i][2*k] + product_over_time[i][2*k+1])/sqrt2;
                temp_fv[(1<<j)+k] = (product_over_time[i][2*k] - product_over_time[i][2*k+1])/sqrt2;
            }
        }
        product_over_time[i] = temp_fv;
#else
        FixedMatrix<double, (1<<(par_jmax+1)), (1<<(par_jmax+1))> temp_fm;
        const double sqrt2(sqrt(2));
        for (int j=par_jmax; j>=0; --j)
        {
            // transform rows
            for (unsigned int k=0; k< (1<< j); ++k)
            {
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    temp_fm.set_entry(k, l,         (product_over_time[i].get_entry(2*k, l) + product_over_time[i].get_entry(2*k+1, l))/sqrt2);
                    temp_fm.set_entry( (1<<j)+k, l, (product_over_time[i].get_entry(2*k, l) - product_over_time[i].get_entry(2*k+1, l))/sqrt2);
                }
            }
            // transform columns
            for (unsigned int k=0; k< (1<< j); ++k)
            {
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    temp_fm.set_entry(k, l,         (product_over_time[i].get_entry(l, 2*k) + product_over_time[i].get_entry(l, 2*k+1))/sqrt2);
                    temp_fm.set_entry( (1<<j)+k, l, (product_over_time[i].get_entry(l, 2*k) - product_over_time[i].get_entry(l, 2*k+1))/sqrt2);
                }
            }
        }
        product_over_time[i] = temp_fm;
#endif
    }

#if 0 // test SparseMatrix
    /* Test SparseMatrix.add */
    SparseMatrix<double> A(3,4);
    A.set_entry(0,0,1);    A.set_entry(0,1,0);    A.set_entry(0,2,0);    A.set_entry(0,3,0);
    A.set_entry(1,0,0);    A.set_entry(1,1,6);    A.set_entry(1,2,7);    A.set_entry(1,3,8);
    A.set_entry(2,0,4);    A.set_entry(2,1,3);    A.set_entry(2,2,2);    A.set_entry(2,3,1);
    cout << "A = " << endl;    A.print(cout,2,4);

    SparseMatrix<double> B(3,4);
    B.set_entry(0,0,10);    B.set_entry(0,1,10);    B.set_entry(0,2,10);    B.set_entry(0,3,10);
    B.set_entry(1,0,0);    B.set_entry(1,1,11);    B.set_entry(1,2,11);    B.set_entry(1,3,0);
    B.set_entry(2,0,12);    B.set_entry(2,1,12);    B.set_entry(2,2,12);    B.set_entry(2,3,12);
    cout << "B = " << endl;    B.print(cout,2,4);
    A.add(1,B);
    cout << "B = " << endl;    B.print(cout,2,4);
    cout << "A = " << endl;    A.print(cout,2,4);
#endif


} // end of main
// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Ulrich Friedrich                                                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PARABOLIC_TOOLS_H
#define _WAVELETTL_PARABOLIC_TOOLS_H

#include <algebra/fixed_vector.h>
#include <algebra/fixed_matrix.h>
#include <parabolic/aff_lin_par_eq.h>
#include <utils/MersenneTwister.h>
#include <iostream>

/*
 * The following makros are introduced in solve_PDE. 
 * Since this file contains useful routines that dont rely on these makros, the following lines provide some default values.
 */
#ifndef _HAAR_JMAX
#define _HAAR_JMAX 4
#endif

#ifndef _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0
#define _PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 0
#endif

#ifndef _FORWARD_PROBLEM_NO
#define _FORWARD_PROBLEM_NO 1
#endif

#ifndef _SECOND_SETUP_COUPLING_MATRIX_W
#define _SECOND_SETUP_COUPLING_MATRIX_W 0
#endif

#ifndef _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL
#define _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL 1
#endif

#ifndef _COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL
#define _COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL 1
#endif

using namespace MathTL;
namespace WaveletTL
{

    /*
     * plot Haar generator coefficients on a certain level to cout
     * dont forget: they are actually scaled with 2^(level*dim)/2)
     * 2D: GeneratorCoeff(i,j) = ith gen in x and jth gen in y direction !! 
     * A plot needs to switch axes and direction in y!
     * NOTE: MATLAB output would need something different: 
     * Entries in a matrix code the x direction from left to right and 
     * the y direction from top to bottom, i.e., A(i,j) resembles i pixels right and j pixels down
     */
    template <unsigned int NUMOFHAARGENERATORS>
    void
    plot_haar_gen_coeffs (FixedVector<double, NUMOFHAARGENERATORS> & plotme, unsigned int level);

    template <unsigned int NUMOFHAARGENERATORS>
    void
    plot_haar_gen_coeffs (FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> & plotme, unsigned int level);

    /*
     * integrate a function given by its wavelet coefficients against a Haar generator
     * \int_0^1 haargen_ij(x)u(x) dx
     * gen_index determines the used generator. same meaning as in INTEGRATE
     * basically a specialized version of integrate
     */
    template <class TBASIS, unsigned int DIM>
    double
    haar_gen_coeff(const TBASIS* basis,
                   const InfiniteVector<double, typename TBASIS::Index> & coeffs,
                   const unsigned int* gen_index,
                   const int haarmaxlevel);

    /*
     * integrate a wavelet against a Haar generator
     * \int_0^1 haargen_ij(x,y) psi_lambda(x,y) dxdy
     * haargen_ij (x,y)= eta_i(y)*eta_j(x) where gen_index[0]=j,gen_index[1]=i
     * (think of the generators arranged in a matrix)
     */
    template <class IBASIS, unsigned int DIM>
    double
    integrate(const WaveletTL::TensorBasis<IBASIS, DIM>* basis,
              const typename WaveletTL::TensorBasis<IBASIS, DIM>::Index lambda,
              const unsigned int * gen_index,
              const int level);

    /*
     * Integrate a function, given by a function pointer against a Haar generator
     */
    template <unsigned int DIM>
    double
    integrate(const Function<DIM,double>* fkt,
              const unsigned int primal_polynomial_degree,
              const unsigned int * gen_index,
              const int level);

    /*
     * Yields the Haar generator coefficients of a function that is given by its wavelet coefficients w.r.t. basis.
     *
     * background: we are interested in the haar-wavelet coefficients  of a function u given by wavelet coefficients coeffs.
     * by evaluating the function on a grid we get the Haar generator coefficients on some level (the dyadic resolution of the grid)
     * In a second step (not done by this method) the Haar wavelet transform yields the result we are really interested in.
     *
     * output = Haar generator coefficients of u
     * 2D: output(i,j) = ith gen in x and jth gen in y direction !! 
     * A plot needs to switch axes and direction in y! -> plot_haar_gen_coeffs
     */    
    template <class TBASIS, unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const TBASIS* basis,
                     const InfiniteVector<double, typename TBASIS::Index>& coeffs,
#if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result
#else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result
#endif
                    );

    template <class TBASIS, unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const TBASIS* basis,
                     const InfiniteVector<double, typename TBASIS::Index>& coeffs,
#if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result,
#else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result,
#endif
                     const unsigned int current_dim,
                     unsigned int* gen_index);
    
    /*
     * compute Haar generator coeficients of a function given by a Function<DIM,double> pointer
     */
    template <unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const Function<DIM,double>* fkt,
                     const unsigned int primal_polynomial_degree,
#if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result
#else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result
#endif
                    );

    template <unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const Function<DIM,double>* fkt,
                     const unsigned int primal_polynomial_degree,
#if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result,
#else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result,
#endif
                     const unsigned int current_dim,
                     unsigned int* gen_index);

    /*
     * Transform given Haar generator coefficients into Haar wavelet coefficients.
     * 1D and 2D case.
     *
     * 2D: generatormatrix(i,j) contains coeff corresponding to the ith gen in x direction \times jth generator in y direction
     * 
     * Example:
     * u = characteristic function on (0.375,0.5)x(0.125,0.25)
     *   = 0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 + 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     * has gen_coeffs
     *   = 0 0 0 0 0 0 0 0
     *     0 0 0 + 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     * has wav_coeffs
     *   = + + - 0 0 - 0 0
     *     + + 0 0 0 0 0 0
     *     + 0 - 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 - 0 0 0 + 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     *     0 0 0 0 0 0 0 0
     * 
     * Haar Wavelet coeffs are stored in a Matrix in the following way (for maxlevel = 2)
     *  0  1  4  5 16 17 18 19
     *  2  3  6  7 20 21 22 23
     *  8  9 12 13 24 25 26 27
     * 10 11 14 15 28 29 30 31
     * 32 33 34 35 48 49 50 51
     * 36 37 38 39 52 53 54 55
     * 40 41 42 43 56 57 58 59
     * 44 45 46 47 60 61 62 63
     */
    template < unsigned int NUMOFHAARGENERATORS>
    void haar_wavelet_transform(const FixedVector<double, NUMOFHAARGENERATORS>& gen_coeffs, FixedVector<double, NUMOFHAARGENERATORS>& wav_coeffs);

    template < unsigned int NUMOFHAARGENERATORS>
    void haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>& gen_coeffs, FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> & wav_coeffs);

    /*
     * Transform given Haar wavelet coefficients into Haar generator coefficients.
     * 1D and 2D case.
     */
    template < unsigned int NUMOFHAARGENERATORS>
    void inverse_haar_wavelet_transform(const FixedVector<double, NUMOFHAARGENERATORS>& wav_coeffs, FixedVector<double, NUMOFHAARGENERATORS>& gen_coeffs);

    template < unsigned int NUMOFHAARGENERATORS>
    void inverse_haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>&  wav_coeffs, FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> & gen_coeffs);
    
    /*
     * Transform Haar Wavelet coefficients from FixedVector/FixedMatrix to InfiniteVector format and vice versa
     */
    template<unsigned int SIZE>
    void transform_fixedToIV(const FixedVector<double, SIZE> & fixed, InfiniteVector<double, int> & iv);
    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void transform_fixedToIV(const FixedMatrix<double, ROWSIZE, COLUMNSIZE> & fixed, InfiniteVector<double, int> & iv);

    template<unsigned int SIZE>
    void transform_IVTofixed(const InfiniteVector<double, int> & iv, FixedVector<double, SIZE> & fixed);
    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void transform_IVTofixed(const InfiniteVector<double, int> & iv, FixedMatrix<double, ROWSIZE, COLUMNSIZE> & fixed);

    /* 
     * Use an AffLinParEq_Precomputed to solve the parabolic equation
     * assumes that
     * parabolic.forward_solution_[0] == parabolic.u0
     * parabolic.backward_solution_[0] = emptyInfiniteVector
     */
    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void solve_parabolic_problem(AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS> & parabolic,
                                  const ROWMethod<InfiniteVector<double, typename PROBLEM::Index> >& method,
                                  const bool time_direction,
                                  const double increment_tolerance,
                                  const double tolerance);
    
    /* 
     * Use an aff_lin_par_eq_qtbasis to solve the parabolic equation
     * assumes that
     * parabolic.forward_solution_[0] == parabolic.u0
     * parabolic.backward_solution_[0] = emptyInfiniteVector
     */
    
//    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT, class IVP>
//    void dummy(AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> & parabolic_problem);
//                                  const ROWMethod<InfiniteVector<double, int>, IVP >& method,
//                                  const bool time_direction,
//                                  const double increment_tolerance,
//                                  const double tolerance);
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT, class IVP>
    void solve_parabolic_problem(AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> & parabolic_problem,
                                  const ROWMethod<InfiniteVector<double, int>, IVP >& method,
                                  const bool time_direction,
                                  const double increment_tolerance,
                                  const double tolerance);

    /*
     * compute -h*u
     * where
     * where at each time step u,h are given by the wavelet coefficients forward_/backward_solution[N-t]
     * the minus comes from the fact that h(t)= (-) v(T-t). This has to be taken into account when using w_update!
     *
     * Sadly it seems too expensive to reuse the Gramian entries already present in the cache! Reason:
     * Strategy to multiply two functions f, g given by their coeffs (in 2 dims)
     *  f = \sum c_lambda \psi_lambda; g =  \sum d_mu \psi_mu
     *  Haar-Generator coeff of product:
     *  f*-coeff == \sum c_lam d_mu \int_support_of_generator \psi_lam\psi_mu
     * The integral is decomposabel into one dimensional integrals and the values of those are alrady present in the cache!
     * Nice! However for each of the possible pairs of c_lam's and d_mu's (complexity #c_lam * #d_mu ) we have to add the Haar-Gen-Coeffs of the integral \psi_lam\psi_mu,
     * e.g. a Matrix in 2 dims.
     * Complexity: #c_lam's * #d_mu's * ONEDIMHAARCOUNT^2
     * BAD!
     * 
     * Alternative:
     * Evaluate f,g seperately. Then generate Haar-Gen-Coeffs
     * Complexity:
     * Evaluate = (#c_lam's + #d_mu's)*ONEDIMHAARCOUNT^2
     * Compute Haar-Gen_coeffs of f*g, given individual Haar-Gen-Coeffs:
     * Complexity: ONEDIMHAARCOUNT^2
     * NICE!
     * 
     * output (in Haar wavelet coefficients) is stored in w_update
     */
    template <class TBASIS, unsigned int NUMOFHAARGENERATORS, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> > & backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update);
    
    
    
    /* for debugging a version that also returns haar_generator coeffs of u and h.
     * output:
     * dont get confused about where t=0 is stored. In this case we have:
     * w_update[i] = u_haar_gen[i]*h_haar_gen[i]
     *
     * apply changes to compute_update_w also here!*/
#if 0
    template <class TBASIS, unsigned int NUMOFHAARGENERATORS, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update,
            Array1D<InfiniteVector<double, int> > & u_haar_gen,
            Array1D<InfiniteVector<double, int> > & h_haar_gen);
#endif
    template <class TBASIS, unsigned int NUMOFHAARGENERATORS, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update,
            std::ofstream& plotstream);
    
    
    /*!
     * write *this to a binary file
     */
    /*
    template <class C, class I>
    void infinite_vector_writeToFile(const char *filename);
     */

    /*!
     * open a file written by writeToFile
     */
    /*
    template <class C, class I>
    void infinite_vector_readFromFile(const char *filename);
     */

    /*
     * Realize the soft shrinkage:
     * apply SHRINK (w + w_update, alpha) and store output in w_next
     * w = old parameter coeffs (w.r.t. Haar wavelet basis)
     * w_update = update coeffs (w.r.t. HWB)
     * SHRINK iterates over the argument coeff vector: coeffs with absolute value
     * below alpha are set to zero, others are decreased (or increased) by alpha
     * 
     * It may be possible to speed up the code using map::transform, although 
     * this may lead to entries with value 0
     */
    void shrinkage_iteration(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double alpha,
            Array1D<InfiniteVector<double, int> >& w_next);
    
    /*
     * if (plotstream != NULL) -> store output in the stream
     * w_update is computed internally.
     * To this end: u and h are sampled with granularity PRECISE_EVALUATE_GRANULARITY in each coordinate.
     * value should be higher than ONEDIMHAARCOUNT
     * 
     * FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & w_next:
     * w_next[time][patch]=sampled version of W, i.e., a matrix
     */
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void shrinkage_iteration(AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> & parabolic_problem,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & w,
            const double alpha,
            FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & w_next,
            std::ofstream * plotstream);

    /*
     * A possible stopping criterion for the inverse iteration.
     * returns true if
     * (sum_i=0^N (||w_update[i]||_2)^2) / sum_i=1^N (||w[i]||_2)^2 <= stopping_tolerance^2
     * where N = number_of_timesteps.
     * if w_update == 0 return true
     * else if case w == 0 return false
     * else normal behaviour
     *
     * ell1 variant:
     * (sum_i=0^N ||w_update[i]||_1) / (sum_i=1^N ||w[i]||_) <= stopping_tolerance
     */
    bool relative_change_criterion_ell2(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double stopping_tolerance,
            double relative_error);
    bool relative_change_criterion_ell1(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double stopping_tolerance,
            double relative_error);

    /*
     * Create a noise coefficient vector.
     * Entries are initialized with normal distributed numbers with zero mean and standard deviation sigma.
     * After initialization the vector is renormalized to the norm delta
     * 
     * the result is densely populated. Non the less we une InfiniteVector. This is not optimal!
     * My reason:
     * Parabolic solver is child of AbstractIVP<InfiniteVector<double, ...>
     * and the computed noise is added to the initial value, i.e., AbstractIVP::u0
     * 
     * More efficient would be: noise as a Vector or FixedVector of length degrees_of_freedom.
     * + same data structure in the AbstractIVP
     * + compatibility with the adaptive code, e.g., CDD1_SOLVE
     * 
     * First noise method: AffLinParEq_Precomputed (with TBasis and precomputed 2d matrices)
     * Second noise method: AffLinParEq_qtbasis (with QTBasis; no precomputations, only onedim integrals are cached)
     */
    template <class TBASIS>
    void create_noise(const double sigma,
                      const double delta,
                      const TBASIS* basis,
                      Array1D<InfiniteVector<double, typename TBASIS::Index> > &noise_coeffs);
    
    template <class QTBASIS, unsigned int NUMOFTIMESTEPSPLUSONE>
    void create_noise(const double sigma,
                      const double delta,
                      const QTBASIS* basis,
                      FixedArray1D<InfiniteVector<double, int>, NUMOFTIMESTEPSPLUSONE > &noise_coeffs);
    
    //template <unsigned int haar_jmax>
    void initialize_coupling_matrix(Array1D<InfiniteVector<double,int> >& wtrue,
                                    const int first_setup_coupling_matrix_w,
                                    const int d);
    
    template < unsigned int ONEDIMHAARCOUNT, unsigned int NUMOFTIMESTEPSPLUSONE >
    void initialize_coupling_matrix(FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPSPLUSONE> &w_true,
                                    const int setup_number);

    /*
     * Write the first lines of the logfile. Explainations of variables, parameter values, initialization with 0 vectors of correct length where needed.
     */
    void initialize_logstream(std::ofstream& logstream, 
            const int compute_noise, 
            const unsigned int first_setup_coupling_matrix_w, 
            const unsigned int max_inverse_iterations, 
            const double stopping_tolerance, 
            const unsigned int spatial_jmax, 
            const unsigned int w_limit_number, 
            const double alpha, 
            const double delta, 
            const double utrue_coeffs_norm, 
            const double uobserverd_coeffs_norm,
            const unsigned int number_of_timesteps);
    
    /*
     * write uexact to files (Matlab output and coefficients)
     * code moved from main to increase readability
     */
    template <class PROBLEM, class CTPROBLEM, unsigned int NUMOFTIMESTEPS>
    void store_uexact(const AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS> &parabolic, 
            const CTPROBLEM& ctgramian, 
            const char *true_forward_solution_filename_start, 
            const char *true_forward_solution_filename_coeffs_end, 
            const char *true_forward_solution_filename_matlab_end, 
            const int resolution);
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void store_uexact(const AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> &parabolic_problem, 
            const char* true_forward_solution_filename_start, 
            const char* true_forward_solution_filename_coeffs_end, 
            const char* true_forward_solution_filename_matlab_end, 
            const int resolution);
    
    /*
     * update_logstream is called in every iteration of the inverse solver.
     * It writes information about computation times and errors, vector norms to the provided filestream.
     * It does not open or close the filestream.
     */
    /*
    template <class TBASIS, unsigned int NUMBER_OF_TIMESTEPS>
    void update_logstream(std::ofstream& logstream, 
            const unsigned int iteration_count, 
            const double total_time, 
            const double assemble_time, 
            const double forward_time, 
            const double backward_time, 
            const double update_time, 
            const double shrinkage_time, 
            const Array1D<InfiniteVector<double, int> >& w, 
            const Array1D<InfiniteVector<double, int> >& wlimit, 
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution, 
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& utrue_coeffs);
    */
    /*
    template <class TBASIS, unsigned int NUMBER_OF_TIMESTEPS>
    void log_times(std::ofstream& logstream, 
            const unsigned int iteration_count, 
            const double total_time, 
            const double assemble_time, 
            const double forward_time, 
            const double backward_time, 
            const double update_time, 
            const double shrinkage_time);
    */
    
    /* 
     * This method computes norms for the two input vectors simultaneous. 
     * This should be faster than calling the individual vector norm routines.
     * 
     * output:
     * diff1 = l_1-norm of x-y
     * diff2 = l_2-norm of x-y
     * x1 = l_1-norm of x
     * x2 = l_2-norm of x
     */
    template <class C, class I>
    void compute_norms(const InfiniteVector<C,I> &x,const InfiniteVector<C,I> &y, double& diff1, double& diff2, double& x1, double& x2);
    
    
    
    template <unsigned int ONEDIMHAARCOUNT>
    void compute_norms(const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > &x,
            const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > &y, 
            double& diff1, 
            double& diff2, 
            double& x1, 
            double& x2);
    
    /* similar. only compute ell2 norms*/
    template <class C, class I>
    void compute_norms2(const InfiniteVector<C,I> &x,const InfiniteVector<C,I> &y, double& diff2, double& x2);
    
    /*
     * Compute norms related to the computed coupling coefficient W. These are
     * 0,1,2 norm of W. 1,2 norm of W-Wlimit.
     * Store output in logstream
     */
    template <unsigned int NUMBER_OF_TIMESTEPS>
    void log_reconstruction_errors(std::ofstream& logstream,
            const Array1D<InfiniteVector<double, int> >& w, 
            const Array1D<InfiniteVector<double, int> >& wlimit);
    
    template <unsigned int NUMBEROFTIMESTEPSPLUSONE, unsigned int ONEDIMHAARCOUNT >
    void log_reconstruction_errors(std::ofstream& logstream,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBEROFTIMESTEPSPLUSONE>& w, 
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBEROFTIMESTEPSPLUSONE>& wlimit);
    
    
    template <class TBASIS, unsigned int NUMBER_OF_TIMESTEPS>
    void log_solution_errors(std::ofstream& logstream,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution, 
            //const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            //const Array1D<InfiniteVector<double,typename TBASIS::Index> >& utrue_coeffs);
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& uobserved_coeffs);
    
    template <unsigned int NUMBER_OF_TIMESTEPS_PLUS_ONE>
    void log_solution_errors(std::ofstream& logstream,
            const FixedArray1D<InfiniteVector<double,int>, NUMBER_OF_TIMESTEPS_PLUS_ONE >& forward_solution, 
            const FixedArray1D<InfiniteVector<double,int>, NUMBER_OF_TIMESTEPS_PLUS_ONE >& uobserved_coeffs);
    
    
    /*
     * Plot additional plots to plotstream.
     * Does not open or close plotstream.
     * 
     * w, wtrue, wlimit
     */
    template<unsigned int NUMBER_OF_TIMESTEPS, unsigned int NUMOFWAVELETSPERCOORDINATE>
    void plot_solutions(std::ofstream& plotstream, 
            const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_true,
            const Array1D<InfiniteVector<double, int> >& w_limit);
    
    template<unsigned int NUMBER_OF_TIMESTEPS_PLUS_ONE, unsigned int ONEDIMHAARCOUNT>
    void plot_solutions(std::ofstream& plotstream, 
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w_true,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w_limit);
}


#include "parabolic_tools.cpp"

#endif

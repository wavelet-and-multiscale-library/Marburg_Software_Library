// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_AFF_LIN_PAR_EQ_H
#define _WAVELETTL_AFF_LIN_PAR_EQ_H

#include <algebra/infinite_vector.h>
#include <algebra/fixed_vector.h>
#include <numerics/ivp.h>
#include <numerics/w_method.h>
#include <utils/function.h>

#include <galerkin/cached_tproblem.h>
#include <galerkin/cached_qtproblem.h>

#include <galerkin/infinite_preconditioner.h>
#include <adaptive/apply.h>
#include <adaptive/apply_tensor.h>
#include <adaptive/cdd1.h>
#include <io/infinite_vector_io.h>

using MathTL::InfiniteVector;
using MathTL::AbstractIVP;
using MathTL::Function;
using MathTL::WMethodPreprocessRHSHelper;

namespace WaveletTL
{
    /*!
     * ROW stage equation helper class for linear parabolic equations (see below).
     * The helper class models the preconditioned stage equation
     *
     *       D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
     *
     * (after multiplication with the Gramian, see implementation)
     *
     * The class has the minimal signature to be used within the APPLY routine
     * or within adaptive solvers like CDD1.
    */
    template <class CACHEDTPROBLEM>
    class AffLinParEq_Precomputed_StageEquationHelper
    : public FullyDiagonalEnergyNormPreconditioner<typename CACHEDTPROBLEM::Index>
    {
        public:
        /*!
         * constructor from alpha, T, the Gramian and y
         */
        AffLinParEq_Precomputed_StageEquationHelper
                (const double alpha,
                const CACHEDTPROBLEM* T,
                const CACHEDTPROBLEM* G,
                const InfiniteVector<double, typename CACHEDTPROBLEM::Index>& y);

        /*!
         * make wavelet basis type accessible
         */
        typedef typename CACHEDTPROBLEM::WaveletBasis WaveletBasis;

        /*!
         * wavelet index class
         */
        typedef typename CACHEDTPROBLEM::Index Index;

        /*!
         * Type of wavelet level
         */
        typedef typename Index::level_type index_lt;

        /*!
         * read access to the basis
         */
        const WaveletBasis& basis() const { return T->basis(); }

        /*!
         * space dimension of the problem
         */
        static const int space_dimension = CACHEDTPROBLEM::space_dimension;

        /*!
         * locality of the operator
         */
        static bool local_operator()
        {
            return true; //CACHEDTPROBLEM::local_operator();
        }

        /*!
         * (half) order t of the operator
         */
        double operator_order() const
        {
            return T->operator_order();
        }

        /*!
         * evaluate the diagonal preconditioner D defined by T
         */
        double D(const Index& lambda) const
        {
            return sqrt(a(lambda, lambda));
        }

        /*!
         * evaluate the (unpreconditioned) bilinear form a
         */
        double a(const Index& lambda,
                 const Index& nu) const
        {
            return alpha * G->a(lambda, nu) + T->a(lambda, nu);
        }

        /*!
         * estimate the spectral norm ||alpha*D^{-2}-A|| from above
         */
        double norm_A() const
        {
            return T->norm_A(); // dirty!
        }

        /*!
         * estimate the spectral norm ||(alpha*D^{-2}-A)^{-1}|| from above
         */
        double norm_Ainv() const
        {
            return T->norm_Ainv(); // dirty!
        }

        /*!
         * estimate compressibility exponent s^*
         */
    //    double s_star() const {
    //      return T->s_star();
    //    }

        /*!
         * alpha_k is used differently for tensorbasis
         */
        double alphak(const unsigned int k) const {
            return T->alphak(10); // argument does not influence the value
        }

        /*!
         * evaluate the (unpreconditioned) right-hand side f
         */
        double f(const Index& lambda) const 
        {
            return y.get_coefficient(lambda);
        }

        /*!
         * approximate the wavelet coefficient set of the preconditioned right-hand side F
         * within a prescribed \ell_2 error tolerance
         */
        void RHS(const double eta,
                 InfiniteVector<double, Index>& coeffs) const
        {
            // dirty
            coeffs = y_scaled;
        }

        /*!
         * compute (or estimate) ||F||_2
         */
        double F_norm() const
        {
            return l2_norm(y_scaled);
        }

        /*
         * Called by add_compresed_column
         * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= range && ||nu|| <= maxlevel)
         * Hack? : works with non dynamic vector w of size = degrees_of_freedom
         *
         * "strategy" argument ist needed for compatibility with compression.h .
         * There is only 1 strategy for TBasis at the moment
         */
        void add_ball(const Index& lambda,
                     //InfiniteVector<double, Index>& w,
                     Vector<double>& w,
                     const int radius,
                     const double factor,
                     const int maxlevel,
                     const CompressionStrategy strategy = tensor_simple,
                     const bool precond = false) const;


        /*! Set the rhs vector y. compute y_scaled for RHS method*/
        void set_y(const InfiniteVector<double, Index> ynew);
        /*! Set the parameter alpha*/
        void set_alpha (const double alpha);

      //protected:
        double alpha;
        const CACHEDTPROBLEM* T;
        const CACHEDTPROBLEM* G;
        InfiniteVector<double, typename CACHEDTPROBLEM::Index> y;
        InfiniteVector<double, typename CACHEDTPROBLEM::Index> y_scaled;
    };
    
    /*!
     * This class will be used for a parameter
     * identification problem. Indeed a simplified version of the problem considered in project DA360/12-1.
     * A parabolic problem is considered in the usual time forward way and also in a time reversed way.
     * 
     * The parabolic equation considered is of the form
     *  u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
     *  u(0) = u_0
     * where -A: \ell_{2,D} -> \ell_{2,D^{-1}} is a positive isomorphism
     * and f: (0,T] -> \ell_{2,D^{-1}}, D being a diagonal preconditioner (see below).
     * An equation of this form arises when equivalently reformulating
     * a problem of the analogous form
     *  v'(t) = Bv(t) + g(t) =: G(t,v(t)), 0 < t <= T
     *  v(0) = v_0
     * with isomorphism -B: H -> H' and right-hand side g: (0,T] -> H'
     * by means of a biorthogonal wavelet basis Psi=\{\psi_\lambda\},
     * setting F(t,v):=<G(t,v^\top\Psi),\tilde\Psi> for all v in \ell_{2,D}.
     * 
     * The (unpreconditioned) stiffness matrix
     *  -A = (a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}
     * is assumed to be of the form
     *  A = (d <grad psi_lambda, grad psi_mu> + < w psi_lambda, psi_mu>)_{lambda,mu}
     * where w is arbitrary in time and space whereas d is assumed constant. w is given
     * by its coefficients with respect to the Haar wavelet basis.
     * 
     * The problem can be solved in a normal time forward way. In this case f=1/2 is used
     * as a source term. Or the problem can be solved with reversed time direction. 
     * In this case v(T)=0 and f=\hat v is used where \hat v is the solution of the 
     * time forward case that has been previously computed.
     * 
     * Precomputed matrices are loaded, but matrix-vector products are still done in an adaptive way
     * (ie. up to tolerance eta, and the compression strategy is as in the non precomputed case)
     * 
     * At each time step the operator matrix W has to be assembled. They are stored, since they are
     * used again for the time reversed case.
     * 
     * The solutions of the time forward and time backward parabolic problem are then used to update
     * the parameter w. To avoid unnecessary disk operations it may be beneficial to reuse already loaded matrices.
     * 
     * Possibly similar coefficients are active in the new w, thus the same component matrices
     * <h_eta psi_lambda, psi_mu> have to be used for the assembly of the new W.
     * Strategies for loading and storing the building matrices into the memory:
     * 
     * 1. On demand: load as many matrices as needed for the assembly of W(t_k) at time t_k. Then forget them.
     * pro: least memory consumption
     * con: highest number of disk operations
     * 2. Cached: as 1, but do not forget the matrices that have already been loaded
     * pro: less disk operations
     * con: higher memory consumption. especially bad, if many different (possibly all
     *   allowed) matrices have to be loaded as the inverse iteration goes.
     *   Worst case: all allowed matrices are loaded.
     *   This determines how large the matrices may be (memory restriction)
     * 3. All: load all matrices possible for possible assemblies of W into the memory
     * pro: fixed structure for matrix storage-> faster access to the matrices
     * con: if there are matrices that are never used in the inverse iteration, then this strategy wastes time and memory.
     * fact: same memory consumption as the worst case in the 2nd strategy.
     * 
     * check:
     * possibilities: store matrices externally or internally.
     * external storage: then we need a constructor for the helper class that works 
     * with references to matrices and the 3 strategies have to implemented elsewhere.
     * internal storage: possible to load the matrices by their filename. Memory strategy is dealt with here. I tend to this option atm ...
     * 
     * Matrices are assumed to be unpreconditioned.
     */
    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    class AffLinParEq_Precomputed
    : public AbstractIVP<InfiniteVector<double, typename PROBLEM::Index> >,
      public WMethodPreprocessRHSHelper<InfiniteVector<double, typename PROBLEM::Index> >
  
    {
    public:
    /*!
     * for convenience: the index class of the basis of the elliptic problem
     * eg PROBLEM == TensorEquation<Basis1d,dim,Basis>
     */
    typedef typename PROBLEM::Index Index;

    /*!
     * space dimension of the problem
     */
    static const int space_dimension = PROBLEM::space_dimension;

    /*!
     * Models the parabolic solver for the linearized parabolic problem.
     * Time forward and backward case are put together in this class in order to reuse as many matrices as possible.      
     * problem = underlying problem, eg TensorEquation<Basis1d,dim,Basis>, needed for the local CompressedProblems
     *           only the basis is needed. so we could replace "problem" with
     *           #if _DIMENSION == 1
     *               Exact_Sol1D<0> will_abort_if_called;
     *           #else
     *               Exact_Sol2D<0> will_abort_if_called;
     *           #endif
     *           AbortBVP<dim> messy_bvp(&will_abort_if_called);
     *           TensorEquation<Basis1d,dim,Basis> unused_problem (&messy_bvp, bc, false);
     *           unused_problem.set_jmax(spatial_jmax,false);
     * u0 = initial value. coeffs w.r.t. primal basis (<u,\tilde \psi_\lambda>), unpreconditioned
     * //f = rhs per time step. coeffs w.r.t. dual basis, unpreconditioned
     * d = diffusion coefficient
     * w = problem coefficients w.r.t. Haar basis
     * time_discretization = vector of time points
     * f_filename = filename of InfiniteVector corresponding to coeffs of f=1/2 w.r.t. dual basis, unpreconditioned, eg
     *            = "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_ffff.iv"
     * haar_gramian_filename_start _end = every char before and after the number eta in the names of the matrices, eg
     *      start = /import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_
     *      end   = _gramian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd.bin
     * laplacian_filename = name of the file containing the matrix of the laplacian, eg
     *            = /import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd.bin
     * //output_path = folder where output is written to,
     * par_jmax = maximal level for the Haar wavelet basis of w. Needed to determine size of matrix storage
     * spatial_jmax = maximal level for the spatial basis {psi_lambda}_lambda
     * tolerance - tolerance of APPLY used in transforming the coeffs of the forward solution into coeffs w.r.t. dual basis. The preprocessed coeffs are afterwards coarsened with tol=1e-13
     */
    AffLinParEq_Precomputed(PROBLEM * problem,
                            const InfiniteVector<double,Index>& u0,
                            //const FixedVector<InfiniteVector<double,Index>, NUMOFTIMESTEPS>& f,
                            const double d,
                            const Array1D<InfiniteVector < double, int > >& w,
                            //const FixedVector <InfiniteVector < double, int > , (NUMOFTIMESTEPS+1) >& w,
                            const Array1D<double>& time_discretization,
                            //const FixedVector<double, (NUMOFTIMESTEPS+1)>& time_discretization,
                            const char* f_filename = "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_ffff.iv",
                            //const char* true_forward_solution_filename_start = "/import/shared/friedrich/source/precomputed/ydata/u_problem_13_wnum_2_primbs_d_dt_3_3_bc_ffff_jmax_5_tstep_10_",
                            //const char* true_forward_solution_filename_end = ".iv",
                            const char* haar_gramian_filename_start = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_",
                            const char* haar_gramian_filename_end = "_gramian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd",
                            const char* laplacian_filename = "/import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd",
                            //const SparseMatrix<double>* laplacian_matrix,
                            //const bool matrix_precond,
                            //const char* output_path,
                            const int par_jmax = 5,
                            const int spatial_jmax = 8,
                            const double tolerance = 1e-6
                            );
 
    /*!
     * as above, but with RHS f given as a (time dependent) function
     * (gcc 2.95 requires implementation of this constructor already in the *.h file!?)
     */
    /*
    AffLinParEq_Precomputed(const PROBLEM* eq,
                            const InfiniteVector<double,Index>& u0,
                            MathTL::Function<PROBLEM::space_dimension,double>* f = 0,
                            const double d,
                            const FixedVector<InfiniteVector<double,Index>, NUMOFTIMESTEPS>& w,
                            const FixedVector<double, NUMOFTIMESTEPS>& time_discretization,
                            const Vector<SparseMatrix<double> *>& matrix_storage,
                            const SparseMatrix<double>* laplacian_matrix,
                            //const bool matrix_precond,
                            const char* output_path,
                            const int par_jmax = 5,
                            const int sol_jmax = 8)
    : problem(eq), f_vector_(), f_function_(f), d_(d), w_(w), time_discretization_(time_discretization), raw_cache(matrix_storage), laplacian_matrix_(laplacian),
    //matrix_precond_(matrix_precond),
    output_path_(output_path), par_jmax_(par_jmax), sol_jmax_(sol_jmax)
    {
        AbstractIVP<InfiniteVector<double,typename PROBLEM::Index> >::u0 = u0;
// TODO: Load and store the laplacian matrix!
    }
    */

    /*
     * store solution.
     * time_direction == true = forward case:
     * write coeffs to forward_solution_ and forward_solution_preprocessed_ at position current_timestep_.
     * Preprocessing is done with the given tolerance.
     * otherwise = backward case:
     * write to backward_solution_[NUMOFTIMESTEPS-current_timestep_]
     *
     * tolerance is used for call to APPLY(identity_) for preprocessing the forward solution
     */
    void set_solution(const InfiniteVector<double, Index>& u, const double& tolerance);
      
    /*
     * Background: The adjoint problem (time backward case) needs the observed data that is stored
     * internally as coeffs w.r.t. the dual basis (==preprocessed)
     * This preprocessed data is used often, so we store it in true_forward_solution_preprocessed_
     * Dont forget that we have to deal with noisy data.
     * 
     * Side effect: parabolic_forward_solution_ is used as a temporary storage!
     * 
     * uobserved_coeffs should contain coeffs of ydata + noise of magnitude delta
     * tolerance = tolerance for the call of APPLY(*identity_)
     */
    void set_true_forward_solution_preprocessed(const Array1D<InfiniteVector<double, Index> > &uobserved_coeffs,
                                                  const double tolerance);

    /*!
     * evaluate the right-hand side F(t,v)=Av+f(t) up to a prescribed tolerance
     * v vector w.r.t. the primal basis
     * result vector w.r.t. dual basis
     *
     * caution: in general stage equations are solved at intermediate time points.
     * t does not necessarily correspond to a time point in time_discretization_.
     * We (might) have to temporal interpolate the operator A!
     * It is checked whether t is close to t_k or t_{k+1} (no interpolation necessary)
     * or not (interpolate!)
     *
     * output 
     * time forward case:  Av+f_
     * time backward case: Av+(forward_solution_preprocessed-true_solution_preprocessed)
     * 
     */
    void evaluate_f(const double t,
                    const InfiniteVector<double,Index>& v,
                    const double tolerance,
                    InfiniteVector<double,Index>& result) const;

    /*!
     * evaluate F_t(t,v)=f'(t) up to a prescribed tolerance
     * result vector w.r.t. dual basis
     *
     * t is always a point of time_discretization_ !!!!
     *
     * wrong:
     * for the time forward case: constant f_ means 0 time derivative
     * time backward case: the solution of the time forward case has
     * to be interpolated to give values for the intermediate time points
     * (t does not necessarily correspond to a point of time_discretization_)
     */
    void evaluate_ft(const double t,
                     const InfiniteVector<double,Index>& v,
                     const double tolerance,
                     InfiniteVector<double,Index>& result) const;

    /*!
     * solve (alpha*I-A)x=y up to a prescribed tolerance
     */
    void solve_ROW_stage_equation(const double t,
                                  const InfiniteVector<double,Index>& v,
                                  const double alpha,
                                  const InfiniteVector<double,Index>& y,
                                  const double tolerance,
                                  InfiniteVector<double,Index>& result) const;

    /*!
     * preprocess the right-hand side share from the previous steps
     * (cf. w_method.h for details)
     * Transforms a primal coefficient vector, i.e. (<u,\tilde\psi_\lambda>) into a dual one, ie (<u,\psi_\lambda>)
     * assumes that the first matrix in raw_storage is the gramian matrix
     */
    void preprocess_rhs_share(InfiniteVector<double,Index>& wbeforeandafter,
                              const double tolerance) const;

    /*
     * return the number of the time discretization point closest to t
     * Numbering starts at tk closest to 0
     */
    //int get_timestep_number(double t) const;

    /*
     * compute the matrices <w\psi_lambda,\psi mu> where w is given by the Haar wavelet coefficients new_w
     *
     * The memory usage strategy of our choice is implemented in this method.
     */
    void assemble_W(const Array1D<InfiniteVector<double,int> >& new_w);      

    inline void set_current_timestep(const unsigned int tk) { current_timestep_ = time_direction_? tk: NUMOFTIMESTEPS-tk; };

    inline void set_time_direction(const bool new_time_direction_) { time_direction_ = new_time_direction_; };

    //protected:
    //! the underlying (uncached) problem
    // PROBLEM = TensorEquation<Basis1d,dim,Basis>
    PROBLEM * problem_;
    
    //const CompressedProblemFromMatrix<PROBLEM >* scaled_laplacian; // ensure theese problems use unpreconditioned matrices, otherwise we would be dependent on problem.D()
    const CompressedProblemFromMatrix<PROBLEM >* identity_;  // ensure theese problems use unpreconditioned matrices

    // to load/store the right filenames:
    //const string par_basisname; // = "haar";
    //const string sol_basisname; // = "primbs";

    // initial value for the time forward case
    // const InfiniteVector<double,Index>& u0; // initial value is stored in forward_solution[0]

    // computed solutions
// TODO InVec*
    Array1D<InfiniteVector<double,Index> > forward_solution_, forward_solution_preprocessed_, backward_solution_, true_forward_solution_preprocessed_;
    //FixedVector<InfiniteVector<double,Index> , NUMOFTIMESTEPS+1> forward_solution_, forward_solution_preprocessed_, backward_solution_; // coeffs should be unpreconditioned

    InfiniteVector<double,Index> f_; // f= 1/2, used for the time forward case, f=forward_solution for the time backward case

    // diffusion coefficient
    const double d_;

    // coefficients of the parameter w
    //FixedVector<InfiniteVector<double,int>, NUMOFTIMESTEPS>& w_; // Haar wavelets use integer indices. coeffs do not need to be stored since the assembled matrices/problems are stored

    // time discretization
    const Array1D<double> time_discretization_;

    bool time_direction_; // true = normal, false = reversed.
    unsigned int current_timestep_; // this solution has already been computed

    // cache for the basic <h_eta psi_lambda, psi_mu> matrices (mutable to overcome the constness of add_column())
    // int = number of the basis function h_eta
    // if matrices would be loaded on demand, we could use:
    //typedef std::map<int, SparseMatrix<double> > MatrixCache;
    // to store the loaded matrices.
    //typedef Vector<SparseMatrix<double> * > MatrixCache;
    typedef Array1D<SparseMatrix<double> > MatrixCache;
    mutable MatrixCache raw_cache_; // haar_gramian matrices
    MatrixCache assembled_matrices_; // matrices of the assembled_problems_

    // matrix of the laplacian (mutable to overcome the constness of add_column()), scales with d_
    mutable SparseMatrix<double> scaled_laplacian_matrix_;

    const char* haar_gramian_filename_start_, *haar_gramian_filename_end_; // and many more paths ...

    const unsigned int spatial_jmax_;
    //const unsigned int par_jmax_, sol_jmax_; // maximal levels
    //const unsigned int par_dof_; // degrees of freedom

 // Macht es Sinn die Matrizen oder gleich die CachedTProblems_Precomputed zu speichern??
    // cache for the assembled matrices d*<grad psi_lambda, grad psi_mu> + <w(t) psi_lambda, psi_mu> (mutable to overcome the constness of add_column())
    //typedef FixedVector<SparseMatrix<double>* , NUMOFTIMESTEPS> AssemblyStorage;
    //mutable AssemblyStorage assembled_matrices;
    //typedef FixedVector< CompressedProblemFromMatrix<PROBLEM >*, (NUMOFTIMESTEPS+1)> ProblemCache;
    typedef Array1D<CompressedProblemFromMatrix<PROBLEM > > ProblemCache;
    mutable ProblemCache assembled_problems_;

    //const double tolerance_; // coeffs of the time forward solution are transformed into coeffs w.r.t. dual basis with this tolerance AND afterwards the result is compressed with tolerance 1e-13
    };
  
    /*!
     * class AffLinParEq_qtbasis
     * 
     * similar to AffLinParEq_Precomputed, but uses qtbasis and does not load matrices from disk
     * 
     * This class will be used for a parameter
     * identification problem. Indeed a simplified version of the problem considered in project DA360/12-1.
     * A parabolic problem is considered in the usual time forward way and also in a time reversed way.
     * The parabolic equation considered is of the form
     *    u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
     *    u(0) = u_0
     * where -A: \ell_{2,D} -> \ell_{2,D^{-1}} is a positive isomorphism
     * and f: (0,T] -> \ell_{2,D^{-1}}, D being a diagonal preconditioner (see below).
     * An equation of this form arises when equivalently reformulating
     * a problem of the analogous form
     *    v'(t) = Bv(t) + g(t) =: G(t,v(t)), 0 < t <= T
     *    v(0) = v_0
     * with isomorphism -B: H -> H' and right-hand side g: (0,T] -> H'
     * by means of a biorthogonal wavelet basis Psi=\{\psi_\lambda\},
     * setting F(t,v):=<G(t,v^\top\Psi),\tilde\Psi> for all v in \ell_{2,D}.
     * The (unpreconditioned) stiffness matrix
     *    -A = (a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}
     * is assumed to be of the form
     *    A = (d <grad psi_lambda, grad psi_mu> + < w psi_lambda, psi_mu>)_{lambda,mu}
     * where w is arbitrary in time and space whereas d is assumed constant. w[t_k] is given
     * by its coefficients with respect to the Haar wavelet basis.
     * 
     * The problem can be solved in a normal time forward way. In this case f=1/2 is used
     * as a source term. Or the problem can be solved with reversed time direction. 
     * In this case v(T)=0 and f=\hat v is used where \hat v is the solution of the 
     * time forward case that has been previously computed.
     * 
     * At the kth time step an elliptic problem using the (scaled) Laplacian and the matrix w[t_k] has to be solved. 
     * To this end this class is a child of CachedQTProblem, but with a different system matrix for each time step. 
     * By this, we achieve that the caches for the onedimensional integrals can be reused for all time steps. 
     * Note, that the Laplacian has only a constant coefficient. 
     * Therefore we only need a single HaarWavelet (i.e., the first Generator) for the computation of the derivatives of the spatial wavelets.
     * 
     * The solutions of the time forward and time backward parabolic problem are then used to update
     * the parameter w. To this end the Haar wavelet coefficients of the product of the time forward and time backward solution need to be computed.
     * This is realized by sampling both functions on a grid with resolution >= the number of Haar wavelets.
     * 
     */

// CLEANUP
    
//    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT>
//    class Dummy
//    : public CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>
//     //AbstractCachedIVP<InfiniteVector<double, int> > // diese Vererbung schmiert immer ab!
//    {
//        public:
//        Dummy(QTBASIS* qtbasis,
//                              //const double d,
//                              const char* onehalf_filename, 
//                              const double normA, // from CQTProblem: estimate for \|dI+W\|
//                              const double normAinv // from CQTProblem:  estimate for \|(dI+W)^{-1}\|
//                              )
//            : CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT,1>(qtbasis, 
//                      Array1D<FixedMatrix<double, 1> > (), 
//                      Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > (),
//                      NULL, onehalf_filename, normA, normAinv)
//        {
//            
//        }
//        
//        void evaluate_f(const double t,
//			    const InfiniteVector<double, int>& v,
//			    const double tolerance,
//			    InfiniteVector<double, int>& result)
//        {
//            
//        }
//        void evaluate_ft(const double t,
//                                 const InfiniteVector<double, int>& v,
//                                 const double tolerance,
//                                 InfiniteVector<double, int>& result)
//        {
//
//        }
//        void solve_ROW_stage_equation(const double t,
//                                              const InfiniteVector<double, int>& v,
//                                              const double alpha,
//                                              const InfiniteVector<double, int>& y,
//                                              const double tolerancs,
//                                              InfiniteVector<double, int>& result)
//        {
//
//        }
//    };
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    class AffLinParEq_qtbasis
    : public CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>,
     public WMethodPreprocessRHSHelper<InfiniteVector<double, int> >
     //public AbstractCachedIVP<InfiniteVector<double, int> > // diese Vererbung schmiert immer ab!
    {
    public:
    /*
     * wavelet index class
     */
    typedef typename QTBASIS::Index Index;
    
    /*!
     * space dimension of the problem
     */
    //static const int space_dimension = QTBASIS::space_dimension;
    static const unsigned int DIM = QTBASIS::space_dimension;

    typedef std::pair<FixedArray1D<double,ONEDIMHAARCOUNT>, FixedArray1D<double,1> > entries;
    typedef typename CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>::Block Block;
    typedef typename CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>::Column Column;
    typedef typename CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>::ColumnCache ColumnCache;
      
    /*!
     * Models the parabolic solver for the linearized parabolic problem.
     * Time forward and backward case are put together in this class and it is a child of CachedQTProblem.
     * Therefore we may use the same cache for the onedimensional integrals in every time step.
     * u0 = initial value. coeffs w.r.t. primal basis (<u,\tilde \psi_\lambda>), unpreconditioned
     * d = diffusion coefficient
     * w = for every time step: for every patch of the geometry: Matrix of Haar-generator coeffs (pixel values) of W
     * time_discretization = vector of time points
     * f_filename = filename of InfiniteVector corresponding to coeffs of f=1/2 w.r.t. dual basis, unpreconditioned, eg
     *            = "/import/shared/friedrich/source/precomputed/functions/onehalf_primbs_d_dt_3_3_bc_ffff.iv"
     * haar_gramian_filename_start _end = every char before and after the number eta in the names of the matrices, eg
     *      start = /import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/haar_wav_
     *      end   = _gramian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd.bin
     * laplacian_filename = name of the file containing the matrix of the laplacian, eg
     *            = /import/shared/friedrich/source/precomputed/stiffness_matrices/unpreconditioned/laplacian_primbs_d_dt_3_3_bc_ffff_jmax_9_npcd.bin
     * //output_path = folder where output is written to,
     * par_jmax = maximal level for the Haar wavelet basis of w. Needed to determine size of matrix storage
     * spatial_jmax = maximal level for the spatial basis {psi_lambda}_lambda
     * tolerance - tolerance of APPLY used in transforming the coeffs of the forward solution into coeffs w.r.t. dual basis. The preprocessed coeffs are afterwards coarsened with tol=1e-13
     * 
     * f = rhs per time step. Should be the constant function 1/2 or the pointer NULL
     *  if it is 1/2: coeffs w.r.t. dual basis, unpreconditioned, are computed and stored in rhs_filename
     *  if it is NULL: rhs_filename is loaded from disk
     * normA, normAinv: some norm estimates that are used in CDD1_SOLVE of the stage equations
     * 
     * uexact != NULL => compute coefficients of the function f==1/2 and store them to filename
     * onehalf_filename != NULL => compute coefficients of the function f==1/2 and store them to filename
     * onehalf_filename == NULL => load the coefficients of the function f==1/2 from filename
     */
    AffLinParEq_qtbasis(QTBASIS* qtbasis,
//                              Function<QTBASIS::space_dimension,double>* uexact_function,
//                              const char* u0_filename,
                              const double d,
//                              const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1>& w,
//                              const FixedArray1D<double, NUMOFTIMESTEPS +1 >& time_discretization,
//                              const Function<QTBASIS::space_dimension>* onehalf_function, // from CQTProblem: rhs as a function
                              const char* onehalf_filename, 
//                              const char* onehalf_precond_gramian_filename, 
//                              const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename,
                              const double normA, // from CQTProblem: estimate for \|dI+W\|
                              const double normAinv // from CQTProblem:  estimate for \|(dI+W)^{-1}\|
//                              const double tolerance_APPLY
                              )
            : CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT,1>(qtbasis, 
                      //Array1D<FixedMatrix<double, 1> > (), 
                      //Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > (),
                      //NULL, onehalf_filename, 
                      normA, normAinv),
      d_(d)//, time_discretization_(time_discretization), 
//            current_timestep_(0),
//            problem_mode_(0)
    {
//            //haar_gramian_filename_start_(haar_gramian_filename_start), 
//            //haar_gramian_filename_end_(haar_gramian_filename_end),
//            //spatial_jmax_(spatial_jmax),
//        time_direction_ = true;
//        
//        if (onehalf_function != 0)
//        {
//            compute_rhs(onehalf_filename, onehalf_precond_gramian_filename, onehalf_mode_filename);
//        }
//        else
//        {
//            load_rhs(onehalf_filename, onehalf_precond_gramian_filename, onehalf_mode_filename);
//            
//        }
//        InfiniteVector<double,int> temp_u0;
//        //Vector<double> temp_u0;
//        
//        if (uexact_function != NULL)
//        {
//            cout << "AffLinParEq_qtbasis:: use CDD1 to compute u0T (u0 w.r.t. primal basis) ..." << endl;
//            uexact_function->set_time(0);
//            qtbasis_->expand(uexact_function, false, temp_u0);
//            current_rhs_ = &temp_u0;
//            current_rhs_l2_norm_ = l2_norm(temp_u0);
//            //another_rhs_ = temp_u0;
//            //current_rhs_ = &another_rhs_;
//            problem_mode_ = 2;
//            CDD1_SOLVE(*this, tolerance_APPLY, temp_u0, qtbasis_->get_jmax(), tensor_simple);
//            temp_u0.scale(this,-1);
//            temp_u0.compress(1e-13);
//            cout << "AffLinParEq_qtbasis:: write u0T to file " << u0_filename << endl;
//            writeIVToFile(temp_u0,u0_filename);
//            cout << "AffLinParEq_qtbasis:: done" << endl;
//        }
//        else
//        {
//            cout << "AffLinParEq_qtbasis:: load u0T from file " << u0_filename << endl;
//            readIVFromFile(temp_u0,u0_filename);
//            cout << "AffLinParEq_qtbasis:: done" << endl;
//        }
//        //AbstractIVP<InfiniteVector<double,int> >::u0 = temp_u0; // is used only for the time forward problem
//        this->u0 = temp_u0;
//        
//        //cout << "u0 = " << endl << u0 << endl;
//        set_solution(temp_u0, tolerance_APPLY);
    }
    
     AffLinParEq_qtbasis(QTBASIS * qtbasis,
             Function<QTBASIS::space_dimension,double>* uexact_function,
             const char* u0_filename,
             const double d,
             const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1>& w,
             const FixedArray1D<double, NUMOFTIMESTEPS +1 >& time_discretization,
             const Function<QTBASIS::space_dimension>* f, // from CQTProblem: rhs as a function. If NULL: load data from disk using the filenames provided. Otherwise: compute everything and store it to the given filenames
             const char* onehalf_filename, // filename of IV that contains the coeffs of onehalf_function
             const char* onehalf_precond_gramian_filename, 
             const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename,
             const double normA, // from CQTProblem: estimate for \|dI+W\|
             const double normAinv, // from CQTProblem:  estimate for \|(dI+W)^{-1}\|
             const double tolerance_APPLY = 1e-6
             );
//             : CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT,1>(qtbasis, 
//                      Array1D<FixedMatrix<double, 1> > (), 
//                      Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > (),
//                      NULL, onehalf_filename, normA, normAinv),
//      d_(d)//, time_discretization_(time_discretization), 
////            current_timestep_(0),
////            problem_mode_(0)
//    {
//         
//    }
      
    /*
     * store solution.
     * time_direction == true = forward case:
     * write coeffs to forward_solution_ and forward_solution_preprocessed_ at position current_timestep_.
     * Preprocessing is done with the given tolerance.
     * otherwise = backward case:
     * write to backward_solution_[NUMOFTIMESTEPS-current_timestep_]
     *
     * APPLY_TENSOR_tolerance is used for call to APPLY_TENSOR(identity_) for preprocessing the forward solution
     */
    void set_solution(const InfiniteVector<double, int>& u, const double& APPLY_TENSOR_tolerance);
      
    /*
     * Background: The adjoint problem (time backward case) needs the observed data that is stored
     * internally as coeffs w.r.t. the dual basis (==preprocessed)
     * This preprocessed data is used often, so we store it in true_forward_solution_preprocessed_
     * Dont forget that we have to deal with noisy data.
     * 
     * uobserved_coeffs should contain coeffs of ydata + noise of magnitude delta
     * tolerance = tolerance for the call of APPLY(*identity_)
     */
    void set_true_forward_solution_preprocessed(const FixedArray1D<InfiniteVector<double, int>, NUMOFTIMESTEPS+1 > &uobserved_coeffs,
                                                const double APPLY_TENSOR_tolerance);
    /*!
     * evaluate the right-hand side F(t,v)=Av+f(t) up to a prescribed tolerance
     * v vector w.r.t. the primal basis
     * result vector w.r.t. dual basis
     *
     * caution: in general stage equations are solved at intermediate time points.
     * t does not necessarily correspond to a time point in time_discretization_.
     * We (might) have to temporal interpolate the operator A!
     * It is checked whether t is close to t_k or t_{k+1} (no interpolation necessary)
     * or not (interpolate!)
     *
     * output 
     * time forward case:  Av+f_
     * time backward case: Av+(forward_solution_preprocessed-true_solution_preprocessed)
     * 
     */
    
//    void evaluate_f(const double t,
//			    const InfiniteVector<double, int>& v,
//			    const double tolerance,
//			    InfiniteVector<double, int>& result)
//        {
//            
//        }
//        void evaluate_ft(const double t,
//                                 const InfiniteVector<double, int>& v,
//                                 const double tolerance,
//                                 InfiniteVector<double, int>& result)
//        {
//
//        }
//        void solve_ROW_stage_equation(const double t,
//                                              const InfiniteVector<double, int>& v,
//                                              const double alpha,
//                                              const InfiniteVector<double, int>& y,
//                                              const double tolerancs,
//                                              InfiniteVector<double, int>& result)
//        {
//
//        }
    void evaluate_f(const double t,
                    const InfiniteVector<double,int>& v,
                    const double tolerance,
                    InfiniteVector<double,int>& result);
    
    /*!
     * evaluate F_t(t,v)=f'(t) up to a prescribed tolerance
     * result vector w.r.t. dual basis
     *
     * t is always a point of time_discretization_ !!!!
     *
     * wrong:
     * for the time forward case: constant f_ means 0 time derivative
     * time backward case: the solution of the time forward case has
     * to be interpolated to give values for the intermediate time points
     * (t does not necessarily correspond to a point of time_discretization_)
     */
    void evaluate_ft(const double t,
                     const InfiniteVector<double,int>& v,
                     const double tolerance,
                     InfiniteVector<double,int>& result);

    /*!
     * solve (alpha*I-A)x=y up to a prescribed tolerance
     */
    void solve_ROW_stage_equation(const double t,
                                  const InfiniteVector<double,int>& v,
                                  const double alpha,
                                  const InfiniteVector<double,int>& y,
                                  const double tolerance,
                                  InfiniteVector<double,int>& result);

    /*!
     * preprocess the right-hand side share from the previous steps
     * (cf. w_method.h for details)
     * Transforms a primal coefficient vector, i.e. (<u,\tilde\psi_\lambda>) into a dual one, ie (<u,\psi_\lambda>)
     */
    void preprocess_rhs_share(InfiniteVector<double,int>& wbeforeandafter,
                              const double tolerance) const;

    /*
     * return the number of the time discretization point closest to t
     * Numbering starts at tk closest to 0
     */
    //int get_timestep_number(double t) const;

    /*
     * Evaluate the bilinear form a of the underlying problem
     * depend on problem_mode_
     * if problem_mode_ 
     *        == 0: *this behaves like T = d*laplace+W
     *        == 1: *this behaves exactly as \alpha G - T
     *        == 2: *this behaves exactly as the Identity problem 
     * Caching: compute and store all onedim integrals needed for the entries a(mu2,nu) with
     * mu2 having the same sublevel as mu. 
     * 
     * We assume that the operator is local
     * 
     */
    double a(const int munum, const int nunum);
    
    /*
     * taken from cached_qtproblem
     * specialized for the parabolic setting and the time discretization
     * depending on problem_mode_ this method behaved different
     * By the assumption that the diffusion coefficient is constant, we have
     * DER_ONEDIMHAARCOUNT == 1 and GRAM_ONEDIMHAARCOUNT == ONEDIMHAARCOUNT
     * in cached_qtproblem
     * In fact, this is the reason, why these two template arguments have been introduced.
     * Input: 
     * patch (lambda's patch)
     * specialized_intinfo (needs to include info from both wavelets! not only one wavelet and (j,p) of the other! get_LMR_info is NOT sufficient)
     * integralshares (one dimensional integrals; lambda_i is reflected to mu's patch if necessary)
     * 
     * output: a(mu, lambda)
     * 
     * specialized_intinfo :
     * this now contains information about both wavelets!
     * eg. "LR" means "nu is extended left and mu is extended right. They are both boundary wavelets/generators"
     * The output of get_LMR_info is different: "LR" means:
     * "nu is extended left and the patch of mu is extended right, e.g., there exists some mu that are extended right, but there also exist some that are NOT extended right and that intersect nu.
     * nu is a boundary wavelets/generators. No information about the type of the other wavelet/generator is included.
     */  
    double compute_sum_over_patches (const int nu_p, 
            const MultiIndex<unsigned int, DIM> intinfo, 
            const FixedArray1D<entries,DIM> integralshares) const;
    /*
     * compose the entries a(mu,nu).
     * For i=0,...,DIM-2: use the onedim-integral values given by block_it.
     * For i=DIM-1: take all entries in the Block starting at block_it[DIM-1] (the remaining dimension)
     * blocksize = size of the block with first element block_it[DIM-1]
     * mu_adapted_info_it: analog to block_it
     * add factor times those entries to w[number],...,w[number+block_it.size()-1]
     * The patches intersecting the support of mu are computed using nu_p and the intinfo provided by mu_adapted_intinfo1 and _2
     * if (precond == true) devide w[number] by D(number)
     * 
     * ATTENTION: header is for arbitrary dimensions!
     * Internally we assume that agencoefs and qgencoeffs are matrices, iff DIM = 2
     * Otherwise we assume ONEDIMHAARCOUNT == 1 and compute the Poisson equation
     */
    void compose_wavelets (Vector<double>& w,
            const unsigned int start,
            const unsigned int blocksize,
            const double factor,
            const unsigned int nu_p,
            const FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it,
            const FixedArray1D< typename Block::iterator, DIM> block_it,
            const bool precond) const;
    
    /*
     * compute square root of the diagonal entries of the gramian
     */        
    void compute_D_gramian();
    
    /*
     * used in constructor and whenever an update to the matrix W has been computed
     * 
     * - stores Wnew
     * - recomputes whatever is necessary
     * - assumes that W_new[time][patch] is a Haar-Wavelet-Coefficient Matrix, i.e.,
     *   application of inverse_haar_wavelet_trafo * 2^{(DIM*level)/2} yields pixel values
     */
    void set_W(const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & W_new);
    //void assemble_W(const Array1D<InfiniteVector<double,int> >& new_w);
    
    

    inline void set_current_timestep(const unsigned int tk) { current_timestep_ = time_direction_? tk: NUMOFTIMESTEPS-tk; };

    inline void set_time_direction(const bool new_time_direction_) { time_direction_ = new_time_direction_; };

        
    /*
     * Evaluate the diagonal preconditioner D, i.e., the square root of the diagonal of A.
     * Values are taken from the cache if possible. 
     * If not, the whole column a_{mu,nu} with |mu|=|nu| is computed and stored.
     * 
     * Maybe this is too much effort!
     * Alternative: a separate cache for the diagonal of A.
     * 
     * diag is included only for compatibility with 
     * FullyDiagonalEnergyNormPreconditioner<int>
     * 
     * enables the use of InfiniteVector . scale (this,-1) <-- useful
     */
    inline double D(const int lambdanum) const
    {
        if ( (problem_mode_ == 0) || (problem_mode_ == 1) )
        {
            return diagonal_cache_mode_[problem_mode_][current_timestep_][lambdanum];
        }
        else
        {
            return diagonal_cache_gramian_[lambdanum];
        }
    }
    
    inline double diag(const int& lambdanum) const
    {
        if ( (problem_mode_ == 0) || (problem_mode_ == 1) )
        {
            return diagonal_cache_mode_[problem_mode_][current_timestep_][lambdanum];
        }
        else
        {
            return diagonal_cache_gramian_[lambdanum];
        }
    }
    
    /* For compatibility with CDD1 */
    inline double F_norm() const
    {
        return current_rhs_l2_norm_;
        //sqrt(onehalf_precond_norm_sqr_mode_[problem_mode_]);
    }
        
    /*
     * set another_rhs_ to the given value.
     */
    inline void set_another_rhs_(InfiniteVector<double, int>* new_rhs)
    {
        //another_rhs_ = new_rhs;
        //current_rhs_ = &another_rhs_;
        current_rhs_ = new_rhs;
        current_rhs_l2_norm_ = l2_norm(*new_rhs);
    }
            
    /*
     * set up the fcoeffs_, f_precond_norm_sqr_
     * store f_coeffs_precond_unsorted_ to disc
     */
    void compute_rhs(const char* onehalf_filename, const char* onehalf_precond_gramian_filename, const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename);

    /*
     * set up the fcoeffs_, f_precond_norm_sqr_
     * to this end: load f_coeffs_precond_unsorted_ from disc
     */
    void load_rhs(const char* onehalf_filename, const char* onehalf_precond_gramian_filename, const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename);
        
    /*
     * approximate the wavelet coefficient set of the preconditioned right-hand side F
     * within a prescribed \ell_2 error tolerance.
     * uses only entries from fcoeffs.
     */
    void RHS(const double eta,
             InfiniteVector<double, int>& coeffs) const;

    /*
     * access to the preconditioned RHS : f(lambda)
     */
    inline double f(const unsigned int lambdanum) const
    {
        // this method should only be called by setup_righthand_side via CDD1_SOLVE->GALERKIN
        //return another_rhs_.get_coefficient(lambdanum);
        return current_rhs_->get_coefficient(lambdanum);
    }
    
    /*
     * Called by APPLY_TENSOR
     * copied from cached_qtproblem
     * 
     * computes
     * w += factor * (stiffness matrix entries in column lambda with ||nu-lambda|| <= radius && ||nu|| <= jmax)
     */
    void add_ball(const unsigned int& lambdanum,
                Vector<double>& w,
                const int radius,
                const double factor,
                const int maxlevel = 99,
                const CompressionStrategy strategy = tensor_simple,
                const bool precond = true);
    
    /*
     * called by add_ball for DIM >= 3
     * copied from cached_qtproblem
     * 
     * It may happen that the minimal levels for the patches are not the same. 
     * Therefore it may happen, that not all patches allow for a level with center_j[0],...,center_j[fixed_dimensions]
     * If this is the case it holds that set_is_restriction == true and j0_restriction_of_patches is a list of all patches that allow the current restrictions by center_j
     */
    void add_leveldisc_recurse(
            const unsigned int& lambdanum,
            Vector<double>& w,
            const MultiIndex<int,DIM> center_j,
            const int fixed_dimensions,
            const unsigned int radius,
            const double factor,
            const list<int> & j0_restriction_of_patches,
            const bool set_is_restriction,
            const bool precond);
    
    /*
     * w += factor * A_levelnum e_lambda
     * A_levelnum = matrix where all entries a_(mu,nu) with nu != lambda and |mu| != levelnum are 0
     * e_lambda = lambda_th unit vector
     */
    void add_level(const unsigned int& lambdanum,
            Vector<double>& w,
            const unsigned int levelnum,
            const double factor,
            const bool precond);
        
    /*! if problem_mode_ 
     *        == 0: *this behaves like T = d*laplace+W
     *        == 1: *this behaves exactly as \alpha G - T
     *        == 2: *this behaves exactly as the Identity problem 
     */
    inline
    void
    set_problem_mode (unsigned int mode)
    {
        problem_mode_ = mode;
    }
    
    /*
     * works only in 2D
     * f_haar_gen_coeffs [patch] = Matrix with pixelvalues of function f given by f_qtbasis_coeffs:
     * This is the same as SampledMapping does, with one difference: one dimensional components of the involved wavelet evaluations are stored in a cache
     * 
     * observe: output are pixel values. NOT Haar-generator coeffs! (Scaling is different)
     */
    void
    cached_sampled_mapping(const InfiniteVector<double,int> & f_qtbasis_coeffs, 
            Array1D < FixedMatrix<double, PRECISE_EVALUATE_GRANULARITY, PRECISE_EVALUATE_GRANULARITY> >& f_haar_gen_coeffs);
    
    /*
     * compute one dimensional integrals 
     * \int g_\eta \psi_{\lam_i} \eta =0,...,PRECISE_EVALUATE_GRANULARITY-1
     * typically PRECISE_EVALUATE_GRANULARITY > ONEDIMHAARCOUNT
     */
    void precise_evaluate (const int lami_j,
            const int lami_e,
            const int lami_k,
            const unsigned int lami_basisnum,
            FixedArray1D<double, PRECISE_EVALUATE_GRANULARITY>& precise_lambda_integrals);
    
    // computed solutions
    mutable FixedArray1D<InfiniteVector<double,int>, NUMOFTIMESTEPS+1 > forward_solution_, forward_solution_preprocessed_, backward_solution_, true_forward_solution_preprocessed_; // coeffs should be unpreconditioned
    
    InfiniteVector<double,int> u0_; // initial value
    protected:
     
    /*! if problem_mode_ 
     *        == 0: *this behaves like T = d*laplace+W
     *        == 1: *this behaves exactly as \alpha G - T
     *        == 2: *this behaves exactly as the Identity problem 
     */
    mutable unsigned int problem_mode_; // mutable to avoid const of evaluate_f (Abstact IVP)
    mutable double mode_one_alpha_; // mutable to avoid const of solve_ROW_stage_equation
    //InfiniteVector<double,int> another_rhs_;
    //Vector<double> another_rhs_;
    mutable InfiniteVector<double,int>* current_rhs_; // This rhs is used in P.rhs and P.f(int); // mutable to avoid const of evaluate_f (Abstact IVP)
    //Vector<double>* current_rhs_; // This rhs is used in P.rhs and P.f(int);
    //double another_rhs_l2_norm_;
    mutable double current_rhs_l2_norm_; // mutable to avoid const of evaluate_f (Abstact IVP)
    mutable Array1D<std::pair<int,double> > * current_rhs_sorted_;
    
    //InfiniteVector<double,int>* mode_one_rhs_; // pointer to the RHS y of the stage_equations
      
    // store the diagonal entries 
    FixedArray1D< FixedArray1D <Vector<double>, NUMOFTIMESTEPS+1 >, 2> diagonal_cache_mode_;
    Vector<double> diagonal_cache_gramian_;
    
    //Vector<double> diagonal_cache_mode_1_;
    //Vector<double> diagonal_cache_mode_2_;
    
    // right-hand-site f = 1/2
    const Function<DIM> * onehalf_function_;
    InfiniteVector<double,int> onehalf_coeffs_;
        
    // \int f \psi_lambda coefficient vector
    // right-hand side coefficients on a fine level, sorted by modulus
    
    FixedArray1D< FixedArray1D<Array1D<std::pair<int,double> > , NUMOFTIMESTEPS+1>, 2> onehalf_coeffs_precond_sorted_mode_;
    FixedArray1D< FixedArray1D<Vector<double>, NUMOFTIMESTEPS+1>, 2> onehalf_coeffs_precond_unsorted_mode_;
    Array1D<std::pair<int,double> > onehalf_coeffs_precond_sorted_gramian_;
    Vector<double> onehalf_coeffs_precond_unsorted_gramian_;
        
    // (squared) \ell_2 norm of the precomputed and preconditioned(!) right-hand side
    FixedArray1D<FixedArray1D<double, NUMOFTIMESTEPS+1>, 2> onehalf_precond_norm_sqr_mode_;
    double onehalf_precond_norm_sqr_gramian_;
        
    //FixedArray1D< double, 3> f_precond_norm_sqr_mode_;
    //double f_precond_norm_sqr_mode_1_;
    //double f_precond_norm_sqr_mode_2_;
      

    //! the underlying QTBasis
    QTBASIS * qtbasis_;
     
    //const CompressedProblemFromMatrix<PROBLEM >* scaled_laplacian; // ensure theese problems use unpreconditioned matrices, otherwise we would be dependent on problem.D()
    //const CompressedProblemFromMatrix<PROBLEM >* identity_;  // ensure theese problems use unpreconditioned matrices

    // to load/store the right filenames:
    //const string par_basisname; // = "haar";
    //const string sol_basisname; // = "primbs";

    // initial value for the time forward case
    // const InfiniteVector<double,Index>& u0; // initial value is stored in forward_solution[0]


    
    //InfiniteVector<double,int> rhs_coeffs_unpreconditioned_; // f= 1/2, used for the time forward case, f=forward_solution for the time backward case

      // diffusion coefficient
    const double d_;

      // coefficients of the parameter w
      //FixedVector<InfiniteVector<double,int>, NUMOFTIMESTEPS>& w_; // Haar wavelets use integer indices. coeffs do not need to be stored since the assembled matrices/problems are stored

      // time discretization
    const FixedArray1D<double, NUMOFTIMESTEPS+1> time_discretization_;

    bool time_direction_; // true = normal, false = reversed.
    mutable unsigned int current_timestep_; // this solution has already been computed

      // cache for the basic <h_eta psi_lambda, psi_mu> matrices (mutable to overcome the constness of add_column())
      // int = number of the basis function h_eta
      // if matrices would be loaded on demand, we could use:
      //typedef std::map<int, SparseMatrix<double> > MatrixCache;
      // to store the loaded matrices.
      //typedef Vector<SparseMatrix<double> * > MatrixCache;
//      typedef Array1D<SparseMatrix<double> > MatrixCache;
//      mutable MatrixCache raw_cache_; // haar_gramian matrices
//      MatrixCache assembled_matrices_; // matrices of the assembled_problems_

    FixedArray1D< Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs_;
      
    // matrix of the laplacian (mutable to overcome the constness of add_column()), scales with d_
//      mutable SparseMatrix<double> scaled_laplacian_matrix_;

//      const char* haar_gramian_filename_start_, *haar_gramian_filename_end_; // and many more paths ...

//      const unsigned int spatial_jmax_;
      //const unsigned int par_jmax_, sol_jmax_; // maximal levels
      //const unsigned int par_dof_; // degrees of freedom

 // Macht es Sinn die Matrizen oder gleich die CachedTProblems_Precomputed zu speichern??
    // cache for the assembled matrices d*<grad psi_lambda, grad psi_mu> + <w(t) psi_lambda, psi_mu> (mutable to overcome the constness of add_column())
    //typedef FixedVector<SparseMatrix<double>* , NUMOFTIMESTEPS> AssemblyStorage;
    //mutable AssemblyStorage assembled_matrices;
    //typedef FixedVector< CompressedProblemFromMatrix<PROBLEM >*, (NUMOFTIMESTEPS+1)> ProblemCache;
//      typedef Array1D<CompressedProblemFromMatrix<PROBLEM > > ProblemCache;
//      mutable ProblemCache assembled_problems_;

    //const double tolerance_; // coeffs of the time forward solution are transformed into coeffs w.r.t. dual basis with this tolerance AND afterwards the result is compressed with tolerance 1e-13
    
    typedef FixedArray1D<double,PRECISE_EVALUATE_GRANULARITY> Entries;
    typedef std::map<int, Entries> OneBasisCache;
    typedef FixedArray1D<OneBasisCache,4> AllBasesCache; // there are 4 types of 1d bases
    AllBasesCache precise_eval_cache_;
    };
  
}
#include "aff_lin_par_eq.cpp"

#endif

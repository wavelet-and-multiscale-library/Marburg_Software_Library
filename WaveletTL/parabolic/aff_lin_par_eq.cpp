
// #include "aff_lin_par_eq.h"
namespace WaveletTL
{
    template <class CACHEDPROBLEM>
    AffLinParEq_Precomputed_StageEquationHelper<CACHEDPROBLEM>
    ::AffLinParEq_Precomputed_StageEquationHelper(const double a,
                                        const CACHEDPROBLEM* elliptic,
                                        const CACHEDPROBLEM* gramian,
                                        const InfiniteVector<double, typename CACHEDPROBLEM::Index>& z)
    : alpha(a), T(elliptic), G(gramian), y(z), y_scaled(z)
    {
// CLEANUP
        //cout << "AffLinParEq_Precomputed_StageEquationHelper" << endl;
        //cout << "z = " << endl << z << endl;
        y_scaled.scale(this,-1);
    }


    template <class CACHEDPROBLEM>
    void
    AffLinParEq_Precomputed_StageEquationHelper<CACHEDPROBLEM>
    ::add_ball(const Index& lambda,
               //InfiniteVector<double, Index>& w,
               Vector<double>& w,
               const int radius,
               const double factor,
               const int maxlevel,
               const CompressionStrategy strategy,
               const bool precond) const
    {
        // clearly caching the gramian and A at the same t15ime would be more efficient
#if 0
        InfiniteVector<double,Index> g;
        Vector<double> g_full(w.size());15
        //G->add_ball(lambda, g_full, radius, alpha*factor/T->D(lambda), maxlevel, strategy);
        G->add_ball(lambda, g_full, radius, alpha*factor/D(lambda), maxlevel, strategy, false);
        // hack: copy full vector g into sparse representation
        for (unsigned int i = 0; i < g_full.size(); i++)
        {
            if (g_full[i] != 0.)
            {
                g.set_coefficient(basis().full_collection[i], g_full[i]);
            }
        }
        g.scale(this, -1);
        T->add_ball(lambda, w, radius, factor, maxlevel, strategy, false);
#else
        double d1 = precond?D(lambda):1.0;
        InfiniteVector<double,Index> g;
        Vector<double> g_full(w.size());
        //G->add_ball(lambda, g_full, radius, alpha*factor/T->D(lambda), maxlevel, strategy);
        G->add_ball(lambda, g_full, radius, alpha*factor/d1, maxlevel, strategy, false); // do not use preconditioning with diagonal of G
        T->add_ball(lambda, g_full, radius, factor/d1, maxlevel, strategy, false); // do not use preconditioning with diagonal of T
        // hack: copy full vector g into sparse representation
        for (unsigned int i = 0; i < g_full.size(); i++)
        {
            if (g_full[i] != 0.)
            {
                g.set_coefficient(basis().get_wavelet(i), g_full[i]);
            }
        }
        g.scale(this, -1);

        Vector<double> g_full_new(w.size());
        for (typename InfiniteVector<double,Index>::const_iterator it(g.begin());
               it != g.end(); ++it) {
          g_full_new[it.index().number()] = *it;
        }
        //w.add(g);
        w += g_full_new;
#endif
    }


    template <class CACHEDPROBLEM>
    void
    AffLinParEq_Precomputed_StageEquationHelper<CACHEDPROBLEM>
    ::set_y(const InfiniteVector<double, Index> y_new)
    {
        y = y_new;
        y_scaled = y_new;
        y_scaled.scale(this,-1);
    }

    template <class CACHEDPROBLEM>
    void
    AffLinParEq_Precomputed_StageEquationHelper<CACHEDPROBLEM>
    ::set_alpha (const double alpha_new)
    {
        alpha = alpha_new;
    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::AffLinParEq_Precomputed(PROBLEM* problem,
                              const InfiniteVector<double,Index>& u0,
                              const double d,
                              const Array1D<InfiniteVector<double, int> >& w,
                              //const FixedVector<InfiniteVector<double, int>, (NUMOFTIMESTEPS+1)>& w,
                              const Array1D<double>& time_discretization,
                              //const FixedVector<double, (NUMOFTIMESTEPS+1)>& time_discretization,
                              const char* f_filename,
                              const char* haar_gramian_filename_start,
                              const char* haar_gramian_filename_end,
                              const char* laplacian_filename,
                              const int par_jmax,
                              const int spatial_jmax,
                              const double tolerance)
    : problem_(problem), d_(d), time_discretization_(time_discretization), haar_gramian_filename_start_(haar_gramian_filename_start), haar_gramian_filename_end_(haar_gramian_filename_end),
      spatial_jmax_(spatial_jmax), current_timestep_(0)
    {
        time_direction_ = true;
        AbstractIVP<InfiniteVector<double,typename PROBLEM::Index> >::u0 = u0; // is used only for the time forward problem
        
        //set_constant_f !!
// TODO check whether these coeffs need to be preconditioned or w.r.t. the dual basis
                
        // f_ uses Index. Could be changed to int. Maybe more efficient!
        readIVFromFile (& problem->basis(),f_,f_filename);
        
// TODO check meaning of tolerance
        
        forward_solution_.resize(NUMOFTIMESTEPS+1);        
        forward_solution_preprocessed_.resize(NUMOFTIMESTEPS+1);
        true_forward_solution_preprocessed_.resize(NUMOFTIMESTEPS+1);
        backward_solution_.resize(NUMOFTIMESTEPS+1);
#if _DIMENSION == 1
        raw_cache_.resize((1<<(par_jmax+1)));
#else
        raw_cache_.resize(1<<( (par_jmax+1)<<1) );
#endif

        assembled_matrices_.resize(NUMOFTIMESTEPS+1);
        assembled_problems_.resize(NUMOFTIMESTEPS+1);
        //cout << "try to create CompressedProblemFromMatrix" << endl;
        
        //CompressedProblemFromMatrix<PROBLEM> temp_problem (dummy,&matrix,5.5,35.0, false);
        char filename [250];        
        //haar_wav_0_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
        sprintf(filename, "%s%d%s",haar_gramian_filename_start_, 0, haar_gramian_filename_end_);
        cout << "aff_lin_par_eq:: Constructor:: load gramian from file " << filename << endl;
        raw_cache_[0].matlab_input(filename);

        cout << "aff_lin_par_eq:: Constructor:: try to create CompressedProblemFromMatrix" << endl;
        //PROBLEM* dummy;
        //CompressedProblemFromMatrix<PROBLEM> temp_problem (dummy,raw_cache_[0],5.5,35.0, false);
        identity_ = new CompressedProblemFromMatrix<PROBLEM> (problem_,&raw_cache_[0],5.5,35.0, false);
        //cout << " raw_cache_[0] " << endl << raw_cache_[0] << endl;
        //cout << identity_->entries_cache->get_entry(0,0);
        //cout << problem_->basis().get_wavelet(0);
        /* 
        // remove const from identity_ definition
        cout << "aff_lin_par_eq::constructor:: load f_ from file = " << endl << f_ << endl;
        ConstantFunction<2> constant_rhs(Vector<double>(1, "0.5"));
        identity_->set_f(&constant_rhs);
        CDD1_SOLVE(*identity_, tolerance, f_, 6, tensor_simple);
        f_.scale(identity_,-1); 
        cout << "aff_lin_par_eq::constructor:: compute f_ directly = " << endl << f_ << endl;
        */

        for (unsigned int i=0; i<= NUMOFTIMESTEPS;++i)
        {
            //PROBLEM* dummy;
            assembled_problems_[i].problem = problem;
            assembled_problems_[i].normA = 5.5;
            assembled_problems_[i].normAinv = 35.0;
            assembled_problems_[i].preconditioned = false;
            //= new CompressedProblemFromMatrix<PROBLEM> (dummy,&matrix,5.5,35.0, false);
        }

        //cout << "u0 = " << endl << u0 << endl;
        set_solution(u0, tolerance);

        //cout << "aff_lin_par_eq:: Constructor:: load laplacian from file " << laplacian_filename << endl;
        scaled_laplacian_matrix_.matlab_input(laplacian_filename);
        scaled_laplacian_matrix_.scale(d);

        assemble_W(w);
    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::set_solution(const InfiniteVector<double, Index>& u, const double& tolerance)
    {
        if (time_direction_ == true)
        {
            forward_solution_[current_timestep_] = u;
            InfiniteVector<double,Index> help(u);
    // TODO check whether preconditioning is useful here. maybe it just slows things down!
    // TODO true, false??
            help.scale(identity_,1); // should be 1 anyways ... although it isn't ...
            //APPLY_TENSOR(*assembled_problems[0], help, tolerance, wbeforeandafter, spatial_jmax, tensor_simple, true); // false?
            
            APPLY_TENSOR(*identity_, help, tolerance, forward_solution_preprocessed_[current_timestep_], spatial_jmax_, tensor_simple, true); // false?
            forward_solution_preprocessed_[current_timestep_].scale(identity_,1);
            forward_solution_preprocessed_[current_timestep_].compress(1e-13);
// CLEANUP
            //cout << "aff_lin_par_eq::set_solution:: forward_solution_[" << current_timestep_ << "] = " << endl << forward_solution_[current_timestep_] << endl;
            //cout << "aff_lin_par_eq::set_solution:: forward_solution_preprocessed_[" << current_timestep_ << "] = " << endl << forward_solution_preprocessed_[current_timestep_] << endl;
            //cout << "aff_lin_par_eq::set_solution:: true_forward_solution_preprocessed_[" << current_timestep_ << "] = " << endl << true_forward_solution_preprocessed_[current_timestep_] << endl;
/*
            char filename[250];
            InfiniteVector<double,Index> temp_coeffs, temp_output;
            sprintf(filename,"%s%d%s","/import/shared/friedrich/source/precomputed/ydata/u_problem_13_wnum_2_primbs_d_dt_3_3_bc_ff_jmax_3_tstep_10_", current_timestep_, ".iv");
            temp_coeffs.readFromFile(filename);
            cout << "aff_lin_par_eq::set_solution:: true_forward_solution_ = " << endl << temp_coeffs << endl;
// TODO check whether preconditioning is useful here. maybe it just slows things down!
            temp_coeffs.scale(identity_,1); // should be 1 anyways ... although it isn't ...
            cout << "aff_lin_par_eq::set_solution:: after scaling = " << endl << temp_coeffs << endl;
            //APPLY_TENSOR(*assembled_problems[0], help, tolerance, wbeforeandafter, spatial_jmax, tensor_simple, true); // false
            APPLY_TENSOR(*identity_, temp_coeffs, tolerance, temp_output, spatial_jmax_, tensor_simple, true);
            cout << "aff_lin_par_eq::set_solution:: after apply = " << endl << temp_output << endl;
            temp_output.scale(identity_,1);
            cout << "aff_lin_par_eq::set_solution:: after scaling = " << endl << temp_output << endl;
            temp_output.compress(1e-13);
            cout << "aff_lin_par_eq::set_solution:: recomputed the preprocessing: temp_output = " << endl << temp_output << endl;
 * */
        }
        else
        {
            backward_solution_[NUMOFTIMESTEPS-current_timestep_] = u;
        }
    }
    
    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void 
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::set_true_forward_solution_preprocessed(const Array1D<InfiniteVector<double, Index> > &uobserved_coeffs,
                                             const double tolerance)
    {
        for (unsigned int i=0; i<= NUMOFTIMESTEPS;++i)
        {
            forward_solution_[i] = uobserved_coeffs[i];
            //forward_solution_[i] = utrue_coeffs[i];
            //forward_solution_[i].add(noise_coeffs[i]);
// TODO check whether preconditioning is useful here. maybe it just slows things down!
            forward_solution_[i].scale(identity_,1); // should be 1 anyways ... although it isn't ...
            //APPLY_TENSOR(*assembled_problems[0], help, tolerance, wbeforeandafter, spatial_jmax, tensor_simple, true); // false
            APPLY_TENSOR(*identity_, forward_solution_[i], tolerance, true_forward_solution_preprocessed_[i], spatial_jmax_, tensor_simple, true);
            true_forward_solution_preprocessed_[i].scale(identity_,1);
            true_forward_solution_preprocessed_[i].compress(1e-13);
        }
    }
    

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::evaluate_f(const double t,
	         const InfiniteVector<double,Index>& v,
	         const double tolerance,
	         InfiniteVector<double,Index>& result) const
    {
// CLEANUP
        /*
        cout << "evaluate_f called with t = " << t << " time_discretization_[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << endl;
        set <Index> temp_set;
        v.support(temp_set);
        cout << "alpe::eval_f: support of v" << endl;
        for ( typename set<Index>::const_iterator it(temp_set.begin());
                it != temp_set.end(); ++it)
            cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
        temp_set.clear();
*/
        result.clear();
        // compute Av
        // interpolate if necessary since t may not belong to time_discretization_
        int next_timestep;
        double relative_t, h;
        if (time_direction_)
        {
            next_timestep = current_timestep_ + 1;
            relative_t = t;
            h = time_discretization_[next_timestep] - time_discretization_[current_timestep_];
        }
        else
        {
            next_timestep = current_timestep_ - 1;
            relative_t = 1-t;
            h = time_discretization_[current_timestep_] - time_discretization_[next_timestep];
        }


        InfiniteVector<double,Index> input_v(v);
        if (abs(time_discretization_[current_timestep_] - relative_t) < h/10)
        { // we are at relative_t=t[i] or at least very close
// CLEANUP
            //cout << "t_k[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << ", t_k[" << next_timestep << "] = " << time_discretization_[next_timestep] << ", relative_t = " << relative_t << "=> use CURRENT" << endl;
            //cout << "alpe_eval_f input_v " << endl << input_v << endl;
            //cout << "aff_lin_par_eq::evaluate_f:: input_v = " << endl << input_v << endl;
            input_v.scale(&assembled_problems_[current_timestep_], 1); // w = Dv
            //cout << "aff_lin_par_eq::evaluate_f:: after scaling: input_v = " << endl << input_v << endl;
            //cout << "alpe_eval_f assembled_problems_["<< current_timestep_ << "].entries_cache " << endl << (*assembled_problems_[current_timestep_].entries_cache) << endl;
            //cout << "alpe_eval_f input_v_scaled " << endl << input_v << endl;
            InfiniteVector<double,Index> temp;
            //APPLY_TENSOR(*elliptic, w, tolerance, temp, jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
            APPLY_TENSOR(assembled_problems_[current_timestep_], input_v, tolerance, temp, spatial_jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
            //cout << "aff_lin_par_eq::evaluate_f:: APPLY => temp = " << endl << temp << endl;
            temp.scale(&assembled_problems_[current_timestep_], 1);
            temp.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            //cout << "aff_lin_par_eq::evaluate_f:: after scaling: temp = " << endl << temp << endl;

            //     // multiply with inverse primal gramian (i.e., switch from dual to primal basis)
            //     G.set_rhs(temp);
            //     CDD1_SOLVE(GC, tolerance, result, jmax_);
            result = temp;
            //cout << "aff_lin_par_eq::evaluate_f:: result = " << endl << result << endl;
            //cout << "aff_lin_par_eq::evaluate_f:: forward_solution_preprocessed_[" << current_timestep_ << "] = " << endl << forward_solution_preprocessed_[current_timestep_] << endl;

            //result.scale(&assembled_problems_[current_timestep_], 1);
            //result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            //cout << "alpe::eval_f A1 " << endl << result << endl;
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[current_timestep_]);
                result.subtract(true_forward_solution_preprocessed_[current_timestep_]);
            }

            //cout << "aff_lin_par_eq::evaluate_f:: result = " << endl << result << endl;

            //cout << "alpe::eval_f A2 " << endl << result << endl;
/*
            result.support(temp_set);
            cout << "alpe::eval_f: support of result A" << endl;
            for ( typename set< Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            temp_set.clear();
 */
        }
        else if (abs(time_discretization_[next_timestep] - relative_t) < h/10)
        { // we are at relative_t=t[i+1] or at least very close
// CLEANUP
            //cout << "t_k[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << ", t_k[" << next_timestep << "] = " << time_discretization_[next_timestep] << ", relative_t = " << relative_t << "=> use NEXT" << endl;

            input_v.scale(&assembled_problems_[next_timestep], 1); // w = Dv
            APPLY_TENSOR(assembled_problems_[next_timestep], input_v, tolerance, result, spatial_jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
            result.scale(&assembled_problems_[next_timestep], 1);
            result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[next_timestep]);
                result.subtract(true_forward_solution_preprocessed_[next_timestep]);
            }
/*
            result.support(temp_set);
            cout << "alpe::eval_f: support of result B" << endl;
            for ( typename set< Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            temp_set.clear();
 */
        }
        else
        { // we are at t[i] < relative_t < t[i+1]; interpolate!
            // t = c*t[i] + (1-c)*t[i+1]; then c == (t-t[i]) / (t[i+1]-t[i])
// CLEANUP
            //cout << "t_k[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << ", t_k[" << next_timestep << "] = " << time_discretization_[next_timestep] << ", relative_t = " << relative_t << "=> INTERPOLATE" << endl;

            double c = (relative_t-time_discretization_[current_timestep_])/h * (time_direction_ ? 1:-1);
            input_v.scale(&assembled_problems_[current_timestep_], 1); // w = Dv
            APPLY_TENSOR(assembled_problems_[current_timestep_], input_v, tolerance, result, spatial_jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
            result.scale(&assembled_problems_[current_timestep_], 1);
            result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[current_timestep_]);
                result.subtract(true_forward_solution_preprocessed_[current_timestep_]);
            }
/*
            result.support(temp_set);
            cout << "alpe::eval_f: support of result C" << endl;
            for ( typename set< Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            temp_set.clear();
*/
// TODO: the next line was disabled, yet the program worked?!            
            input_v = v;
            InfiniteVector<double, Index> output_v;
            input_v.scale(&assembled_problems_[next_timestep], 1); // w = Dv
            APPLY_TENSOR(assembled_problems_[next_timestep], input_v, tolerance, output_v, spatial_jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
            output_v.scale(&assembled_problems_[next_timestep], 1);
            output_v.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            if (!time_direction_)
            {
                output_v.add(forward_solution_preprocessed_[next_timestep]);
                output_v.subtract(true_forward_solution_preprocessed_[next_timestep]);
            }
            result.scale(c);
            result.add(1-c,output_v);
/*
            result.support(temp_set);
            cout << "alpe::eval_f: support of result D" << endl;
            for ( typename set< Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            temp_set.clear();
 */
        }

        // add constant driving term
        if (time_direction_)
        {
            //cout << "aff_lin_par_eq::evaluate_f:: f_ = " << endl << f_ << endl;
            result.add(f_);
            //cout << "aff_lin_par_eq::evaluate_f:: result = " << endl << result << endl;
            //cout << "alpe::eval_f Beitrag von f_" << endl << f_ << endl;
/*
            f_.support(temp_set);
            cout << "alpe::eval_f: support of result f_" << endl;
            for ( typename set< Index>::const_iterator it(temp_set.begin());
                    it != temp_set.end(); ++it)
                cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
            temp_set.clear();
 */
        }
        //cout << "alpe::eval_f B " << endl << result << endl;
/*
        result.support(temp_set);
        cout << "alpe::eval_f: support of result E" << endl;
        for ( typename set< Index>::const_iterator it(temp_set.begin());
                it != temp_set.end(); ++it)
            cout << "*it = " << *it << "(*it).number() = " << (*it).number() << endl;
        temp_set.clear();
 */
    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::evaluate_ft(const double t,
                  const InfiniteVector<double,Index>& v,
		  const double tolerance,
		  InfiniteVector<double,Index>& result) const
    {
        if (abs(time_discretization_[current_timestep_] - (time_direction_?t:(1-t))) > 1e-6)
        {
            cout << "evaluate_ft:: time step does not seem to lie in the time discretization. time_discretization_[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
            abort();
        }

        result.clear();
        // time forward case: f=1/2 => f'=0 (we are done)
        // old: backward direction: f consists of the solution of the time-forward case
        // new: backward direction: f consists of the solution of the time-forward case - the observed data
        if (!time_direction_) // reverse direction: rhs consists of the solution of the time-forward case - ydata
        {
            switch (current_timestep_)
            {
                case 0:
                    result = forward_solution_preprocessed_[0];
                    result.subtract(true_forward_solution_preprocessed_[0]);
                    result.subtract(forward_solution_preprocessed_[1]); // assume that time_discretization_ has more than one element
                    result.add(true_forward_solution_preprocessed_[1]);
                    result.scale(1./(time_discretization_[1]-time_discretization_[0])); // always: time_discretization_[0] == 0?
                    break;
                case NUMOFTIMESTEPS: // time_discretization_.length()
                    result = forward_solution_preprocessed_[NUMOFTIMESTEPS-1];
                    result.subtract(true_forward_solution_preprocessed_[NUMOFTIMESTEPS-1]);
                    result.subtract(forward_solution_preprocessed_[NUMOFTIMESTEPS]);
                    result.add(true_forward_solution_preprocessed_[NUMOFTIMESTEPS]);
                    result.scale(1./(time_discretization_[NUMOFTIMESTEPS]-time_discretization_[NUMOFTIMESTEPS-1]));
                    break;
                default:
                    result = forward_solution_preprocessed_[current_timestep_-1];
                    result.subtract(true_forward_solution_preprocessed_[current_timestep_-1]);
                    result.subtract(forward_solution_preprocessed_[current_timestep_+1]);
                    result.add(true_forward_solution_preprocessed_[current_timestep_+1]);
                    result.scale(1./(time_discretization_[current_timestep_+1]-time_discretization_[current_timestep_-1]));
                    break;
            }
        }
    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::solve_ROW_stage_equation(const double t,
			       const InfiniteVector<double,Index>& v,
			       const double alpha,
			       const InfiniteVector<double,Index>& y,
			       const double tolerance,
			       InfiniteVector<double,Index>& result) const
    {
// CLEANUP
        if (time_direction_) // reverse direction: rhs consists of the solution of the time-forward case
        {
            if (abs(time_discretization_[current_timestep_] - t) > 1e-6)
            {
                //double unused = time_discretization_[current_timestep_];
                //cout << unused << endl;
                cout << "solve_ROW_stage_equation:: time step does not seem to lie in the time discretization. time_discretization_[current_timestep_] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
                abort();
            }
        }
        else
        {
            // assume that t is a time step from time_discretization_
            if (abs(time_discretization_[current_timestep_] - (1-t)) > 1e-6)
            {
                cout << "solve_ROW_stage_equation:: time step does not seem to lie in the time discretization. time_discretization_[current_timestep_] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
                abort();
            }
        }
// TODO y precond? what does the helper need, what does this method assume? -> write to comments!
        AffLinParEq_Precomputed_StageEquationHelper<CompressedProblemFromMatrix<PROBLEM > > helper(alpha, &assembled_problems_[current_timestep_], identity_, y);
        /*
        cout << "alpha = " << alpha << endl;
        cout << "identity_.entries_cache" << endl << *identity_->entries_cache << endl;
        cout << "raw_cache_[0] (should be identity)"<< endl << raw_cache_[0] << endl;
        cout << "scaled_laplacian_matrix_" << endl << scaled_laplacian_matrix_ << endl;
        cout << "assembled_problems_[" << current_timestep_<< "].entries_cache (should be composed out of the above)" << endl << *assembled_problems_[current_timestep_].entries_cache << endl;
        cout << "helper.a(first_gen,first_gen) = " << helper.a(identity_->basis().first_generator(),identity_->basis().first_generator()) << endl;
        cout << "helper.matrix= " << endl;
        for (unsigned int i=0; i<18;++i)
        {
            for (unsigned int j=0; j<18;++j)
            {
                cout << helper.a(identity_->basis().get_wavelet(i),identity_->basis().get_wavelet(j)) << "  ";
            }
            cout << endl;
        }

        cout << "helper.norm_A() " << helper.norm_A() << endl;
        cout << "helper.norm_Ainv() " << helper.norm_Ainv() << endl;
        cout << "helper.alphak(10) " << helper.alphak(10) << endl;
        cout << "helper.f" << endl;
        for (unsigned int i=0; i<18;++i)
        {
            cout << helper.f(identity_->basis().get_wavelet(i)) << " " << endl;
        }
        InfiniteVector<double,Index> temp_iv;
        helper.RHS(1e-5,temp_iv);
        cout << "helper.RHS(1e-5)" << endl << temp_iv << endl;
        cout << "helper.F_norm = " << helper.F_norm() << endl;

        cout << "aborting ..." << endl;
        abort();
        */

        CDD1_SOLVE(helper, tolerance, result, spatial_jmax_); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
        result.scale(&helper, -1); // Dx -> x

    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::preprocess_rhs_share(InfiniteVector<double,Index>& wbeforeandafter,
                           const double tolerance) const
    {
        InfiniteVector<double,Index> help(wbeforeandafter);
// TODO check whether preconditioning is useful here. maybe it just slows things down!
        help.scale(identity_,1); // should be 1 anyways ... although it isn't ...
        APPLY_TENSOR(*identity_, help, tolerance, wbeforeandafter, spatial_jmax_, tensor_simple, true); // or "false" and comment the lines before and after
        wbeforeandafter.scale(identity_,1);

        /*
        Vector<double> help(par_dof);
        for (typename InfiniteVector<double,Index>::const_iterator it(wbeforeandafter.begin()); it != wbeforeandafter.end(); ++it)
        {
            help[it.index().number()] = *it;
        }
        help = raw_cache[0].apply(help);
// TODO SPEEDUP
        
        //unsigned int k=0;
        //++k in die schleife
        
        set<Index> Lambda_IndexSet;
        for (Index lambda(problem.basis().first_generator(eq.basis().j0())), lambda_end(problem.basis().last_wavelet(sol_jmax));; ++lambda)
        {
            Lambda_IndexSet.insert(lambda);
            if (lambda == lambda_end) break;
            // if (k == sol_dof) break;
        }
        unsigned int id = 0;
        wbeforeandafter.clear();
        for (typename set<Index>::const_iterator it = Lambda_IndexSet.begin(), itend = Lambda_IndexSet.end(); it != itend; ++it, ++id)
        {
            wbeforeandafter.set_coefficient(*it, xk[id]);
        }
        */
    }

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void
    AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS>
    ::assemble_W(const Array1D<InfiniteVector<double,int> >& new_w)
    {
        // assemble the matrix
        // we aim for A = -delta-<W\Psi,\Psi> // look out for the scale(-1) below, as we actually solve u'=(delta+W)u+f
        //SparseMatrix<double> W (laplacian_);
        
        for (unsigned int t_step=0; t_step <= NUMOFTIMESTEPS; ++t_step)
        {

            //assembled_problems_[t_step].entries_cache = scaled_laplacian_matrix_;
            assembled_matrices_[t_step] = scaled_laplacian_matrix_;
            //SparseMatrix<double> matrix(scaled_laplacian_matrix_);
// CLEANUP
            //cout << "assemble_W::  t_step = " << t_step << endl << " new_w[t_step] = " << endl << new_w[t_step] << endl;
            for (typename InfiniteVector<double,int>::const_iterator it(new_w[t_step].begin()), it_end(new_w[t_step].end()); it != it_end; ++it)
            {
                //cout << "assemble_W:: try to add *it = " << (*it) << " it.index() = " << it.index() << endl;
                //assembled_matrices[t_step] -= (*it) * raw_cache[it.index()]);
                if ((raw_cache_[it.index()].row_dimension() == 1) && (raw_cache_[it.index()].get_entry(0,0) == 0))
                //if (raw_cache_[it.index()].row_dimension() == 1)
                //if (raw_cache_[it.index()] == 0)
                {  // load matrix from file
                    char matrix_filename[250];
                    //haar_wav_2_gramian_primbs_d_dt_3_3_bc_tf_jmax_3_npcd
                    sprintf(matrix_filename, "%s%d%s",haar_gramian_filename_start_, it.index(), haar_gramian_filename_end_);
                    //cout << "assemble_W:: load and store. t_step = " << t_step << " it.index() = " << it.index() << endl << "   filename = " << matrix_filename << endl;
                    //cout << "assemble_W:: pre matlab_input: this is stored in memory: raw_cache_["<< it.index() << "] =" << endl << (raw_cache_[it.index()]) << endl;
                    //cout << "assemble_W:: A =" << endl << raw_cache_[it.index()].row_dimension() << endl;
                    raw_cache_[it.index()].matlab_input(matrix_filename);
                    //cout << "B " << raw_cache_[it.index()].row_dimension() << endl;
                    //cout << "assemble_W:: this is stored in memory: raw_cache_["<< it.index() << "] = " << (raw_cache_[it.index()]) << endl;
                }
                /*
                else
                {
                    //cout << "matrix corresponding to it.index() = " << it.index() << " seems to be already loaded to raw_cache" << endl;                    
                    //cout << "assemble_W:: this is stored in memory: raw_cache_["<< it.index() << "] =" << endl << (raw_cache_[it.index()]) << endl;
                    //cout << "assemble_W:: A =" << endl << raw_cache_[it.index()].row_dimension() << endl;
                }
                 */
                //cout << "C " << raw_cache_[it.index()].row_dimension() << endl;
                //cout << "assemble_W:: this is stored in memory: raw_cache_["<< it.index() << "] = " << (raw_cache_[it.index()]) << endl;
                //cout << "matrix.rd,cd,raw_cache["<<it.index()<<"].rd,cd = " << matrix.row_dimension() << " " << matrix.column_dimension() << " " << raw_cache_[it.index()].row_dimension() << " " << raw_cache_[it.index()].column_dimension() << endl;
                //assembled_problems_[t_step].entries_cache.add(-(*it),raw_cache_[it.index()]);
                assembled_matrices_[t_step].add((*it),raw_cache_[it.index()]);
            }
            assembled_problems_[t_step].entries_cache = &assembled_matrices_[t_step];
// TODO either include a "problem" in constructor, or write a CompressedProblem that does not need an underlying PROBLEM
// TODO check the meaning of 5.5 and 35.0
        }
    }    
    
/* ----------------------------  */
    
#if 1   
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::AffLinParEq_qtbasis(QTBASIS* qtbasis,
                              //const InfiniteVector<double,int>& u0,
                              Function<QTBASIS::space_dimension,double>* uexact_function,
                              const char* u0_filename,
                              const double d,
                              //const FixedArray1D<InfiniteVector<double, int>, NUMOFTIMESTEPS+1 >& w,
                              const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1>& w,
                              //const Array1D<double>& time_discretization,
                              const FixedArray1D<double, NUMOFTIMESTEPS +1 >& time_discretization,
                              //const char* f_filename,
                              //const char* haar_gramian_filename_start,
                              //const char* haar_gramian_filename_end,
                              //const char* laplacian_filename,
                              //const int par_jmax,
                              //const int spatial_jmax,
                              const Function<QTBASIS::space_dimension>* onehalf_function, // from CQTProblem: rhs as a function
                              const char* onehalf_filename, 
                              const char* onehalf_precond_gramian_filename, 
                              const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename,
                              const double normA, // from CQTProblem: estimate for \|dI+W\|
                              const double normAinv, // from CQTProblem:  estimate for \|(dI+W)^{-1}\|
                              const double tolerance_APPLY
                              )
    : CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT,1>(qtbasis, 
                      //Array1D<FixedMatrix<double, 1> > (), 
                      //Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > (),
                      //NULL, onehalf_filename, 
                      normA, normAinv),
      d_(d), time_discretization_(time_discretization), 
            current_timestep_(0),
            problem_mode_(0),
            qtbasis_(qtbasis),
            onehalf_function_(onehalf_function)
    {
            //haar_gramian_filename_start_(haar_gramian_filename_start), 
            //haar_gramian_filename_end_(haar_gramian_filename_end),
            //spatial_jmax_(spatial_jmax),
        time_direction_ = true;
        
        compute_D_gramian();
        set_W(w);
        
        if (onehalf_function_ != 0)
        {
            compute_rhs(onehalf_filename, onehalf_precond_gramian_filename, onehalf_mode_filename);
        }
        else
        {
            load_rhs(onehalf_filename, onehalf_precond_gramian_filename, onehalf_mode_filename);
            
        }
        u0_.clear();
        
        //Vector<double> temp_u0;
        
        if (uexact_function != NULL)
        {
            cout << "AffLinParEq_qtbasis:: use CDD1 to compute u0T (u0 w.r.t. primal basis) ..." << endl;
            uexact_function->set_time(0);
            InfiniteVector<double,int> u0T;
            qtbasis_->expand(uexact_function, false, u0T);
            current_rhs_ = &u0T;
            current_rhs_l2_norm_ = l2_norm(u0T);
            
            Array1D<std::pair<int,double> > another_rhs_sorted;
            another_rhs_sorted.resize(0); // clear eventual old values
            another_rhs_sorted.resize(u0T.size());
            unsigned int id = 0;
            for (typename InfiniteVector<double,int>::const_iterator it(u0T.begin()), itend(u0T.end());
                    it != itend; ++it, ++id)
            {
                another_rhs_sorted[id] = std::pair<int,double>(it.index(), *it);
            }
            sort(another_rhs_sorted.begin(), another_rhs_sorted.end(), typename InfiniteVector<double,int>::decreasing_order());
            current_rhs_sorted_ = &another_rhs_sorted;

            problem_mode_ = 2;
            CDD1_SOLVE(*this, tolerance_APPLY, u0_, qtbasis_->get_jmax(), tensor_simple);
            u0_.scale(this,-1);
            u0_.compress(1e-13);
            //cout << "after call to CDD1_SOLVE: u0_ = " << endl << u0_ << endl;
            cout << "AffLinParEq_qtbasis:: write u0T to file " << u0_filename << endl;
            writeIVToFile(u0_,u0_filename);
            cout << "AffLinParEq_qtbasis:: done" << endl;
        }
        else
        {
            cout << "AffLinParEq_qtbasis:: load u0T from file " << u0_filename << endl;
            readIVFromFile(u0_,u0_filename);
            cout << "AffLinParEq_qtbasis:: done" << endl;
        }
        set_solution(u0_, tolerance_APPLY);
    }
#endif
    
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::set_solution(const InfiniteVector<double, int>& u, const double& APPLY_TENSOR_tolerance)
    {
        if (time_direction_ == true)
        {
            unsigned int oldmode = problem_mode_;
            problem_mode_ = 2;
            forward_solution_[current_timestep_] = u;
            InfiniteVector<double,int> help2(u);
            APPLY_TENSOR(*this, help2, APPLY_TENSOR_tolerance, forward_solution_preprocessed_[current_timestep_], qtbasis_->get_jmax(), tensor_simple, false);
            forward_solution_preprocessed_[current_timestep_].compress(1e-13);
            problem_mode_ = oldmode;
        }
        else
        {
            backward_solution_[NUMOFTIMESTEPS-current_timestep_] = u;
        }
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::set_true_forward_solution_preprocessed(const FixedArray1D<InfiniteVector<double, int>, NUMOFTIMESTEPS+1 > &uobserved_coeffs,
                                             const double APPLY_TENSOR_tolerance)
    {
        unsigned int old_mode = problem_mode_;
        problem_mode_ = 2;
        for (unsigned int i=0; i<= NUMOFTIMESTEPS;++i)
        {
            //forward_solution_[i] = uobserved_coeffs[i];
            //forward_solution_[i].scale(this,1); // should be 1 anyways ... although it isn't ...
            APPLY_TENSOR(*this, uobserved_coeffs[i], APPLY_TENSOR_tolerance, true_forward_solution_preprocessed_[i], qtbasis_->get_jmax(), tensor_simple, false);
            //true_forward_solution_preprocessed_[i].scale(this,1);
            true_forward_solution_preprocessed_[i].compress(1e-13);
        }
        problem_mode_ = old_mode;
    }

#if 1    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::evaluate_f(const double t,
	         const InfiniteVector<double,int>& v,
	         const double tolerance,
	         InfiniteVector<double,int>& result)
    {
        unsigned int oldmode = problem_mode_;
        problem_mode_ = 0;
        result.clear();
        // compute Av
        // interpolate if necessary since t may not belong to time_discretization_
        int next_timestep;
        double relative_t, h;
        if (time_direction_)
        {
            next_timestep = current_timestep_ + 1;
            relative_t = t;
            cout << " time_discretization_ = " << time_discretization_ << " current_timestep_ = " << current_timestep_ << "; next_timestep = " << next_timestep << endl;
            h = time_discretization_[next_timestep] - time_discretization_[current_timestep_];
        }
        else
        {
            next_timestep = current_timestep_ - 1;
            relative_t = 1-t;
            h = time_discretization_[current_timestep_] - time_discretization_[next_timestep];
        }
        
        if (abs(time_discretization_[current_timestep_] - relative_t) < h/10)
        { // we are at relative_t=t[i] or at least very close
//            input_v.scale(this, 1); // w = Dv
//            APPLY_TENSOR(*this, input_v, tolerance, result, qtbasis_->get_jmax(), tensor_simple); // yields -D^{-1}AD^{-1}w
//            result.scale(this, 1);
            
            InfiniteVector<double,int> temp2;
            APPLY_TENSOR(*this, v, tolerance, result, qtbasis_->get_jmax(), tensor_simple,false); // yields -D^{-1}AD^{-1}w
            
            result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            //result = temp;
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[current_timestep_]);
                result.subtract(true_forward_solution_preprocessed_[current_timestep_]);
            }
        }
        else if (abs(time_discretization_[next_timestep] - relative_t) < h/10)
        { // we are at relative_t=t[i+1] or at least very close
            unsigned int old_timestep(current_timestep_);
            current_timestep_ = next_timestep;
            InfiniteVector<double,int> input_v(v);
            input_v.scale(this, 1); // w = Dv
            APPLY_TENSOR(*this, input_v, tolerance, result, qtbasis_->get_jmax(), tensor_simple); // yields -D^{-1}AD^{-1}w
            result.scale(this, 1);
            
            InfiniteVector<double,int> temp2;
            APPLY_TENSOR(*this, v, tolerance, temp2, qtbasis_->get_jmax(), tensor_simple,false); // yields -D^{-1}AD^{-1}w
            cout << "zum Vergleich:\n result = " << result << "\n temp2 = " << temp2 << "\n l2_norm(result) = " << l2_norm(result) << "; l2_norm(temp2) = " << l2_norm(temp2) << endl;
            
            result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            current_timestep_ = old_timestep;
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[next_timestep]);
                result.subtract(true_forward_solution_preprocessed_[next_timestep]);
            }
        }
        else
        { // we are at t[i] < relative_t < t[i+1]; interpolate!
            // t = c*t[i] + (1-c)*t[i+1]; then c == (t-t[i]) / (t[i+1]-t[i])
            double c = (relative_t-time_discretization_[current_timestep_])/h * (time_direction_ ? 1:-1);
            InfiniteVector<double,int> input_v(v);
            input_v.scale(this, 1); // w = Dv
            APPLY_TENSOR(*this, input_v, tolerance, result, qtbasis_->get_jmax(), tensor_simple); // yields -D^{-1}AD^{-1}w
            result.scale(this, 1);
            
            InfiniteVector<double,int> temp2;
            APPLY_TENSOR(*this, v, tolerance, temp2, qtbasis_->get_jmax(), tensor_simple,false); // yields -D^{-1}AD^{-1}w
            cout << "zum Vergleich:\n result = " << result << "\n temp2 = " << temp2 << "\n l2_norm(result) = " << l2_norm(result) << "; l2_norm(temp2) = " << l2_norm(temp2) << endl;
            
            result.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            if (!time_direction_)
            {
                result.add(forward_solution_preprocessed_[current_timestep_]);
                result.subtract(true_forward_solution_preprocessed_[current_timestep_]);
            }
            InfiniteVector<double, int> output_v;
            unsigned int old_timestep(current_timestep_);
            current_timestep_ = next_timestep;
            input_v = v;
            input_v.scale(this, 1); // w = Dv
            APPLY_TENSOR(*this, input_v, tolerance, output_v, qtbasis_->get_jmax(), tensor_simple); // yields -D^{-1}AD^{-1}w
            output_v.scale(this, 1);
            
            APPLY_TENSOR(*this, v, tolerance, temp2, qtbasis_->get_jmax(), tensor_simple,false); // yields -D^{-1}AD^{-1}w
            cout << "zum Vergleich:\n output_v = " << output_v << "\n temp2 = " << temp2 << "\n l2_norm(output_v) = " << l2_norm(output_v) << "; l2_norm(temp2) = " << l2_norm(temp2) << endl;
            
            output_v.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
            current_timestep_ = old_timestep;
            if (!time_direction_)
            {
                output_v.add(forward_solution_preprocessed_[next_timestep]);
                output_v.subtract(true_forward_solution_preprocessed_[next_timestep]);
            }
            result.scale(c);
            result.add(1-c,output_v);
        }
        // add constant driving term
        if (time_direction_)
        {
            result.add(onehalf_coeffs_);
        }
        problem_mode_ = oldmode;
    }
#endif
#if 1    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::evaluate_ft(const double t,
                  const InfiniteVector<double,int>& v,
		  const double tolerance,
		  InfiniteVector<double,int>& result)
    {
        if (abs(time_discretization_[current_timestep_] - (time_direction_?t:(1-t))) > 1e-6)
        {
            cout << "evaluate_ft:: time step does not seem to lie in the time discretization. time_discretization_[" << current_timestep_ << "] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
            abort();
        }
        result.clear();
        // time forward case: f=1/2 => f'=0 (we are done)
        // backward direction: f consists of the solution of the time-forward case - the observed data
        if (!time_direction_) // reverse direction: rhs consists of the solution of the time-forward case - ydata
        {
            switch (current_timestep_)
            {
                case 0:
                    result = forward_solution_preprocessed_[0];
                    result.subtract(true_forward_solution_preprocessed_[0]);
                    result.subtract(forward_solution_preprocessed_[1]); // assume that time_discretization_ has more than one element
                    result.add(true_forward_solution_preprocessed_[1]);
                    result.scale(1./(time_discretization_[1]-time_discretization_[0])); // always: time_discretization_[0] == 0?
                    break;
                case NUMOFTIMESTEPS: // time_discretization_.length()
                    result = forward_solution_preprocessed_[NUMOFTIMESTEPS-1];
                    result.subtract(true_forward_solution_preprocessed_[NUMOFTIMESTEPS-1]);
                    result.subtract(forward_solution_preprocessed_[NUMOFTIMESTEPS]);
                    result.add(true_forward_solution_preprocessed_[NUMOFTIMESTEPS]);
                    result.scale(1./(time_discretization_[NUMOFTIMESTEPS]-time_discretization_[NUMOFTIMESTEPS-1]));
                    break;
                default:
                    result = forward_solution_preprocessed_[current_timestep_-1];
                    result.subtract(true_forward_solution_preprocessed_[current_timestep_-1]);
                    result.subtract(forward_solution_preprocessed_[current_timestep_+1]);
                    result.add(true_forward_solution_preprocessed_[current_timestep_+1]);
                    result.scale(1./(time_discretization_[current_timestep_+1]-time_discretization_[current_timestep_-1]));
                    break;
            }
        }
    }
#endif
#if 1
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::solve_ROW_stage_equation(const double t,
			       const InfiniteVector<double,int>& v,
			       const double alpha,
			       const InfiniteVector<double,int>& y,
			       const double tolerance,
			       InfiniteVector<double,int>& result)
    {
// CLEANUP
        if (time_direction_) // reverse direction: rhs consists of the solution of the time-forward case
        {
            if (abs(time_discretization_[current_timestep_] - t) > 1e-6)
            {
                //double unused = time_discretization_[current_timestep_];
                //cout << unused << endl;
                cout << "solve_ROW_stage_equation:: time step does not seem to lie in the time discretization. time_discretization_[current_timestep_] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
                abort();
            }
        }
        else
        {
            // assume that t is a time step from time_discretization_
            if (abs(time_discretization_[current_timestep_] - (1-t)) > 1e-6)
            {
                cout << "solve_ROW_stage_equation:: time step does not seem to lie in the time discretization. time_discretization_[current_timestep_] = " << time_discretization_[current_timestep_] << " t = " << t << endl;
                abort();
            }
        }
        //AffLinParEq_Precomputed_StageEquationHelper<CompressedProblemFromMatrix<PROBLEM > > helper(alpha, &assembled_problems_[current_timestep_], identity_, y);
        unsigned int old_mode(problem_mode_);
        problem_mode_ = 1;
        //current_rhs_ = &y;
        current_rhs_l2_norm_ = l2_norm(y);
        InfiniteVector<double,int> another_rhs = y;
        Array1D<std::pair<int,double> > another_rhs_sorted;
        another_rhs_sorted.resize(0); // clear eventual old values
        another_rhs_sorted.resize(y.size());
        unsigned int id = 0;
        for (typename InfiniteVector<double,int>::const_iterator it(y.begin()), itend(y.end());
                it != itend; ++it, ++id)
        {
            another_rhs_sorted[id] = std::pair<int,double>(it.index(), *it);
        }
        sort(another_rhs_sorted.begin(), another_rhs_sorted.end(), typename InfiniteVector<double,int>::decreasing_order());
                
        current_rhs_ = &another_rhs;
        current_rhs_sorted_ = &another_rhs_sorted;
        mode_one_alpha_ = alpha;
        
        CDD1_SOLVE(*this, tolerance, result, qtbasis_->get_jmax(),tensor_simple); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
        result.scale(this, -1); // Dx -> x
        problem_mode_ = old_mode;
    }
#endif
#if 1    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::preprocess_rhs_share(InfiniteVector<double,int>& wbeforeandafter,
                           const double tolerance) const
    {
        InfiniteVector<double,int> help(wbeforeandafter);
        
        unsigned int old_mode(problem_mode_);
        problem_mode_ = 2;
        
        InfiniteVector<double,int> temp_iv;
        APPLY_TENSOR(*this, help, tolerance, temp_iv, qtbasis_->get_jmax(), tensor_simple, false); // or "false" and comment the lines before and after
        
        
        help.scale(this,1); // should be 1 anyways ... although it isn't ...
        APPLY_TENSOR(*this, help, tolerance, wbeforeandafter, qtbasis_->get_jmax(), tensor_simple, true); // or "false" and comment the lines before and after
        wbeforeandafter.scale(this,1);
        
        cout << "whats the difference?"<< endl << "temp_iv = " << temp_iv <<"\n wbeforeandafter = " << wbeforeandafter << "; l2_norm(temp_iv) = " << l2_norm(temp_iv) << "; l2_norm(wbeforeandafter) = " << l2_norm(wbeforeandafter) << endl;
        problem_mode_ = old_mode;
    }
#endif
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    double
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::a(const int munum,
            const int nunum)
    {
        
        //für festes nu, berechne alle mu2 mit selbem Level wie mu die nu überlappen
        double r = 0;
        Index mu(qtbasis_->get_wavelet(munum));
        Index nu (qtbasis_->get_wavelet(nunum));

        const MultiIndex<int, DIM> nu_j(nu.j()), nu_e(nu.e()), nu_k(nu.k()),mu_j(mu.j()), mu_e(mu.e()), mu_k(mu.k());
        const unsigned int nu_p(nu.p()), mu_p(mu.p());

        MultiIndex<unsigned int, DIM> intinfo;
        //bool temp_b(qtbasis_->get_LMR_info(nunum,mu_j,mu_p, intinfo));
        bool temp_b(qtbasis_->get_LMR_info(nunum,mu_p, intinfo));
        if (!temp_b) 
            return 0;
        FixedArray1D<entries,DIM> integralshares; // store the value of the one dimensional integrals mu_i nu_i against Haar generators
        int kmingeni, kmaxgeni, kminwavi, kmaxwavi;
        unsigned int nui_basisnum, mui_basisnum;
        bool mu_min_type_i;
        for (unsigned int i=0; i<DIM; i++)
        {
            nui_basisnum = (((qtbasis_->get_bc()[nu_p][2*i])?0:2) + ((qtbasis_->get_bc()[nu_p][2*i+1])?0:1));
            mui_basisnum = (((qtbasis_->get_bc()[mu_p][2*i])?0:2) + ((qtbasis_->get_bc()[mu_p][2*i+1])?0:1));
            mu_min_type_i = (mu_j[i]==qtbasis_->j0()[mu_p][i]) ? true:false;
            ColumnCache* cachePointer;
            Column* cachePointerGen;
            switch (intinfo[i])
            {
                case 0:
                case 4:
                case 8: // LL MM RR -> type 1
                    cachePointer = & this->typeIcache_[4*nui_basisnum + mui_basisnum];
                    if (mu_min_type_i)
                        cachePointerGen = & this->typeIcachedGens_[4*nui_basisnum + mui_basisnum];
                    break;
                case 1:
                case 7: 
                case 2:
                case 6: // LM RM LR RL -> type 2
                    cachePointer = &this->typeIIcache_[4*(nui_basisnum-1) + mui_basisnum];
                    if (mu_min_type_i)
                        cachePointerGen = & this->typeIIcachedGens_[4*(nui_basisnum-1) + mui_basisnum];
                    break;
                case 3: // ML -> type 3             
                    // EXTREME CAUTION:: mui_basisnum - (intinfo[i] == 5)?0:1 produces garbage!
                    assert (mui_basisnum != 1);
                    cachePointer = &this->typeIIIcache_[4*nui_basisnum + mui_basisnum - 1];
                    if (mu_min_type_i)
                        cachePointerGen = & this->typeIIIcachedGens_[4*nui_basisnum + mui_basisnum - 1];
                    break;
                case 5: // MR -> type 3                
                    assert (mui_basisnum != 2);
                    cachePointer = &this->typeIIIcache_[4*nui_basisnum + ((mui_basisnum == 1)?0:3)];
                    if (mu_min_type_i)
                        cachePointerGen = & this->typeIIIcachedGens_[4*nui_basisnum + ((mui_basisnum == 1)?0:3)];
                    break;
                default:
                    abort();
            }
            // search for column 'mu_i'
            typename QTBASIS::IntervalBasis::Index nui(nu_j[i],nu_e[i],nu_k[i],qtbasis_->get_bases_infact()[nui_basisnum]);
            unsigned int nui_num = nui.number();
            typename ColumnCache::iterator col_lb(cachePointer->lower_bound(nui_num));
            typename ColumnCache::iterator col_it(col_lb);
            if (col_lb == cachePointer->end() ||
                cachePointer->key_comp()(nui_num, col_lb->first))
            {
                // the nui-th column has not been requested so far
                // insert a new column and continue with the blocks
                typedef typename ColumnCache::value_type value_type;
                col_it = cachePointer->insert(col_lb, value_type(nui_num, Column()));
            }
            // col_it points to the column of \nu_i
            Column& col(col_it->second);
            // check whether the level 'mu_i' belongs to has already been calculated
            int mui_levelnum (mu_j[i] - qtbasis_->j0()[mu_p][i]);
            typename Column::iterator lb(col.lower_bound(mui_levelnum));
            typename Column::iterator it(lb);
            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                // no entries have ever been computed for this column and this level
                // insert a new block, 
                // and then compute the whole level block 
                //      (\int psi_mui psi_nui)_mui, (\int psi_mui' psi_nui')_mui
                // for all mui that intersect nui. 
                // if mui is on the minimal level we also need to add a new Block() to
                // the generator integral cache, i.e., we need to 
                // integrate against all generators on the lowest level
                typedef typename Column::value_type value_type;
                it = col.insert(lb, value_type(mui_levelnum, Block()));
                Block& block(it->second);
                bool gen_intersection_i, wav_intersection_i;
                qtbasis_->get_onedim_intersections(intinfo[i],
                        nu_j[i],
                        nu_e[i],
                        nu_k[i],
                        nui_basisnum,
                        (mui_levelnum == 0),
                        mu_j[i],
                        mui_basisnum,
                        kmingeni,
                        kmaxgeni,
                        kminwavi,
                        kmaxwavi,
                        gen_intersection_i,
                        wav_intersection_i);
// cleanup                
                if ((nu_j[i] == 3) && (nu_e[i] == 0 && ((nu_k[i] == 3) && ((mu_j[i] == 3) && ((nui_basisnum == 0) && (mui_basisnum == 3) && (intinfo[i] == 3))))))
                {
                    cout << "// CACHED_QTPROBLEM:: //" << endl;
                    cout << "intinfoi = " << intinfo[i] << endl;
                    cout << "kmingeni = " << kmingeni << "; "
                        << "kmaxgeni = " << kmaxgeni << "; "
                        << "kminwavi = " << kminwavi << "; "
                        << "kmaxwavi = " << kmaxwavi << "; "
                        << "gen_intersection_i = " << gen_intersection_i << "; "
                        << "wav_intersection_i = " << wav_intersection_i << endl;
                }
                if ((nu_j[i] == 3) && (nu_e[i] == 0 && ((nu_k[i] == 3) && ((mu_j[i] == 3) && ((nui_basisnum == 0) && (mui_basisnum == 3) && (intinfo[i] == 5))))))
                {
                    cout << "// CACHED_QTPROBLEM:: //" << endl;
                    cout << "intinfoi = " << intinfo[i] << endl;
                    cout << "kmingeni = " << kmingeni << "; "
                        << "kmaxgeni = " << kmaxgeni << "; "
                        << "kminwavi = " << kminwavi << "; "
                        << "kmaxwavi = " << kmaxwavi << "; "
                        << "gen_intersection_i = " << gen_intersection_i << "; "
                        << "wav_intersection_i = " << wav_intersection_i << endl;
                }
                FixedArray1D<double,ONEDIMHAARCOUNT> gram;
                FixedArray1D<double,1> der;
                if (mu_min_type_i)
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    typename Column::iterator it2(lb2);
                    if (lb2 == cachePointerGen->end() ||
                        cachePointerGen->key_comp()(nui_num, lb2->first))
                    {
                        // no entries have ever been computed for this column and this level
                        // insert a new block, 
                        // and then compute the whole level block 
                        //      (\int psi_nui psi_mui)_nui, (\int psi_nui' psi_mui')_nui
                        // for all nui that intersect mui. 
                        // if mui is on the minimal level we also need to add a new Block() to
                        // the generator integral cache, i.e., we need to 
                        // integrate against all generators on the lowest level
                        typedef typename Column::value_type value_type;
                        it2 = cachePointerGen->insert(lb2, value_type(nui_num, Block()));
                        Block& block2(it2->second);
                        if (gen_intersection_i)
                        {
                            
                            for (int kgen = kmingeni; kgen <= kmaxgeni; ++kgen)
                            {
                                if ( (kgen == mu_k[i]) && (mu_e[i] == 0))
                                {
                                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                            mu_j[i], 0, kgen, mui_basisnum,
                                            gram,
                                            der);
                                    integralshares[i] = make_pair(gram,der);
                                    typedef typename Block::value_type value_type_block;
                                    block2.insert(block2.end(), value_type_block(kgen, integralshares[i]));
                                    temp_b = false;
//                                    cout << "inserted new generator Blocks. entry for mu is nontrivial" << endl;
                                }
                                else
                                {
                                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                            mu_j[i], 0, kgen, mui_basisnum,
                                            gram,
                                            der);
                                    
                                    typedef typename Block::value_type value_type_block;
                                    block2.insert(block2.end(), value_type_block(kgen, make_pair(gram,der)));
                                }
                            }
// CLEANUP                            
//                            if (temp_b)
//                            {
//                                cout << "inserted new generator Blocks. entry for mu is trivial. method should return 0" << endl;
//                            }
                        }
                    }
                    else
                    {
                        cout << "just inserted a wavelet block, but gen block was already existent?!!" << endl;
                        abort();
                    }
                }
                // code invariant: there exists a Block() (maybe empty) corresponding to nu and mu
                
                // if (no intersections) => current block (wavCache) remains empty
                // if we are on the lowest possible level for mui, we add an empty Block to the genCache
                // advantage: genCache[nuinum] exists and contains a block. 
                // code invariant: to access an entry: load the block and check if there is an entry (simple! same whether there is an entry or not!)
                // it would be possible to simply leave the genCache as it is.
                // advantage: nothing to do at this point
                // access to an entry: check whether there is a column in genCache. If not: value 0. If there is: get Block and take value from it.
                // This leads to 2 checks per call to a() instead of one. 
                // This is cheaper than method 1 above if nu does not intersect at all with basis functions mu.
                // However, most nus will have some intersection? In this case this variant leads to 1 additional check
                // so ... hopefully this makes the average access time to the cache faster
// cleanup                
                assert (wav_intersection_i || (!gen_intersection_i));
                if (!wav_intersection_i) return 0;
                
                // wav_intersection_i == true guarantees kminwavi <=kmaxwavi and that the values are meaningful
                for (int kwav = kminwavi; kwav <= kmaxwavi; ++kwav)
                {
                    if ( (kwav == mu_k[i]) && (mu_e[i] == 1))
                    {
                        // nu reflected?  = !(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) )
                        compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), 
                                nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                mu_j[i], 1, kwav, mui_basisnum,
                                gram,
                                der);
                        integralshares[i] = make_pair (gram,der);
                        typedef typename Block::value_type value_type_block;
                        block.insert(block.end(), value_type_block(kwav, integralshares[i]));
                        temp_b = false;
//                        cout << "inserted new wavelet Blocks. entry for mu is nontrivial" << endl;
                    }
                    else
                    {
                        compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                mu_j[i], 1, kwav, mui_basisnum,
                                gram,
                                der);
                        
                        typedef typename Block::value_type value_type_block;
                        block.insert(block.end(), value_type_block(kwav, make_pair (gram,der) ));
                    }
                }
                if (temp_b) 
                {
                    return 0;
                }
                temp_b = true;
            }
            else // column and levelblock already exist
            {
                if (mu_e[i] == 0)
                {
// CLEANUP
                    assert (mu_min_type_i);
                    typename Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    Block& block(lb2->second);
                    typename Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    typename Block::iterator block_it(block_lb);
                    // intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        return 0;
                    }
                    else 
                    {
                        integralshares[i] = block_it->second;
                    }
                }
                else
                {
                    Block& block(it->second);
                    typename Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    typename Block::iterator block_it(block_lb);
                    // level exists, but in row 'mui_k' no entry is available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        return 0;
                    }
                    else 
                    {
                        integralshares[i] = block_it->second;
                    }
                }
            }
        } // end of loop over dim
        // sum over all patches
        // include information about mu into LMR info:
        for (unsigned int i=0; i<DIM; ++i)
        {
            // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
            switch (intinfo[i])
            {
                case 0:
                    // LL: check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < qtbasis_->get_numofbw() ) ) ) )
                    {
                        // nu is not extended left
                        intinfo[i] = 4; // MM
                    }
                    break;
                case 2:
                    // LR: check whether mu is a right boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaRmax(  qtbasis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) ) ))
                    {
                        intinfo[i] = 1; // LM
                    }
                    break;
                case 3:
                    // ML : check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < qtbasis_->get_numofbw()) ) ) )
                    {
                        // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                        abort();
                    }
                    break;
                case 5:
                    // MR: check whether mu is a right boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) ) ))
                    {
                        // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                        cout << "cached_qtproblem::a ; i = " << i << "; intinfo = " << intinfo << endl;
                        cout << "nu_j = " << nu_j << "; nu_e = " << nu_e << "; nu_k = " << nu_k << "; nu_p = " << nu_p << endl;
                        cout << "mu_j = " << mu_j << "; mu_e = " << mu_e << "; mu_k = " << mu_k << "; mu_p = " << mu_p << endl;
                        cout << "qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) ) = " << qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) ) << endl;
                        cout << "(1<<mu_j[i])- 1 - qtbasis_->get_numofbw() = " << (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() << endl;
                        MultiIndex<unsigned int, DIM> intinfo2;
                        
                        //temp_b = (qtbasis_->get_LMR_info(nunum,mu_j,mu_p, intinfo2));
                        temp_b = (qtbasis_->get_LMR_info(nunum,mu_p, intinfo2));
                        cout << "nunum = " << nunum << "; mu_j = " << mu_j << "; mu_p = " << mu_p << "; intinfo2 = " << intinfo2 << "; temp_b = " << temp_b << endl;
                        temp_b = (qtbasis_->get_LMR_info(nunum,mu_p, intinfo2));
                        //a(munum,nunum);
                        abort();
                    }
                    break;
                case 6:
                    // RL : check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < qtbasis_->get_numofbw()) ) ) )
                    {
                        // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                        intinfo[i] = 7; // RM
                    }
                    break;
                case 8:
                    // RR: check whether mu is a right boundary gen/wav
// CLEANUP
                    /*
                    if (munum == 151)
                        if (nunum == 80)
                        {
                            cout << "qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) ) = " << qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) ) << endl;
                            cout << "qtbasis_->get_numofbw() = " << qtbasis_->get_numofbw() << endl;
                            cout << "(1<<mu_j[i])- 1 - qtbasis_->get_numofbw() = " << (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() << endl;
                        }
                     */
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) ) ))
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
        return compute_sum_over_patches (nu_p, intinfo, integralshares);
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    double
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::compute_sum_over_patches (const int nu_p, 
                const MultiIndex<unsigned int, DIM> intinfo, 
                const FixedArray1D<entries,DIM> integralshares) const
    {
        double result = 0;
        double temp_d1, temp_d2;
        unsigned int geometry_type;
        int centerpatchnumber;
        MultiIndex<bool, DIM> orientation;
        qtbasis_->get_intersection_geometry(nu_p, intinfo, geometry_type, centerpatchnumber, orientation);
        if (DIM == 2)
        {
            int north, east, northeast;
            if (geometry_type == 3)
            {
                if (centerpatchnumber == -1)
                {
                    north = qtbasis_->get_neighbours(nu_p,0);
                    east = qtbasis_->get_neighbours(nu_p,2);
                    geometry_type = 8;
                }
                else
                {
                    north = qtbasis_->get_neighbours(centerpatchnumber,3);
                    east = qtbasis_->get_neighbours(centerpatchnumber,1);
                    if (north != -1)
                    {
                        northeast = qtbasis_->get_neighbours(north,1);
                        if (northeast == -1)
                        {
                            geometry_type = 10;
                        }
                        else
                        {
                            if (east == -1)
                            {
                                geometry_type = 9;
                            }
                        }
                    }
                    else
                    {
                        northeast = qtbasis_->get_neighbours(east,3);
                        geometry_type = 11;
                    }
                }
            }
            switch (geometry_type)
            {
                case 0: // 1 patch
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // FixedArray1D< Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs;
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * integralshares[0].first[ONEDIMHAARCOUNT-1-x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y)) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // FixedArray1D< Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs;
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y)) * integralshares[0].first[ONEDIMHAARCOUNT-1-x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += W_times_gram[y] * integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += W_times_gram[y] * integralshares[1].first[ONEDIMHAARCOUNT-1-y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 1: // 2 patches: centerpatch and the patch right of it
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    const int rightneighbour(qtbasis_->get_neighbours(centerpatchnumber,1));
                    assert (rightneighbour != -1); // assert existence
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 2*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][rightneighbour] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // FixedArray1D< Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs;
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][rightneighbour] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][rightneighbour] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][rightneighbour] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += W_times_gram[y] * integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += W_times_gram[y] * integralshares[1].first[ONEDIMHAARCOUNT-1-y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 2 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 2: // 2 patches: centerpatch and the one above
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    const int upperneighbour(qtbasis_->get_neighbours(centerpatchnumber,3));
                    assert (upperneighbour != -1); // assert existence
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 2*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[1]) // <-- 1 !!
                            {
                                for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                                {
                                    for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[x] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][upperneighbour] (x,ONEDIMHAARCOUNT-1-y) ) * integralshares[1].first[y];
                                    }
                                }
                            }
                            else
                            {
                                for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                                {
                                    for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[x] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,ONEDIMHAARCOUNT-1-y) + Wgencoeffs_[current_timestep_][upperneighbour] (x,y) ) * integralshares[1].first[y];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[1]) // <-- 1 !!
                            {
                                for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                                {
                                    for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[x] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][upperneighbour] (x,ONEDIMHAARCOUNT-1-y) ) * integralshares[1].first[y];
                                    }
                                }
                            }
                            else
                            {
                                for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                                {
                                    for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W_times_gram[x] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,ONEDIMHAARCOUNT-1-y) + Wgencoeffs_[current_timestep_][upperneighbour] (x,y) ) * integralshares[1].first[y];
                                    }
                                }
                            }
                        }
                        if (orientation[0])
                        {
                            for (unsigned int x=0; x<ONEDIMHAARCOUNT; ++x)
                            {
                                result += W_times_gram[x] * integralshares[0].first[x];
                            }
                        }
                        else
                        {
                            for (unsigned int x=0; x<ONEDIMHAARCOUNT; ++x)
                            {
                                result += W_times_gram[x] * integralshares[0].first[ONEDIMHAARCOUNT-1-x];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 2 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 3: // 4 patches in a square. centerpatch is the lower left one. 
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W1_times_gram[eta] = 0;
                        W2_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 4*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y] )* integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y] )* integralshares[1].first[y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 4 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 8: // 3 patches, L-shaped support. Southwest is missing, i.e., centerpatch == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W1_times_gram[eta] = 0;
                        W2_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 3*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = Wgencoeffs_[current_timestep_][east] (x,y) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y] )* integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y] )* integralshares[1].first[y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 3 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 9: // 3 patches, L-shaped support. Southeast is missing, i.e., east == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W1_times_gram[eta] = 0;
                        W2_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 3*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y] )* integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y] )* integralshares[1].first[y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 3 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 10: // 3 patches, L-shaped support. Northeast is missing, i.e., northeast == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W1_times_gram[eta] = 0;
                        W2_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 3*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y)  ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y)  ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y] )* integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y] )* integralshares[1].first[y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 3 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;
                case 11: // 3 patches, L-shaped support. North is missing, i.e., north == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        W1_times_gram[eta] = 0;
                        W2_times_gram[eta] = 0;
                    }
                    temp_d1 = 0;
                    temp_d2 = 0;
                    if (problem_mode_ == 2)
                    {
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result = 4*temp_d1 * temp_d2;
                    }
                    else 
                    {
                        if (problem_mode_ == 0)
                        {
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        else
                        {
// CLEANUP
                            assert (problem_mode_ == 1);
                            if (orientation[0])
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                            else
                            {
                                for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                {
                                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                    {
                                        // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                        W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * integralshares[0].first[x];
                                        W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * integralshares[0].first[x];
                                    }
                                }
                            }
                        }
                        if (orientation[1])
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y] )* integralshares[1].first[y];
                            }
                        }
                        else
                        {
                            for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                            {
                                result += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y] )* integralshares[1].first[y];
                            }
                        }
                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                        {
                            temp_d1 += integralshares[0].first[x];
                            temp_d2 += integralshares[1].first[x];
                        }
                        result += 3 * d_ * (temp_d1 * integralshares[1].second[0]
                                        + temp_d2 * integralshares[0].second[0] );
                    }
                    return result;
                }
                    break;    
                default:
                    cout << "geometry is not implemented" << endl;
                    abort();
                    break;
            } // end of switch(geometry))
        } // end of DIM == 2
        else
        {
            // DIM == 3. Only Poisson equation implemented!
            cout << "warning: calling aff_lin_par_eq::compute_sum_over_patches with DIM >2.\nParameter reconstruction for the 3D case is not implemented.\nNeed to implement 3D Gram problem with arbitrary matrix W first." << endl;
            abort();
            return result;
        }
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::compose_wavelets(Vector<double>& w,
                const unsigned int start,
                const unsigned int blocksize,
                const double factor,
                const unsigned int nu_p,
                const FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it,
                const FixedArray1D< typename Block::iterator, DIM> block_it,
                const bool precond) const
                //const MultiIndex<unsigned int, DIM> intinfo, 
                //const FixedArray1D<entries,DIM> integralshares) const
    {
        assert ((DIM == 2) || (DIM == 3));
        assert (blocksize > 0);
        // The first OR the last wavelet in the block may have a different intinfo than the others.
        // Therefore, in this method, the geometry has to be computed exactly once or twice.
        unsigned int geometry_type;
        int centerpatchnumber;
        MultiIndex<bool, DIM> orientation;
        MultiIndex<unsigned int, DIM> intinfo;
        
        unsigned int current_offset(0);
        typename Block::iterator current_shares(block_it[DIM-1]);
        typename Array1D<unsigned int>::const_iterator current_last_intinfo (mu_adapted_intinfo_it[DIM-1]);
        bool recompute_geometry(false);
        unsigned int temp_i;
        double temp_d, temp_d1, temp_d2; //, temp_d3;
        while (current_offset < blocksize)
        {
            for (unsigned int i=0; i<DIM-1; ++i)
            {
                intinfo[i] = *mu_adapted_intinfo_it[i];
            }
            intinfo[DIM-1] = *current_last_intinfo;
            qtbasis_->get_intersection_geometry(nu_p, intinfo, geometry_type, centerpatchnumber, orientation);

            if (DIM == 2)
            {
                int north, east, northeast;
                if (geometry_type == 3)
                {
                    if (centerpatchnumber == -1)
                    {
                        north = qtbasis_->get_neighbours(nu_p,0);
                        east = qtbasis_->get_neighbours(nu_p,2);
                        geometry_type = 8;
                    }
                    else
                    {
                        north = qtbasis_->get_neighbours(centerpatchnumber,3);
                        east = qtbasis_->get_neighbours(centerpatchnumber,1);
                        if (north != -1)
                        {
                            northeast = qtbasis_->get_neighbours(north,1);
                            if (northeast == -1)
                            {
                                geometry_type = 10;
                            }
                            else
                            {
                                if (east == -1)
                                {
                                    geometry_type = 9;
                                }
                            }
                        }
                        else
                        {
                            northeast = qtbasis_->get_neighbours(east,3);
                            geometry_type = 11;
                        }
                    }
                }
                
                if (problem_mode_ == 2)
                {
                    unsigned int geometry_factor;
                    temp_d1 = 0;
                    for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                    {
                        temp_d1 += (block_it[0]->second).first[x];
                    }
                    switch (geometry_type)
                    {
                        case 0: // 1 patch
                            geometry_factor = 1;
                            break;
                        case 1: // 2 patches
                        case 2:
                            geometry_factor = 2;
                            break;
                        case 3: // 4 patches
                            geometry_factor = 4;
                            break;
                        case 8: // 3 patches
                        case 9:
                        case 10:
                        case 11:
                            geometry_factor = 3;
                            break;
                        default:
                            cout << "aff_lin_par_eq::compose_wavelets: this geometry is not implemented" << endl;
                            abort();
                            break;
                    } // end of switch(geometry))
                    while(current_offset < blocksize)
                    {
                        temp_d2 = 0;
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            temp_d2 += (current_shares->second).first[y];
                        }
                        w[start+current_offset] += precond? (factor*geometry_factor*temp_d1*temp_d2/ diagonal_cache_gramian_[start+current_offset]) : (factor*geometry_factor*temp_d1*temp_d2);
                        temp_i = *current_last_intinfo;
                        ++current_last_intinfo;
                        ++current_shares;
                        ++current_offset;
                        if (current_offset != blocksize)
                        {
                            if (temp_i != *current_last_intinfo)
                            {
                                assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                recompute_geometry = true;
                                break;
                            }
                        }
                    }
                }
                else
                { // problem_mode_ == 0 or 1
                    switch (geometry_type)
                    {
                        case 0: // 1 patch
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // FixedArray1D< Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs;
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * (block_it[0]->second).first[ONEDIMHAARCOUNT-1-x];
                                        }
                                    }
                                }
                            }
                            else
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y)) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y)) * (block_it[0]->second).first[ONEDIMHAARCOUNT-1-x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += W_times_gram[y] * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += W_times_gram[y] * (current_shares->second).first[ONEDIMHAARCOUNT-1-y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
//                                diagonal_cache_mode_0_[start+current_offset]
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 1: // 2 patches: centerpatch and the patch right of it
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            const int east(qtbasis_->get_neighbours(centerpatchnumber,1));
                            assert (east != -1); // assert existence
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // FixedArray1D< Array1D<FixedMatrix<double, _ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS+1> Wgencoeffs;
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += W_times_gram[y] * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += W_times_gram[y] * (current_shares->second).first[ONEDIMHAARCOUNT-1-y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 2 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 2: // 2 patches: centerpatch and the one above
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            const int north(qtbasis_->get_neighbours(centerpatchnumber,3));
                            assert (north != -1); // assert existence
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = Wgencoeffs_[current_timestep_][north] (x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) )* (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) )* (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) )* (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) )* (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 2 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);                                
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 3: // 4 patches in a square. centerpatch is the lower left one
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 2);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 4 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 8: // 3 patches, L-shaped support. Southwest is missing, i.e., centerpatch == -1
                        {
                            
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][east] (x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 3 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 9: // 3 patches, L-shaped support. Southeast is missing, i.e., east == -1
                            {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 3 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 10: // 3 patches, L-shaped support. Northeast is missing, i.e., northeast == -1
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = Wgencoeffs_[current_timestep_][north] (x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][north] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 3 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        case 11: // 3 patches, L-shaped support. North is missing, i.e., north == -1
                        {
                            FixedArray1D<double,ONEDIMHAARCOUNT> W1_times_gram, W2_times_gram;
                            for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                            {
                                W1_times_gram[eta] = 0;
                                W2_times_gram[eta] = 0;
                            }
                            temp_d1 = 0;
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                temp_d1 += (block_it[0]->second).first[x];
                            }
                            if (problem_mode_ == 0)
                            {
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] =  Wgencoeffs_[current_timestep_][northeast] (x,y) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            else
                            {
    // CLEANUP
                                assert (problem_mode_ == 1);
                                if (orientation[0])
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (x,y) + Wgencoeffs_[current_timestep_][east] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][northeast] (ONEDIMHAARCOUNT-1-x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                                else
                                {
                                    for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                                    {
                                        for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                        {
                                            // Wgencoeffs_[timestep][centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                            W1_times_gram[y] = (2*mode_one_alpha_ + Wgencoeffs_[current_timestep_][centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + Wgencoeffs_[current_timestep_][east] (x,y) ) * (block_it[0]->second).first[x];
                                            W2_times_gram[y] = (mode_one_alpha_ + Wgencoeffs_[current_timestep_][northeast] (x,y) ) * (block_it[0]->second).first[x];
                                        }
                                    }
                                }
                            }
                            while(current_offset < blocksize)
                            {
                                temp_d = 0;
                                if (orientation[1])
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[y] + W2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).first[y];
                                    }
                                }
                                else
                                {
                                    for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                    {
                                        temp_d += (W1_times_gram[ONEDIMHAARCOUNT-1-y] + W2_times_gram[y]) * (current_shares->second).first[y];
                                    }
                                }
                                temp_d2 = 0;
                                for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                                {
                                    temp_d2 += (current_shares->second).first[y];
                                }
                                temp_d += 3 * d_ * (temp_d1 * (current_shares->second).second[0] + temp_d2 * (block_it[0]->second).second[0]);
                                w[start+current_offset] += precond? (factor*temp_d/diagonal_cache_mode_[problem_mode_][current_timestep_][start+current_offset]) : (factor*temp_d);
                                temp_i = *current_last_intinfo;
                                ++current_last_intinfo;
                                ++current_shares;
                                ++current_offset;
                                if (current_offset != blocksize)
                                {
                                    if (temp_i != *current_last_intinfo)
                                    {
                                        assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                        recompute_geometry = true;
                                        break;
                                    }
                                }
                            }
                        }
                            break;
                        default:
                            cout << "this geometry is not implemented" << endl;
                            abort();
                            break;
                    } // end of switch(geometry))
                } // end of else problem_mode_ == 0 or == 2
            } // end of DIM == 2
            else
            {
                // DIM == 3. Only Poisson equation implemented!
                cout << "warning: calling aff_lin_par_eq::compose_wavelets with DIM >2.\nParameter reconstruction for the 3D case is not implemented.\nNeed to implement 3D Gram problem with arbitrary matrix W first." << endl;
                abort();
            }
        } // end of while (current_offset < blocksize)
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::compute_D_gramian()
    {
//        unsigned int oldmode = problem_mode_;
//        unsigned int oldtime = current_timestep_;
//        for (int mode=0; mode<2; ++mode)
//        {
//            problem_mode_ = mode;
//            for (int t=0; t<NUMOFTIMESTEPS+1; ++t)
//            {
//                current_timestep_ = t;
//                diagonal_cache_mode_[mode][current_timestep_].resize(qtbasis_->degrees_of_freedom());
//                for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
//                {
//                    this->a(i,i);
//                    CachedQTProblem<QTBASIS, ONEDIMHAARCOUNT, 1>::a(i,i);
//                    
//                    diagonal_cache_mode_[mode][current_timestep_][i] = sqrt(this->a(i,i));
//                }
//            }
//        }
        diagonal_cache_gramian_.resize(qtbasis_->degrees_of_freedom());
        problem_mode_ = 2;
        for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
        {
            diagonal_cache_gramian_[i] = sqrt(a(i,i));
        }
//        problem_mode_ = oldmode;
//        current_timestep_ = oldtime;
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void 
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::set_W(const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & W_new)
    {
        Wgencoeffs_ = W_new;
        // implement inverse haar-wav-trafo+scaling as in the commentary
        cout << "aff_lin_par_eq::set_W: Caution: matrix norms are not recomputed!" << endl;
//        unsigned int oldmode = problem_mode_;
//        unsigned int oldtime = current_timestep_;
        InfiniteVector<double,int> fhelp;
        int id;
        double temp_d;
        for (int mode=0; mode<2; ++mode)
        {
            problem_mode_ = mode;
            for (int t=0; t<NUMOFTIMESTEPS+1; ++t)
            {
                fhelp.clear();
                current_timestep_ = t;
                diagonal_cache_mode_[mode][current_timestep_].resize(qtbasis_->degrees_of_freedom());
                for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
                {
                    //temp_d = sqrt(CachedQTProblem<QTBASIS, 1, ONEDIMHAARCOUNT>::a(i,i));
                    temp_d = sqrt(this->a(i,i));
                    diagonal_cache_mode_[mode][current_timestep_][i] = temp_d;
                    temp_d = onehalf_coeffs_.get_coefficient(i)/temp_d;
                    if (fabs(temp_d)>1e-15)
                    {
                        fhelp.set_coefficient(i, temp_d);
                        onehalf_coeffs_precond_unsorted_mode_[mode][t][i] = temp_d;
                        onehalf_precond_norm_sqr_mode_[mode][t] += temp_d*temp_d;
                    }
                }
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(0); // clear eventual old values
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(fhelp.size());
                id = 0;
                for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                        it != itend; ++it, ++id)
                {
                    onehalf_coeffs_precond_sorted_mode_[mode][t][id] = std::pair<int,double>(it.index(), *it);
                }
                sort(onehalf_coeffs_precond_sorted_mode_[mode][t].begin(), onehalf_coeffs_precond_sorted_mode_[mode][t].end(), typename InfiniteVector<double,int>::decreasing_order());
            }
        }
//        problem_mode_ = oldmode;
//        current_timestep_ = oldtime;
    }
    
    
    
    
    
    
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::compute_rhs(const char* onehalf_filename, const char* onehalf_precond_gramian_filename, const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename)
    {
        // precompute the right-hand side on a fine level
        unsigned int oldtime = current_timestep_;
        InfiniteVector<double,int> fhelp;
        double coeff, temp_d;
        unsigned int id;
        fhelp.clear();
        onehalf_coeffs_precond_unsorted_gramian_.resize(0);
        onehalf_coeffs_precond_unsorted_gramian_.resize(qtbasis_->degrees_of_freedom());
        onehalf_precond_norm_sqr_gramian_ = 0;
        for (int mode = 0; mode < 2; ++ mode)
        { 
            for (unsigned int i = 0; i <= NUMOFTIMESTEPS; ++i)
            {
                onehalf_coeffs_precond_unsorted_mode_[mode][i].resize(0);
                onehalf_coeffs_precond_unsorted_mode_[mode][i].resize(qtbasis_->degrees_of_freedom());
                onehalf_precond_norm_sqr_mode_[mode][i] = 0;
// CLEANUP        
                for (unsigned int j = 0; j< qtbasis_->degrees_of_freedom();++j)
                {
                    assert (onehalf_coeffs_precond_unsorted_mode_[mode][i][j] == 0);
                }
            }
        }
        
        for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
        {
            temp_d = qtbasis_->integrate(onehalf_function_,i);
            onehalf_coeffs_.set_coefficient(i,temp_d);
            coeff = temp_d / diagonal_cache_gramian_[i];
            if (fabs(coeff)>1e-15)
            {
                fhelp.set_coefficient(i, coeff);
                onehalf_coeffs_precond_unsorted_gramian_[i] = coeff;
                onehalf_precond_norm_sqr_gramian_ += coeff*coeff;
            }
        }
        writeIVToFile(onehalf_coeffs_,onehalf_filename);
        std::ofstream ofs(onehalf_precond_gramian_filename,std::ofstream::binary);
        writeVectorToFile(onehalf_coeffs_precond_unsorted_gramian_, ofs);
        onehalf_coeffs_precond_sorted_gramian_.resize(0); // clear eventual old values
        onehalf_coeffs_precond_sorted_gramian_.resize(fhelp.size());
        id = 0;
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                it != itend; ++it, ++id)
        {
            onehalf_coeffs_precond_sorted_gramian_[id] = std::pair<int,double>(it.index(), *it);
        }
        sort(onehalf_coeffs_precond_sorted_gramian_.begin(), onehalf_coeffs_precond_sorted_gramian_.end(), typename InfiniteVector<double,int>::decreasing_order());
        for (int t = 0; t < NUMOFTIMESTEPS+1; ++t)
        {
            current_timestep_ = t;
            for (int mode = 0; mode < 2; ++mode)
            {
                fhelp.clear();
                for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
                {
                    temp_d = onehalf_coeffs_.get_coefficient(i);
                    coeff = temp_d / diagonal_cache_mode_[mode][t][i];
                    if (fabs(coeff)>1e-15)
                    {
                        fhelp.set_coefficient(i, coeff);
                        onehalf_coeffs_precond_unsorted_mode_[mode][t][i] = coeff;
                        onehalf_precond_norm_sqr_mode_[mode][t] += coeff*coeff;
                    }
                }
                // store on disc
                std::ofstream ofs(onehalf_mode_filename[mode][t],std::ofstream::binary);
                writeVectorToFile(onehalf_coeffs_precond_unsorted_mode_[mode][t], ofs);
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(0); // clear eventual old values
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(fhelp.size());
                id = 0;
                for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                        it != itend; ++it, ++id)
                {
                    onehalf_coeffs_precond_sorted_mode_[mode][t][id] = std::pair<int,double>(it.index(), *it);
                }
                sort(onehalf_coeffs_precond_sorted_mode_[mode][t].begin(), onehalf_coeffs_precond_sorted_mode_[mode][t].end(), typename InfiniteVector<double,int>::decreasing_order());
            }
        }
        current_timestep_ = oldtime;
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::load_rhs(const char* onehalf_filename, const char* onehalf_precond_gramian_filename, const FixedArray1D< FixedArray1D<char*, NUMOFTIMESTEPS+1>,2> onehalf_mode_filename)
    {
        unsigned int oldtime = current_timestep_;
        InfiniteVector<double,int> fhelp;
        double coeff;
        unsigned int id;
        readIVFromFile(onehalf_coeffs_,onehalf_filename);
        readVectorFromFile(onehalf_coeffs_precond_unsorted_gramian_,onehalf_precond_gramian_filename);
        assert (qtbasis_->degrees_of_freedom() == onehalf_coeffs_precond_unsorted_gramian_.size());
        fhelp.clear();
        onehalf_precond_norm_sqr_gramian_ = 0;
        for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
        {
            //coeff = qtbasis_->integrate(f_,i)/D(i);
            coeff = onehalf_coeffs_precond_unsorted_gramian_[i];
            if (coeff != 0)
            {
                fhelp.set_coefficient(i, coeff);
                onehalf_precond_norm_sqr_gramian_ += coeff*coeff;
            }
        }
        onehalf_coeffs_precond_sorted_gramian_.resize(0); // clear eventual old values
        onehalf_coeffs_precond_sorted_gramian_.resize(fhelp.size());
        id = 0;
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                it != itend; ++it, ++id)
        {
            onehalf_coeffs_precond_sorted_gramian_[id] = std::pair<int,double>(it.index(), *it);
        }
        sort(onehalf_coeffs_precond_sorted_gramian_.begin(), onehalf_coeffs_precond_sorted_gramian_.end(), typename InfiniteVector<double,int>::decreasing_order());
        for (int t = 0; t < NUMOFTIMESTEPS+1; ++t)
        {
            current_timestep_ = t;
            for (int mode = 0; mode < 2; ++mode)
            {
                readVectorFromFile(onehalf_coeffs_precond_unsorted_mode_[mode][t], onehalf_mode_filename[mode][t]);
                assert (qtbasis_->degrees_of_freedom() == onehalf_coeffs_precond_unsorted_mode_[mode][t].size());
                onehalf_precond_norm_sqr_mode_[mode][t] = 0;
                fhelp.clear();
                for (unsigned int i = 0; i< qtbasis_->degrees_of_freedom();i++)
                {
                    //coeff = qtbasis_->integrate(f_,i)/D(i);
                    coeff = onehalf_coeffs_precond_unsorted_mode_[mode][t][i];
                    if (coeff != 0)
                    {
                        fhelp.set_coefficient(i, coeff);
                        onehalf_precond_norm_sqr_mode_[mode][t] += coeff*coeff;
                    }
                }
                // sort the coefficients into fcoeffs
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(0); // clear eventual old values
                onehalf_coeffs_precond_sorted_mode_[mode][t].resize(fhelp.size());
                id = 0;
                for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                        it != itend; ++it, ++id)
                {
                    onehalf_coeffs_precond_sorted_mode_[mode][t][id] = std::pair<int,double>(it.index(), *it);
                }
                sort(onehalf_coeffs_precond_sorted_mode_[mode][t].begin(), onehalf_coeffs_precond_sorted_mode_[mode][t].end(), typename InfiniteVector<double,int>::decreasing_order());
            }
        }
        current_timestep_ = oldtime;
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::RHS(const double eta,
            InfiniteVector<double, int>& coeffs) const
    {
        coeffs.clear();
        if (eta < 1e-8)
        {
            coeffs = *current_rhs_;
//            if (problem_mode_ == 2)
//            {
//                for (unsigned int i=0; i< onehalf_coeffs_precond_unsorted_gramian_.size(); ++i)
//                {
//                    coeffs.set_coefficient(i, onehalf_coeffs_precond_unsorted_gramian_[i]);
//                }
//                cout << "wait a bit" << endl;
//            }
//            else
//            {
//                for (unsigned int i=0; i< onehalf_coeffs_precond_unsorted_mode_.size(); ++i)
//                {
//                    coeffs.set_coefficient(i, onehalf_coeffs_precond_unsorted_mode_[problem_mode_][current_timestep_][i]);
//                }
//                cout << "wait a bit" << endl;
//            }
            return;
        }
        double coarsenorm(0);
        double bound;
//        typename Array1D<std::pair<int, double> >::const_iterator it, itend;
        typename Array1D<std::pair<int, double> >::const_iterator it(current_rhs_sorted_->begin()), 
                itend(current_rhs_sorted_->end());
//        if (problem_mode_ == 2)
//        {
//            bound = (onehalf_precond_norm_sqr_gramian_ - eta*eta);
//            it = onehalf_coeffs_precond_sorted_gramian_.begin();
//            itend = onehalf_coeffs_precond_sorted_gramian_.end();
//        }
//        else
//        {
//            bound = (onehalf_precond_norm_sqr_mode_[problem_mode_][current_timestep_] - eta*eta);
//            it = onehalf_coeffs_precond_sorted_mode_[problem_mode_][current_timestep_].begin();
//            itend = onehalf_coeffs_precond_sorted_mode_[problem_mode_][current_timestep_].end();
//        }
        bound = current_rhs_l2_norm_*current_rhs_l2_norm_ - eta*eta;
        while ((it != itend) && (coarsenorm < bound))
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
            ++it;
        }
    }

    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::add_ball(const unsigned int& lambdanum,
                                                   Vector<double>& w,
                                                   const int radius,
                                                   const double factor,
                                                   const int maxlevel,
                                                   const CompressionStrategy strategy,
                                                   const bool precond)
    {
        /*
         * For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
         * (Attention: possibly different minimal levels on different patches!)
         * For higher dimensions this is done recursivly.
         */
        Index temp_lam(qtbasis_->get_wavelet(lambdanum));
        MultiIndex<int,DIM> lam_j(temp_lam.j());
        double d1 = precond ? D(lambdanum) 
                            : 1.0;
        MultiIndex<int,DIM> temp_j;
        unsigned int temp_p;
        if (DIM == 1)
        {
            abort();
            //confer cached_qtproblem version
        }
        else if (DIM == 2)
        {
            add_leveldisc_recurse(lambdanum,w,lam_j,-1,radius,factor/d1,list<int>(),false,precond);
        }
        else // dim > 2. iteration over all levels in 'range' is done recursivly for arbitrary dimensions
        {
            // iterate over all possible values for the first component of the level j
            // for each j[0]: add levels with valid values for j[1],...,j[DIM-1]
            // this can be done recursively
            // for the last two dimensions instead of doing a final recursion, proceed as for the case DIM=2 (above in this method)
            
            int xstart, xend;
            xstart = lam_j[0]-radius;
            bool set_is_restriction(false);
            // qtbasis_->j0()[patch_with_minimal_j0_[0].first][0] == lowest possible value for j[0]
            if (xstart <= qtbasis_->j0()[this->patch_with_minimal_j0_[0].front()][0])
            {
                xstart = qtbasis_->j0()[this->patch_with_minimal_j0_[0].front()][0];
                set_is_restriction =(this->patch_with_minimal_j0_[0].size() != qtbasis_->get_nop()); // check whether all patches allow for the lowest possible value for j[0]
            }
            temp_j = lam_j;
            assert (radius >= 0);
            assert (qtbasis_->get_jmax() - multi_degree(lam_j) >= 0);
            xend = lam_j[0] + std::min (radius, (int) (qtbasis_->get_jmax() - multi_degree(lam_j)));
            
            temp_j[0] = xstart;
            add_leveldisc_recurse(lambdanum,w,temp_j,0,radius-abs(xstart-lam_j[0]),factor/d1,this->patch_with_minimal_j0_[0],set_is_restriction,precond);
                
            for (unsigned int x = xstart+1; x<= xend; ++x)
            {
                // first component of the current level is set to "x"
                // add all possible levels with this value
                // ASSUMPTION: minimal levels differ at most by 1 => for x>xstart the value x is possible for each patch
                temp_j[0] = x;
                add_leveldisc_recurse(lambdanum,w,temp_j,0,radius-abs(x-lam_j[0]),factor/d1,list<int>(),false,precond);
            }
        }
    }
    
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS,  QTBASIS,  PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT >
    ::add_leveldisc_recurse(
                const unsigned int& lambdanum,
                Vector<double>& w,
                const MultiIndex<int,DIM> center_j,
                const int fixed_dimensions,
                const unsigned int radius,
                const double factor,
                const list<int> & j0_restriction_of_patches,
                const bool set_is_restriction,
                const bool precond)
    {
        assert (DIM-fixed_dimensions>2);
        MultiIndex<int,DIM> temp_j;
        int xstart, xend;
        list<int> intersection;
        
        if (DIM-fixed_dimensions > 3)
        {
            xstart = center_j[fixed_dimensions+1]-radius;
            bool restricted(set_is_restriction);
            // qtbasis_->j0()[patch_with_minimal_j0_[0].first][0] == lowest possible value for j[0]
            if (xstart <= qtbasis_->j0()[this->patch_with_minimal_j0_[fixed_dimensions+1].front()][fixed_dimensions+1])
            {
                xstart = qtbasis_->j0()[this->patch_with_minimal_j0_[fixed_dimensions+1].front()][fixed_dimensions+1];
                if (this->patch_with_minimal_j0_[fixed_dimensions+1].size() != qtbasis_->get_nop())
                {
                    if (restricted)
                    {
                        set_intersection(this->patch_with_minimal_j0_[fixed_dimensions+1].begin(), this->patch_with_minimal_j0_[fixed_dimensions+1].end(),
                                j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                back_inserter(intersection));
                        if (intersection.empty())
                        {
                            return;
                        }
                    }
                    else
                    {
                        restricted = true;
                        intersection = this->patch_with_minimal_j0_[fixed_dimensions+1];
                    }
                }
            }
            temp_j = center_j;
            assert (radius >= 0);
            assert (qtbasis_->get_jmax() - multi_degree(center_j) >= 0);
            xend = center_j[fixed_dimensions+1] + std::min ((int) radius, (int) (qtbasis_->get_jmax()) - (int)multi_degree(center_j));
            
            temp_j[fixed_dimensions+1] = xstart;
            add_leveldisc_recurse(lambdanum,w,temp_j,fixed_dimensions+1,radius-abs(xstart-center_j[fixed_dimensions+1]),factor,  intersection,restricted,precond);
            for (unsigned int x = xstart+1; x<= xend; ++x)
            {
                // first component of the current level is set to "x"
                // add all possible levels with this value
                // ASSUMPTION: minimal levels differ at most by 1 => for x>xstart the value x is possible for each patch
                temp_j[fixed_dimensions+1] = x;
                add_leveldisc_recurse(lambdanum,w,temp_j,fixed_dimensions+1,radius-abs(x-center_j[fixed_dimensions+1]),factor,  j0_restriction_of_patches,set_is_restriction,precond);
            }
        }
        else
        {
            // there are only 2 dimensions left
            // if set_is_restriction == false: the code to iterate over the last 2 dimensions looks the same as in add_ball for DIM==2
            // else: the current levelline has holes, i.e., only patches in the restriction set may be active
            
            // The ball can be described of levellines consisting of levels with the same multidegree
            // The first is determined with the distance of lambda.j and (minimal) j0. The last with qtbasis_->get_jmax()
            // The first level in a levelline is determined with minx = max(j0[DIM-2], lambda.j[DIM-2]-radius)
            // The last level in a levelline is determined with miny = max(j0[DIM-1], lambda.j[DIM-1]-radius)
            
            MultiIndex<int, 2> min_j, max_j;
            unsigned int temp_p;
            qtbasis_->get_level(0, temp_j,temp_p);
            
            //index_lt j0(this->basis().j0());
            int lambdaline = multi_degree(center_j);
            int lowestline = multi_degree(temp_j);
            int dist2j0=lambdaline-lowestline;
            int dist2maxlevel=qtbasis_->get_jmax()-lambdaline;
            int ystart,temp_i;
            bool is_minimal_x(false), is_minimal_y(false); // is_minimal_x == true means that the first level of the current levelline (in the current ball) hits j0[DIM-2]
            list<int> intersection_x, intersection_y; // intersection_x contains (if needed) the intersection of j0_restriction_of_patches and patch_with_minimal_j0_[DIM-2]  
            
            //bool j0_restricts_x(set_is_restriction), j0_restrics_y(set_is_restriction); // false means that all patches yield valid levels
            //List<int>* list_x; // if j0_restricts_x == true this points to the list with the valid patches
            //list<int>* list_y;
            //unsigned int first_levelnum, last_levelnum; // first and last level on the current levelline
            
            
            // iterate the levellines. offset relative to center_j's levelline
            for (int offset = -std::min(dist2j0,(int)radius); offset < std::min(dist2maxlevel,(int)radius)+1; offset++)
            {
                // iterate over the levels on the levelline
                // ignoring restrictions by j0 for the moment, we have:
                
                xstart = center_j[DIM-2]-radius+(int)ceil((radius+offset)/2.0); //x coordinate of the first level
                xend = center_j[DIM-2]+(int)floor((radius+offset)/2.0); // same for the last
                ystart = center_j[DIM-1]+(int)floor((radius+offset)/2.0); // and for the second dimension
                // we walk (make steps) through the levelline
                // The level (xstart,ystart) in the current levelline hast steps=0.
                // restrictions by j0[DIM-2] or j0[DIM-1] mean:
                // The first level on the current levelline may have steps >0.
                
                // ignoring the restriction by set_is_restriction for the moment, we have
                
                temp_i = qtbasis_->j0()[this->patch_with_minimal_j0_[DIM-2].front()][DIM-2];
                if (xstart <= temp_i)
                {
                    min_j[0]= temp_i;
                    min_j[1]= ystart-(temp_i-xstart);
                    is_minimal_x = true; // restrictions by patch_with_minimal_j0_[DIM-2] apply
                }
                else
                {
                    min_j[0] = xstart;
                    min_j[1] = ystart;
                }
                temp_i = min(xend-xstart,ystart-qtbasis_->j0()[this->patch_with_minimal_j0_[DIM-1].front()][DIM-1]);
                max_j[0]= xstart+temp_i;
                max_j[1]= ystart-temp_i;

                if (ystart-qtbasis_->j0()[this->patch_with_minimal_j0_[DIM-1].front()][DIM-1] <= xend - xstart )
                {
                    is_minimal_y = true; // restrictions by patch_with_minimal_j0_[DIM-1] apply
                }
                
                // how many j's are legal between xmin and xmax?
                temp_i = max_j[0] - min_j[0];
                if (temp_i < 0)
                {
                    return;
                }
// cleanup
                assert (temp_i == (min_j[1]-max_j[1]));
                // incorporate set_is_restriction
                temp_j = center_j;
                temp_j[DIM-2] = min_j[0];
                temp_j[DIM-1] = min_j[1];          

                if (set_is_restriction)
                {
                    // restrictions by j0_restriction_of_patches apply
                    if (is_minimal_x)
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(this->patch_with_minimal_j0_[DIM-2].begin(), this->patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
                                set_intersection(this->patch_with_minimal_j0_[DIM-1].begin(), this->patch_with_minimal_j0_[DIM-1].end(),
                                        intersection_x.begin(), intersection_x.end(),
                                        back_inserter(intersection));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(this->patch_with_minimal_j0_[DIM-2].begin(), this->patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
                                set_intersection(this->patch_with_minimal_j0_[DIM-1].begin(), this->patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in intersection_x and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            qtbasis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(this->patch_with_minimal_j0_[DIM-2].begin(), this->patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection_x.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(this->patch_with_minimal_j0_[DIM-2].begin(), this->patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in intersection_x and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<=temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            qtbasis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                            }
                        }
                    }
                    else
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] do not apply                        
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(this->patch_with_minimal_j0_[DIM-1].begin(), this->patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection_y.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection_y and with j = min_j
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(this->patch_with_minimal_j0_[DIM-1].begin(), this->patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in j0_restriction_of_patches and with j = min_j, steps = 0,...,temp_i-1
                                for (int i=0; i<temp_i;++i)
                                {
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            qtbasis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply

                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (int i=0; i<=temp_i;++i)
                            {
                                for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                            }
// CLEANUP
                            assert ((temp_j[DIM-2]-1) == max_j[0]);
                            assert ((temp_j[DIM-1]+1) == max_j[1]);
                        }
                    }
                }
                else
                {
                    // restrictions by j0_restriction_of_patches do not apply
                    if (is_minimal_x)
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(this->patch_with_minimal_j0_[DIM-2].begin(), this->patch_with_minimal_j0_[DIM-2].end(),
                                        this->patch_with_minimal_j0_[DIM-1].begin(), this->patch_with_minimal_j0_[DIM-1].end(),
                                        back_inserter(intersection));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
//                                cout << "start: temp_j = " << temp_j << endl;
                                for (list<int>::const_iterator it(this->patch_with_minimal_j0_[DIM-2].begin()), itend(this->patch_with_minimal_j0_[DIM-2].end()); it!= itend; ++it)
                                {
//                                    cout << "start: adding level temp_j = " << temp_j << "; p = " << *it << endl;
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
//                                    cout << "center: temp_j = " << temp_j << endl;
                                    for (unsigned int p=0; p<qtbasis_->get_nop();++p)
                                    {
//                                        cout << "center: adding level temp_j = " << temp_j << "; p = " << p << endl;
                                        add_level(lambdanum,
                                            w,
                                            qtbasis_->get_levelnum(temp_j,p),
                                            factor,
                                            precond);
                                    }
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
//                                cout << "end: temp_j = " << temp_j << endl;
//                                cout << "min_j = " << min_j << endl;
//                                cout << "max_j = " << max_j << endl;
                                
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(this->patch_with_minimal_j0_[DIM-1].begin()), itend(this->patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
//                                    cout << "end: adding level temp_j = " << temp_j << "; p = " << *it << endl;
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            } 
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply
                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (list<int>::const_iterator it(this->patch_with_minimal_j0_[DIM-2].begin()), itend(this->patch_with_minimal_j0_[DIM-2].end()); it!= itend; ++it)
                            {
                                add_level(lambdanum,
                                    w,
                                    qtbasis_->get_levelnum(temp_j,*it),
                                    factor,
                                    precond);
                            }
                            for (int i=1; i<=temp_i;++i)
                            {
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                                for (unsigned int p=0; p<qtbasis_->get_nop();++p)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,p),
                                        factor,
                                        precond);
                                }
                            }
// CLEANUP
                            assert (temp_j[DIM-2] == max_j[0]);
                            assert (temp_j[DIM-1] == max_j[1]);
                        }
                    }
                    else
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] do not apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                // iterate over all levels with p\in intersection_y and with j = min_j
                                for (list<int>::const_iterator it(this->patch_with_minimal_j0_[DIM-1].begin()), itend(this->patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                // iterate over all levels with p\in j0_restriction_of_patches and with j = min_j, steps = 0,...,temp_i-1
                                for (int i=0; i<temp_i;++i)
                                {
                                    for (unsigned int p=0; p<qtbasis_->get_nop();++p)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            qtbasis_->get_levelnum(temp_j,p),
                                            factor,
                                            precond);
                                    }
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(this->patch_with_minimal_j0_[DIM-1].begin()), itend(this->patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply

                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (int i=0; i<=temp_i;++i)
                            {
                                for (unsigned int p=0; p<qtbasis_->get_nop();++p)
                                {
                                    add_level(lambdanum,
                                        w,
                                        qtbasis_->get_levelnum(temp_j,p),
                                        factor,
                                        precond);
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                            }
// CLEANUP
                            assert ((temp_j[DIM-2]-1) == max_j[0]);
                            assert ((temp_j[DIM-1]+1) == max_j[1]);
                        }
                    }
                }
            }
        }
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT>
    ::add_level(const unsigned int& lambdanum,
                Vector<double>& w,
                const unsigned int levelnum,
                const double factor,
                const bool precond)
    {
        // for fixed lambda: compute entries a(mu_2,lambda) for all mu_2 on level levelnum that intersect lambda
        // part 1: for each dimension: extract 1d Block a(mu2_i, lam_i) of 1d integrals with |mu2_i| = |levelnum_j[i]|
        // as in a() these blocks may or may not already be in the cache. So maybe they have to be computed here
        // part 2: compose all relevant entries of w from the extracted 1d Blocks

        // Part 1: extract 1d Blocks from the cache (similar to a() ))
        Index lambda (qtbasis_->get_wavelet(lambdanum));
        const MultiIndex<int, DIM> lambda_j(lambda.j()), lambda_e(lambda.e()), lambda_k(lambda.k());
        const unsigned int lambda_p(lambda.p());
        MultiIndex<unsigned int, DIM> intinfo;
        MultiIndex<int, DIM> mu_j;
        unsigned int mu_p;
        qtbasis_->get_level(levelnum,mu_j,mu_p);
        //bool temp_b(qtbasis_->get_LMR_info(lambdanum,mu_j,mu_p, intinfo));
        bool temp_b(qtbasis_->get_LMR_info(lambdanum,mu_p, intinfo));
        if (!temp_b) 
            return;
        FixedArray1D<int,DIM> kmingen, kmaxgen, kminwav, kmaxwav;
        unsigned int lami_basisnum, mui_basisnum;
        FixedArray1D<bool,DIM> mu_min_type;
        FixedArray1D< Block*, DIM> waveletBlock;
        FixedArray1D< Block*, DIM> generatorBlock;
        for (unsigned int i=0; i<DIM; i++)
        {
            lami_basisnum = (((qtbasis_->get_bc()[lambda_p][2*i])?0:2) + ((qtbasis_->get_bc()[lambda_p][2*i+1])?0:1));
            mui_basisnum = (((qtbasis_->get_bc()[mu_p][2*i])?0:2) + ((qtbasis_->get_bc()[mu_p][2*i+1])?0:1));
            mu_min_type[i] = (mu_j[i]== qtbasis_->j0()[mu_p][i]) ? true:false;
            ColumnCache* cachePointer;
            Column* cachePointerGen;
            switch (intinfo[i])
            {
                case 0:
                case 4:
                case 8: // LL MM RR -> type 1
//                    typedef std::map<int, Column> ColumnCache;
//                    typedef FixedArray1D<ColumnCache,16> typeIcache;
//                    typeIcache typeIcache_;
//                    cachePointer = & typeIcache_[4*lami_basisnum + mui_basisnum];
//                    (this->typeIcache_[0])
                    cachePointer = &this->typeIcache_[4*lami_basisnum + mui_basisnum];
                    if (mu_min_type[i])
                        cachePointerGen = &this->typeIcachedGens_[4*lami_basisnum + mui_basisnum];
                    break;
                case 1:
                case 7: 
                case 2:
                case 6: // LM RM LR RL -> type 2
                    cachePointer = &this->typeIIcache_[4*(lami_basisnum-1) + mui_basisnum];
                    if (mu_min_type[i])
                        cachePointerGen = & this->typeIIcachedGens_[4*(lami_basisnum-1) + mui_basisnum];
                    break;
                case 3: // ML -> type 3
                    // EXTREME CAUTION:: mui_basisnum - (intinfo[i] == 5)?0:1 produces garbage!
                    assert (mui_basisnum != 1);
                    cachePointer = & this->typeIIIcache_[4*lami_basisnum + mui_basisnum - 1];
                    if (mu_min_type[i])
                        cachePointerGen = & this->typeIIIcachedGens_[4*lami_basisnum + mui_basisnum - 1];
                    break;
                case 5: // MR -> type 3
                    assert (mui_basisnum != 2);
                    cachePointer = & this->typeIIIcache_[4*lami_basisnum + ((mui_basisnum == 1)?0:3)];
                    if (mu_min_type[i])
                        cachePointerGen = & this->typeIIIcachedGens_[4*lami_basisnum + ((mui_basisnum == 1)?0:3)];
                    break;
                default:
                    abort();
            }
            // search for column 'lami'
            typename QTBASIS::IntervalBasis::Index lami(lambda_j[i],lambda_e[i],lambda_k[i],qtbasis_->get_bases_infact()[lami_basisnum]);
            unsigned int lami_num = lami.number();
            typename ColumnCache::iterator col_lb(cachePointer->lower_bound(lami_num));
            typename ColumnCache::iterator col_it(col_lb);
            if (col_lb == cachePointer->end() ||
                cachePointer->key_comp()(lami_num, col_lb->first))
            {
                // the lami-th column has not been requested so far
                // insert a new column and continue with the blocks
                typedef typename ColumnCache::value_type value_type;
                col_it = cachePointer->insert(col_lb, value_type(lami_num, Column()));
            }
            // col_it points to the column of lami
            Column& col(col_it->second);
            // check whether the level 'mu_i' belongs to has already been calculated
            int mui_levelnum (mu_j[i] - qtbasis_->j0()[mu_p][i]);
            typename Column::iterator lb(col.lower_bound(mui_levelnum));
            typename Column::iterator it(lb);
            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                // no entries have ever been computed for this column and this level
                // insert a new block, 
                // and then compute the whole level block 
                //      (\int psi_mui psi_lami)_mui, (\int psi_mui' psi_lami')_mui
                // for all mui that intersect lami. 
                // if mui is on the minimal level we also need to add a new Block() to
                // the generator integral cache, i.e., we need to 
                // integrate against all generators on the lowest level
                typedef typename Column::value_type value_type;
                it = col.insert(lb, value_type(mui_levelnum, Block()));
                waveletBlock[i] = & it->second;
                bool gen_intersection_i, wav_intersection_i;
                qtbasis_->get_onedim_intersections(intinfo[i],
                        lambda_j[i],
                        lambda_e[i],
                        lambda_k[i],
                        lami_basisnum,
                        (mui_levelnum == 0),
                        mu_j[i],
                        mui_basisnum,
                        kmingen[i],
                        kmaxgen[i],
                        kminwav[i],
                        kmaxwav[i],
                        gen_intersection_i,
                        wav_intersection_i);
                FixedArray1D<double,ONEDIMHAARCOUNT> gram;
                FixedArray1D<double,1> der;
                if (mu_min_type[i])
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(lami_num));
                    typename Column::iterator it2(lb2);
                    if (lb2 == cachePointerGen->end() ||
                        cachePointerGen->key_comp()(lami_num, lb2->first))
                    {
                        // no entries have ever been computed for this column and this level
                        // insert a new block, 
                        // and then compute the whole level block 
                        //      (\int psi_nui psi_mui)_nui, (\int psi_nui' psi_mui')_nui
                        // for all nui that intersect mui. 
                        // if mui is on the minimal level we also need to add a new Block() to
                        // the generator integral cache, i.e., we need to 
                        // integrate against all generators on the lowest level
                        typedef typename Column::value_type value_type;
                        it2 = cachePointerGen->insert(lb2, value_type(lami_num, Block()));
                        generatorBlock[i] =  & it2->second;
                        //Block& block2(it2->second);
                        if (gen_intersection_i)
                        {
                            for (int kgen = kmingen[i]; kgen <= kmaxgen[i]; ++kgen)
                            {
                                compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), lambda_j[i], lambda_e[i], lambda_k[i], lami_basisnum,
                                        mu_j[i], 0, kgen, mui_basisnum,
                                        gram,
                                        der);
                                typedef typename Block::value_type value_type_block;
                                generatorBlock[i]->insert(generatorBlock[i]->end(), value_type_block(kgen, make_pair(gram,der)));
                            }
                        }
                    }
                    else
                    {
                        cout << "just inserted a wavelet block, but gen block was already existent?!!" << endl;
                        abort();
                    }
                }
                // code invariant: there exists a Block() (maybe empty) corresponding to lami and mui
                
                // if (no intersections) => current block (wavCache) remains empty
                // if we are on the lowest possible level for mui, we add an empty Block to the genCache
                // advantage: genCache[laminum] exists and contains a block. 
                // code invariant: to access an entry: load the block and check if there is an entry (simple! same whether there is an entry or not!)
                // it would be possible to simply leave the genCache as it is.
                // advantage: nothing to do at this point
                // access to an entry: check whether there is a column in genCache. If not: value 0. If there is: get Block and take value from it.
                // This leads to 2 checks per call to a() instead of one. 
                // This is cheaper than method 1 above if lambda does not intersect at all with basis functions mu.
                // However, most lambdas will have some intersection? In this case this variant leads to 1 additional check
                // so ... hopefully this makes the average access time to the cache faster
                if (!wav_intersection_i) 
                    return;
                // wav_intersection_i == true guarantees kminwavi <=kmaxwavi and that the values are meaningful
                for (int kwav = kminwav[i]; kwav <= kmaxwav[i]; ++kwav)
                {
                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), lambda_j[i], lambda_e[i], lambda_k[i], lami_basisnum,
                            mu_j[i], 1, kwav, mui_basisnum,
                            gram,
                            der);

                    typedef typename Block::value_type value_type_block;
                    waveletBlock[i]->insert(waveletBlock[i]->end(), value_type_block(kwav, make_pair (gram,der) ));
                }
            }
            else // column and levelblock already exist
            {
                waveletBlock[i] = & it->second;
                if (waveletBlock[i]->size() == 0)
                    return;
                if (mu_min_type[i])
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(lami_num));
                    generatorBlock[i] = & lb2->second;
                }
            }
        } // end of loop over dim
        // part 2: compose all relevant entries of w from the extracted 1d Blocks
        // part 2a: include information about mu into LMR info:
        FixedArray1D<Array1D<unsigned int> ,DIM> mu_gen_adapted_intinfo;
        FixedArray1D<Array1D<unsigned int> ,DIM> mu_wav_adapted_intinfo;
        for (unsigned int i=0; i<DIM; ++i)
        {
            kminwav[i] = waveletBlock[i]->begin()->first;
            kmaxwav[i] = waveletBlock[i]->rbegin()->first;
            assert (waveletBlock[i]->size() == (kmaxwav[i]-kminwav[i]+1));
            if (mu_min_type[i] && (generatorBlock[i]->size() != 0))
            {
                // loop over all intersecting generators in current dimension
                mu_gen_adapted_intinfo[i].resize(generatorBlock[i]->size());
                kmingen[i] = generatorBlock[i]->begin()->first;
                kmaxgen[i] = generatorBlock[i]->rbegin()->first;
                assert (generatorBlock[i]->size() == (kmaxgen[i]-kmingen[i]+1));
                typename Block::iterator block_lb(generatorBlock[i]->lower_bound(kmingen[i]));
                if (block_lb == generatorBlock[i]->end() ||
                        generatorBlock[i]->key_comp()(kmingen[i], block_lb->first))
                {
                    cout << "error! 1d cache contains generator block for current level but valid entry kmin_gen[i] is missing!" << endl;
                    abort();
                }
                for (int k = kmingen[i], n(0); k<= kmaxgen[i]; ++k, ++n)
                {
                    // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                    mu_gen_adapted_intinfo[i][n] = intinfo[i];
                    switch (intinfo[i])
                    {
                        case 0:
                            // LL: check whether mu is a left boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )
                            {
                                // nu is not extended left
                                mu_gen_adapted_intinfo[i][n] = 4; // MM
                            }
                            break;
                        case 2:
                            // LR: check whether mu is a right boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaRmax(  qtbasis_->get_j0(mu_p,i) )) )
                              
                            {
                                mu_gen_adapted_intinfo[i][n] = 1; // LM
                            }
                            break;
                        case 3:
                            // ML : check whether mu is a left boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                            {
                                // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                                abort();
                            }
                            break;
                        case 5:
                            // MR: check whether mu is a right boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) )) )
                            {
                                // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                                abort();
                            }
                            break;
                        case 6:
                            // RL : check whether mu is a left boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                            {
                                // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                                mu_gen_adapted_intinfo[i][n] = 7; // RM
                            }
                            break;
                        case 8:
                            // RR: check whether mu is a right boundary gen/wav
                            if (!(k == qtbasis_-> get_bases(mu_p,i)->DeltaRmax( qtbasis_->get_j0(mu_p,i) )) ) 
                            {
                                // nu is not extended right
                                mu_gen_adapted_intinfo[i][n] = 4; // MM
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
                } // end of loop over generatorBlock
            }
            // loop over all intersecting wavelets in current dimension
            mu_wav_adapted_intinfo[i].resize(kmaxwav[i]-kminwav[i]+1);
            typename Block::iterator block_lb(waveletBlock[i]->lower_bound(kminwav[i]));
            if (block_lb == waveletBlock[i]->end() ||
                    waveletBlock[i]->key_comp()(kminwav[i], block_lb->first))
            {
                cout << "error! 1d cache contains wavelet block for current level but valid entry kmin_wav[i] is missing!" << endl;
                abort();
            }
            for ( int k = kminwav[i], n(0); k<= kmaxwav[i]; ++k, ++n)
            {
                // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                mu_wav_adapted_intinfo[i][n] = intinfo[i];
                switch (intinfo[i])
                {
                    case 0:
                        // LL: check whether mu is a left boundary wav
                        if (!(k < qtbasis_->get_numofbw() ) ) 
                        {
                            // mui is not extended left
                            mu_wav_adapted_intinfo[i][n] = 4; // MM
                        }
                        break;
                    case 2:
                        // LR: check whether mui is a right boundary wav
                        if (!( (k > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) ) )
                        {
                            mu_wav_adapted_intinfo[i][n] = 1; // LM
                        }
                        break;
                    case 3:
                        // ML : check whether mu is a left boundary wav
                        if (!( k < qtbasis_->get_numofbw()) )
                        {
                            // mui is not extended left. There is no intersection of mu and nu! This should not happen at this point
                            cout << "lambdanum = " << lambdanum 
                                    << "; lambda = " << *qtbasis_->get_wavelet(lambdanum) 
                                    << "; levelnum = " << levelnum 
                                    << "; mu_j = " << mu_j << "; mu_p = " << mu_p << endl;
                            cout << "i = " << i << ";k = " << k << "; (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() = " << (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() << endl;
                            cout << "mu_wav_adapted_intinfo = " << mu_wav_adapted_intinfo 
                                    << "; mu_gen_adapted_intinfo = " << mu_gen_adapted_intinfo
                                    << "; intinfo = " << intinfo << endl;
                                    cout << "kminwav = " << kminwav << "; kmingen = " << kmingen << endl;
                            assert(false);
                        }
                        break;
                    case 5:
                        // MR: check whether mu is a right boundary wav
                        if (!(k > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) )
                        {
                            // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                            
                            cout << "lambdanum = " << lambdanum 
                                    << "; lambda = " << *qtbasis_->get_wavelet(lambdanum) 
                                    << "; levelnum = " << levelnum 
                                    << "; mu_j = " << mu_j << "; mu_p = " << mu_p << endl;
                            cout << "i = " << i << ";k = " << k << "; (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() = " << (1<<mu_j[i])- 1 - qtbasis_->get_numofbw() << endl;
                            cout << "mu_wav_adapted_intinfo = " << mu_wav_adapted_intinfo 
                                    << "; mu_gen_adapted_intinfo = " << mu_gen_adapted_intinfo
                                    << "; intinfo = " << intinfo << endl;
                                    cout << "kminwav = " << kminwav << "; kmingen = " << kmingen << endl;
                            assert(false);
                        }
                        break;
                    case 6:
                        // RL : check whether mu is a left boundary gen/wav
                        if (!(k < qtbasis_->get_numofbw()) )
                        {
                            // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                            mu_wav_adapted_intinfo[i][n] = 7; // RM
                        }
                        break;
                    case 8:
                        // RR: check whether mu is a right boundary gen/wav
                        if (!(k > (1<<mu_j[i])- 1 - qtbasis_->get_numofbw()) )
                        {
                            // nu is not extended right
                            mu_wav_adapted_intinfo[i][n] = 4; // MM
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
            } // end of loop over waveletBlock
        }
        // part 2b: use mu-adapted LMR-info and 1d integrals in wavelet/generatorBlocks to efficiently compute entries of w
        FixedArray1D< typename Block::iterator, DIM> block_it; //, block_it_begin;
        FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it;
        typedef typename Index::type_type type_type;
        type_type current_type;
        for (unsigned int i=0;i<DIM;i++)
        {
            if (mu_min_type[i] == true)
            {
                // current_type[i] is coded as an integer, 0 = gen, 1 = wav
                current_type[i] = 0;
            }
            else
            {
                current_type[i] = 1;
            }
        }
        unsigned int number;
        Index temp_ind(qtbasis_->first_wavelet(qtbasis_->get_levelnum(mu_j,mu_p)));
        number = temp_ind.number(); // first wavelet with current_j,current_p
        FixedArray1D<int,DIM> jump_before, jump_after;
        unsigned int cjd; // current jump dimension == dimension (array index) where we are increasing k
        int blocksize;
        bool skip_this_type(false);
        for (unsigned int i=0; i<DIM; ++i)
        {
            if (mu_min_type[i] && (generatorBlock[i]->size() == 0))
            {
                skip_this_type=true;
                break;
            }
        }
        bool done = false;
        while (!done)
        {
            //compose the information in intersections into QTBasis::Index information and keep track of the number of the indices
            // for the current type e:
            // we iterate over all possible k
            // On the lowest dimension this is simple: every increase in k[DIM-1] increses the number of the wavelet by 1.
            // In other dimensions we need to keep track how the number changes, if we increase k[cjd] and set k[j] to the lowest possible value.
            // this is stored in jump_before and jump_after
            if (skip_this_type)
            {
                done = true; // == increase type e
                blocksize = 1;
                for (unsigned int i=0; i<DIM; ++i)
                {
                    if (current_type[i] == 0)
                    {
                        blocksize *= qtbasis_->bases()[mu_p][i]->Deltasize(mu_j[i]);
                    }
                    else
                    {
                        blocksize *= qtbasis_->bases()[mu_p][i]->Nablasize(mu_j[i]);
                    }
                }
                number += blocksize;
            }
            else
            {
                cjd=0;
                if (DIM == 1)
                {
                    if (current_type[0] == 0)
                    {
                        number += kmingen[0] - qtbasis_->bases()[mu_p][0]->DeltaLmin();
                        blocksize = kmaxgen[0] - kmingen[0] +1;
                    } 
                    else
                    {
                        number += kminwav[0] - qtbasis_->bases()[mu_p][0]->Nablamin();
                        blocksize = kmaxwav[0] - kminwav[0] +1;
                    }
                    for (unsigned int n = number; n < number + blocksize; ++n)
                    {
    // CLEANUP                        
                        assert (qtbasis_->get_wavelet(n)->p() == mu_p);
                        assert (qtbasis_->get_wavelet(n)->j() == mu_j);
                        assert (qtbasis_->get_wavelet(n)->e() == current_type);
                        assert (qtbasis_->get_wavelet(n)->k()[0] == ((current_type[0] == 0)?kmingen[0]:kminwav[0])+n-number);
                        assert (qtbasis_->get_wavelet(n)->number() == n);
                        cout << "sorry! qtbasis::add_ball not completely implemented for 1d case" << endl;
                        abort();
                    }
                    assert ( number + ((current_type[0] == 0)? (qtbasis_->bases()[mu_p][0]->DeltaRmax(mu_j[0]) - kmaxgen[0]) : (qtbasis_->bases()[mu_p][0]->Nablamax(mu_j[0]) - kmaxwav[0]) )
                            ==
                            ((qtbasis_->bases()[mu_p][0]->Deltasize(mu_j[0])) + ((current_type[0] == 1) ? (qtbasis_->bases()[mu_p][0]->Nablasize(mu_j[0])):0)) );
                    number = ((qtbasis_->bases()[mu_p][0]->Deltasize(mu_j[0])) + ((current_type[0] == 1) ? (qtbasis_->bases()[mu_p][0]->Nablasize(mu_j[0])):0));
                    // if e[0] == 0: number now points to the first wavelet with the type e=1 but with the same (p,j). we expect that ++it does not increase p
                    // if e[0] == 1: this patch is finished and we expect that ++it increases p. Iteration should stop now
                    done = true; // i.e.,  increase type e at the end of the while loop
                }
                else
                {
                    int basf(0); // "blocks added so far" (on the current type e) via modulo calculus we can deduce which block we have to add next, i.e., where we have to increase k
                    while (!done)
                    {
                        if (cjd == 0)
                        {
                            // nothing has been done so far for this type
                            // compute all jump_numbers that are needed for the current type e
                            int temp_int(1);
                            for (int i = DIM-1; i >= 0; i--)
                            {
                                if (current_type[i] == 0)
                                { 
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kmingen[i] - qtbasis_->bases()[mu_p][i]->DeltaLmin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= qtbasis_->bases()[mu_p][i]->DeltaRmax(mu_j[i]) - kmaxgen[i];
                                    temp_int *= qtbasis_->bases()[mu_p][i]->Deltasize(mu_j[i]);
                                    block_it[i] = generatorBlock[i]->begin();
                                    mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                                }
                                else
                                {
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kminwav[i] - qtbasis_->bases()[mu_p][i]->Nablamin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= qtbasis_->bases()[mu_p][i]->Nablamax(mu_j[i]) - kmaxwav[i];
                                    temp_int *= qtbasis_->bases()[mu_p][i]->Nablasize(mu_j[i]);
                                    block_it[i] = waveletBlock[i]->begin();
                                    mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                                }
                            }
                            for (int i = DIM-2; i >= 0; --i)
                            {
                                jump_before[i] += jump_before[i+1];
                                jump_after[i] += jump_after[i+1];
                            }
                            if (current_type[DIM-1] == 0)
                            {
                                blocksize = kmaxgen[DIM-1] - kmingen[DIM-1] +1;
                            }
                            else
                            {
                                blocksize = kmaxwav[DIM-1] - kminwav[DIM-1] +1;
                            }
                            basf = 0; // no blocks with the current type were added so far
                        }
                        // "jump_before" for all dimensions (current_jump_dim, current_jump_dim+1,...,dim)
                        number += jump_before[cjd];
                        // add wavelets
                        // block_it points to the first_wavelet we want to add
    // CLEANUP                        
                        for (unsigned int n = number; n < number + blocksize; ++n)
                        {
                            assert (qtbasis_->get_wavelet(n)->p() == mu_p);
                            assert (qtbasis_->get_wavelet(n)->j() == mu_j);
                            assert (qtbasis_->get_wavelet(n)->e() == current_type);
                            assert (qtbasis_->get_wavelet(n)->k()[DIM-1] == ( (current_type[DIM-1] == 0)? kmingen[DIM-1]:kminwav[DIM-1]) +n-number); // this only checks the last entry of k
                            assert (qtbasis_->get_wavelet(n)->number() == n);
                        }
                        compose_wavelets(w, number, blocksize, factor, lambda_p, 
                                mu_adapted_intinfo_it,
                                block_it,
                                precond);
                        ++basf; // we have added a block!
                        number += blocksize;
                        // "increase k"
                        // find the first position (from the back) such that k[i] < kmax[i]
                        // then increase number by jump_after[j]. 
                        // Effect: the new number is the number of the first wavelet with the new value of k[i] and k[j]=kmin[j], j=i+1,..,DIM-1
                        // have we arrived at the very last possible k? -> done = true
                        // compute the dimension in which we want to increase k
                        int temp_i(basf); //, temp_count;
                        for (int i(DIM-1); i >= 0; --i)
                        {
                            if (i == 0) //we have already added the last possible block for this type, i.e., k[0] = maxk[0]
                            {
                                done = true; // apparently we have just added the last possible k for the current type vector e. We need to increase it
                                cjd = 0;
                                break;
                            }
                            if (current_type[i] == 0)
                            {
                                block_it[i] = generatorBlock[i]->begin();
                                mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                            }
                            else
                            {
                                block_it[i] = waveletBlock[i]->begin();
                                mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                            }
                            // use cjd as temporary int variable
                            cjd = (current_type[i-1] == 0)? (kmaxgen[i-1] - kmingen[i-1] +1):(kmaxwav[i-1]-kminwav[i-1]+1); // blocksize in this dimension
                            if ((temp_i % cjd) != 0)
                            {
                                cjd = i; // == "we can increase k[cjd]"
                                ++block_it[i-1];
                                ++mu_adapted_intinfo_it[i-1];
                                break;
                            }
                            temp_i = temp_i / cjd;
                        }
                        number += jump_after[cjd]; // number is the number of the first wavelet with the new value for k[cjd] and k[j] = Nablamin/Deltamin for j>cjd
                    }
                } // end of if (DIM > 1)
            } // end of   if (skip_this_type) {} else {}
            // increase the type e
            // number should be the number of the first wavelet with the new type
            // "small loop":
            // try to increase currenttype
            // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
            
            skip_this_type = false;
            int i(DIM-1);
            for (; i >= 0; i--)
            {
                // find first position on level j0
                if (mu_min_type[i])
                {
                    if (current_type[i] == 1)
                    {
                        current_type[i]=0;
                        if (generatorBlock[i]->size() == 0)
                        {
                            skip_this_type = true;
                        }
                        else
                        {
                            block_it[i] = generatorBlock[i]->begin();
                            mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                        }
                    }
                    else
                    {
                        //block_it_begin[i] = waveletBlock[i]->begin(); // value for i=0 is not used
                        current_type[i]=1;
                        block_it[i] = waveletBlock[i]->begin();
                        mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                        done = false;
                        break;
                    }
                }
            }
            for (int j=0; j<i; ++j)
            {
                if ((current_type[j] == 0) && (generatorBlock[j]->size() == 0))
                {
                    skip_this_type=true;
                    break;
                }
            }
            // done == true means that all components with currentlevel[i]=j0[i] were wavelets, i.e.,
            // we have arrived at the type vector (1,1,...,1) and the iteration is finished
        } // end of while(true), i.e., all intersecting wavelets have been composed and added to the cache
    }
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void
    AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT>
    ::cached_sampled_mapping(const InfiniteVector<double,int> & f_qtbasis_coeffs, 
            Array1D < FixedMatrix<double, PRECISE_EVALUATE_GRANULARITY, PRECISE_EVALUATE_GRANULARITY> >& f_haar_gen_coeffs)
    {
        
        assert (PRECISE_EVALUATE_GRANULARITY >= ONEDIMHAARCOUNT);
        //typename QTBASIS::Index Index;
        Index lambda;
        
        FixedArray1D<FixedArray1D<double, PRECISE_EVALUATE_GRANULARITY>, DIM > precise_onedim_integrals;
        MultiIndex<int, DIM> lambda_j, lambda_e, lambda_k;
        unsigned int lambda_p;
        MultiIndex<unsigned int, DIM> intinfo;
//        bool temp_b;
//        FixedArray1D<int,DIM> kmingen, kmaxgen, kminwav, kmaxwav;
        unsigned int lami_basisnum,lami_num;
//        FixedArray1D<bool,DIM> mu_min_type;
//        FixedArray1D< Block*, DIM> waveletBlock;
//        FixedArray1D< Block*, DIM> generatorBlock;
        for (typename InfiniteVector<double,int>::const_iterator it(f_qtbasis_coeffs.begin()), itend(f_qtbasis_coeffs.end()); it!=itend; ++it)
        {
            lambda = qtbasis_->get_wavelet(it.index());
            lambda_j = lambda.j();
            lambda_e = lambda.e();
            lambda_k = lambda.k();
            lambda_p = lambda.p();
            qtbasis_->get_LMR_info(it.index(),lambda_p, intinfo);
            for (unsigned int i=0; i<DIM; i++)
            {
                lami_basisnum = (((qtbasis_->get_bc()[lambda_p][2*i])?0:2) + ((qtbasis_->get_bc()[lambda_p][2*i+1])?0:1));
                typename QTBASIS::IntervalBasis::Index lami(lambda_j[i],lambda_e[i],lambda_k[i],qtbasis_->get_bases_infact()[lami_basisnum]);
                lami_num = lami.number();
                // search for column 'lami'
                //Entries entry;
                typename OneBasisCache::iterator cache_lb(precise_eval_cache_[lami_basisnum].lower_bound(lami_num));
                //typename OneBasisCache::iterator cache_it(cache_lb);
                
                if (cache_lb == precise_eval_cache_[lami_basisnum].end() ||
                    precise_eval_cache_[lami_basisnum].key_comp()(lami_num, cache_lb->first))
                {
                    // the lami-th entry has not been requested so far
                    precise_evaluate (lambda_j[i],lambda_e[i],lambda_k[i],lami_basisnum,precise_onedim_integrals[i]);
                    precise_eval_cache_[lami_basisnum][lami_num] = precise_onedim_integrals[i];
                }
                else
                {
                    precise_onedim_integrals[i] = precise_eval_cache_[lami_basisnum][lami_num];
                }
                
            } // end of loop over dim

            
            unsigned int geometry_type, min_x,max_x,min_y,max_y;
            int centerpatch;
            int north, east, northeast;
            MultiIndex<bool, DIM> orientation;
            double temp_d;
            qtbasis_->get_intersection_geometry(lambda_p, intinfo, geometry_type, centerpatch, orientation);
            if (geometry_type == 3)
            {
                if (centerpatch == -1)
                {
                    north = qtbasis_->get_neighbours(lambda_p,0);
                    east = qtbasis_->get_neighbours(lambda_p,2);
                    geometry_type = 8;
                }
                else
                {
                    north = qtbasis_->get_neighbours(centerpatch,3);
                    east = qtbasis_->get_neighbours(centerpatch,1);
                    if (north != -1)
                    {
                        northeast = qtbasis_->get_neighbours(north,1);
                        if (northeast == -1)
                        {
                            geometry_type = 10;
                        }
                        else
                        {
                            if (east == -1)
                            {
                                geometry_type = 9;
                            }
                        }
                    }
                    else
                    {
                        northeast = qtbasis_->get_neighbours(east,3);
                        geometry_type = 11;
                    }
                }
            }
            for (unsigned int i = 0; i<PRECISE_EVALUATE_GRANULARITY; ++ i)
            {
                if (precise_onedim_integrals[0][i] != 0)
                {
                    min_x = i;
                    break;
                }
            }
            for (unsigned int i = 0; i<PRECISE_EVALUATE_GRANULARITY; ++ i)
            {
                if (precise_onedim_integrals[1][i] != 0)
                {
                    min_y = i;
                    break;
                }
            }
            for (int i = PRECISE_EVALUATE_GRANULARITY-1; i>=0; -- i)
            {
                if (precise_onedim_integrals[0][i] != 0)
                {
                    max_x = i;
                    break;
                }
            }
            for (int i = PRECISE_EVALUATE_GRANULARITY-1; i>=0; -- i)
            {
                if (precise_onedim_integrals[1][i] != 0)
                {
                    max_y = i;
                    break;
                }
            }
            assert (min_x <= max_x);
            assert (min_y <= max_y);
            switch (geometry_type)
            {
                case 0:
                    //Array1D < FixedMatrix<double, PRECISE_EVALUATE_GRANULARITY, PRECISE_EVALUATE_GRANULARITY> >& f_haar_gen_coeffs
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 1:
                    if (orientation[0])
                    {
                        north = qtbasis_->get_neighbours(lambda_p,1); // "north" is just reused. meaning here: neighbour in x direction
                    }
                    else 
                    {
                        north = centerpatch; // "north" is just reused. meaning here: neighbour in x direction
                    }
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[north](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 2:
                    if (orientation[1])
                    {
                        north = qtbasis_->get_neighbours(lambda_p,3); // "north" is just reused. meaning here: neighbour in x direction
                    }
                    else 
                    {
                        north = centerpatch; // "north" is just reused. meaning here: neighbour in x direction
                    }
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[north](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 3:
                    if (orientation[0])
                    {
                        if (orientation[1])
                        {
                            // centerpatch = lambda_p;
                            for (int i = min_x; i<= max_x; ++ i)
                            {
                                temp_d = *it*precise_onedim_integrals[0][i];
                                for (unsigned int j = min_y; j<= max_y; ++ j)
                                {
                                    f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[east](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[north](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[northeast](PRECISE_EVALUATE_GRANULARITY-1-i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                }
                            }
                        }
                        else
                        {
                            // north = lambda_p;
                            for (int i = min_x; i<= max_x; ++ i)
                            {
                                temp_d = *it*precise_onedim_integrals[0][i];
                                for (unsigned int j = min_y; j<= max_y; ++ j)
                                {
                                    f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[northeast](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[centerpatch](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[east](PRECISE_EVALUATE_GRANULARITY-1-i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                }
                            }
                        }
                    }
                    else
                    {
                        if (orientation[1])
                        {
                            // east = lambda_p;
                            for (int i = min_x; i<= max_x; ++ i)
                            {
                                temp_d = *it*precise_onedim_integrals[0][i];
                                for (unsigned int j = min_y; j<= max_y; ++ j)
                                {
                                    f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[centerpatch](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[northeast](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[north](PRECISE_EVALUATE_GRANULARITY-1-i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                }
                            }
                        }
                        else
                        {
                            // northeast = lambda_p;
                            for (int i = min_x; i<= max_x; ++ i)
                            {
                                temp_d = *it*precise_onedim_integrals[0][i];
                                for (unsigned int j = min_y; j<= max_y; ++ j)
                                {
                                    f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[north](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[east](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                    f_haar_gen_coeffs[centerpatch](PRECISE_EVALUATE_GRANULARITY-1-i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                                }
                            }
                        }
                    }
                    break;
                case 8:
                    // centerpatch == southeast is missing; lambda is the northeast
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[north](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[east](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 9:
                    // east is missing; lambda is in the north
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[northeast](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[centerpatch](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 10:
                    // northeast is missing; lambda is the centerpatch
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[east](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[north](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                case 11:
                    // north is missing; lambda is the east
                    for (int i = min_x; i<= max_x; ++ i)
                    {
                        temp_d = *it*precise_onedim_integrals[0][i];
                        for (unsigned int j = min_y; j<= max_y; ++ j)
                        {
                            f_haar_gen_coeffs[lambda_p](i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[centerpatch](PRECISE_EVALUATE_GRANULARITY-1-i,j)+=temp_d*precise_onedim_integrals[1][j];
                            f_haar_gen_coeffs[northeast](i,PRECISE_EVALUATE_GRANULARITY-1-j)+=temp_d*precise_onedim_integrals[1][j];
                        }
                    }
                    break;
                default:
                    cout << "Strange geometry type = " << geometry_type << endl;
            }
        }
    };
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void 
    AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT>
    ::precise_evaluate (const int nui_j,
            const int nui_e,
            const int nui_k,
            const unsigned int nui_basisnum,
            FixedArray1D<double, PRECISE_EVALUATE_GRANULARITY>& precise_lambda_integrals)
    {
        for (unsigned int eta = 0; eta < PRECISE_EVALUATE_GRANULARITY; ++eta)
        {
            precise_lambda_integrals[eta] = 0;
        }
        int nui_k1,nui_k2, haar_scale, temp_i;
        this->basis_->get_bases_infact()[nui_basisnum] -> support(nui_j,nui_e,nui_k,nui_k1,nui_k2);
        unsigned int scale_nui(nui_j+nui_e);
        int scale(scale_nui);
        haar_scale = log2(PRECISE_EVALUATE_GRANULARITY); // assert that wavelet resolution is finer than Haar resolution
        if ( scale < haar_scale)
        {
            scale = haar_scale;
        }
        nui_k1 = nui_k1<<(scale-scale_nui);
        nui_k2 = nui_k2<<(scale-scale_nui);
        
        const unsigned int N_Gauss = this->basis()->get_bases_infact()[nui_basisnum]->primal_polynomial_degree();
        const double h = ldexp(1.0, -scale);
        Array1D<double> gauss_points, nui_val;
        double gauss_weight;
        gauss_points.resize(N_Gauss*(nui_k2-nui_k1));

        // Set up Gauss points and weights for a composite quadrature formula:
        for (int patch = nui_k1, id = 0; patch < nui_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
            for (unsigned int n = 0; n < N_Gauss; n++, id++)
                gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

        // - compute point values of the integrands
        evaluate(*qtbasis_->get_bases_infact()[nui_basisnum], 0, nui_j, nui_e, nui_k, gauss_points, nui_val);
        // - add all integral shares
        for (int patch = nui_k1, id = 0; patch < nui_k2; patch++)
        {
            // number of the active haar-generator patch, i.e., eta is the biggest integer with
            // 2^{-haar_scale} eta \leq 2^{-scale} k_1
            temp_i = patch / (1<<(scale-haar_scale));
            for (unsigned int n = 0; n < N_Gauss; n++, id++) 
            {
                gauss_weight = GaussWeights[N_Gauss-1][n] * h;
                precise_lambda_integrals[temp_i] += nui_val[id]  * gauss_weight;
            }
        }
        for (unsigned int eta=0; eta < PRECISE_EVALUATE_GRANULARITY; ++eta)
        {
            if (abs(precise_lambda_integrals[eta]) < 1e-15)
                precise_lambda_integrals[eta] = 0;
        }
    };
    
}

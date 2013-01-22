// implementation for lin_par_eq_tensor.h

namespace WaveletTL
{
    template <class CACHEDTPROBLEM>
    LinParEqTenROWStageEquationHelper<CACHEDTPROBLEM>
    ::LinParEqTenROWStageEquationHelper(const double a,
                                        const CACHEDTPROBLEM* elliptic,
                                        const CACHEDTPROBLEM* gramian,
                                        const InfiniteVector<double, typename CACHEDTPROBLEM::Index>& z)
    : alpha(a), T(elliptic), G(gramian), y(z), y_scaled(z)
    {
        y_scaled.scale(this,-1);
    }


    template <class CACHEDPROBLEM>
    void
    LinParEqTenROWStageEquationHelper<CACHEDPROBLEM>
    ::add_ball(const Index& lambda,
               //InfiniteVector<double, Index>& w,
               Vector<double>& w,
               const int radius,
               const double factor,
               const int maxlevel,
               const CompressionStrategy strategy,
               const bool precond) const
 {
        // clearly caching the gramian and A at the same time would be more efficient
#if 0
        InfiniteVector<double,Index> g;
        Vector<double> g_full(w.size());
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
                g.set_coefficient(basis().full_collection[i], g_full[i]);
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
    LinParEqTenROWStageEquationHelper<CACHEDPROBLEM>
    ::set_y(const InfiniteVector<double, Index> y_new)
    {
        y = y_new;
        y_scaled = y_new;
        y_scaled.scale(this,-1);
    }

    template <class CACHEDPROBLEM>
    void
    LinParEqTenROWStageEquationHelper<CACHEDPROBLEM>
    ::set_alpha (const double alpha_new)
    {
        alpha = alpha_new;
    }

    template <class CACHEDTPROBLEM>
    LinearParabolicEquationTensor<CACHEDTPROBLEM>
    ::LinearParabolicEquationTensor(const CACHEDTPROBLEM* ellipt,
                                    CACHEDTPROBLEM* identi,
                                    const InfiniteVector<double,typename CACHEDTPROBLEM::Index>& initial,
                                    const InfiniteVector<double,typename CACHEDTPROBLEM::Index>& f,
                                    const int jmax)
    : elliptic(ellipt), identity(identi),
      constant_f_(f), f_(0), jmax_(jmax)
    {
        AbstractIVP<InfiniteVector<double,typename CACHEDTPROBLEM::Index> >::u0 = initial;
    }

    template <class CACHEDTPROBLEM>
    void
    LinearParabolicEquationTensor<CACHEDTPROBLEM>
    ::evaluate_f(const double t,
	         const InfiniteVector<double,Index>& v,
	         const double tolerance,
	         InfiniteVector<double,Index>& result) const
    {
// CLEANUP
        cout << "evaluate_f called with t = " << t << endl;

        result.clear();
        InfiniteVector<double,Index> w(v), temp;
        w.scale(elliptic, 1); // w = Dv
        APPLY(*elliptic, w, tolerance, temp, jmax_, tensor_simple); // yields -D^{-1}AD^{-1}w
        temp.scale(elliptic, 1);
        temp.scale(-1.0); // result = -D(-D^{-1}AD^{-1}Dv) = Av
        //     // multiply with inverse primal gramian (i.e., switch from dual to primal basis)
        //     G.set_rhs(temp);
        //     CDD1_SOLVE(GC, tolerance, result, jmax_);
        result = temp;
//CLEANUP
        //SampledMapping<2> s(evaluate(elliptic->basis(), result, true, 5));
        //std::ofstream ustream;
        //ustream.open("eval_f_apply.m");
        //s.matlab_output(ustream);
        //ustream.close();


        // add constant driving term (if present)
        if (!constant_f_.empty())
        {
            result.add(constant_f_);
//CLEANUP
            //ustream.open("eval_f_constant_f.m");
            //SampledMapping<2> s2(evaluate(elliptic->basis(), constant_f_, true, 5));
            //s2.matlab_output(ustream);
            //ustream.close();
        }
        // add time-dependent driving term (if present)
        if (f_ != 0)
        {
            f_->set_time(t);
            w.clear();
#if _EXPANSIONTYPE_F_
            // hier stand verwirrrenderweise: elliptic->basis().expand(f_, true, jmax_, w); // expand in the dual (!) basis
            // Zum berechnen der Entwicklungskoeffizienten bezüglich der primalen Basis:
            identity->set_f(f_);
            CDD1_SOLVE(*identity, tolerance, w, jmax_, tensor_simple);
            w.scale(identity,-1);
#else
            //identity->set_f(f_);
            elliptic->basis().expand(f_, false, jmax_, w); // expand in the primal basis
#endif
            //reactivate this line...?
            //w.compress(1e-14);
//CLEANUP
            //ustream.open("eval_f_w.m");
            //SampledMapping<2> s3(evaluate(elliptic->basis(), w, true, 5));
            //s3.matlab_output(ustream);
            //ustream.close();

            result.add(w);
        }
//CLEANUP
        
        //ustream.open("eval_f_result.m");
        //SampledMapping<2> s_4(evaluate(elliptic->basis(), result, true, 5));
        //s_4.matlab_output(ustream);
        //ustream.close();
    }

    template <class CACHEDTPROBLEM>
    void
    LinearParabolicEquationTensor<CACHEDTPROBLEM>
    ::evaluate_ft(const double t,
                  const InfiniteVector<double,Index>& v,
		  const double tolerance,
		  InfiniteVector<double,Index>& result) const
    {
// CLEANUP
        cout << "evaluate_ft called with t = " << t << endl;

        result.clear(); // from the constant driving term
        // approximate derivative of time-dependent driving term (if present)
        if (f_ != 0)
        {
#if _EXPANSIONTYPE_FT_ // compute <f',\tilde psi_lambda>
            const double h = 1e-6;
            InfiniteVector<double,Index> fhelp;
            f_->set_time(t);
            //elliptic->basis().expand(f_, false, jmax_, fhelp); // expand in the primal basis
            // verwirrenderweise: elliptic->basis().expand(f_, true, jmax_, fhelp); // expand in the dual (!) basis
            // Berechnung der dualen Entwicklung:
            identity->set_f(f_);
            CDD1_SOLVE(*identity, tolerance, fhelp, jmax_, tensor_simple); // => fhelp = (<f,\tilde\psi_\lambda>_\lambda)
            fhelp.scale(identity,-1);
            fhelp.compress(1e-14);

            f_->set_time(t+h);
            //elliptic->basis().expand(f_, false, jmax_, result);
            // verwirrenderweise: elliptic->basis().expand(f_, true, jmax_, result);
            // Berechnung der dualen Entwicklung:
            identity->set_f(f_);
            CDD1_SOLVE(*identity, tolerance, result, jmax_, tensor_simple);
            result.scale(identity,-1);
            result.compress(1e-14);
            result.add(-1., fhelp);
            result.scale(1./h);
#else // compute <f'psi_lambda>
            const double h = 1e-6;
            InfiniteVector<double,Index> fhelp;
            f_->set_time(t);
            elliptic->basis().expand(f_, false, jmax_, fhelp); // expand in the primal basis
            fhelp.compress(1e-14);
            f_->set_time(t+h);
            elliptic->basis().expand(f_, false, jmax_, result);
            result.compress(1e-14);
            result.add(-1., fhelp);
            result.scale(1./h);
#endif
        }
    }

    template <class CACHEDTPROBLEM>
    void
    LinearParabolicEquationTensor<CACHEDTPROBLEM>
    ::solve_ROW_stage_equation(const double t,
			       const InfiniteVector<double,Index>& v,
			       const double alpha,
			       const InfiniteVector<double,Index>& y,
			       const double tolerance,
			       InfiniteVector<double,Index>& result) const
    {
// CLEANUP
        cout << "solve_ROW_stage_equation called with t = " << t << endl;
        
//      // multiply everything with the Gramian
//      InfiniteVector<double,Index> Gy;
//      APPLY(GC, y, tolerance, Gy, jmax_, St04a);
//      LinParEqROWStageEquationHelper<CACHEDTPROBLEM> helper(alpha, elliptic, GC, Gy);
        LinParEqTenROWStageEquationHelper<CACHEDTPROBLEM> helper(alpha, elliptic, identity, y);
        CDD1_SOLVE(helper, tolerance, result, jmax_); // D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y
        result.scale(&helper, -1); // Dx -> x
    }

}
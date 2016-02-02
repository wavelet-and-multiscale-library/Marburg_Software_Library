// implementation for cdd1.h
namespace WaveletTL
{
  
    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // start with zero guess
        CDD1_SOLVE(P, epsilon,
                InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>(),
                u_epsilon, jmax, strategy);
    }
    
    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            InfiniteVector<double, int>& u_epsilon,
            const int jmax,
            const CompressionStrategy strategy)
    {
        CDD1_SOLVE(P, epsilon, InfiniteVector<double, int>(), u_epsilon,
                1.0/P.norm_Ainv(), P.norm_A(),
                jmax, strategy);
    }
    
    template <class PROBLEM>
    void CDD1_SOLVE_LOGGED(std::ofstream& logstream,
            PROBLEM& P, const double epsilon,
            InfiniteVector<double, int>& u_epsilon,
            const int jmax,
            const CompressionStrategy strategy)
    {
        CDD1_SOLVE_LOGGED(logstream,
                P, epsilon, InfiniteVector<double, int>(), u_epsilon,
                1.0/P.norm_Ainv(), P.norm_A(),
                jmax, strategy);
    }

    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // start with pessimistic parameters c1, c2
        CDD1_SOLVE(P, epsilon, guess, u_epsilon,
                1.0/P.norm_Ainv(), P.norm_A(),
                jmax, strategy);
    }

    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            const InfiniteVector<double, int>& guess,
            InfiniteVector<double, int>& u_epsilon,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // start with pessimistic parameters c1, c2
        CDD1_SOLVE(P, epsilon, guess, u_epsilon,
                1.0/P.norm_Ainv(), P.norm_A(),
                jmax, strategy);
    }
    
    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
            const double c1,
            const double c2,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // INIT, cf. [BB+] 
        CDD1Parameters params;
        params.c1 = c1;
        params.c2 = c2;
        params.kappa = params.c2/params.c1;
        params.gamma = 0.8;
        params.F = P.F_norm();
        // determination of q=q1=q2=q3,q4 according to [CDD1, (7.23)ff]
        params.q4 = 1. / (20. * params.kappa);
        const double A = params.c1 / (20. * (3. + params.c1 / params.c2));
        const double B = params.c2 * (0.1 - params.q4 * sqrt(params.kappa));
        const double C = params.q4 / (1 / params.c2 + 6. * (params.gamma + 1.) / (params.gamma * params.c1));
        params.q1 = params.q2 = params.q3 = std::min(A, std::min(B, C));
        params.q0 = sqrt(params.kappa) + params.q3/params.c2;
        params.theta = sqrt(1 - params.c1 * params.gamma * params.gamma / (4. * params.c2));
        params.theta_bar = 1 - 1. / (6. * params.kappa);
        params.K = (unsigned int) floor(log(20 * params.kappa) / fabs(log(params.theta))) + 1;
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE parameters:" << endl;
        cout << "c1=" << params.c1 << ", c2=" << params.c2 << ", kappa=" << params.kappa << endl;
        cout << "gamma=" << params.gamma << endl;
        cout << "F=" << params.F << endl;
        cout << "q0=" << params.q0 << ", q1=" << params.q1 << ", q2=" << params.q2
                << ", q3=" << params.q3 << ", q4=" << params.q4 << endl;
        cout << "theta=" << params.theta << ", theta_bar=" << params.theta_bar << endl;
        cout << "K=" << params.K << endl;
#endif
        typedef typename PROBLEM::WaveletBasis::Index Index;
        set<Index> Lambda, Lambda_hat;
        u_epsilon = guess;
        u_epsilon.support(Lambda);
        double delta = params.F;
        InfiniteVector<double,Index> v_hat, r_hat, u_bar, F;
        P.RHS(2*params.q2*epsilon, F);
        while (delta > epsilon) { // sqrt(params.c1)*epsilon) { // check the additional factor c1^{1/2} in [BB+] !?
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "CDD1_SOLVE: delta=" << delta << endl;
#endif
            NPROG(P, params, F, Lambda, u_epsilon, delta, v_hat, Lambda_hat, r_hat, u_bar, jmax, strategy);
            if (l2_norm(r_hat)+(params.q1+params.q2+(1+1./params.kappa)*params.q3)*delta <= params.c1*epsilon)
            {
                u_epsilon.swap(u_bar);
                break;
            }
            else
            {
                u_epsilon.swap(v_hat);
                Lambda.swap(Lambda_hat);
            }
            //       delta *= 0.5; // original
            delta *= 0.1; // tuned
        }
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE: done!" << endl;
#endif
    }

    template <class PROBLEM>
    void CDD1_SOLVE(PROBLEM& P, const double epsilon,
            const InfiniteVector<double, int>& guess,
            InfiniteVector<double, int>& u_epsilon,
            const double c1,
            const double c2,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // INIT, cf. [BB+] 
        CDD1Parameters params;
        params.c1 = c1;
        params.c2 = c2;
        params.kappa = params.c2/params.c1;
        params.gamma = 0.8;
        params.F = P.F_norm();
        // determination of q=q1=q2=q3,q4 according to [CDD1, (7.23)ff]
        params.q4 = 1. / (20. * params.kappa);
        const double A = params.c1 / (20. * (3. + params.c1 / params.c2));
        const double B = params.c2 * (0.1 - params.q4 * sqrt(params.kappa));
        const double C = params.q4 / (1 / params.c2 + 6. * (params.gamma + 1.) / (params.gamma * params.c1));
        params.q1 = params.q2 = params.q3 = std::min(A, std::min(B, C));
        params.q0 = sqrt(params.kappa) + params.q3/params.c2;
        params.theta = sqrt(1 - params.c1 * params.gamma * params.gamma / (4. * params.c2));
        params.theta_bar = 1 - 1. / (6. * params.kappa);
        params.K = (unsigned int) floor(log(20 * params.kappa) / fabs(log(params.theta))) + 1;
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE parameters:" << endl;
        cout << "c1=" << params.c1 << ", c2=" << params.c2 << ", kappa=" << params.kappa << endl;
        cout << "gamma=" << params.gamma << endl;
        cout << "F=" << params.F << endl;
        cout << "q0=" << params.q0 << ", q1=" << params.q1 << ", q2=" << params.q2
                << ", q3=" << params.q3 << ", q4=" << params.q4 << endl;
        cout << "theta=" << params.theta << ", theta_bar=" << params.theta_bar << endl;
        cout << "K=" << params.K << endl;
#endif
        typedef typename PROBLEM::WaveletBasis::Index Index;
        set<int> Lambda, Lambda_hat;
        u_epsilon = guess;
        u_epsilon.support(Lambda);
        double delta = params.F;
        InfiniteVector<double,int> v_hat, r_hat, u_bar, F;
        cout << "CDD1:: 2*params.q2*epsilon = " << (2*params.q2*epsilon) << endl;
        P.RHS(2*params.q2*epsilon, F);
        while (delta > epsilon) { // sqrt(params.c1)*epsilon) { // check the additional factor c1^{1/2} in [BB+] !?
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "CDD1_SOLVE: delta=" << delta << endl;
#endif
            NPROG(P, params, F, Lambda, u_epsilon, delta, v_hat, Lambda_hat, r_hat, u_bar, jmax, strategy);
            if (l2_norm(r_hat)+(params.q1+params.q2+(1+1./params.kappa)*params.q3)*delta <= params.c1*epsilon)
            {
                u_epsilon.swap(u_bar);
                break;
            }
            else
            {
                u_epsilon.swap(v_hat);
                Lambda.swap(Lambda_hat);
            }
            //       delta *= 0.5; // original
            delta *= 0.1; // tuned
        }
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE: done!" << endl;
#endif
    }
    
    template <class PROBLEM>
    void CDD1_SOLVE_LOGGED(std::ofstream& logstream, 
            PROBLEM& P, const double epsilon,
            const InfiniteVector<double, int>& guess,
            InfiniteVector<double, int>& u_epsilon,
            const double c1,
            const double c2,
            const int jmax,
            const CompressionStrategy strategy)
    {
        // INIT, cf. [BB+] 
        CDD1Parameters params;
        params.c1 = c1;
        params.c2 = c2;
        params.kappa = params.c2/params.c1;
        params.gamma = 0.8;
        params.F = P.F_norm();
        // determination of q=q1=q2=q3,q4 according to [CDD1, (7.23)ff]
        params.q4 = 1. / (20. * params.kappa);
        const double A = params.c1 / (20. * (3. + params.c1 / params.c2));
        const double B = params.c2 * (0.1 - params.q4 * sqrt(params.kappa));
        const double C = params.q4 / (1 / params.c2 + 6. * (params.gamma + 1.) / (params.gamma * params.c1));
        params.q1 = params.q2 = params.q3 = std::min(A, std::min(B, C));
        params.q0 = sqrt(params.kappa) + params.q3/params.c2;
        params.theta = sqrt(1 - params.c1 * params.gamma * params.gamma / (4. * params.c2));
        params.theta_bar = 1 - 1. / (6. * params.kappa);
        params.K = (unsigned int) floor(log(20 * params.kappa) / fabs(log(params.theta))) + 1;
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE parameters:" << endl;
        cout << "c1=" << params.c1 << ", c2=" << params.c2 << ", kappa=" << params.kappa << endl;
        cout << "gamma=" << params.gamma << endl;
        cout << "F=" << params.F << endl;
        cout << "q0=" << params.q0 << ", q1=" << params.q1 << ", q2=" << params.q2
                << ", q3=" << params.q3 << ", q4=" << params.q4 << endl;
        cout << "theta=" << params.theta << ", theta_bar=" << params.theta_bar << endl;
        cout << "K=" << params.K << endl;
#endif
        typedef typename PROBLEM::WaveletBasis::Index Index;
        set<int> Lambda, Lambda_hat;
        u_epsilon = guess;
        u_epsilon.support(Lambda);
        double delta = params.F;
        InfiniteVector<double,int> v_hat, r_hat, u_bar, F;
        cout << "CDD1:: 2*params.q2*epsilon = " << (2*params.q2*epsilon) << endl;
        P.RHS(2*params.q2*epsilon, F);
        double res, Dres;
        InfiniteVector<double,int> Dr_hat;
        
        while (delta > epsilon) { // sqrt(params.c1)*epsilon) { // check the additional factor c1^{1/2} in [BB+] !?
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "CDD1_SOLVE: delta=" << delta << endl;
#endif
            NPROG(P, params, F, Lambda, u_epsilon, delta, v_hat, Lambda_hat, r_hat, u_bar, jmax, strategy);
            res = l2_norm(r_hat);
            Dr_hat = r_hat;
            Dr_hat.scale(&P,1);
            Dres = l2_norm(Dr_hat);
            logstream << "erg = [erg ; " << v_hat.size() << " " << res << " " << Dres << "];" << endl;
            cout << "CDD1_SOLVE_LOGGED: active coeffs = " << v_hat.size() << "; current residual: res = " << res << "; scaled with diagonal = Dres = " << Dres << endl;
            
            if (res+(params.q1+params.q2+(1+1./params.kappa)*params.q3)*delta <= params.c1*epsilon)
            {
                u_epsilon.swap(u_bar);
                break;
            }
            else
            {
                u_epsilon.swap(v_hat);
                Lambda.swap(Lambda_hat);
            }
            //       delta *= 0.5; // original
            delta *= 0.1; // tuned
        }
        
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "CDD1_SOLVE: done!" << endl;
#endif
    }
    
    template <class PROBLEM>
    void NPROG(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
            const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
            const double delta,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v_hat,
            set<typename PROBLEM::WaveletBasis::Index>& Lambda_hat,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r_hat,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_Lambda_k,
            const int jmax,
            const CompressionStrategy strategy)
    {
        typedef typename PROBLEM::WaveletBasis::Index Index;
        set<Index> Lambda_k(Lambda), Lambda_kplus1;
        unsigned int k = 0;
        GALERKIN(P, params, F, Lambda_k, v, delta, params.q3*delta/params.c2, u_Lambda_k, jmax, strategy);
        while (true) 
        {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "NPROG: k=" << k << " (K=" << params.K << ")" << endl;
#endif
            NGROW(P, params, F, Lambda_k, u_Lambda_k, params.q1*delta, params.q2*delta, Lambda_kplus1, r_hat, jmax, strategy);
            if (l2_norm(r_hat) <= params.c1*delta/20. || k == params.K || Lambda_k.size() == Lambda_kplus1.size()) {
                //  	u_Lambda_k.COARSE(2.*delta/5., v_hat);
                v_hat = u_Lambda_k;
                v_hat.support(Lambda_hat);
                break;
            }
            GALERKIN(P, params, F, Lambda_kplus1, u_Lambda_k, params.q0*delta, params.q3*delta/params.c2, v_hat, jmax, strategy);
            u_Lambda_k.swap(v_hat);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
            cout << "u_Lambda_k = " << endl << u_Lambda_k << endl;
#endif
            Lambda_k.swap(Lambda_kplus1);
            k++;
        }
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "NPROG: done" << endl;
#endif
    }
    
    template <class PROBLEM>
    void NPROG(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, int>& F,
            const set<int>& Lambda,
            const InfiniteVector<double, int>& v,
            const double delta,
            InfiniteVector<double, int>& v_hat,
            set<int>& Lambda_hat,
            InfiniteVector<double, int>& r_hat,
            InfiniteVector<double, int>& u_Lambda_k,
            const int jmax,
            const CompressionStrategy strategy)
    {
        set<int> Lambda_k(Lambda), Lambda_kplus1;
        unsigned int k = 0;
        GALERKIN(P, params, F, Lambda_k, v, delta, params.q3*delta/params.c2, u_Lambda_k, jmax, strategy);
        while (true) 
        {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "NPROG: k=" << k << " (K=" << params.K << ")" << endl;
#endif
            NGROW(P, params, F, Lambda_k, u_Lambda_k, params.q1*delta, params.q2*delta, Lambda_kplus1, r_hat, jmax, strategy);
            if (l2_norm(r_hat) <= params.c1*delta/20. || k == params.K || Lambda_k.size() == Lambda_kplus1.size()) {
                //  	u_Lambda_k.COARSE(2.*delta/5., v_hat);
                v_hat = u_Lambda_k;
                v_hat.support(Lambda_hat);
                break;
            }
            GALERKIN(P, params, F, Lambda_kplus1, u_Lambda_k, params.q0*delta, params.q3*delta/params.c2, v_hat, jmax, strategy);
            u_Lambda_k.swap(v_hat);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
            cout << "u_Lambda_k = " << endl << u_Lambda_k << endl;
#endif
            Lambda_k.swap(Lambda_kplus1);
            k++;
        }
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "NPROG: done" << endl;
#endif
    }
    
    template <class PROBLEM>
    void GALERKIN(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
            const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
            const double delta,
            const double eta,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_bar,
            const int jmax,
            const CompressionStrategy strategy)
    {
        typedef typename PROBLEM::WaveletBasis::Index Index;
#if 0
        // original GALERKIN version from [CDD1],[BB+]
        cout << "GALERKIN called..." << endl;
        InfiniteVector<double,Index> r;
        u_bar = v;
        double mydelta = delta;
        cout << "GALERKIN, internal residuals: " << endl;
        while (true) 
        {
            INRESIDUAL(P, params, F, Lambda, u_bar, params.c1*eta/6., params.c1*eta/6., r, jmax, strategy);
            const double inresidual_norm = l2_norm(r);
            cout << inresidual_norm << " ";
            cout.flush();
            if (eta >= std::min(params.theta_bar*mydelta, inresidual_norm/params.c1+eta/3.)) 
            {
                cout << "... GALERKIN done, norm of internal residual: " << inresidual_norm << endl;
                break;
            }
            //       u_bar += 1/params.c2 * r; // original [CDD1] relaxation parameter
            u_bar += 2./(params.c2+1./params.c1) * r; // optimal relaxation parameter
            mydelta *= params.theta_bar;
        }
#else
        // conjugate gradient version (no theory for its complexity available yet, but very fast)
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "GALERKIN called..." << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "... with Lambda=" << endl;
        for (typename set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                it != itend; ++it) 
        {
            cout << *it << endl;
        }
#endif
        u_bar.clear();
        if (Lambda.size() > 0) 
        {
            // setup A_Lambda and f_Lambda
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(P, Lambda, A_Lambda);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
            cout << "... GALERKIN: A_Lambda=" << endl << A_Lambda;
#endif
            Vector<double> F_Lambda;
            setup_righthand_side(P, Lambda, F_Lambda);
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN: stiffness matrix and right-hand side set up, iterating ..." << endl;
#endif
            // setup initial approximation xk
            Vector<double> xk(Lambda.size());
            unsigned int id = 0;
            for (typename set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                    it != itend; ++it, ++id)
                xk[id] = v.get_coefficient(*it);
            unsigned int iterations = 0;
            //       CG(A_Lambda, F_Lambda, xk, eta, 150, iterations);
            CG(A_Lambda, F_Lambda, xk, 1e-15, 250, iterations);
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN done, " << iterations << " CG iterations needed" << endl;
#endif
            id = 0;
            for (typename set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                    it != itend; ++it, ++id)
                u_bar.set_coefficient(*it, xk[id]);
        } 
        else 
        {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN done, no CG iteration needed" << endl;
#endif
        }
#endif
    }
    
    template <class PROBLEM>
    void GALERKIN(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, int>& F,
            const set<int>& Lambda,
            const InfiniteVector<double, int>& v,
            const double delta,
            const double eta,
            InfiniteVector<double, int>& u_bar,
            const int jmax,
            const CompressionStrategy strategy)
    {
#if 0
        // original GALERKIN version from [CDD1],[BB+]
        cout << "GALERKIN called..." << endl;
        InfiniteVector<double,int> r;
        u_bar = v;
        double mydelta = delta;
        cout << "GALERKIN, internal residuals: " << endl;
        while (true) 
        {
            INRESIDUAL(P, params, F, Lambda, u_bar, params.c1*eta/6., params.c1*eta/6., r, jmax, strategy);
            const double inresidual_norm = l2_norm(r);
            cout << inresidual_norm << " ";
            cout.flush();
            if (eta >= std::min(params.theta_bar*mydelta, inresidual_norm/params.c1+eta/3.)) 
            {
                cout << "... GALERKIN done, norm of internal residual: " << inresidual_norm << endl;
                break;
            }
            //       u_bar += 1/params.c2 * r; // original [CDD1] relaxation parameter
            u_bar += 2./(params.c2+1./params.c1) * r; // optimal relaxation parameter
            mydelta *= params.theta_bar;
        }
#else
        // conjugate gradient version (no theory for its complexity available yet, but very fast)
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "GALERKIN called..." << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "... with Lambda=" << endl;
        for (typename set<int>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                it != itend; ++it) 
        {
            cout << *it << endl;
        }
#endif
        u_bar.clear();
        if (Lambda.size() > 0) 
        {
            // setup A_Lambda and f_Lambda
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(P, Lambda, A_Lambda);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
            cout << "... GALERKIN: A_Lambda=" << endl << A_Lambda;
#endif
            Vector<double> F_Lambda;
            setup_righthand_side(P, Lambda, F_Lambda);
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN: stiffness matrix and right-hand side set up, iterating ..." << endl;
#endif
            // setup initial approximation xk
            Vector<double> xk(Lambda.size());
            unsigned int id = 0;
            for (typename set<int>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                    it != itend; ++it, ++id)
                xk[id] = v.get_coefficient(*it);
            unsigned int iterations = 0;
            //       CG(A_Lambda, F_Lambda, xk, eta, 150, iterations);
            CG(A_Lambda, F_Lambda, xk, 1e-15, 250, iterations);
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN done, " << iterations << " CG iterations needed" << endl;
#endif
            id = 0;
            for (typename set<int>::const_iterator it = Lambda.begin(), itend = Lambda.end();
                    it != itend; ++it, ++id)
                u_bar.set_coefficient(*it, xk[id]);
        } 
        else 
        {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
            cout << "... GALERKIN done, no CG iteration needed" << endl;
#endif
        }
#endif
    }
    

    template <class PROBLEM>
    void NGROW(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
            const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_bar,
            const double xi1,
            const double xi2,
            set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
            const int jmax,
            const CompressionStrategy strategy)
    {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "NGROW called..." << endl;
#endif
        typedef typename PROBLEM::WaveletBasis::Index Index;
        set<Index> Lambda_c;
        NRESIDUAL(P, params, F, Lambda, u_bar, xi1, xi2, r, Lambda_c, jmax, strategy);
        const double residual_norm = l2_norm(r);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "* NGROW: current residual is " << endl << r << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "* NGROW: current residual norm is " << residual_norm << endl;
#endif
        if (residual_norm > 1000)
        {
            cout << "* NGROW: residual norm is greater than 1000, some debugging output:" << endl;
            cout << "u_bar=" << endl << u_bar << endl;
            cout << "xi1=" << xi1 << ", xi2=" << xi2 << endl;
#ifdef ENABLE_NGROW_THROW_EXCEPTION
            std::stringstream ss;
            ss << "* NGROW: residual norm is greater than 1000!\nSome debugging output:\n"
               << "u_bar=\n" << u_bar << "\nxi1=" << xi1 << ", xi2=" << xi2 << std::endl;
            throw std::runtime_error(ss.str());
#else
            abort();
#endif
        }
        InfiniteVector<double,Index> pr;
        r.COARSE(sqrt(1-params.gamma*params.gamma)*residual_norm, pr);
        pr.support(Lambda_c);
        Lambda_tilde.clear();
        
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        const unsigned int Lambdasize = Lambda.size();
        const unsigned int Lambdacsize = Lambda_c.size();
        cout << "* NGROW: size of old index set Lambda is " << Lambdasize << ", size of increment set Lambda_c is " << Lambdacsize << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "* NGROW: old index set Lambda is..." << endl;
        for (typename set<typename PROBLEM::WaveletBasis::Index>::const_iterator it(Lambda.begin());
                it != Lambda.end(); ++it)
            cout << *it << endl;
        cout << "* NGROW: increment set Lambda_c is..." << endl;
        for (typename set<typename PROBLEM::WaveletBasis::Index>::const_iterator it(Lambda_c.begin());
                it != Lambda_c.end(); ++it)
            cout << *it << endl;
#endif
        std::set_union(Lambda.begin(), Lambda.end(),
                Lambda_c.begin(), Lambda_c.end(),
                inserter(Lambda_tilde, Lambda_tilde.end()));
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "... NGROW done, size of new index set: " << Lambda_tilde.size() << endl;
#endif
    }

    template <class PROBLEM>
    void NGROW(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, int>& F,
            const set<int>& Lambda,
            const InfiniteVector<double, int>& u_bar,
            const double xi1,
            const double xi2,
            set<int>& Lambda_tilde,
            InfiniteVector<double, int>& r,
            const int jmax,
            const CompressionStrategy strategy)
    {
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "NGROW called..." << endl;
#endif
        set<int> Lambda_c;
        NRESIDUAL(P, params, F, Lambda, u_bar, xi1, xi2, r, Lambda_c, jmax, strategy);
        const double residual_norm = l2_norm(r);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "* NGROW: current residual is " << endl << r << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "* NGROW: current residual norm is " << residual_norm << endl;
#endif
        if (residual_norm > 1000)
        {
            cout << "* NGROW: residual norm is greater than 1000, some debugging output:" << endl;
            cout << "u_bar=" << endl << u_bar << endl;
            cout << "xi1=" << xi1 << ", xi2=" << xi2 << endl;
#ifdef ENABLE_NGROW_THROW_EXCEPTION
            std::stringstream ss;
            ss << "* NGROW: residual norm is greater than 1000!\nSome debugging output:\n"
               << "u_bar=\n" << u_bar << "\nxi1=" << xi1 << ", xi2=" << xi2 << std::endl;
            throw std::runtime_error(ss.str());
#else
            abort();
#endif
        }
        InfiniteVector<double,int> pr;
        r.COARSE(sqrt(1-params.gamma*params.gamma)*residual_norm, pr);
        pr.support(Lambda_c);
        Lambda_tilde.clear();
        
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        const unsigned int Lambdasize = Lambda.size();
        const unsigned int Lambdacsize = Lambda_c.size();
        cout << "* NGROW: size of old index set Lambda is " << Lambdasize << ", size of increment set Lambda_c is " << Lambdacsize << endl;
#endif
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "* NGROW: old index set Lambda is..." << endl;
        for (typename set<int>::const_iterator it(Lambda.begin());
                it != Lambda.end(); ++it)
            cout << *it << endl;
        cout << "* NGROW: increment set Lambda_c is..." << endl;
        for (typename set<int>::const_iterator it(Lambda_c.begin());
                it != Lambda_c.end(); ++it)
            cout << *it << endl;
#endif
        std::set_union(Lambda.begin(), Lambda.end(),
                Lambda_c.begin(), Lambda_c.end(),
                inserter(Lambda_tilde, Lambda_tilde.end()));
#if _WAVELETTL_CDD1_VERBOSITY >= 1
        cout << "... NGROW done, size of new index set: " << Lambda_tilde.size() << endl;
#endif
    }
    
    template <class PROBLEM>
    void INRESIDUAL(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
            const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
            const double eta1,
            const double eta2,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
            const int jmax,
            const CompressionStrategy strategy)
    {
        typedef typename PROBLEM::WaveletBasis::Index Index;
        InfiniteVector<double,Index> w, g(F);
        // TODO: speed up the following two lines
#if _WAVELETTL_USE_TBASIS == 1
        APPLY_TENSOR(P,v,eta1,w,jmax,strategy);
#else
        APPLY(P, v, eta1, w, jmax, strategy);
#endif
        w.clip(Lambda);
        g.clip(Lambda);
        g.COARSE(eta2, r);
        r -= w;
    }
    
    template <class PROBLEM>
    void INRESIDUAL(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, int>& F,
            const set<int>& Lambda,
            const InfiniteVector<double, int>& v,
            const double eta1,
            const double eta2,
            InfiniteVector<double, int>& r,
            const int jmax,
            const CompressionStrategy strategy)
    {
        InfiniteVector<double,int> w, g(F);
        // TODO: speed up the following two lines
#if _WAVELETTL_USE_TBASIS == 1
        APPLY_TENSOR(P,v,eta1,w,jmax,strategy);
#else
        APPLY(P, v, eta1, w, jmax, strategy);
#endif
        w.clip(Lambda);
        g.clip(Lambda);
        g.COARSE(eta2, r);
        r -= w;
    }
  
    template <class PROBLEM>
    void NRESIDUAL(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& F,
            const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
            const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
            const double eta1,
            const double eta2,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& r,
            set<typename PROBLEM::WaveletBasis::Index>& Lambda_tilde,
            const int jmax,
            const CompressionStrategy strategy)
    {
        typedef typename PROBLEM::WaveletBasis::Index Index;
        InfiniteVector<double,Index> w;
#if _WAVELETTL_USE_TBASIS == 1
        APPLY_TENSOR(P,v,eta1,w,jmax,strategy);
#else
        APPLY(P, v, eta1, w, jmax, strategy);
#endif
        F.COARSE(eta2, r);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "NRESIDUAL: v = "<<endl << v << endl<<"NRESIDUAL: w = "<< endl << w << endl<<"NRESIDUAL: F.coarse("<<eta2<<")= "<<endl <<r<<endl;
#endif
        r -= w;
        r.support(Lambda_tilde);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "NRESIDUAL: r = "<<endl << r << endl;
        cout << "NRESIDUAL: Lambda_tilde is..." << endl;
        for (typename set<typename PROBLEM::WaveletBasis::Index>::const_iterator it(Lambda_tilde.begin());
                it != Lambda_tilde.end(); ++it)
            cout << *it << endl;
#endif
    }
    
    template <class PROBLEM>
    void NRESIDUAL(PROBLEM& P, const CDD1Parameters& params,
            const InfiniteVector<double, int>& F,
            const set<int>& Lambda,
            const InfiniteVector<double, int>& v,
            const double eta1,
            const double eta2,
            InfiniteVector<double, int>& r,
            set<int>& Lambda_tilde,
            const int jmax,
            const CompressionStrategy strategy)
    {
        InfiniteVector<double, int> w;
#if _WAVELETTL_USE_TBASIS == 1
        APPLY_TENSOR(P,v,eta1,w,jmax,strategy, true);
#else
        APPLY(P, v, eta1, w, jmax, strategy);
#endif
        F.COARSE(eta2, r);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "NRESIDUAL: v = "<<endl << v << endl<<"NRESIDUAL: w = "<< endl << w << endl<<"NRESIDUAL: F.coarse("<<eta2<<")= "<<endl <<r<<endl;
#endif
        r -= w;
        r.support(Lambda_tilde);
#if _WAVELETTL_CDD1_VERBOSITY >= 2
        cout << "NRESIDUAL: r = "<<endl << r << endl;
        cout << "NRESIDUAL: Lambda_tilde is..." << endl;
        for (typename set<int>::const_iterator it(Lambda_tilde.begin());
                it != Lambda_tilde.end(); ++it)
            cout << *it << endl;
#endif
    }
}

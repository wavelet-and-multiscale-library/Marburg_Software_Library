// -*- c++ -*-

// +--------------------------------------------------------------------+
// | stevenson_AWGM.cpp, Copyright (c) 2018                             |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of WaveletTL - the Wavelet Template Library.     |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+


// implementation for stevenson_AWGM.h


namespace WaveletTL
{


template <class PROBLEM>
void AWGM_SOLVE(const PROBLEM& P, const double epsilon,
               InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
               const int jmax,
               MathTL::AbstractConvergenceLogger& logger,
               const double alpha,
               const double omega,
               const double gamma,
               const double theta,
               const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
               const CompressionStrategy strategy)
{
    // Start with nu_{-1} = ||F||_2
    AWGM_SOLVE(P, epsilon, u_epsilon, jmax, P.F_norm(), logger, alpha, omega, gamma, theta, guess, strategy);
}



template <class PROBLEM>
void AWGM_SOLVE(const PROBLEM& P, const double epsilon,
               InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
               const int jmax,
               const double nu_neg1,
               MathTL::AbstractConvergenceLogger& logger,
               const double alpha,
               const double omega,
               const double gamma,
               const double theta,
               const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
               const CompressionStrategy strategy)
{
    unsigned int k = 0; // iteration counter
    double nu = nu_neg1;

    typedef typename PROBLEM::WaveletBasis::Index Index;

    set<Index> Lambda, supp_r_coarse;
    InfiniteVector<double, Index> r, r_help, g;

    u_epsilon = guess;
    u_epsilon.support(Lambda);

    logger.startClock();

    while(true)
    {
        logger.checkAbortConditions();

        unsigned int res_loop_counter = 0;
        RES(P, u_epsilon, theta*nu*omega/(1-omega), omega, epsilon, jmax, r, nu, res_loop_counter, strategy, false);

        cout << "GHS_SOLVE: k=" << k << ", nu=" << nu << endl;
        cout << "       (epsilon=" << epsilon << "), support size: " << u_epsilon.size() << endl;
        cout << "       Number of loops needed in RES: " << res_loop_counter << endl;

        logger.logConvergenceData(u_epsilon.size(), nu);

        if(nu <= epsilon) break;

        double norm_r = l2_norm(r);
        r_help = r;
        r_help.clip(Lambda);
        r -= r_help;
        r.COARSE(sqrt(1-alpha*alpha)*norm_r, r_help);
        r_help.support(supp_r_coarse);
        Lambda.insert(supp_r_coarse.begin(), supp_r_coarse.end());

        P.RHS(gamma*nu, g);
        g.clip(Lambda);
        GALSOLVE(P, Lambda, g, u_epsilon, (1+gamma)*nu, gamma*nu);

        ++k;
    }

    cout << "GHS_SOLVE: done!" << endl;
}



template <class PROBLEM>
void GALSOLVE(const PROBLEM& P, const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
              const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& g_Lambda,
              InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& w_Lambda,
              const double delta,
              const double epsilon)
{
    typedef typename PROBLEM::WaveletBasis::Index Index;

    // setup A_Lambda
    SparseMatrix<double> A_Lambda;
    setup_stiffness_matrix(P, Lambda, A_Lambda);

    // setup right-hand side
    Vector<double> g(Lambda.size());
    unsigned int id = 0;
    for (typename set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end()); it != itend; ++it, ++id)
        g[id] = g_Lambda.get_coefficient(*it);

    // setup initial approximation xk
    Vector<double> xk(Lambda.size());
    id = 0;
    for (typename set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end()); it != itend; ++it, ++id)
        xk[id] = w_Lambda.get_coefficient(*it);

#if _WAVELETTL_GHS_VERBOSITY >= 2
    cout << "       GALSOLVE: stiffness matrix and right-hand side set up, iterating ..." << endl;
#endif

    unsigned int iterations = 0;
    if(!MathTL::CG(A_Lambda, g, xk, epsilon/delta, 250, iterations))
        cout << "GALSOLVE: CG could not reach tolerance within 250 iterations!" << endl;

    id = 0;
    w_Lambda.clear();

    for (typename set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end()); it != itend; ++it, ++id)
            w_Lambda.set_coefficient(*it, xk[id]);

#if _WAVELETTL_GHS_VERBOSITY >= 2
    cout << "       ... GALSOLVE done, " << iterations << " CG iterations needed" << endl;
#endif
}


}

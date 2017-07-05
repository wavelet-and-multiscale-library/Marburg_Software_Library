// implementation for tframe_equation.h



namespace WaveletTL
{
    //PERFORMANCE : wozu compute_rhs? welches jmax in den Konstruktoren? (wenn überhaupt)
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(EllipticBVP<DIM>* bvp,
                                                           const FixedArray1D<bool,2*DIM>& bc,
                                                           const bool precompute)
    : bvp_(bvp), frame_(bc), normA(0.0), normAinv(0.0)
    {
        if (precompute == true)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
            compute_rhs();
        }
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const EllipticBVP<DIM>* bvp,
                                                     const FixedArray1D<int,2*DIM>& bc,
                                                     const bool precompute)
    : bvp_(bvp), frame_(bc), normA(0.0), normAinv(0.0)
    {
        if (precompute == true)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0); // for a first quick hack
            compute_rhs();
        }
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const TensorFrameEquation& eq)
    : bvp_(eq.bvp_), frame_(eq.frame_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
    {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
    }

// TODO PERFORMANCE:: use setup_full_collection entries
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::compute_rhs()
    {
        //cout << "TensorFrameEquation(): precompute right-hand side..." << endl;
        typedef typename QuarkletFrame::Index Index;
        // precompute the right-hand side on a fine level
        InfiniteVector<double,Index> fhelp;
        fnorm_sqr = 0;
        for (unsigned int i = 0; (int)i< frame_.degrees_of_freedom();i++)
        {
            const double coeff = f(frame_.get_quarklet(i)) / D(frame_.get_quarklet(i));
            if (fabs(coeff)>1e-15)
            {
                fhelp.set_coefficient(frame_.get_quarklet(i), coeff);
                fnorm_sqr += coeff*coeff;
            }
        }
        //cout << "... done, sort the entries in modulus..." << endl;
        // sort the coefficients into fcoeffs
        fcoeffs.resize(0); // clear eventual old values
        fcoeffs.resize(fhelp.size());
        unsigned int id(0);
        for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                it != itend; ++it, ++id)
        {
            fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
//            cout << it.index() << ", " << *it << endl;
        }
        sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
        //cout << "... done, all integrals for right-hand side computed!" << endl;
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    inline
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::D(const typename QuarkletFrame::Index& lambda) const
    {
//        return sqrt(a(lambda, lambda));
        double hspreconditioner(0), l2preconditioner(1);
        for (int i=0; i<space_dimension; i++){
            hspreconditioner+=pow(1+lambda.p()[i],8)*ldexp(1.0,2*lambda.j()[i]);
            l2preconditioner*=pow(1+lambda.p()[i],2);
        }
        double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
        return preconditioner;
        
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    inline
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::a(const typename QuarkletFrame::Index& lambda,
                                          const typename QuarkletFrame::Index& mu) const
    {
        return a(lambda, mu, IFRAME::primal_polynomial_degree()*IFRAME::primal_polynomial_degree());
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::a(const typename QuarkletFrame::Index& lambda,
                                              const typename QuarkletFrame::Index& mu,
                                              const unsigned int p) const
    {
        // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
        double r = 0;
        // first decide whether the supports of psi_lambda and psi_mu intersect
        typedef typename QuarkletFrame::Support Support;
        Support supp;
        if (intersect_supports(frame_, lambda, mu, supp))
        {
//            cout << "support in a-routine for index " << lambda << ", " << mu << endl;
//    cout << "j[0]: " << supp.j[0] << ", j[1]: " << supp.j[1] << endl;
//    
//            cout << "Vergleichswerte Cube: [" << supp.a[0] << ", " << supp.b[0] << "] x [" << supp.a[1] << ", " << supp.b[1] <<"]" << endl;
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = (p+1)/2+(multi_degree(lambda.p())+multi_degree(mu.p())+1)/2;
            //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
            // FixedArray1D<double,DIM> h; // granularity for the quadrature
            double hi; // granularity for the quadrature
            FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
            for (unsigned int i = 0; i < DIM; i++) {
                hi = ldexp(1.0, -supp.j[i]);
                gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
                                = hi*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
                                = hi*GaussWeights[N_Gauss-1][n];
                    }
            }
            // compute point values of the integrand (where we use that it is a tensor product)
            FixedArray1D<Array1D<double>,DIM> psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
                                              psi_mu_values,         // -"-, for psi_mu
                                              psi_lambda_der_values, // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]
                                              psi_mu_der_values;     // -"-, for psi_mu
            for (unsigned int i = 0; i < DIM; i++) {
                evaluate(*frame_.frames()[i], 0,
                                                lambda.p()[i],
                                                lambda.j()[i],
                                                lambda.e()[i],
                                                lambda.k()[i],
                         gauss_points[i], psi_lambda_values[i]);
//                cout << "Gauss Points Cube: " << gauss_points[i] << endl;
//                cout << psi_lambda_values[i] << endl;
                evaluate(*frame_.frames()[i], 1,
                                                lambda.p()[i],
                                                lambda.j()[i],
                                                lambda.e()[i],
                                                lambda.k()[i],
                         gauss_points[i], psi_lambda_der_values[i]);
//                cout << psi_lambda_der_values[0] << endl;
                evaluate(*frame_.frames()[i], 0,
                                                mu.p()[i],
                                                mu.j()[i],
                                                mu.e()[i],
                                                mu.k()[i],
                         gauss_points[i], psi_mu_values[i]);
//                cout << psi_mu_values[i] << endl;
                evaluate(*frame_.frames()[i], 1,
                                                mu.p()[i],
                                                mu.j()[i],
                                                mu.e()[i],
                                                mu.k()[i],
                         gauss_points[i], psi_mu_der_values[i]);
//                cout << psi_mu_der_values[i] << endl;
            }
            // iterate over all points and sum up the integral shares
            int index[DIM]; // current multiindex for the point values
            for (unsigned int i = 0; i < DIM; i++)
                index[i] = 0;
            Point<DIM> x;
            const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
            const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
            double grad_psi_lambda[DIM], grad_psi_mu[DIM], weights;
            if (bvp_->constant_coefficients())
            {
                while (true)
                {
                    for (unsigned int i = 0; i < DIM; i++)
                        x[i] = gauss_points[i][index[i]];
                    // product of current Gauss weights
                    weights = 1.0;
                    for (unsigned int i = 0; i < DIM; i++)
                        weights *= gauss_weights[i][index[i]];
                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
                    for (unsigned int i = 0; i < DIM; i++)
                    {
                        grad_psi_lambda[i] = 1.0;
                        grad_psi_mu[i] = 1.0;
                        for (unsigned int s = 0; s < DIM; s++) 
                        {
                            if (i == s)
                            {
                                grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
                                grad_psi_mu[i]     *= psi_mu_der_values[i][index[i]];
                            } else
                            {
                                grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
                                grad_psi_mu[i] *= psi_mu_values[s][index[s]];
                            }
                        }
                    }
                    double share = 0;
                    for (unsigned int i = 0; i < DIM; i++)
                        share += grad_psi_lambda[i]*grad_psi_mu[i];
                    r += ax * weights * share;
                    // compute the share q(x)psi_lambda(x)psi_mu(x)
                    share = qx * weights;
                    for (unsigned int i = 0; i < DIM; i++)
                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    for (unsigned int i = 0; i < DIM; i++)
                    {
                        if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1)
                        {
                            index[i] = 0;
                            exit = (i == DIM-1);
                        } else
                        {
                            index[i]++;
                            break;
                        }
                    }
                    if (exit) break;
                }
            } else // coefficients are not constant:
            {
                while (true) {
                    for (unsigned int i = 0; i < DIM; i++)
                        x[i] = gauss_points[i][index[i]];
                    // product of current Gauss weights
                    weights = 1.0;
                    for (unsigned int i = 0; i < DIM; i++)
                        weights *= gauss_weights[i][index[i]];
                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
                    for (unsigned int i = 0; i < DIM; i++) {
                        grad_psi_lambda[i] = 1.0;
                        grad_psi_mu[i] = 1.0;
                        for (unsigned int s = 0; s < DIM; s++) {
                            if (i == s) {
                                grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
                                grad_psi_mu[i]     *= psi_mu_der_values[i][index[i]];
                            } else {
                                grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
                                grad_psi_mu[i] *= psi_mu_values[s][index[s]];
                            }
                        }
                    }
                    double share = 0;
                    for (unsigned int i = 0; i < DIM; i++)
                        share += grad_psi_lambda[i]*grad_psi_mu[i];
                    r += bvp_->a(x) * weights * share;
                    // compute the share q(x)psi_lambda(x)psi_mu(x)
                    share = bvp_->q(x) * weights;
                    for (unsigned int i = 0; i < DIM; i++)
                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    for (unsigned int i = 0; i < DIM; i++) {
                        if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
                            index[i] = 0;
                            exit = (i == DIM-1);
                        } else {
                            index[i]++;
                            break;
                        }
                    }
                    if (exit) break;
                }
            }
        }
        return r;
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::f(const typename QuarkletFrame::Index& lambda) const
    {

        // f(v) = \int_0^1 g(t)v(t) dt
#if 0
        double r = 1.;
        for (unsigned int i = 0; i < DIM; i++)
            r *= evaluate(*frame_.frames()[i], 0,
                          typename IFRAME::Index(lambda.j()[i],
                                                 lambda.e()[i],
                                                 lambda.k()[i],
                                                 frame_.frames()[i]),
                          0.5);
        return r;

#endif
#if 1
        double r = 0;
        // first compute supp(psi_lambda)
        typename QuarkletFrame::Support supp;
        support(frame_, lambda, supp);
        // setup Gauss points and weights for a composite quadrature formula:
        const int N_Gauss = 5;
        // FixedArray1D<double,DIM> h; // = ldexp(1.0, -supp.j); // granularity for the quadrature
        double hi; // granularity for the quadrature
        FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
        for (unsigned int i = 0; i < DIM; i++) {
            hi = ldexp(1.0,-supp.j[i]);
            gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
            gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
            for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
                for (int n = 0; n < N_Gauss; n++) {
                    gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
                            = hi*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                    gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
                            = hi*GaussWeights[N_Gauss-1][n];
                }
        }
        // compute the point values of the integrand (where we use that it is a tensor product)
        for (unsigned int i = 0; i < DIM; i++)
            evaluate(*frame_.frames()[i], 0,
                     typename IFRAME::Index(lambda.p()[i],
                                            lambda.j()[i],
                                            lambda.e()[i],
                                            lambda.k()[i],
                                            frame_.frames()[i]),
                     gauss_points[i], v_values[i]);
        // iterate over all points and sum up the integral shares
        int index[DIM]; // current multiindex for the point values
        for (unsigned int i = 0; i < DIM; i++)
            index[i] = 0;
        Point<DIM> x;
        while (true) {
            for (unsigned int i = 0; i < DIM; i++)
                x[i] = gauss_points[i][index[i]];
            double share = bvp_->f(x);
            for (unsigned int i = 0; i < DIM; i++)
                share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
            r += share;
            // "++index"
            bool exit = false;
            for (unsigned int i = 0; i < DIM; i++) {
                if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
                    index[i] = 0;
                    exit = (i == DIM-1);
                } else {
                    index[i]++;
                    break;
                }
            }
            if (exit) break;
        }
        return r;
#endif
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::RHS(const double eta,
                                                InfiniteVector<double, typename QuarkletFrame::Index>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(fnorm_sqr - eta*eta);
        typedef typename QuarkletFrame::Index Index;
        typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
        while (it != fcoeffs.end() && coarsenorm < bound)
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
            ++it;
        }
    }
    
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::RHS(const double eta,
                                                InfiniteVector<double, int>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(fnorm_sqr - eta*eta);
        typedef typename QuarkletFrame::Index Index;
        typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
        while (it != fcoeffs.end() && coarsenorm < bound)
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first.number(), it->second);
            ++it;
        }
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::set_bvp(const EllipticBVP<DIM>* bvp)
    {
        bvp_ = bvp;
        compute_rhs();
    }


    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::set_f(const Function<DIM>* fnew)
    {
        bvp_->set_f(fnew);
        compute_rhs();
    }

//    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
//    double
//    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::s_star() const
//    {
//        return 1000.1; // a big number, since s_star == \infty
        /*
        // notation from [St04a]
        const double t = operator_order();
        const int n = DIM;
        const int dT = QuarkletFrame::primal_vanishing_moments();
        const double gamma = QuarkletFrame::primal_regularity();
        return (n == 1
                ? t+dT // [St04a], Th. 2.3 for n=1
                : std::min((t+dT)/(double)n, (gamma-t)/(n-1.))); // [St04a, Th. 2.3]
         */
//    }

    // PERFORMANCE :: use setup_full_collection
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::norm_A() const
    {
        if (normA == 0.0) {
            typedef typename QuarkletFrame::Index Index;
            typedef typename Index::polynomial_type polynomial_type;
            int offsetj, offsetp;
            switch (space_dimension) {
                case 1:
                    offsetj = 2, offsetp = std::min((int)frame_.get_pmax(),2);
                    break;
                case 2:
                    offsetj = 1, offsetp = std::min((int)frame_.get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
            std::set<Index> Lambda;
            
            // cout << "tframe_equation.norm_A :: last quarklet = " << (frame_.last_quarklet(multi_degree(frame_.j0())+offset)) << endl;
            polynomial_type p;
            for (Index lambda = frame_.first_generator(frame_.j0(), p) ;;) {
                Lambda.insert(lambda);
                if (lambda == frame_.last_quarklet(multi_degree(frame_.j0())+offsetj, p)){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = frame_.first_generator(frame_.j0(), p);
                }
                else
                    ++lambda;
                    
            }
//            typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
//            for(int zaehler=0;it!=itend;it++){
//                ++zaehler;
//                cout << zaehler << ", " << *it << endl;
//            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
#if 1
            double help;
            unsigned int iterations;
            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            normAinv = 1./help;
#else
            Vector<double> xk(Lambda.size(), false);
            xk = 1;
            unsigned int iterations;
            normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
        }
        return normA;
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::norm_Ainv() const
    {
        if (normAinv == 0.0) {
            typedef typename QuarkletFrame::Index Index;
            typedef typename Index::polynomial_type polynomial_type;
            std::set<Index> Lambda;
            polynomial_type p;
            //const int j0 = frame().j0();
            //const int jmax = j0+3;
            int offsetj, offsetp;
            switch (space_dimension) {
                case 1:
                    offsetj = 2, offsetp = std::min((int)frame_.get_pmax(),2);
                    break;
                case 2:
                    offsetj = 1, offsetp = std::min((int)frame_.get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
            for (Index lambda = frame_.first_generator() ;;) {
                Lambda.insert(lambda);
                
                if (lambda == frame_.last_quarklet(multi_degree(frame_.j0())+offsetj, p) ){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = frame_.first_generator(frame_.j0(), p);
                }
                else
                    ++lambda;
                    
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
#if 1
            double help;
            unsigned int iterations;
            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            normAinv = 1./help;
#else
            Vector<double> xk(Lambda.size(), false);
            xk = 1;
            unsigned int iterations;
            normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
        }
        return normAinv;
    }
}

namespace WaveletTL
{
    template <unsigned int NUMOFHAARGENERATORS>
    void
    plot_haar_gen_coeffs (FixedVector<double, NUMOFHAARGENERATORS> & plotme, unsigned int level)
    {
        cout << "[ ";
        for (unsigned int i=0; i< NUMOFHAARGENERATORS; ++i)
        {
            cout << twotothejhalf(level)*plotme[i] << " ";
        }
        cout << "]" << endl;
    }

    template <unsigned int NUMOFHAARGENERATORS>
    void
    plot_haar_gen_coeffs (FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> & plotme, unsigned int level)
    {
        for (unsigned int i=0; i< NUMOFHAARGENERATORS; ++i)
        {
            for (unsigned int j=0; j< NUMOFHAARGENERATORS; ++j)
            {
                cout << (1<<level)*plotme.get_entry(j,NUMOFHAARGENERATORS -1 - i) << " ";
            }
            cout << endl;
        }
    }

    /*
     * it is possible to compute the product of two functions given by their
     * coefficient vectors w.r.t. the Haar basis. This is more complicated than
     * multiplying the functions, if they are given by their generator coefficients.
     * It would look like this:
     * return the coefficient vector of the product of the functions corresponding to u*v
     */
    //void haar_product(const Vector<double>& u, const Vector<double>& v, const unsigned int dim, Vector<double>& w)
    //{
        /*
         * u = sum ci*hi, v = sum di*hi. u*v=w=sum xi hi.
         * xk = int hk*u*v = sum ci*dj*int hk*hi*hj
         * 3 cases for the integral:
         * a) highest level occours once:
         * => wavelet on highest level is constant on whole support of the remaining two
         * they are part of an l2 orthogonal basis!
         * b) highest level occours twice:
         * product of the two wavelets on the highest level is constant on support of the third
         * the third is a wavelet. Integrated against a constant function => result = 0
         * c) highest level occours 3 times:
         * 1D: level > 0: on this level only one wavelet with intersecting support exists
         * => we need to consider int xk xk xk which equals 0
         * level = 0. nontrivial integrals consist of 1 gen and 2 wav or 3 gen
         * 2D: the 1D integrals are != 0 iff 3 gens or 1 gen and 2 wavelets are integrated. This results in 4 combinations (ggg, gww, wgw, wwg) for each dimension, resulting in 16 nontrivial integrals for 2D. observe that only 6 of them do not contain a 2D generator.
         * => level = 0 : 16 combinations ( integrals: 8 for h0, 8 for h1) , otherwise 6 types (2 integrals for each h_k)
         */

    //}

    template <class TBASIS, unsigned int DIM>
    double
    haar_gen_coeff(const TBASIS* basis,
                   const InfiniteVector<double, typename TBASIS::Index> & coeffs,
                   const unsigned int* gen_index,
                   const int haarmaxlevel)
    {
        /*
        cout << "haar_gen_coeff:: called with gen_index = ";
        if (DIM == 1)
        {
            cout << "[ " << gen_index[0] << " ]" << endl;
        }
        else
        {
            cout << "[ " << gen_index[0] << ", " << gen_index[1] << " ]" << endl;
        }
         */
        double result = 0;
        //for (InfiniteVector<double, typename TBASIS::Index > ::const_iterator it(coeffs.begin()), itend(coeffs.end()); it != itend; ++it)

        typedef typename TBASIS::Index Index;

        for (typename InfiniteVector<double, Index > ::const_iterator it(coeffs.begin()), itend(coeffs.end()); it != itend; ++it)
        {
            //cout << "index = " << it.index() << " (*it) = " << (*it) << " haarmaxlevel = " << haarmaxlevel << " integrate = " << integrate(basis, it.index(), gen_index, haarmaxlevel+1) << endl;
            result+= (*it)*  integrate(basis, it.index(), gen_index, haarmaxlevel+1);
        }

        return result;
    };


    template <class IBASIS, unsigned int DIM>
    double
    integrate(const TensorBasis<IBASIS, DIM>* basis,
              const typename TensorBasis<IBASIS, DIM>::Index lambda,
              const unsigned int * gen_index,
              const int level)
    {
     // \int_Omega haar_gen(x)psi_lambda(x) dx
        double r = 0;
        // first decide whether the supports of psi_lambda and haar_gen intersect
        typedef typename TensorBasis<IBASIS, DIM>::Support Support;
        Support supp_lambda, supp_haar, supp_inter; // 2^{-j_}<a_,b_> = 2^{-j_1}[a_1,b_1]x...x2^{-j_n}[a_n,b_n]
        basis->support(lambda, supp_lambda);
        // compute support intersection (compare intersecting support in tbasis_support.h

        bool finished = false;
        for (unsigned int i=0;i<DIM;i++)
        {
            supp_inter.j[i] = std::max(supp_lambda.j[i], level);
            supp_haar.a[i]=gen_index[i];
            supp_haar.b[i]=gen_index[i]+1;
            supp_haar.j[i]=level; //propably unused
            if (supp_lambda.j[i] > level) {
                const int adjust = 1<<(supp_lambda.j[i]-level);
                supp_haar.a[i] *= adjust;
                supp_haar.b[i] *= adjust;
                finished = finished || ( (lambda.e()[i] == 1) && (supp_haar.a[i] <= supp_lambda.a[i]) && (supp_lambda.b[i] <= supp_haar.b[i])); // we assume one vanishing moment for the wavelets
            } else {
                const int adjust = 1<<(level-supp_lambda.j[i]);
                supp_lambda.a[i] *= adjust;
                supp_lambda.b[i] *= adjust;
            }
            supp_inter.a[i] = std::max(supp_lambda.a[i],supp_haar.a[i]);
            supp_inter.b[i] = std::min(supp_lambda.b[i],supp_haar.b[i]);
            //if (supp.a[i] >= supp.b[i])
              //  return false;
            finished = finished || (supp_inter.a[i] >= supp_inter.b[i]);
        }

        if (!finished)
        {
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = (IBASIS::primal_polynomial_degree()*IBASIS::primal_polynomial_degree()+1)/2;

            //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
            // FixedArray1D<double,DIM> h; // granularity for the quadrature
            double hi; // granularity for the quadrature
            FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
            for (unsigned int i = 0; i < DIM; i++) {
                hi = ldexp(1.0, -supp_inter.j[i]);
                gauss_points[i].resize(N_Gauss*(supp_inter.b[i]-supp_inter.a[i]));
                gauss_weights[i].resize(N_Gauss*(supp_inter.b[i]-supp_inter.a[i]));
                for (int patch = supp_inter.a[i]; patch < supp_inter.b[i]; patch++)
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(patch-supp_inter.a[i])*N_Gauss+n]
                                = hi*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(patch-supp_inter.a[i])*N_Gauss+n]
                                = hi*GaussWeights[N_Gauss-1][n];
                    }
            }
            // compute point values of the integrand (where we use that it is a tensor product)
            FixedArray1D<Array1D<double>,DIM> psi_lambda_values;     // values of the components of psi_lambda at gauss_points[i]
            for (unsigned int i = 0; i < DIM; i++)
            {
/*
                cout << (*(basis->bases()[i])).get_s0();
                typedef typename IBASIS::Index LocalIndex;
                LocalIndex temp_index(lambda.j()[i],
                                      lambda.e()[i],
                                      lambda.k()[i],
                                      basis->bases()[i]);
*/
                evaluate(*(basis->bases()[i]), 0,
                         typename IBASIS::Index(lambda.j()[i],
                                                lambda.e()[i],
                                                lambda.k()[i],
                                                basis->bases()[i]),
                         gauss_points[i], psi_lambda_values[i]);
            }
            // iterate over all points and sum up the integral shares
            int index[DIM]; // current multiindex for the point values
            for (unsigned int i = 0; i < DIM; i++)
                index[i] = 0;
            Point<DIM> x;
            //const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
            //const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
            double grad_psi_lambda[DIM], grad_psi_mu[DIM], weights;
            //if (bvp_->constant_coefficients())
            double share = 0;
            while (true)
            {
                for (unsigned int i = 0; i < DIM; i++)
                    x[i] = gauss_points[i][index[i]];
                // product of current Gauss weights
                weights = 1.0;
                for (unsigned int i = 0; i < DIM; i++)
                    weights *= gauss_weights[i][index[i]];
                // compute the share q(x)psi_lambda(x)psi_mu(x)
                share = twotothejhalf(level*DIM) * weights; // here enters the value of the Haar generator
                for (unsigned int i = 0; i < DIM; i++)
                    share *= psi_lambda_values[i][index[i]];
                r += share;
                // "++index"
                bool exit = false;
                for (unsigned int i = 0; i < DIM; i++)
                {
                    if (index[i] == N_Gauss*(supp_inter.b[i]-supp_inter.a[i])-1)
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
        }
        return r;
    };

    template <unsigned int DIM>
    double
    integrate(const Function<DIM,double>* fkt,
              const unsigned int primal_polynomial_degree,
              const unsigned int * gen_index,
              const int level)
    {
        // \int_Omega haar_gen(x)fkt(x) dx
        double r = 0;
        // setup Gauss points and weights for a composite quadrature formula:
        //const int N_Gauss = (IBASIS::primal_polynomial_degree()*IBASIS::primal_polynomial_degree()+1)/2;
        const int N_Gauss = (primal_polynomial_degree*primal_polynomial_degree+1)/2;
        //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
        // FixedArray1D<double,DIM> h; // granularity for the quadrature
        double hi; // granularity for the quadrature
        FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
        for (unsigned int i = 0; i < DIM; i++) {
            hi = ldexp(1.0, -level);
            gauss_points[i].resize(N_Gauss);
            gauss_weights[i].resize(N_Gauss);
            for (int n = 0; n < N_Gauss; n++) {
                gauss_points[i][n]
                        = hi*(2*gen_index[i]+1+GaussPoints[N_Gauss-1][n])/2.;
                gauss_weights[i][n]
                        = hi*GaussWeights[N_Gauss-1][n];
            }
        }
// CLEANUP
        //cout << "integrate:: gauss_points[0] = " << gauss_points[0] << endl;
        //cout << "integrate:: gauss_weights[0] = " << gauss_weights[0] << endl;

        // iterate over all points and sum up the integral shares
        int index[DIM]; // current multiindex for the point values
        for (unsigned int i = 0; i < DIM; i++)
            index[i] = 0;
        Point<DIM> x;
        //const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
        //const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
        double weights;
        //if (bvp_->constant_coefficients())
        double share = 0;
        while (true)
        {
            for (unsigned int i = 0; i < DIM; i++)
                x[i] = gauss_points[i][index[i]];
            // product of current Gauss weights
            weights = 1.0;
            for (unsigned int i = 0; i < DIM; i++)
                weights *= gauss_weights[i][index[i]];
// CLEANUP
            //cout << "fkt->value("<<x<<") = " << (fkt->value(x)) << endl;

            r += twotothejhalf(level*DIM) * weights * (fkt->value(x)); // here enters the value of the Haar generator
            // "++index"
            bool exit = false;
            for (unsigned int i = 0; i < DIM; i++)
            {
                if (index[i] == N_Gauss-1)
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

        return r;
    };

    template <class TBASIS, unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const TBASIS* basis,
                     const InfiniteVector<double, typename TBASIS::Index>& coeffs,
    #if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result
    #else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result
    #endif
                    )
    {
        //result.clear();
        //unsigned int* gen_index = new unsigned int[DIM];
        unsigned int gen_index[DIM];
        /*
        for (unsigned int i=0; i< DIM; ++i)
        {
            gen_index[i] = 0;
        }
         */
        precise_evaluate<TBASIS, NUMOFHAARWAVELETS, DIM>(basis, coeffs, result, 0, gen_index);
    }

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
                     unsigned int* gen_index)
    {
        // assume that NUMOFHAARWAVELETS is the number of wavelets in one direction
        if (current_dim < DIM)
        {
            for (unsigned int i=0; i< NUMOFHAARWAVELETS; ++i)
            {
                gen_index[current_dim] = i;
                precise_evaluate<TBASIS, NUMOFHAARWAVELETS, DIM>(basis, coeffs, result, current_dim+1, gen_index);
            }
        }
        else
        {
            int level = log2(NUMOFHAARWAVELETS) -1;
    #if _DIMENSION == 1
            //cout << "precise evaluate:: result[" << gen_index[0] << "] = " << haar_gen_coeff<TBASIS, DIM>(basis,coeffs,gen_index, level) << endl;
            result[gen_index[0]] = haar_gen_coeff<TBASIS, DIM>(basis,coeffs,gen_index, level);
    #else
            result.set_entry(gen_index[0], gen_index[1] ,haar_gen_coeff<TBASIS, DIM> (basis,coeffs,gen_index, level));
    #endif
        }
    };

    template <unsigned int NUMOFHAARWAVELETS, unsigned int DIM>
    void
    precise_evaluate(const Function<DIM,double>* fkt,
                     const unsigned int primal_polynomial_degree,
#if _DIMENSION == 1
                     FixedVector<double, NUMOFHAARWAVELETS> & result
#else
                     FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> & result
#endif
                    )
    {
        unsigned int gen_index[DIM];
        precise_evaluate<NUMOFHAARWAVELETS, DIM>(fkt, primal_polynomial_degree, result, 0, gen_index);
    };

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
                     unsigned int* gen_index)
    {
        // assume that NUMOFHAARWAVELETS is the number of wavelets in one direction
        int level = log2(NUMOFHAARWAVELETS) -1;
        if (current_dim < DIM)
        {
// DEBUG:
            if ((1<< (level+1)) != NUMOFHAARWAVELETS)
            {
                cout << "got ya! parabolic_tools::precise_evaluate Line ~= 367." << endl;
                abort(); // hier und in der anderen precise evaluate aendern...
            }
            for (unsigned int i=0; i< NUMOFHAARWAVELETS; ++i)
            {
                gen_index[current_dim] = i;
                precise_evaluate<NUMOFHAARWAVELETS, DIM>(fkt, primal_polynomial_degree, result, current_dim+1, gen_index);
            }
        }
        else
        {
    #if _DIMENSION == 1
            //cout << "precise evaluate:: result[" << gen_index[0] << "] = " << integrate<DIM>(fkt, primal_polynomial_degree, gen_index, level+1) << endl;
            result[gen_index[0]] = integrate<DIM>(fkt, primal_polynomial_degree, gen_index, level+1);
    #else
            result.set_entry(gen_index[0], gen_index[1] ,integrate<DIM>(fkt, primal_polynomial_degree, gen_index, level+1));
    #endif
        }
    };

    template < unsigned int NUMOFHAARGENERATORS>
    void
    haar_wavelet_transform(const FixedVector<double, NUMOFHAARGENERATORS>& gen_coeffs, FixedVector<double, NUMOFHAARGENERATORS>& wav_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARGENERATORS)-1;
        FixedVector<double, NUMOFHAARGENERATORS> temp_fv;
        wav_coeffs = gen_coeffs;
        for (int j=haar_jmax; j>=0; --j)
        {
            for (unsigned int k=0; k< (1<<j); ++k)
            {
                temp_fv[k]        = (wav_coeffs[2*k] + wav_coeffs[2*k+1])/sqrt2;
                temp_fv[(1<<j)+k] = (wav_coeffs[2*k] - wav_coeffs[2*k+1])/sqrt2;
            }
            for (unsigned int k=0; k< (1<<j); ++k)
            {
                wav_coeffs[k]        = temp_fv[k];
                wav_coeffs[(1<<j)+k] = temp_fv[(1<<j)+k];
            }
        }
    };

    template < unsigned int NUMOFHAARGENERATORS>
    void
    haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>&  gen_coeffs, FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>& wav_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARGENERATORS)-1;
        FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>  temp_coeffs;
        wav_coeffs = gen_coeffs;

        for (int j=haar_jmax; j>=0; --j)
        {
            for (unsigned int col = 0; col < (1<<(j+1)); ++col)
            {
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    //temp_coeffs.set_entry(k,        col, (wav_coeffs.get_entry(2*k, col) + wav_coeffs.get_entry(2*k+1, col))/sqrt2);
                    //temp_coeffs.set_entry((1<<j)+k, col, (wav_coeffs.get_entry(2*k, col) - wav_coeffs.get_entry(2*k+1, col))/sqrt2);
                    temp_coeffs.set_entry(k,        col, (wav_coeffs.get_entry(2*k, col) + wav_coeffs.get_entry(2*k+1, col)) );
                    temp_coeffs.set_entry((1<<j)+k, col, (wav_coeffs.get_entry(2*k, col) - wav_coeffs.get_entry(2*k+1, col)) );
                }
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    //wav_coeffs[k]        = temp_fv[k];
                    //wav_coeffs[(1<<j)+k] = temp_fv[(1<<j)+k];
                    wav_coeffs.set_entry(k,        col, temp_coeffs.get_entry(k,        col));
                    wav_coeffs.set_entry((1<<j)+k, col, temp_coeffs.get_entry((1<<j)+k, col));
                }
            }

            for (unsigned int row = 0; row < (1<<(j+1)); ++row)
            {
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    //temp_fv[k]        = (wav_coeffs[2*k] + wav_coeffs[2*k+1])/sqrt2;
                    //temp_fv[(1<<j)+k] = (wav_coeffs[2*k] - wav_coeffs[2*k+1])/sqrt2;

                    temp_coeffs.set_entry(row, l,        (wav_coeffs.get_entry(row, 2*l) + wav_coeffs.get_entry(row, 2*l+1))/2);
                    temp_coeffs.set_entry(row, (1<<j)+l, (wav_coeffs.get_entry(row, 2*l) - wav_coeffs.get_entry(row, 2*l+1))/2);
                }
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    //wav_coeffs[k]        = temp_fv[k];
                    //wav_coeffs[(1<<j)+k] = temp_fv[(1<<j)+k];
                    wav_coeffs.set_entry(row, l,        temp_coeffs.get_entry(row,        l));
                    wav_coeffs.set_entry(row, (1<<j)+l, temp_coeffs.get_entry(row, (1<<j)+l));
                }
            }
        }
    }

    template < unsigned int NUMOFHAARGENERATORS>
    void
    anisotropic_haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>&  gen_coeffs, FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>& wav_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARGENERATORS)-1;
        FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS>  temp_coeffs;
        wav_coeffs = gen_coeffs;
        // transform columns
        for (unsigned int col = 0; col < NUMOFHAARGENERATORS; ++col)
        {
            for (int j=haar_jmax; j>=0; --j)
            {
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    //temp_fv[k]        = (wav_coeffs[2*k] + wav_coeffs[2*k+1])/sqrt2;
                    //temp_fv[(1<<j)+k] = (wav_coeffs[2*k] - wav_coeffs[2*k+1])/sqrt2;

                    temp_coeffs.set_entry(k,        col, (wav_coeffs.get_entry(2*k, col) + wav_coeffs.get_entry(2*k+1, col))/sqrt2);
                    temp_coeffs.set_entry((1<<j)+k, col, (wav_coeffs.get_entry(2*k, col) - wav_coeffs.get_entry(2*k+1, col))/sqrt2);
                }
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    //wav_coeffs[k]        = temp_fv[k];
                    //wav_coeffs[(1<<j)+k] = temp_fv[(1<<j)+k];
                    wav_coeffs.set_entry(k,        col, temp_coeffs.get_entry(k,        col));
                    wav_coeffs.set_entry((1<<j)+k, col, temp_coeffs.get_entry((1<<j)+k, col));
                }
            }
        }
        // transform rows
        for (unsigned int row = 0; row < NUMOFHAARGENERATORS; ++row)
        {
            for (int j=haar_jmax; j>=0; --j)
            {
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    //temp_fv[k]        = (wav_coeffs[2*k] + wav_coeffs[2*k+1])/sqrt2;
                    //temp_fv[(1<<j)+k] = (wav_coeffs[2*k] - wav_coeffs[2*k+1])/sqrt2;

                    temp_coeffs.set_entry(row, l,        (wav_coeffs.get_entry(row, 2*l) + wav_coeffs.get_entry(row, 2*l+1))/sqrt2);
                    temp_coeffs.set_entry(row, (1<<j)+l, (wav_coeffs.get_entry(row, 2*l) - wav_coeffs.get_entry(row, 2*l+1))/sqrt2);
                }
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    //wav_coeffs[k]        = temp_fv[k];
                    //wav_coeffs[(1<<j)+k] = temp_fv[(1<<j)+k];
                    wav_coeffs.set_entry(row, l,        temp_coeffs.get_entry(row,        l));
                    wav_coeffs.set_entry(row, (1<<j)+l, temp_coeffs.get_entry(row, (1<<j)+l));
                }
            }
        }
    };

    template < unsigned int NUMOFHAARWAVELETS>
    void
    inverse_haar_wavelet_transform(const FixedVector<double, NUMOFHAARWAVELETS>& wav_coeffs, FixedVector<double, NUMOFHAARWAVELETS>& gen_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARWAVELETS)-1;
        FixedVector<double, NUMOFHAARWAVELETS> temp_fv;
        gen_coeffs = wav_coeffs;
        for (int j=0; j<=haar_jmax; ++j)
        {
            for (unsigned int k=0; k< (1<<j); ++k)
            {
                temp_fv[2*k]   = (gen_coeffs[k] + gen_coeffs[(1<<j)+k])/sqrt2;
                temp_fv[2*k+1] = (gen_coeffs[k] - gen_coeffs[(1<<j)+k])/sqrt2;
            }
            for (unsigned int k=0; k< (1<<j); ++k)
            {
                gen_coeffs[2*k]   = temp_fv[2*k];
                gen_coeffs[2*k+1] = temp_fv[2*k+1];
            }
        }
        //coarsening
        for (unsigned int i=0; i< (1<<(haar_jmax+1) ); ++i)
        {
            if (abs(gen_coeffs[i]) <= 1e-14)
            {
                gen_coeffs[i] = 0;
            }
        }
    };

    template< unsigned int NUMOFHAARWAVELETS>
    void
    inverse_haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS>& wav_coeffs, FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS>& gen_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARWAVELETS)-1;
        FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> temp_coeffs;
        gen_coeffs = wav_coeffs;
        for (int j=0; j<=haar_jmax; ++j)
        {
            //cout << "j = " << j << endl << "before doing rows: gen_coeffs = " << endl << gen_coeffs << endl;
            for (unsigned int row = 0; row < (1<<(j+1)); ++row)
            {
                //cout << " row = " << row << " j = " << j << " (1<<j)+(1<<j)-1 = " << ((1<<j)+(1<<j)-1) << endl;
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    temp_coeffs.set_entry(row, 2*l,   (gen_coeffs.get_entry(row, l) + gen_coeffs.get_entry(row, (1<<j)+l))/sqrt2);
                    temp_coeffs.set_entry(row, 2*l+1, (gen_coeffs.get_entry(row, l) - gen_coeffs.get_entry(row, (1<<j)+l))/sqrt2);
                }
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    gen_coeffs.set_entry(row, 2*l,   temp_coeffs.get_entry(row, 2*l));
                    gen_coeffs.set_entry(row, 2*l+1, temp_coeffs.get_entry(row, 2*l+1));
                }
            }
            //cout << "after doing rows: gen_coeffs = " << endl << gen_coeffs << endl;
            for (unsigned int col = 0; col < (1<<(j+1)); ++col)
            {
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    //CLEANUP
                    //cout << "gen_coeffs.get_entry("<<k<<", " <<col<<") +- = "<<gen_coeffs.get_entry(k, col) << "+-"
                    temp_coeffs.set_entry(2*k,   col, (gen_coeffs.get_entry(k, col) + gen_coeffs.get_entry((1<<j)+k, col))/sqrt2);
                    temp_coeffs.set_entry(2*k+1, col, (gen_coeffs.get_entry(k, col) - gen_coeffs.get_entry((1<<j)+k, col))/sqrt2);
                }
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    gen_coeffs.set_entry(2*k,   col, temp_coeffs.get_entry(2*k,   col));
                    gen_coeffs.set_entry(2*k+1, col, temp_coeffs.get_entry(2*k+1, col));
                }
            }
            //cout << "after doing columns: gen_coeffs = " << endl << gen_coeffs << endl;
        }
        //coarsening
        for (unsigned int i=0; i< (1<<(haar_jmax+1) ); ++i)
        {
            for (unsigned int j=0; j< (1<<(haar_jmax+1) ); ++j)
            {
                if (abs(gen_coeffs.get_entry(i,j)) <= 1e-14)
                {
                    gen_coeffs.set_entry(i,j,0);
                }
            }
        }
    };
    template < unsigned int NUMOFHAARWAVELETS>
    void
    anisotropic_inverse_haar_wavelet_transform(const FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS>& wav_coeffs, FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS>& gen_coeffs)
    {
        const double sqrt2(sqrt(2));
        const int haar_jmax = log2(NUMOFHAARWAVELETS)-1;
        FixedMatrix<double, NUMOFHAARWAVELETS, NUMOFHAARWAVELETS> temp_coeffs;
        gen_coeffs = wav_coeffs;
        // transform rows
        for (unsigned int row = 0; row < NUMOFHAARWAVELETS; ++row)
        {
            for (int j=0; j<=haar_jmax; ++j)
            {
                //cout << " row = " << row << " j = " << j << " (1<<j)+(1<<j)-1 = " << ((1<<j)+(1<<j)-1) << endl;
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    temp_coeffs.set_entry(row, 2*l,   (gen_coeffs.get_entry(row, l) + gen_coeffs.get_entry(row, (1<<j)+l))/sqrt2);
                    temp_coeffs.set_entry(row, 2*l+1, (gen_coeffs.get_entry(row, l) - gen_coeffs.get_entry(row, (1<<j)+l))/sqrt2);
                }
                for (unsigned int l=0; l< (1<<j); ++l)
                {
                    gen_coeffs.set_entry(row, 2*l,   temp_coeffs.get_entry(row, 2*l));
                    gen_coeffs.set_entry(row, 2*l+1, temp_coeffs.get_entry(row, 2*l+1));
                }
            }
        }
        // transform columns
        for (unsigned int col = 0; col < NUMOFHAARWAVELETS; ++col)
        {
            for (int j=0; j<=haar_jmax; ++j)
            {
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    temp_coeffs.set_entry(2*k,   col, (gen_coeffs.get_entry(k, col) + wav_coeffs.get_entry((1<<j)+k, col))/sqrt2);
                    temp_coeffs.set_entry(2*k+1, col, (gen_coeffs.get_entry(k, col) - wav_coeffs.get_entry((1<<j)+k, col))/sqrt2);
                }
                for (unsigned int k=0; k< (1<<j); ++k)
                {
                    gen_coeffs.set_entry(2*k,   col, temp_coeffs.get_entry(2*k,   col));
                    gen_coeffs.set_entry(2*k+1, col, temp_coeffs.get_entry(2*k+1, col));
                }
            }
        }
    };




    template<unsigned int SIZE>
    void
    transform_fixedToIV(const FixedVector<double, SIZE> & fixed, InfiniteVector<double, int> & iv)
    {
        iv.clear();
        for (unsigned int i=0; i < SIZE; ++i)
        {
            if (fixed[i] != 0)
            {
                iv.set_coefficient(i, fixed[i]);
            }
        }
    };

    template<unsigned int SIZE>
    void
    transform_IVTofixed(const InfiniteVector<double, int> & iv, FixedVector<double, SIZE> & fixed)
    {
        fixed.scale(0.0);
        for (InfiniteVector<double, int>::const_iterator it(iv.begin()), itend(iv.end()); it!=itend; ++it)
        {
//CLEANUP
            //if (!( (it.index() >= 0) && (it.index() < SIZE) ))
            //{
            //    cout << "transform_IVToFixed:: problem detected. it.index() = " << it.index() << " SIZE = " << SIZE << " iv = " << endl << iv << endl;
            //}
            assert ( (it.index() >= 0) && (it.index() < SIZE) );
            fixed[it.index()]=*it;
        }
    };

    /*
     * Transform Haar Wavelet coefficients from FixedVector/FixedMatrix to InfiniteVector format and vice versa
     * Haar Wavelet coeffs are stored in a Matrix in the following way (for maxlevel = 2)
     * CORRECT:
     *  0  1  4  5 16 17 18 19
     *  2  3  6  7 20 21 22 23
     *  8  9 12 13 24 25 26 27
     * 10 11 14 15 28 29 30 31
     * 32 33 34 35 48 49 50 51
     * 36 37 38 39 52 53 54 55
     * 40 41 42 43 56 57 58 59
     * 44 45 46 47 60 61 62 63
     * WRONG:
     *  0  2  8 10 32 36 40 44
     *  1  3  9 11 33 37 41 45
     *  4  6 12 14 34 38 42 46
     *  5  7 13 15 35 39 43 47
     * 16 20 24 28 48 52 56 60
     * 17 21 25 29 49 53 57 61
     * 18 22 26 30 50 54 58 62
     * 19 23 27 31 51 55 59 63
     */
    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void
    transform_fixedToIV(const FixedMatrix<double, ROWSIZE, COLUMNSIZE> & fixed, InfiniteVector<double, int> & iv)
    {
        assert(ROWSIZE == COLUMNSIZE);
        int maxlevel = log2(ROWSIZE)-1;
        double temp_d;
        int onthislevel;
        iv.clear();
        for (unsigned int level = 0; level <= maxlevel; ++level)
        {
            // each level above level 0 contains 3*2^(2j) many Wavelets
            if (level == 0 )
            {
                temp_d = fixed.get_entry(0,0); if (abs(temp_d) > 1e-15) { iv.set_coefficient(0, temp_d); }
                // transposed
                temp_d = fixed.get_entry(1,0); if (abs(temp_d) > 1e-15) { iv.set_coefficient(2, temp_d); }
                // transposed
                temp_d = fixed.get_entry(0,1); if (abs(temp_d) > 1e-15) { iv.set_coefficient(1, temp_d); }
                temp_d = fixed.get_entry(1,1); if (abs(temp_d) > 1e-15) { iv.set_coefficient(3, temp_d); }
            }
            else
            {
                onthislevel = (1<<(2*level));
                for (unsigned int i=1; i <= 3; ++i) // the 3 types
                {
                    for (unsigned int l=0; l< (1<<level);++l)
                    {
                        for (unsigned int k=0; k< (1<<level);++k)
                        {
                            //assert ((1<<level)*(i%2) == ((i==2)?0:(1<<level)));
                            // WRONG:
                            //temp_d = fixed.get_entry((1<<level)*(i%2)+ k,((i==1)?0:(1<<level))+l);
                            // CORRECT:
                            temp_d = fixed.get_entry(((i==1)?0:(1<<level))+l,(1<<level)*(i%2)+ k);
                            if (abs(temp_d) > 1e-15)
                            {
                                iv.set_coefficient(i*onthislevel+k+l*(1<<level), temp_d);
                            }
                        }
                    }
                }
            }
        }
    };

    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void
    transform_IVTofixed(const InfiniteVector<double, int> & iv, FixedMatrix<double, ROWSIZE, COLUMNSIZE> & fixed)
    {
        assert(ROWSIZE == COLUMNSIZE);
        int maxlevel = log2(ROWSIZE)-1;
        double temp_d;
        //InfiniteVector<double, int>::const_iterator it(iv.begin()), itend(iv.end());
        int current_number; // = it.index();
        int onthislevel;
        int current_level, current_type;
        fixed.scale(0.0);
        for (InfiniteVector<double, int>::const_iterator it(iv.begin()), itend(iv.end()); it != itend; ++it)
        {
            //determine level of it
            current_number = it.index();
            if (current_number < 4)
            {
                // transposed
                //fixed.set_entry(current_number%2, current_number/2, *it);
                fixed.set_entry(current_number/2, current_number%2, *it);
                
            }
            else
            {
                //double temp_d1 = log(current_number+1)/M_LN2;
                //double temp_d2 = log(current_number+1)/M_LN2/2;
                //int temp_i2 = ceil(log(current_number+1)/M_LN2/2);
                //cout << "temp_d1 = " << temp_d << " temp_d2 = " << temp_i1 << " temp_i2 = " << temp_i2 << endl;
                current_level = ceil(log(current_number+1)/M_LN2/2)-1;
                current_type = current_number / (1<<(2*current_level));
// CLEANUP
                //if (!( (current_type >= 1) && (current_type <= 3) ))
                //{
                //    cout << "Problem detected!" << endl;
                //}
                assert ( (current_type >= 1) && (current_type <= 3) );
                current_number = current_number - current_type * (1<<(2*current_level)); // local number in current block
                // transpose:
                //fixed.set_entry((1<<current_level)*(current_type%2)+ current_number % (1<<current_level) ,((current_type==1)?0:(1<<current_level))+current_number / (1<<current_level), *it);
                fixed.set_entry(((current_type==1)?0:(1<<current_level))+current_number / (1<<current_level),(1<<current_level)*(current_type%2)+ current_number % (1<<current_level), *it);
            }
        }
        /* works, if IV is ordered ... :
        for (unsigned int level = 1; level <= maxlevel; ++level)
        {
            // each level above level 0 contains 3*2^(2j) many Wavelets
            onthislevel = (1<<(2*level));
            for (unsigned int i=1; i <= 3; ++i) // the 3 types
            {
                for (unsigned int l=0; l< (1<<level);++l)
                {
                    for (unsigned int k=0; k< (1<<level);++k)
                    {
                        assert (i*onthislevel+k+l*(1<<level) <= it.index());
                        if (i*onthislevel+k+l*(1<<level) == current_number)
                        {
                            fixed.set_entry((1<<level)*(i%2)+ k,((i==1)?0:(1<<level))+l, *it);
                            ++it;
                            current_number = it.index();
                        }
                    }
                }
            }
        }
         * */
    };

    template <class PROBLEM, unsigned int NUMOFTIMESTEPS>
    void solve_parabolic_problem(AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS> &parabolic,
                                  const ROWMethod<InfiniteVector<double, typename PROBLEM::Index> >& method,
                                  const bool time_direction,
                                  const double increment_tolerance,
                                  const double tolerance)
    {
        typedef typename PROBLEM::Index Index;
        typedef typename PROBLEM::WaveletBasis Basis;
        InfiniteVector<double,Index> temp_iv, result, error_estimate;
        unsigned int number_of_timesteps = parabolic.time_discretization_.size()-1;
        parabolic.set_time_direction(time_direction);
        if (time_direction)
        {
// CLEANUP
            //cout << "parabolic_tools::solve_parabolic_equation: time forward case" << endl;

            temp_iv = parabolic.u0;
            //cout << "parabolic_tools::u0 = " << endl << temp_iv << endl;
        }
        else
        {
            //cout << "parabolic_tools::solve_parabolic_equation: time backward case" << endl;
            temp_iv.clear();
            //cout << "parabolic_tools::u0 = " << endl << temp_iv << endl;
        }
        parabolic.set_current_timestep(0);
        for (unsigned int i=1; i<= number_of_timesteps; ++i)
        {
            //cout << "parabolic_tools::solve_parabolic_equation: increment with i = " << i << endl;

            method.increment(&parabolic, (i-1)/(double)number_of_timesteps, temp_iv, 1.0/number_of_timesteps, result, error_estimate, increment_tolerance);
            /*
            cout << "parabolic_tools:: increment arguments:" << endl
                    << "  (i-1)/(double)number_of_timesteps = " << ((i-1)/(double)number_of_timesteps) << endl
                    << "  temp_iv = " << temp_iv << endl
                    << "  1.0/number_of_timesteps = " << (1.0/number_of_timesteps) << endl
                    << "  result = " << result << endl
                    << "  error_estimate = " << error_estimate << endl
                    << "  increment_tolerance = " << increment_tolerance << endl;
             */
            parabolic.set_current_timestep(i);
            parabolic.set_solution(result,tolerance);
            temp_iv = result;
        }
    };

    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT, class IVP>
    void solve_parabolic_problem(AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> & parabolic_problem,
                                  const ROWMethod<InfiniteVector<double, int>, IVP >& method,
                                  const bool time_direction,
                                  const double increment_tolerance,
                                  const double tolerance)
    {
        //typedef typename PROBLEM::Index Index;
        //typedef typename PROBLEM::WaveletBasis Basis;
        InfiniteVector<double,int> temp_iv, result, error_estimate;
        //unsigned int number_of_timesteps = parabolic.time_discretization_.size()-1;
        parabolic_problem.set_time_direction(time_direction);
        if (time_direction)
        {
// CLEANUP
            //cout << "parabolic_tools::solve_parabolic_equation: time forward case" << endl;
            temp_iv = parabolic_problem.u0_;
            //cout << "parabolic_tools::u0 = " << endl << temp_iv << endl;
        }
        else
        {
            //cout << "parabolic_tools::solve_parabolic_equation: time backward case" << endl;
            temp_iv.clear();
            //cout << "parabolic_tools::u0 = " << endl << temp_iv << endl;
        }
        parabolic_problem.set_current_timestep(0);
        for (int i=0; i< NUMOFTIMESTEPS; ++i)
        {
            //cout << "parabolic_tools::solve_parabolic_equation: increment with i = " << i << endl;
            method.increment(&parabolic_problem, (double)(i)/(double)NUMOFTIMESTEPS, temp_iv, 1.0/NUMOFTIMESTEPS, result, error_estimate, increment_tolerance);
            /*
            cout << "parabolic_tools:: increment arguments:" << endl
                    << "  (i-1)/(double)number_of_timesteps = " << ((i-1)/(double)number_of_timesteps) << endl
                    << "  temp_iv = " << temp_iv << endl
                    << "  1.0/number_of_timesteps = " << (1.0/number_of_timesteps) << endl
                    << "  result = " << result << endl
                    << "  error_estimate = " << error_estimate << endl
                    << "  increment_tolerance = " << increment_tolerance << endl;
             */
            parabolic_problem.set_current_timestep(i+1);
            parabolic_problem.set_solution(result,tolerance);
            temp_iv = result;
        }
    };
    
    template <class TBASIS, unsigned int NUMOFHAARGENERATORS, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update)
    {
        unsigned int number_of_timesteps = forward_solution.size()-1;
        w_update.resize(number_of_timesteps+1);
#if _DIMENSION == 1
        FixedVector<double, NUMOFHAARGENERATORS> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedVector<double, NUMOFHAARGENERATORS> > product_over_time;
#else
        FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> > product_over_time;
#endif
        product_over_time.resize(number_of_timesteps+1);
        // at each timestep: multiply generator coefficients => coefficients of the product
        // dont forget the weights at the Haar generators (they do evaluate to 0 or 1/L_2norm = 2^((level*dim)/2)
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            precise_evaluate<TBASIS, NUMOFHAARGENERATORS, _DIMENSION>(basis,forward_solution[i], u_haar_gen_coeffs);
            precise_evaluate<TBASIS, NUMOFHAARGENERATORS, _DIMENSION>(basis,backward_solution[number_of_timesteps - i], h_haar_gen_coeffs);
            //precise_evaluate<TBASIS, NUMOFHAARGENERATORS, _DIMENSION>(basis,backward_solution[i], v_haar_gen_coeffs);
#if _DIMENSION == 1
            for (unsigned int j = 0; j< NUMOFHAARGENERATORS ;++j)
            {
                // alt : product_over_time[i][j] = (u_haar_gen_coeffs[j]-ydata[i][j]) * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFHAARGENERATORS));
                product_over_time[i][j] = u_haar_gen_coeffs[j] * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFHAARGENERATORS));
            }
#else
            // alt: transform_IVTofixed(ydata[i], temp_coeffs); // we do not need to transform the InfVec in the 1D case, since its entries are indexed as in the FixedVector case!

            for (unsigned int j = 0; j < NUMOFHAARGENERATORS ;++j)
            {
                for (unsigned int k = 0; k < NUMOFHAARGENERATORS ;++k)
                {
                    // alt: product_over_time[i].set_entry(j,k, (u_haar_gen_coeffs.get_entry(j,k) - temp_coeffs.get_entry(j,k))* h_haar_gen_coeffs.get_entry(j,k) * NUMOFHAARGENERATORS );
                    product_over_time[i].set_entry(j,k, u_haar_gen_coeffs.get_entry(j,k) * h_haar_gen_coeffs.get_entry(j,k) * NUMOFHAARGENERATORS );
                }
            }
#endif

#if _COMPUTE_UPDATE_W_VERBOSITY > 0
            cout << "parabolic_tools::compute_update_w:: plotting u_haar_gen_coeffs i = " << i << endl;
            plot_haar_gen_coeffs(u_haar_gen_coeffs,log2(NUMOFHAARGENERATORS));
            //cout << "parabolic_tools::compute_update_w:: plotting ydata" << endl;

#if _DIMENSION == 1
            FixedVector<double, NUMOFHAARGENERATORS> temp_storage;
#else
            FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> temp_storage;
#endif
            //transform_IVTofixed(ydata[i], temp_storage);
            //plot_haar_gen_coeffs(temp_storage,log2(NUMOFHAARGENERATORS));
            cout << "parabolic_tools::compute_update_w:: plotting h_haar_gen_coeffs i = " << i << endl;
            plot_haar_gen_coeffs(h_haar_gen_coeffs,log2(NUMOFHAARGENERATORS));
            cout << "parabolic_tools::compute_update_w:: plotting product_over_time[" << i << "]" << endl;
            plot_haar_gen_coeffs(product_over_time[i],log2(NUMOFHAARGENERATORS));

            /*
            cout << "parabolic_tools::compute_update_w::raw coeff vectors" << endl;
            cout << "parabolic_tools::compute_update_w:: plotting ydata" << endl << ydata[i] << endl;
            cout << "parabolic_tools::compute_update_w::forward_solution[" << i << "] = " << endl << forward_solution[i] << endl;
            cout << "parabolic_tools::compute_update_w::precise eval => u_haar_gen_coeffs = " << endl << u_haar_gen_coeffs << endl;
            cout << "parabolic_tools::compute_update_w::backward_solution[" << (number_of_timesteps - i) << "] = " << endl << backward_solution[number_of_timesteps - i ] << endl;
            cout << "parabolic_tools::compute_update_w::precise eval => v_haar_gen_coeffs = " << endl << v_haar_gen_coeffs << endl;
*/
            //cout << "parabolic_tools::compute_update_w:: product_over_time[" << i << "] = " << endl << product_over_time[i] << endl;
#endif

#if _DIMENSION == 1
            FixedVector<double, NUMOFHAARGENERATORS> temp_haar_wav_coeffs;
#else
            FixedMatrix<double, NUMOFHAARGENERATORS, NUMOFHAARGENERATORS> temp_haar_wav_coeffs;
#endif

            haar_wavelet_transform(product_over_time[i], temp_haar_wav_coeffs);
            //cout << "parabolic_tools::compute_update_w:: HWT => temp_haar_wav_coeffs = " << endl << temp_haar_wav_coeffs << endl;
            transform_fixedToIV(temp_haar_wav_coeffs, w_update[i]);
#if _COMPUTE_UPDATE_W_VERBOSITY > 0
            //cout << "parabolic_tools::compute_update_w:: w_update[" << i << "] = " << endl << w_update[i] << endl;

            transform_IVTofixed(w_update[i], temp_haar_wav_coeffs);
            //cout << "IVTofixed => " << endl << data_haar_wav_coeffs << endl;
            inverse_haar_wavelet_transform(temp_haar_wav_coeffs, temp_storage);
            //cout << "inverse_haar_wavelet_transform => " << endl << data_haar_gen_coeffs << endl;
            //cout << "plot data_haar_gen_coeffs" << endl;

            cout << "parabolic_tools::compute_update_w:: product_over_time["<<i<<"] = " << endl << product_over_time[i] << endl;
            cout << "parabolic_tools::compute_update_w:: haar_wavelet_transform => temp_haar_wav_coeffs = " << endl << temp_haar_wav_coeffs << endl;
            cout << "parabolic_tools::compute_update_w:: inverse_haar_wavelet_transform => temp_storage = " << endl << temp_storage << endl;
            cout << "parabolic_tools::compute_update_w:: w_update (with plot_haar_gen_coeffs)" << endl;
            plot_haar_gen_coeffs(temp_storage, log2(temp_storage.size()));
#endif
        }
    };

    

#if 0
    template <class TBASIS, unsigned int NUMOFWAVELETSPERCOORDINATE, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update,
            Array1D<InfiniteVector<double, int> > & u_haar_gen,
            Array1D<InfiniteVector<double, int> > & h_haar_gen)
    {        
        unsigned int number_of_timesteps = forward_solution.size()-1;
        //int level_plusone = (int)log2(NUMOFWAVELETSPERCOORDINATE);
        w_update.resize(number_of_timesteps+1);
#if _DIMENSION == 1
        FixedVector<double, NUMOFWAVELETSPERCOORDINATE> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedVector<double, NUMOFWAVELETSPERCOORDINATE> > product_over_time;
#else
        FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> > product_over_time;
#endif
        product_over_time.resize(number_of_timesteps+1);
        // at each timestep: multiply generator coefficients => coefficients of the product
        // dont forget the weights at the Haar generators (they do evaluate to 0 or 1/L_2norm = 2^((level*dim)/2)
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,forward_solution[i], u_haar_gen_coeffs);
            transform_fixedToIV(u_haar_gen_coeffs,u_haar_gen[i]);
            precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,backward_solution[number_of_timesteps - i], h_haar_gen_coeffs);
            //transform_fixedToIV(h_haar_gen_coeffs, h_haar_gen[number_of_timesteps - i]);
            transform_fixedToIV(h_haar_gen_coeffs, h_haar_gen[i]);
            //precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,backward_solution[i], h_haar_gen_coeffs);
            //transform_fixedToIV(h_haar_gen_coeffs, h_haar_gen[i]);
            
            //precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,backward_solution[i], v_haar_gen_coeffs);
#if _DIMENSION == 1
            for (unsigned int j = 0; j< NUMOFWAVELETSPERCOORDINATE ;++j)
            {
                // alt : product_over_time[i][j] = (u_haar_gen_coeffs[j]-ydata[i][j]) * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
                product_over_time[i][j] = u_haar_gen_coeffs[j] * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
                //product_over_time[i][j] = u_haar_gen[i][j] * h_haar_gen_coeffs[number_of_timesteps-i][j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
            }
#else
            // alt: transform_IVTofixed(ydata[i], temp_coeffs); // we do not need to transform the InfVec in the 1D case, since its entries are indexed as in the FixedVector case!

            for (unsigned int j = 0; j < NUMOFWAVELETSPERCOORDINATE ;++j)
            {
                for (unsigned int k = 0; k < NUMOFWAVELETSPERCOORDINATE ;++k)
                {
                    // alt: product_over_time[i].set_entry(j,k, (u_haar_gen_coeffs.get_entry(j,k) - temp_coeffs.get_entry(j,k))* h_haar_gen_coeffs.get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE );
                    product_over_time[i].set_entry(j,k, u_haar_gen_coeffs.get_entry(j,k) * h_haar_gen_coeffs.get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE ); // faktor = 2^( (hjmax+1)*DIM/2 )
                    //product_over_time[i].set_entry(j,k, u_haar_gen[i].get_entry(j,k) * h_haar_gen[i].get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE );
                }
            }
#endif

#if _COMPUTE_UPDATE_W_VERBOSITY > 0

            
#if _DIMENSION == 1
            //FixedVector<double, NUMOFWAVELETSPERCOORDINATE> temp_storage;
#else
            //FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> temp_storage;
#endif
            cout << "parabolic_tools::compute_update_w:: plotting u_haar_gen_coeffs i = " << i << endl;
            plot_haar_gen_coeffs(u_haar_gen_coeffs,log2(NUMOFWAVELETSPERCOORDINATE));
            //
            cout << "parabolic_tools::compute_update_w:: plotting u_haar_gen i = " << i << endl;
            cout << u_haar_gen[i] << endl;
            //plot_haar_gen_coeffs(u_haar_gen[i],log2(NUMOFWAVELETSPERCOORDINATE));
            
            cout << "parabolic_tools::compute_update_w:: plotting h_haar_gen_coeffs i = " << i << endl;
            plot_haar_gen_coeffs(h_haar_gen_coeffs,log2(NUMOFWAVELETSPERCOORDINATE));
            //
            cout << "parabolic_tools::compute_update_w:: plotting h_haar_gen i = " << i << endl;
            cout << h_haar_gen << endl;
            //plot_haar_gen_coeffs(h_haar_gen[i],log2(NUMOFWAVELETSPERCOORDINATE));
            
            cout << "parabolic_tools::compute_update_w:: plotting product_over_time[" << i << "]" << endl;
            plot_haar_gen_coeffs(product_over_time[i],log2(NUMOFWAVELETSPERCOORDINATE));

            // /*
            cout << "parabolic_tools::compute_update_w::raw coeff vectors" << endl;
            //cout << "parabolic_tools::compute_update_w:: plotting ydata" << endl << ydata[i] << endl;
            cout << "parabolic_tools::compute_update_w::forward_solution[" << i << "] = " << endl << forward_solution[i] << endl;
            cout << "parabolic_tools::compute_update_w::precise eval => u_haar_gen_coeffs = " << endl << u_haar_gen_coeffs << endl;
            cout << "parabolic_tools::compute_update_w::backward_solution[" << (number_of_timesteps - i) << "] = " << endl << backward_solution[number_of_timesteps - i ] << endl;
            //cout << "parabolic_tools::compute_update_w::precise eval => v_haar_gen_coeffs = " << endl << v_haar_gen_coeffs << endl;

            // */
            //
            cout << "parabolic_tools::compute_update_w:: product_over_time[" << i << "] = " << endl << product_over_time[i] << endl;
#endif

#if _DIMENSION == 1
            FixedVector<double, NUMOFWAVELETSPERCOORDINATE> temp_haar_wav_coeffs;
#else
            FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> temp_haar_wav_coeffs;
#endif

            haar_wavelet_transform(product_over_time[i], temp_haar_wav_coeffs);

            transform_fixedToIV(temp_haar_wav_coeffs, w_update[i]);
#if _COMPUTE_UPDATE_W_VERBOSITY > 1
            // check Haar wavelet trafo:
            //cout << "parabolic_tools::compute_update_w:: HWT => temp_haar_wav_coeffs = " << endl << temp_haar_wav_coeffs << endl;
            //cout << "parabolic_tools::compute_update_w:: w_update[" << i << "] = " << endl << w_update[i] << endl;
            //transform_IVTofixed(w_update[i], temp_haar_wav_coeffs);
            //cout << "IVTofixed => " << endl << temp_haar_wav_coeffs << endl;
            //inverse_haar_wavelet_transform(temp_haar_wav_coeffs, temp_storage);
            //cout << "inverse_haar_wavelet_transform => " << endl;
            //plot_haar_gen_coeffs(temp_storage,log2(NUMOFHAARGENERATORS)); // should be the same as product_over_time
            //cout << "plot data_haar_gen_coeffs" << endl;
            //cout << "parabolic_tools::compute_update_w:: product_over_time["<<i<<"] = " << endl << product_over_time[i] << endl;
            //cout << "parabolic_tools::compute_update_w:: haar_wavelet_transform => temp_haar_wav_coeffs = " << endl << temp_haar_wav_coeffs << endl;
            //cout << "parabolic_tools::compute_update_w:: inverse_haar_wavelet_transform => temp_storage = " << endl << temp_storage << endl;
            //cout << "parabolic_tools::compute_update_w:: w_update" << endl;
            //plot_haar_gen_coeffs(temp_storage, log2(NUMOFWAVELETSPERCOORDINATE)); // should be the same as product_over_time
            //cout << "parabolic_tools::compute_update_w:: product_over_time["<<i<<"] = " << endl;
            //plot_haar_gen_coeffs(product_over_time[i], log2(NUMOFWAVELETSPERCOORDINATE));
#endif
        }
    };
#endif    
    
    template <class TBASIS, unsigned int NUMOFWAVELETSPERCOORDINATE, unsigned int DIM>
    void compute_update_w(const TBASIS* basis,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            Array1D<InfiniteVector<double, int> > & w_update,
            std::ofstream& plotstream)
    {        
        unsigned int number_of_timesteps = forward_solution.size()-1;
        
        Array1D<InfiniteVector<double, int> > u_haar_gen, h_haar_gen;
        u_haar_gen.resize(number_of_timesteps+1);
        h_haar_gen.resize(number_of_timesteps+1);
        
        //int level_plusone = (int)log2(NUMOFWAVELETSPERCOORDINATE);
        w_update.resize(number_of_timesteps+1);
#if _DIMENSION == 1
        FixedVector<double, NUMOFWAVELETSPERCOORDINATE> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedVector<double, NUMOFWAVELETSPERCOORDINATE> > product_over_time;
#else
        FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> u_haar_gen_coeffs, h_haar_gen_coeffs;
        Array1D<FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> > product_over_time;
#endif
        product_over_time.resize(number_of_timesteps+1);
        // at each timestep: multiply generator coefficients => coefficients of the product
        // dont forget the weights at the Haar generators (they do evaluate to 0 or 1/L_2norm = 2^((level*dim)/2)
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,forward_solution[i], u_haar_gen_coeffs);
            //transform_fixedToIV(u_haar_gen_coeffs,u_haar_gen[i]);
            precise_evaluate<TBASIS, NUMOFWAVELETSPERCOORDINATE, _DIMENSION>(basis,backward_solution[number_of_timesteps - i], h_haar_gen_coeffs);
            //transform_fixedToIV(h_haar_gen_coeffs, h_haar_gen[i]);
#if _DIMENSION == 1
            for (unsigned int j = 0; j< NUMOFWAVELETSPERCOORDINATE ;++j)
            {
                // alt : product_over_time[i][j] = (u_haar_gen_coeffs[j]-ydata[i][j]) * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
                product_over_time[i][j] = u_haar_gen_coeffs[j] * h_haar_gen_coeffs[j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
                //product_over_time[i][j] = u_haar_gen[i][j] * h_haar_gen_coeffs[number_of_timesteps-i][j] * twotothejhalf(log2(NUMOFWAVELETSPERCOORDINATE));
            }
#else
            // alt: transform_IVTofixed(ydata[i], temp_coeffs); // we do not need to transform the InfVec in the 1D case, since its entries are indexed as in the FixedVector case!
            for (unsigned int j = 0; j < NUMOFWAVELETSPERCOORDINATE ;++j)
            {
                for (unsigned int k = 0; k < NUMOFWAVELETSPERCOORDINATE ;++k)
                {
                    // alt: product_over_time[i].set_entry(j,k, (u_haar_gen_coeffs.get_entry(j,k) - temp_coeffs.get_entry(j,k))* h_haar_gen_coeffs.get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE );
                    product_over_time[i].set_entry(j,k, u_haar_gen_coeffs.get_entry(j,k) * h_haar_gen_coeffs.get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE ); // faktor = 2^( (hjmax+1)*DIM/2 )
                    //product_over_time[i].set_entry(j,k, u_haar_gen[i].get_entry(j,k) * h_haar_gen[i].get_entry(j,k) * NUMOFWAVELETSPERCOORDINATE );
                }
            }
#endif

            // transform haar generator coefficients (fixed vector) into haar wavelet coefficients (infinite vector)
#if _DIMENSION == 1
            FixedVector<double, NUMOFWAVELETSPERCOORDINATE> temp_haar_wav_coeffs;
#else
            FixedMatrix<double, NUMOFWAVELETSPERCOORDINATE, NUMOFWAVELETSPERCOORDINATE> temp_haar_wav_coeffs;
#endif
            haar_wavelet_transform(product_over_time[i], temp_haar_wav_coeffs);
            transform_fixedToIV(temp_haar_wav_coeffs, w_update[i]);
            
            // plot u,h,w_update
            plotstream << "u_haar_gen_coeffs{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE)*u_haar_gen_coeffs[l];
            }
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*u_haar_gen_coeffs.get_entry(l,k); // MATLAB: y goes from top to bottom 
                }
                plotstream << "; ";
            }
#endif     
            plotstream << "];\n";
            plotstream << "h_haar_gen_coeffs{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                //plotstream << " " << twotothejhalf(haar_jmax+1)*h_haar_gen_coeffs[l];
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE)*h_haar_gen_coeffs[l];
            }
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*h_haar_gen_coeffs.get_entry(l,k); // MATLAB: y goes from top to bottom
                }
                plotstream << "; ";
            }
#endif  
            plotstream << "];\n";
            plotstream << "w_update_plot{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE) * product_over_time[i][l];
            }
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*product_over_time[i].get_entry(l,k); // MATLAB: y goes from top to bottom
                }
                plotstream << "; ";
            }
#endif
            plotstream << "];\n";
        }
    };

/*
    template <class C, class I>
    void infinite_vector_writeToFile(const char *filename)
    {
        std::ofstream bin_file(filename);
        bin_file.is_open();
        if (bin_file.is_open())
        {
            try
            {
                int dim, temp_i;
                double temp_d;
                for (typename InfiniteVector<C,I>::const_iterator it((*this).begin()); it != (*this).end(); ++it)
                {
                    //cout << "writing: it = " << it.index() << " value = " << (*it) << endl;
#if 0
                    int j=it.index().j()[0];
                    int e=it.index().e()[0];
                    int k=it.index().k()[0];
                    double d = (*it);

                    bin_file.write((char*)(&j), sizeof(int));
                    bin_file.write((char*)(&e), sizeof(int));
                    bin_file.write((char*)(&k), sizeof(int));
                    bin_file.write(reinterpret_cast<char*>(&d), sizeof(double));
#else
                    dim = it.index().j().size();
                    for (unsigned int i=0;i<dim;++i)
                    {
                        temp_i = it.index().j()[i];
                        //cout << "j["<<i<<"]="<<temp_i<<endl;
                        bin_file.write((char*)(&temp_i), sizeof(int));
                    }
                    for (unsigned int i=0;i<dim;++i)
                    {
                        temp_i = it.index().e()[i];
                        //cout << "e["<<i<<"]="<<temp_i<<endl;
                        bin_file.write((char*)(&temp_i), sizeof(int));
                    }
                    for (unsigned int i=0;i<dim;++i)
                    {
                        temp_i = it.index().k()[i];
                        //cout << "k["<<i<<"]="<<temp_i<<endl;
                        bin_file.write((char*)(&temp_i), sizeof(int));
                    }
                    temp_d = (*it);
                    //cout << "temp_d="<<temp_d<<endl;
                    bin_file.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
#endif
                }
            }
            catch (...)
            {
                cout << "InfiniteVector::writeToFile: Error while writing to " << filename << endl;
            }
            bin_file.close();
        }
        else
        {
            cout << "InfiniteVector::writeToFile: Could not write to " << filename << endl;
        }
    };
*/

/*
    template <class C, class I>
    void infinite_vector_readFromFile(const char *filename)
    {
        std::ifstream bin_file(filename);
        if (bin_file.is_open())
        {
            (*this).clear(); // going to fill this InfiniteVector with data from file

            //while (!bin_file.eof()) // doesn't work, see comment after .eof()
            while(true)
            {
                MultiIndex<int,_DIMENSION> j,e,k;
                int temp_i;
                double temp_d;
#if 1 // WAVELETTL_USE_TBASIS == 1
                try
                {
                    //unsigned int i = 0;
                    for (unsigned int i=0;i<_DIMENSION;++i)
                    {
                        bin_file.read((char*)(&temp_i), sizeof(int));
                        j[i]=temp_i;
                        //cout << "j["<<i<<"]="<<temp_i << " eof = " << bin_file.eof() <<endl;
                    }
                    if (bin_file.eof()) // the last meaningful .read operation sets the last bit to the last bit of the last variable read. that means that we are NOT at the end of file. So the whole try-block is executed again, resulting in one wrong entry in *this.
                    {
                        break;
                    }
                    for (unsigned int i=0;i<_DIMENSION;++i)
                    {
                        bin_file.read((char*)(&temp_i), sizeof(int));
                        e[i]=temp_i;
                        //cout << "e["<<i<<"]="<<temp_i<< " eof = " << bin_file.eof() <<endl;
                    }
                    for (unsigned int i=0;i<_DIMENSION;++i)
                    {
                        bin_file.read((char*)(&temp_i), sizeof(int));
                        k[i]=temp_i;
                        //cout << "k["<<i<<"]="<<temp_i<< " eof = " << bin_file.eof() <<endl;
                    }
                    bin_file.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    //cout << "d="<<temp_d<<endl;
                    I ind(j,e,k,0); // NOTE: Index should be given a BASIS, here it is 0 instead
                    (*this).set_coefficient(ind, temp_d);
                }
#else
                try
                {
                    int j,e,k;
                    double d;
                    bin_file.read((char*)(&j), sizeof(int));
                    bin_file.read((char*)(&e), sizeof(int));
                    bin_file.read((char*)(&k), sizeof(int));
                    bin_file.read(reinterpret_cast<char*>(&d), sizeof(double));

                    I ind(j, e, k, 0);  // NOTE: Index should be given a BASIS, here it is 0 instead
                    cout << "loaded: ind = " << ind << endl;
                    (*this).set_coefficient(ind, d);
                }
#endif
                catch (...)
                {
                    cout << "InfiniteVector::readFromFile: Read error in " << filename << endl;
                }
            }
            bin_file.close();
        }
        else
        {
            cout << "InfiniteVector::readFromFile: Could not read from " << filename << endl;
        }
    };
*/

    void shrinkage_iteration(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double alpha,
            Array1D<InfiniteVector<double, int> >& w_next)
    {
        assert (w.size() == w_update.size());
        int number_of_timesteps = w.size() -1;
                
// CLEANUP
        /*
        cout << "shrinkage_iteration:: w (post shrinkage)" << endl;
        for (unsigned int i=0; i<= number_of_timesteps;++i)
        {
            cout << "shrinkage_iteration:: w[" << i << "] = " << endl << w[i] << endl;
        }
        cout << "shrinkage_iteration:: w_update" << endl;
        for (unsigned int i=0; i<= number_of_timesteps;++i)
        {
            cout << "shrinkage_iteration:: w_update[" << i << "] = " << endl << w_update[i] << endl;
        }
        */
        w_next.resize(number_of_timesteps+1);
        for (unsigned int i=0; i <= number_of_timesteps; ++i)
        {
            //assert(w_next[i].empty());
            //w_next[i].clear();
            // this can be done in a more clever way: have a look at InfiniteVector.add
            w_next[i] = w[i];
            w_next[i].add(w_update[i]);
            for (InfiniteVector<double, int>::const_iterator it(w_next[i].begin()), itend(w_next[i].end()); it != itend; ++it)
            {
                if (abs(*it) <= alpha)
                {
                    // w_next[i].erase(it.index()); // erase is protected ...
                    w_next[i].add_coefficient(it.index(), - *it); // this leads to a call of erase
                    assert (w_next[i].get_coefficient(it.index()) == 0);
                }
                else
                {
                    (w_next[i]).add_coefficient(it.index(), ((*it > 0)?-alpha:alpha) );
                }
            }
            //cout << "parabolic_tools:shrinkage_iteration: wnext[" << i << "] = " << wnext [i] << endl;
        }
// CLEANUP
        //cout << "shrinkage_iteration:: w_next (post shrinkage)" << endl;
        //for (unsigned int i=0; i<= number_of_timesteps;++i)
        //{
        //    cout << "shrinkage_iteration:: w_next[" << i << "] = " << endl << w_next[i] << endl;
        //}
    };

    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void shrinkage_iteration(AffLinParEq_qtbasis< NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT> & parabolic_problem,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & w,
            const double alpha,
            FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPS +1> & w_next,
            std::ofstream * plotstream)
    {
        // FixedArray1D<InfiniteVector<double,int>, NUMOFTIMESTEPS+1 > forward_solution_, backward_solution_
        Array1D < FixedMatrix<double, PRECISE_EVALUATE_GRANULARITY, PRECISE_EVALUATE_GRANULARITY> > u_haar_gen_coeffs, h_haar_gen_coeffs;
        //FixedArray1D< Array1D< FixedMatrix<double, ONEDIMHAARCOUNT> > > product_over_time;
        FixedMatrix<double, ONEDIMHAARCOUNT> product_time_patch;
        // at each timestep: multiply generator coefficients => coefficients of the product
        // dont forget the weights at the Haar generators (they do evaluate to 0 or 1/L_2norm = 2^((level*dim)/2)
        unsigned int factor = PRECISE_EVALUATE_GRANULARITY / ONEDIMHAARCOUNT;
        double entry = 0;
        FixedMatrix<double, ONEDIMHAARCOUNT> temp_haar_wav_coeffs;
        for (int i=0; i <= NUMOFTIMESTEPS; ++i )
        {            
            parabolic_problem.cached_sampled_mapping(parabolic_problem.forward_solution_[i], u_haar_gen_coeffs);
            parabolic_problem.cached_sampled_mapping(parabolic_problem.backward_solution_[NUMOFTIMESTEPS - i], h_haar_gen_coeffs);
            for (unsigned int p=0; p<u_haar_gen_coeffs.size(); ++p)
            {
                for (unsigned int j = 0; j < ONEDIMHAARCOUNT ;++j)
                {
                    for (unsigned int k = 0; k < ONEDIMHAARCOUNT ;++k)
                    {
                        entry = 0;
                        for (unsigned int l=0; l< factor; ++l)
                        {
                            for (unsigned int m=0; m< factor; ++m)
                            {
                                entry += u_haar_gen_coeffs[p].get_entry(j*factor+l,k*factor+m) * h_haar_gen_coeffs[p].get_entry(j*factor+l,k*factor+m);
                            }
                        }
                        entry /= (factor*factor);
                        product_time_patch.set_entry(j,k, entry);
                    }
                }
                if (ONEDIMHAARCOUNT > 1)
                {
                    product_time_patch.scale(1.0/ONEDIMHAARCOUNT); // gen-coeffs!
                }
                haar_wavelet_transform(product_time_patch, temp_haar_wav_coeffs);
                for (unsigned int j=0; j>ONEDIMHAARCOUNT; ++j)
                {
                    for (unsigned int k=0; k>ONEDIMHAARCOUNT; ++k)
                    {
                        entry = temp_haar_wav_coeffs(j,k) + w[i][p](j,k);
                        if (abs(entry) > alpha)
                        {
                            w_next[i][p](j,k) = (entry > 0)?(entry-alpha):(entry+alpha);
                        }
                    }
                }
            }
        }
    }
    
    
    bool relative_change_criterion_ell2(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double stopping_tolerance,
            double relative_error)
    {
        assert(w.size() == w_update.size());
        double norm_w_sqr(0), norm_wu_sqr(0);
        for (unsigned int i = 0; i < w.size(); ++i)
        {
            norm_wu_sqr += l2_norm_sqr(w_update[i]);
        }
        if (norm_wu_sqr == 0)
        {
            relative_error = 0;
            return true;
        }
        for (unsigned int i = 0; i < w.size(); ++i)
        {
            norm_w_sqr += l2_norm_sqr(w[i]);
        }
        if (norm_w_sqr == 0)
        {
            relative_error = 0;
            return false;
        }
        relative_error = sqrt(norm_wu_sqr/norm_w_sqr);
// cleanup
        cout << "parabolic_tools::relative_change_criterion:: norm_wu_sqr = " << norm_wu_sqr << " norm_w_sqr = " << norm_w_sqr << " relative_error = " << relative_error << " stopping_tolerance = " << stopping_tolerance << endl;
        return (relative_error < stopping_tolerance);
    };

    bool relative_change_criterion_ell1(const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_update,
            const double stopping_tolerance,
            double relative_error)
    {
        assert(w.size() == w_update.size());
        double norm_w(0), norm_wu(0);
        relative_error = 0;
        for (unsigned int i = 0; i < w.size(); ++i)
        {
            norm_wu += l1_norm(w_update[i]);
        }
        if (norm_wu == 0)
        {
            return true;
        }
        for (unsigned int i = 0; i < w.size(); ++i)
        {
            norm_w += l1_norm(w[i]);
        }
        if (norm_w == 0)
        {
            return false;
        }
        relative_error = norm_wu/norm_w;
// cleanup
        cout << "parabolic_tools::relative_change_criterion:: ||w_update||_ell1 = " << (norm_wu/w.size()) << " ||w||_ell1 = " << (norm_w/w.size()) << " relative_error = " << relative_error << " stopping_tolerance = " << stopping_tolerance << endl;
        return (relative_error < stopping_tolerance);
    };

    template <class TBASIS>
    void create_noise(const double sigma,
                      const double delta,
                      const TBASIS* basis,
                      Array1D<InfiniteVector<double, typename TBASIS::Index> > &noise_coeffs)
    {
        // get a random number generator

// CLEANUP
        //cout << "create_noise:: sigma = " << sigma << endl;
        //cout << "create_noise:: delta = " << delta << endl;
        //cout << "create_noise:: noise_coeffs.size() = " << (noise_coeffs.size()) << endl;

        MTRand mtrand1;
        typedef typename TBASIS::Index Index;
        Index start (0,basis), lauf (0,basis);
        Array1D<double> noise_norms_sqr;
        noise_norms_sqr.resize(noise_coeffs.size());
// CLEANUP

        for (unsigned int i=0; i< noise_coeffs.size(); ++i)
        {
            noise_norms_sqr[i] = 0;
            //cout << "create_noise:: noise_norms_sqr[" << i<< "] = " << noise_norms_sqr[i] << endl;
        }

        double temp_double;

        for (unsigned int j=0; j< basis->degrees_of_freedom(); ++j)
        {
            for (unsigned int i=0; i< noise_coeffs.size(); ++i)
            {
                temp_double = mtrand1.randNorm(0.0,sigma);
                noise_coeffs[i].set_coefficient(lauf,temp_double);
                noise_norms_sqr[i]+= temp_double*temp_double;
// CLEANUP
                //cout << "create_noise:: lauf = " << lauf << ", temp_double = " << temp_double << ", noise_norms_sqr[" << i << "] = " << noise_norms_sqr[i] << endl;
            }
            ++lauf;
        }
        // renormalize everything
        temp_double = noise_norms_sqr[0];
// CLEANUP
        //cout << "create_noise:: noise_norms_sqr[0] = " << noise_norms_sqr[0] << endl;
        for (unsigned int i=1; i < noise_coeffs.size();++i)
        {
            temp_double+=noise_norms_sqr[i];
            //cout << "noise_norms_sqr[" << i<< "] = " << noise_norms_sqr[i] << endl;
        }
        //cout << "create_noise:: temp_double = " << temp_double << endl;
        //cout << "create_noise:: delta/sqrt(temp_double) = " << (delta/sqrt(temp_double)) << endl;
        for (unsigned int i=0; i< noise_coeffs.size();++i)
        {
            noise_coeffs[i].scale(delta/sqrt(temp_double));
        }
    };

   template <class QTBASIS, unsigned int NUMOFTIMESTEPSPLUSONE>
    void create_noise(const double sigma,
                      const double delta,
                      const QTBASIS* basis,
                      FixedArray1D<InfiniteVector<double, int>, NUMOFTIMESTEPSPLUSONE > &noise_coeffs)
    {
        // get a random number generator

// CLEANUP
        //cout << "create_noise:: sigma = " << sigma << endl;
        //cout << "create_noise:: delta = " << delta << endl;
        //cout << "create_noise:: noise_coeffs.size() = " << (noise_coeffs.size()) << endl;

        MTRand mtrand1;
        //typedef typename TBASIS::Index Index;
        //Index start (0,basis), lauf (0,basis);
        FixedArray1D<double, NUMOFTIMESTEPSPLUSONE> noise_norms_sqr;
        //noise_norms_sqr.resize(noise_coeffs.size());
// CLEANUP

        for (int i=0; i< NUMOFTIMESTEPSPLUSONE; ++i)
        {
            noise_norms_sqr[i] = 0;
            //cout << "create_noise:: noise_norms_sqr[" << i<< "] = " << noise_norms_sqr[i] << endl;
        }

        double temp_double;

        for (unsigned int j=0; j< basis->degrees_of_freedom(); ++j)
        {
            for (int i=0; i< NUMOFTIMESTEPSPLUSONE; ++i)
            {
                temp_double = mtrand1.randNorm(0.0,sigma);
                noise_coeffs[i].set_coefficient(j,temp_double);
                noise_norms_sqr[i]+= temp_double*temp_double;
// CLEANUP
                //cout << "create_noise:: lauf = " << lauf << ", temp_double = " << temp_double << ", noise_norms_sqr[" << i << "] = " << noise_norms_sqr[i] << endl;
            }
            //++lauf;
        }
        // renormalize everything
        temp_double = noise_norms_sqr[0];
// CLEANUP
        //cout << "create_noise:: noise_norms_sqr[0] = " << noise_norms_sqr[0] << endl;
        for (int i=1; i < NUMOFTIMESTEPSPLUSONE;++i)
        {
            temp_double+=noise_norms_sqr[i];
            //cout << "noise_norms_sqr[" << i<< "] = " << noise_norms_sqr[i] << endl;
        }
        //cout << "create_noise:: temp_double = " << temp_double << endl;
        //cout << "create_noise:: delta/sqrt(temp_double) = " << (delta/sqrt(temp_double)) << endl;
        for (unsigned int i=0; i< noise_coeffs.size();++i)
        {
            noise_coeffs[i].scale(delta/sqrt(temp_double));
        }
    };
    
    void initialize_coupling_matrix(Array1D<InfiniteVector<double, int> >  & wtrue,
                                    const int first_setup_coupling_matrix_w,
                                    const int d)
    {
        if (first_setup_coupling_matrix_w <= 12)
        {
            switch (first_setup_coupling_matrix_w)
            {
                case 0:
                    for (unsigned int i = 0; i < wtrue.size(); ++i)
                    {
                        wtrue[i].clear(); // an empty w vector should result in the heat equation (for d=1)
                    }
                    break;
                case 1:
                    for (unsigned int i = 0; i < wtrue.size(); ++i)
                    {
                        wtrue[i].set_coefficient(0,1); // this adds the Gram matrix, resulting in A= -delta+Id
                    }
                    //parabolic.assemble_W(w);
                    break;
                case 2: // toy example: 1 active coefficient localized in time and space
    #if _DIMENSION == 1
                    wtrue[5].set_coefficient(6,1);
    #else
                    wtrue[5].set_coefficient(58,1);
    #endif
                    break;
                case 3:
    #if _DIMENSION == 1
                    wtrue[0].set_coefficient(6,1);
                    wtrue[1].set_coefficient(6,1);
                    wtrue[2].set_coefficient(6,1);
                    wtrue[3].set_coefficient(6,1);
                    wtrue[4].set_coefficient(6,1);
                    wtrue[5].set_coefficient(6,1);
                    wtrue[6].set_coefficient(6,1);
                    wtrue[7].set_coefficient(6,1);    
    #else
                    wtrue[0].set_coefficient(58,1);
                    wtrue[1].set_coefficient(58,1);
                    wtrue[2].set_coefficient(58,1);
                    wtrue[3].set_coefficient(58,1);
                    wtrue[4].set_coefficient(58,1);
                    wtrue[5].set_coefficient(58,1);
                    wtrue[6].set_coefficient(58,1);
                    wtrue[7].set_coefficient(58,1);
    #endif
                    break;
                case 4: // only one active coef to reconstruct.
    #if _DIMENSION == 1
                    wtrue[3].set_coefficient(6,1);
                    wtrue[4].set_coefficient(6,1);
                    wtrue[5].set_coefficient(6,1);
                    wtrue[6].set_coefficient(6,1);
                    wtrue[7].set_coefficient(6,1);
    #else
                    wtrue[3].set_coefficient(58,1);
                    wtrue[4].set_coefficient(58,1);
                    wtrue[5].set_coefficient(58,1);
                    wtrue[6].set_coefficient(58,1);
                    wtrue[7].set_coefficient(58,1);
    #endif
                    break;
                case 5:                
    #if _DIMENSION == 1
                    // would be one genetrator to reconstruct = 3 coeffs in 1D ... if not for the wrong sign!
                    wtrue[3].set_coefficient(0,0.5); wtrue[3].set_coefficient(1,0.5); wtrue[3].set_coefficient(3,1.0/sqrt(2.0));
                    wtrue[4].set_coefficient(0,0.5); wtrue[4].set_coefficient(1,0.5); wtrue[4].set_coefficient(3,1.0/sqrt(2.0));
                    wtrue[5].set_coefficient(0,0.5); wtrue[5].set_coefficient(1,0.5); wtrue[5].set_coefficient(3,1.0/sqrt(2.0));
                    wtrue[6].set_coefficient(0,0.5); wtrue[6].set_coefficient(1,0.5); wtrue[6].set_coefficient(3,1.0/sqrt(2.0));
                    wtrue[7].set_coefficient(0,0.5); wtrue[7].set_coefficient(1,0.5); wtrue[7].set_coefficient(3,1.0/sqrt(2.0));
    #else
                    // one generator at [0.5,0.75]^2
    #if _NUMBER_OF_TIME_STEPS == 10
                    for (unsigned int i=3; i<=7;++i)
    #else
    #if _NUMBER_OF_TIME_STEPS == 6
                    for (unsigned int i=2; i<=4;++i)
    #else
                    cout << "wnum == 5 requires number_of_timesteps == 6 or 10" << endl;
                    abort();
                    for (unsigned int i=1; i<=1;++i)
    #endif
    #endif
                    {
                        wtrue[i].set_coefficient(0,0.25);
                        wtrue[i].set_coefficient(1,-0.25);
                        wtrue[i].set_coefficient(2,-0.25);
                        wtrue[i].set_coefficient(3,0.25);
                        wtrue[i].set_coefficient(7,0.5);
                        wtrue[i].set_coefficient(11,0.5);
                        wtrue[i].set_coefficient(15,0.5);
                    }
    #endif
                    break;
                case 6: // use Exact_SolXD case 17 for certain time points
                {
    #if _DIMENSION == 1
                    Exact_Sol1D<17> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
    #else
                    Exact_Sol2D<17> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[2]);
                    for (unsigned int i=3; i<= 4; ++i)
                    {
                        wtrue[i] = wtrue[2];
                    }
    #endif
    // CLEANUP
                    // cout << "w_gen_coeffs = " << endl << w_gen_coeffs << endl;
                    // cout << "w_wav_coeffs = " << endl << w_wav_coeffs << endl;
                    //plot_haar_gen_coeffs(w_gen_coeffs, _HAAR_JMAX+1);

                }
                    break;            
                case 7: // use Exact_SolXD case 22 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<18> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<22> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                        wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                    }
    #endif
                }
                    break;
                case 8: // use Exact_SolXD case 18 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<18> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<18> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                        wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                    }
    #endif
                }
                    break;
                case 9: // use Exact_SolXD case 19 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<18> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<19> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        wtrue[i].scale(std::max(0.0, 5.0*(double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                    }
    #endif
                }
                    break;
                case 10: // use Exact_SolXD case 20 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<20> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<20> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                        wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                    }
    #endif
                }
                    break;
                case 11: // use Exact_SolXD case 21 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<21> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<21> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                    }
    #endif
                }
                    break;
                case 12: // use Exact_SolXD case 21 for certain time points
                {
    #if _DIMENSION == 1
                    cout << "angepasst an Rudolfs cutoff!" << endl;
                    Exact_Sol1D<21> wexact; // w as a function
                    FixedVector<double, (1<<(_HAAR_JMAX+1))> w_gen_coeffs,w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);
                    transform_fixedToIV (w_wav_coeffs, wtrue[3]);
                    for (unsigned int i=4; i<= 7; ++i)
                    {
                        wtrue[i] = wtrue[3];
                    }
                    for (unsigned int i=3; i<= 7; ++i)
                    {
                        wtrue[i].scale((std::max(0.0,25.0/16.0*((i-1)/10.0-0.9)*(0.1-(i-1)/10.0) )));
                    }
    #else
                    // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                    Exact_Sol2D<22> wexact;
                    FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                    precise_evaluate(&wexact, d, w_gen_coeffs);
                    haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                    unsigned int start_i, end_i;
                    switch (wtrue.size()-1) // number_of_timesteps
                    {
                        case 4:
                            start_i = 1;
                            end_i = 3;
                            break;
                        case 6:
                            start_i = 2;
                            end_i = 4;
                            break;
                        case 8:
                            start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                            end_i = 6;
                            break;
                        case 10:
                            start_i = 3;
                            end_i = 7;
                            break;
                        default:
                            cout << "cannot handle value for number_of_timesteps" << endl;
                            abort();
                            break;
                    }
                    transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                    for (unsigned int i=start_i+1; i<= end_i; ++i)
                    {
                        wtrue[i] = wtrue[start_i];
                    }
                    // temporal cutoff: 0 <= start_i and >=end_i

                    for (unsigned int i=start_i; i<=end_i; ++i)
                    {
                        wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                        //wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*2.0/(double)(end_i-start_i+2)     ));
                    }
    #endif
                }
                    break;
                default:
                    cout << "parabolic_tools:: first setup for coupling matrix W:: setup number " << first_setup_coupling_matrix_w << " unknown. aborting" << endl;
                    abort();
                    break;
            }
        }
        else //first_setup_coupling_matrix_w matrix > 12
            if (first_setup_coupling_matrix_w <= 28)
            {
                for (unsigned int i=0; i<wtrue.size(); i++)
                {
                    wtrue[i].clear();
                    wtrue[i].set_coefficient(first_setup_coupling_matrix_w-13,1);
                }
            }
            else // first_setup_coupling_matrix > 28
            {
#if _DIMENSION == 1
                cout << "parabolic_tools.cpp::initialize_coupling_matrix case " << first_setup_coupling_matrix_w << endl << "case not implemented!" << endl;
                abort();
#endif
                switch (first_setup_coupling_matrix_w)
                {
                    case 28: // bump: centered in x, up in y on [0,1/2]
                    {
#if _DIMENSION == 1
                        //not tested!
                        Exact_Sol1D<27> wexact;
                        FixedVector<double,(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#else
                        Exact_Sol2D<27> wexact;
                        FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#endif
                        precise_evaluate(&wexact, d, w_gen_coeffs);
                        haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                        transform_fixedToIV (w_wav_coeffs, wtrue[0]);
                        for (unsigned int i=0; i<= (wtrue.size()-1)/2; ++i)
                        {
                            wtrue[i] = wtrue[0];
                        }
                        break;
                    }
                    case 29: // bump: centered in x, down in y on [0,1/2]
                    {
#if _DIMENSION == 1
                        //not tested!
                        Exact_Sol1D<28> wexact;
                        FixedVector<double,(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#else
                        Exact_Sol2D<28> wexact;
                        FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#endif
                        precise_evaluate(&wexact, d, w_gen_coeffs);
                        haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                        transform_fixedToIV (w_wav_coeffs, wtrue[0]);
                        for (unsigned int i=0; i<= (wtrue.size()-1)/2; ++i)
                        {
                            wtrue[i] = wtrue[0];
                        }
                        break;
                    }
                    case 30: // bump: left in x, centered in y on [0,1/2]
                    {
#if _DIMENSION == 1
                        //not tested!
                        Exact_Sol1D<29> wexact;
                        FixedVector<double,(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#else
                        Exact_Sol2D<29> wexact;
                        FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#endif
                        precise_evaluate(&wexact, d, w_gen_coeffs);
                        haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                        transform_fixedToIV (w_wav_coeffs, wtrue[0]);
                        for (unsigned int i=0; i<= (wtrue.size()-1)/2; ++i)
                        {
                            wtrue[i] = wtrue[0];
                        }
                        break;
                    }
                    case 31: // bump: right in x, centered in y on [0,1/2]
                    {
#if _DIMENSION == 1
                        //not tested!
                        Exact_Sol1D<30> wexact;
                        FixedVector<double,(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#else
                        Exact_Sol2D<30> wexact;
                        FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
#endif
                        precise_evaluate(&wexact, d, w_gen_coeffs);
                        haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                        transform_fixedToIV (w_wav_coeffs, wtrue[0]);
                        for (unsigned int i=0; i<= (wtrue.size()-1)/2; ++i)
                        {
                            wtrue[i] = wtrue[0];
                        }
                        break;
                    }
                    case 32:
                    {
#if _DIMENSION == 1
                        abort();
#else
                        // andere Abschneidefunktion, quadratisch symmetrisch in der Zeit!
                        Exact_Sol2D<32> wexact;
                        FixedMatrix<double,(1<<(_HAAR_JMAX+1)),(1<<(_HAAR_JMAX+1))> w_gen_coeffs, w_wav_coeffs;
                        precise_evaluate(&wexact, d, w_gen_coeffs);
                        haar_wavelet_transform(w_gen_coeffs,w_wav_coeffs);

                        unsigned int start_i, end_i;
                        switch (wtrue.size()-1) // number_of_timesteps
                        {
                            case 4:
                                start_i = 1;
                                end_i = 3;
                                break;
                            case 6:
                                start_i = 2;
                                end_i = 4;
                                break;
                            case 8:
                                start_i = 2; // 5 active, 4 inactive positions. A possible alternative would be 3active/6inactive
                                end_i = 6;
                                break;
                            case 10:
                                start_i = 3;
                                end_i = 7;
                                break;
                            default:
                                cout << "cannot handle value for number_of_timesteps" << endl;
                                abort();
                                break;
                        }
                        transform_fixedToIV (w_wav_coeffs, wtrue[start_i]);
                        for (unsigned int i=start_i+1; i<= end_i; ++i)
                        {
                            wtrue[i] = wtrue[start_i];
                        }
                        // temporal cutoff: 0 <= start_i and >=end_i

                        for (unsigned int i=start_i; i<=end_i; ++i)
                        {
                            wtrue[i].scale(std::max(0.0, (double)(i-start_i+1)*(double)(end_i+1-i)*4.0/(double)((end_i-start_i+2)*(end_i-start_i+2))     ));
                        }
#endif
                        break;
                    }
                default:
                
                    cout << "parabolic_tools:: first setup for coupling matrix W:: setup number " << first_setup_coupling_matrix_w << " unknown. aborting" << endl;
                        abort();
                        break;
                } // end of switch
            }
    }
    
    
    template < unsigned int ONEDIMHAARCOUNT, unsigned int NUMOFTIMESTEPSPLUSONE >
    void initialize_coupling_matrix(FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMOFTIMESTEPSPLUSONE> &w_true,
                                    const int setup_number)
    {
        int onequarter = ONEDIMHAARCOUNT/4;
        switch (setup_number)
        {
            case 10: // 1 patch, constant value 0, constant in time
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(1);
                    w_true[t][0] = FixedMatrix<double, ONEDIMHAARCOUNT> (0);
                }
                break;
            case 11: // 1 patch, constant value 1, constant in time
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(1);
                    w_true[t][0] = FixedMatrix<double, ONEDIMHAARCOUNT> (1);
                }
                break;
            case 12: // 1 patch, constant value 1 + 1 bump, constant in time
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(1);
                    w_true[t][0] = FixedMatrix<double, ONEDIMHAARCOUNT> (1);
                    for (int i=onequarter; i< 2*onequarter;++i)
                    {
                        for (int j=onequarter; j<= 2*onequarter;++j)
                        {
                            w_true[t][0].set_entry(i,j,3);
                        }
                    }
                }
                break;
            case 30: // 3 patches, L-domain (northeast is missing), constant value 0
                    // constant value 1 + bump in the upper domain (number 0) + bump across the interface of the lower domains (number 1,2), constant in time
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(3);
                    for (int p = 0; p<3; ++p)
                    {
                        w_true[t][p] = FixedMatrix<double, ONEDIMHAARCOUNT> (0);
                    }
                }
                break;
            case 31: // 3 patches, L-domain (northeast is missing), 
                    // constant value 1 + bump in the upper domain (number 0) + bump across the interface of the lower domains (number 1,2), constant in time
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(3);
                    for (int p = 0; p<3; ++p)
                    {
                        w_true[t][p] = FixedMatrix<double, ONEDIMHAARCOUNT> (1);
                    }
                    for (int i=0; i< onequarter;++i)
                    {
                        for (int j=onequarter; j<= 2*onequarter;++j)
                        {
                            w_true[t][0].set_entry(i+onequarter,j+onequarter,3);
                            w_true[t][1].set_entry(ONEDIMHAARCOUNT-i, j+2*onequarter,3);
                            w_true[t][2].set_entry(i, j+2*onequarter,3);
                        }
                    }
                }
                break;
            case 32: // 3 patches, L-domain (northeast is missing), 
                    // constant value 1 + bump in the upper domain (number 0) + bump across the interface of the lower domains (number 1,2), 
                    // bumps are only active from timestep 2,...,num_of_timesteps-2 
                for (int t=0; t<NUMOFTIMESTEPSPLUSONE;++t)
                {
                    w_true[t].resize(3);
                    for (int p = 0; p<3; ++p)
                    {
                        w_true[t][p] = FixedMatrix<double, ONEDIMHAARCOUNT> (1);
                    }
                }
                assert (NUMOFTIMESTEPSPLUSONE > 1);
                for (int t=2; t<NUMOFTIMESTEPSPLUSONE-2;++t)
                {
                    for (int i=0; i< onequarter;++i)
                    {
                        for (int j=onequarter; j<= 2*onequarter;++j)
                        {
                            w_true[t][0].set_entry(i+onequarter,j+onequarter,3);
                            w_true[t][1].set_entry(ONEDIMHAARCOUNT-i, j+2*onequarter,3);
                            w_true[t][2].set_entry(i, j+2*onequarter,3);
                        }
                    }
                }
                break;    
            default:
                cout << "initialize_coupling_matrix: unknown setup_number!" << endl;
                abort();
                break;
        }
    }
    
    
    
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
                              const double uobserved_coeffs_norm, 
                              const unsigned int number_of_timesteps)
    {
        logstream << "% Matlab Logfile!!\n%";
        logstream << "\n% Unknown Parameter Setting (wnum) = ";
        logstream << first_setup_coupling_matrix_w;
        if (_PROBLEMS_USED_FOR_MORE_THAN_ONLY_U0 == 1)
        {
            logstream << "\n% Initial value u0 is given by _FORWARD_PROBLEM_NO = " << _FORWARD_PROBLEM_NO;
            logstream << "\n% ydata is computed using given W and u0";
        }
        else
        {
            logstream << "\n% Initial value u0 and ydata are given by _FORWARD_PROBLEM_NO = " << _FORWARD_PROBLEM_NO;
        }
        logstream << "\n% Criterion used to stop iteration:";
        if (max_inverse_iterations == 0)
        {
            logstream << "\n% relative l1 error of update <= tolerance = " << stopping_tolerance;
        }
        else
        {
            logstream << "\n% number of iteration limited by " << max_inverse_iterations << " Iterations";
        }
        logstream << "\n% Maximal wavelet level for the forward solver jmax = " << spatial_jmax;
        if (w_limit_number == 0)
        {
            logstream << "\n% There is no approximation to the reconstructed parameter W_limit available. Approximation errors will be computed with respect to (the unknown and in general different) W_true given by wnum = " << first_setup_coupling_matrix_w;
        }
        else
        {
            logstream << "\n% Approximation errors will be computed with the approximation of W_limit after " << w_limit_number << " Iterations.";
        }
        logstream << "\n% Shrinkage parameter alpha = " << alpha;
        if (compute_noise == 1)
        {
            logstream << "\n% A new noise vector is computed and stored on disk";
        }
        else
        {
            logstream << "\n% The noise is loaded from disk";
        }
        logstream << "\n% Noise level delta = " << delta;;
        logstream << "\n% ||y||_ell2 = " << utrue_coeffs_norm;
        logstream << "\n% ||y^delta||_ell2 = " << uobserved_coeffs_norm;
        logstream << "\n% Relative noise level = delta/||y||_ell2 = " << (delta / utrue_coeffs_norm);
        logstream << "\n% W for the inverse iteration is initialized with setting = ";
        logstream << _SECOND_SETUP_COUPLING_MATRIX_W;
        logstream << "\n% time_needed = [ total, assemble, solve forward problem, solve backward problem, compute update, shrinkage]\n%\n";

        logstream << "t = [0 0 0 0 0 0];\n";

        if (_COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL == 1)

        {
            logstream << "we1 = zeros (1," << (number_of_timesteps+1) << "); % ||W-W_true||_ell1 (for each timestep) \n"; // ell 1 norm of w-w_limit
            logstream << "we2 = zeros (1," << (number_of_timesteps+1) << "); % ||W-W_limit||_ell2\n"; // ell 2 norm of w-w_limit
            logstream << "w0 = zeros (1," << (number_of_timesteps+1) << ");% ||W||_ell0\n"; // number of coeffs of w
            logstream << "w1 = zeros (1," << (number_of_timesteps+1) << ");% ||W||_ell1\n"; // ell 1 norm of w
            logstream << "w2 = zeros (1," << (number_of_timesteps+1) << ");% ||W||_ell2\n"; // ell 2 norm of w
        }
        if (_COMPARE_FORWARD_SOLUTION_WITH_TRUE_SOL == 1)
        {
            logstream << "ue2 = zeros (1," << (number_of_timesteps+1) << ");\n";
            logstream << "u0 = zeros (1," << (number_of_timesteps+1) << ");\n";
            logstream << "u2 = zeros (1," << (number_of_timesteps+1) << ");\n";
            //logstream << "v0 = zeros (1," << (number_of_timesteps+1) << ");\n";
            //logstream << "v2 = zeros (1," << (number_of_timesteps+1) << ");\n";
        }
    };
    
    template <class PROBLEM, class CTPROBLEM, unsigned int NUMOFTIMESTEPS>
    void store_uexact(const AffLinParEq_Precomputed<PROBLEM, NUMOFTIMESTEPS> &parabolic, const CTPROBLEM& ctgramian, const char* true_forward_solution_filename_start, const char* true_forward_solution_filename_coeffs_end, const char* true_forward_solution_filename_matlab_end, const int resolution)
    {
        char filename[250];
        std::ofstream filestream;
        for (unsigned int i=0; i<=NUMOFTIMESTEPS;++i)
        {
            std::ostringstream output_filename;
            output_filename << true_forward_solution_filename_start << i << true_forward_solution_filename_matlab_end;
            filestream.open(output_filename.str().c_str());
            SampledMapping<PROBLEM::space_dimension> ui_plot(evaluate(ctgramian.basis(), parabolic.forward_solution_[i], true, resolution));
            sprintf(filename,"%s%d%s",true_forward_solution_filename_start,i,true_forward_solution_filename_coeffs_end);
            cout << "store_uexact:: saving the true_forward_solution[" << i << "] to files" << endl 
                 << true_forward_solution_filename_start << i << true_forward_solution_filename_matlab_end << endl
                 << filename << endl;
            ui_plot.matlab_output(filestream);
            filestream.close();
            cout << "sore_uexact:: write utrue_coeffs[" << i << "] to file = " << endl << filename << endl;
            cout << "current version uses Index class. Maybe inefficient?" << endl;
            writeIVToFile(parabolic.forward_solution_[i],filename);
            
            
            
            /*
            f_.clear();

        InfiniteVector<double,int> temp_iv;
        readIVFromFile (temp_iv,f_filename);
        for (InfiniteVector<double, int>::const_iterator it(temp_iv.begin()), itend(temp_iv.end()); it!=itend; ++it)
        {
            cout << it.index() << endl;
            //cout << it.index().number() <<endl;
            f_.set_coefficient(problem_->basis().get_wavelet(it.index()),*it);
        }
        */
        
        }
    };
    
    template <unsigned int NUMOFTIMESTEPS, class QTBASIS, unsigned int PRECISE_EVALUATE_GRANULARITY, unsigned int ONEDIMHAARCOUNT>
    void store_uexact(const AffLinParEq_qtbasis<NUMOFTIMESTEPS, QTBASIS, PRECISE_EVALUATE_GRANULARITY, ONEDIMHAARCOUNT > &parabolic_problem, 
            const char* true_forward_solution_filename_start, 
            const char* true_forward_solution_filename_coeffs_end, 
            const char* true_forward_solution_filename_matlab_end, 
            const int resolution)
    {
        
        Array1D<SampledMapping<QTBASIS::space_dimension> > ergebnis;
        char filename[250];
        std::ofstream filestream;
        for (unsigned int i=0; i<=NUMOFTIMESTEPS;++i)
        {
            std::ostringstream output_filename;
            output_filename << true_forward_solution_filename_start << i << true_forward_solution_filename_matlab_end;
            //SampledMapping<PROBLEM::space_dimension> ui_plot(evaluate(ctgramian.basis(), parabolic.forward_solution_[i], true, resolution));
            sprintf(filename,"%s%d%s",true_forward_solution_filename_start,i,true_forward_solution_filename_coeffs_end);
            cout << "store_uexact:: saving the true_forward_solution[" << i << "] to files" << endl 
                 << true_forward_solution_filename_start << i << true_forward_solution_filename_matlab_end << endl
                 << filename << endl;
            filestream.open(output_filename.str().c_str());
            ergebnis = parabolic_problem.basis_->sampled_output(parabolic_problem.forward_solution_[i],
                    true,
                    resolution);
            matlab_output(ergebnis, filestream);
            //ui_plot.matlab_output(filestream);
            filestream.close();
            cout << "sore_uexact:: write utrue_coeffs[" << i << "] to file = " << endl << filename << endl;
            //parabolic_problem.forward_solution_[i].writeToFile(filename);
            writeIVToFile(parabolic_problem.forward_solution_[i],filename);
        }
    };
    
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
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& utrue_coeffs)
    {
    }
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
            const double shrinkage_time)
    {
        logstream << "% Iteration " << iteration_count << "\n";
        // write computation time to logfile
        logstream << "t = [t; " << total_time << " " << assemble_time << " " << forward_time << " " << backward_time << " " << update_time << " " << shrinkage_time << "];\n";

    }
    */
    
    template <class C, class I>
    void compute_norms(const InfiniteVector<C,I> &x,const InfiniteVector<C,I> &y, double& diff1, double& diff2, double& x1, double& x2)
    {
        typename InfiniteVector<C,I>::const_iterator itx(x.begin()), itxend(x.end()),
          ity(y.begin()), ityend(y.end());
        diff1=0;diff2=0;x1=0;x2=0;
        while (itx != itxend && ity != ityend)
          {
            if (itx.index() < ity.index())
              {
                diff1 += fabs(*itx);
                diff2 += *itx * *itx;
                x1 += fabs(*itx);
                x2 += *itx * *itx;
                ++itx;
              }
            else
              {
                if (ity.index() < itx.index())
                  {
                    diff1 += fabs(*ity);
                    diff2 += *ity * *ity;
                    ++ity;
                  }
                else
                  {
                    const C value(*itx - *ity);
                    if (value != C(0)) {
                        diff1 += fabs(value);
                        diff2 += value * value;
                    }
                    x1 += fabs(*itx);
                    x2 += *itx * *itx;
                    ++itx;
                    ++ity;
                  }
              }
          }
        while (itx != itxend)
          {
            diff1 += fabs(*itx);
            diff2 += *itx * *itx;
            x1 += fabs(*itx);
            x2 += *itx * *itx;
            ++itx;
          }
        while (ity != ityend)
          {
            diff1 += fabs(*ity);
            diff2 += *ity * *ity;
            ++ity;
          }
        diff2 = sqrt(diff2);
        x2 = sqrt(x2);
    }
    
    
    template <class C, class I, unsigned int ONEDIMHAARCOUNT>
    void compute_norms(const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > &x,
            const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > &y, 
            double& diff1, 
            double& diff2, 
            double& x1, 
            double& x2)
    {
        const unsigned int nop(x[0].size());
        double temp_d;
        diff1=0;diff2=0;x1=0;x2=0;
        for (unsigned int p=0; p<nop;++p)
        {
            for (unsigned int row = 0; row<ONEDIMHAARCOUNT; ++row)
            {
                for (unsigned int col = 0; col<ONEDIMHAARCOUNT; ++col)
                {
                    temp_d = x[p](row,col);
                    x1 += temp_d;
                    x2 += temp_d*temp_d;
                    temp_d = temp_d-y[p](row,col);
                    diff1 += abs (temp_d);
                    diff2 += temp_d*temp_d;
                }
            }
        }
        diff2 = sqrt(diff2);
        x2 = sqrt(x2);
    }
    
    template <class C, class I>
    void compute_norms2(const InfiniteVector<C,I> &x,const InfiniteVector<C,I> &y, double& diff2, double& x2)
    {
        typename InfiniteVector<C,I>::const_iterator itx(x.begin()), itxend(x.end()),
                ity(y.begin()), ityend(y.end());
        diff2=0;x2=0;
        while ((itx != itxend) && (ity != ityend))
        {
            if (itx.index() < ity.index())
            {
                diff2 += *itx * *itx;
                x2 += *itx * *itx;
                ++itx;
            }
            else
            {
                if (ity.index() < itx.index())
                {
                    diff2 += *ity * *ity;
                    ++ity;
                }
                else
                {
                    const C value(*itx - *ity);
                    if (value != C(0)) 
                    {
                        diff2 += value * value;
                    }
                    x2 += *itx * *itx;
                    ++itx;
                    ++ity;
                }
            }
        }
        while (itx != itxend)
        {
            diff2 += *itx * *itx;
            x2 += *itx * *itx;
            ++itx;
        }
        while (ity != ityend)
        {
            diff2 += *ity * *ity;
            ++ity;
        }
        diff2 = sqrt(diff2);
        x2 = sqrt(x2);
    }
    
    template <unsigned int NUMBER_OF_TIMESTEPS>
    void log_reconstruction_errors(std::ofstream& logstream, 
            const Array1D<InfiniteVector<double, int> >& w, 
            const Array1D<InfiniteVector<double, int> >& wlimit)
    {
#if _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL == 1
        //include we1,we2,w0,w1,w2 in logfile
            Array1D<double> we1(NUMBER_OF_TIMESTEPS+1),we2(NUMBER_OF_TIMESTEPS+1),w1(NUMBER_OF_TIMESTEPS+1),w2(NUMBER_OF_TIMESTEPS+1);
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                compute_norms<double, int>(w[i],wlimit[i],we1[i],we2[i],w1[i],w2[i]);
            }
            
            logstream << "we1 = [we1;";
            cout << "update_logstream:: ell_1 error for the estimated parameter w" << endl;
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                //temp_d = l1_norm(w[i]-wlimit[i]);
                cout << "    t=" << i/(double)NUMBER_OF_TIMESTEPS << ": " << we1[i] << endl;
                logstream << " " << we1[i];
            }
            logstream << "];\n";
            //cout << "main:: ell_2 error for the estimated parameter w" << endl;
            logstream << "we2 = [we2;";
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                //temp_d = l2_norm(wlimit[i]-w[i]);
                //cout << "    t=" << i/(double)number_of_timesteps << ": " << temp_d << endl;
                logstream << " " << we2[i];
            }
            logstream << "];\n";
            logstream << "w0 = [w0;";
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                //temp_i = w[i].size();
                logstream << " " << w[i].size();
            }
            logstream << "];\n";
            logstream << "w1 = [w1;";
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                //temp_d = l1_norm(w[i]);
                logstream << " " << w1[i];
            }
            logstream << "];\n";
            logstream << "w2 = [w2;";
            for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
            {
                //temp_d = l2_norm(w[i]);
                logstream << " " << w2[i];
            }
            logstream << "];\n";
#endif // finished writing errors related to w
    }
    
                           
    template <unsigned int NUMBEROFTIMESTEPSPLUSONE, unsigned int ONEDIMHAARCOUNT >
    void log_reconstruction_errors(std::ofstream& logstream,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBEROFTIMESTEPSPLUSONE>& w, 
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBEROFTIMESTEPSPLUSONE>& wlimit)
    {
#if _COMPARE_INVERSE_SOLUTION_WITH_TRUE_SOL == 1
        //include we1,we2,w0,w1,w2 in logfile
        FixedArray1D<double, NUMBEROFTIMESTEPSPLUSONE> we1, we2, w1, w2;
        for (int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            compute_norms(w[i],wlimit[i],we1[i],we2[i],w1[i],w2[i]);
        }

        logstream << "we1 = [we1;";
        cout << "update_logstream:: ell_1 error for the estimated parameter w" << endl;
        for (int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            //temp_d = l1_norm(w[i]-wlimit[i]);
            cout << "    t=" << i/(double)(NUMBEROFTIMESTEPSPLUSONE-1) << ": " << we1[i] << endl;
            logstream << " " << we1[i];
        }
        logstream << "];\n";
        //cout << "main:: ell_2 error for the estimated parameter w" << endl;
        logstream << "we2 = [we2;";
        for (int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            //temp_d = l2_norm(wlimit[i]-w[i]);
            //cout << "    t=" << i/(double)number_of_timesteps << ": " << temp_d << endl;
            logstream << " " << we2[i];
        }
        logstream << "];\n";
        logstream << "w0 = [w0;";
        for (int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            //temp_i = w[i].size();
            logstream << " " << w[i].size();
        }
        logstream << "];\n";
        logstream << "w1 = [w1;";
        for (unsigned int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            //temp_d = l1_norm(w[i]);
            logstream << " " << w1[i];
        }
        logstream << "];\n";
        logstream << "w2 = [w2;";
        for (unsigned int i=0; i< NUMBEROFTIMESTEPSPLUSONE; ++i)
        {
            //temp_d = l2_norm(w[i]);
            logstream << " " << w2[i];
        }
        logstream << "];\n";
#endif // finished writing errors related to w
    }
            
    
    template <class TBASIS, unsigned int NUMBER_OF_TIMESTEPS>
    void log_solution_errors(std::ofstream& logstream, 
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& forward_solution,
            //const Array1D<InfiniteVector<double,typename TBASIS::Index> >& backward_solution,
            //const Array1D<InfiniteVector<double,typename TBASIS::Index> >& utrue_coeffs)
            const Array1D<InfiniteVector<double,typename TBASIS::Index> >& uobserved_coeffs)
    {
        // time forward case:
        Array1D<double> diff2(NUMBER_OF_TIMESTEPS+1),fs2(NUMBER_OF_TIMESTEPS+1);
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            compute_norms2<double, typename TBASIS::Index>(forward_solution[i],uobserved_coeffs[i],diff2[i],fs2[i]);
        }
        logstream << "ue2 = [ue2;";
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            //temp_d = l2_norm(parabolic.forward_solution_[i] - utrue_coeffs[i]);
            //cout << "    t=" << i/(double)number_of_timesteps << ": " << temp_d << endl;
            logstream << " " << diff2[i]; //l2_norm(forward_solution[i] - utrue_coeffs[i]);
        }
        logstream << "];\n";
        logstream << "u0 = [u0;";
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            //temp_i = parabolic.forward_solution_[i].size();
            logstream << " " << forward_solution[i].size();
        }
        logstream << "];\n";
        //logstream << "u1 = [u1;";
        //for (unsigned int i=0; i<= number_of_timesteps; ++i)
        //{
        //    temp_d = l1_norm(parabolic.forward_solution_[i]);
        //    logstream << " " << temp_d;
        //}
        //logstream << "];\n";
        logstream << "u2 = [u2;";
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            //temp_d = l2_norm(parabolic.forward_solution_[i]);
            logstream << " " << fs2[i];//l2_norm(forward_solution[i]);
        }
        logstream << "];\n";

        // time backward case:

        //cout << "main:: ell_2 errors for the time backward case:" << endl;
        /*
        logstream << "v0 = [v0;";
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            //temp_i = parabolic.backward_solution_[i].size();
            logstream << " " << backward_solution[i].size();
        }
        logstream << "];\n";
         */
        /*
        logstream << "v1 = [u1;";
        for (unsigned int i=0; i<= number_of_timesteps; ++i)
        {
            //temp_d = l1_norm(parabolic.backward_solution_[i]);
            logstream << " " << l1_norm(backward_solution[i]);
        }
        logstream << "];\n";
         */
        /*
        logstream << "v2 = [v2;";
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            //temp_d = l2_norm(parabolic.backward_solution_[i]);
            logstream << " " << l2_norm(backward_solution[i]);
        }
        logstream << "];\n";
         */
        // finished writing errors related to u,v
    }
    
    template <unsigned int NUMBER_OF_TIMESTEPS_PLUS_ONE>
    void log_solution_errors(std::ofstream& logstream,
            const FixedArray1D<InfiniteVector<double,int>, NUMBER_OF_TIMESTEPS_PLUS_ONE >& forward_solution, 
            const FixedArray1D<InfiniteVector<double,int>, NUMBER_OF_TIMESTEPS_PLUS_ONE >& uobserved_coeffs)
    {
        // time forward case:
        FixedArray1D<double, NUMBER_OF_TIMESTEPS_PLUS_ONE> diff2,fs2;
        for (unsigned int i=0; i< NUMBER_OF_TIMESTEPS_PLUS_ONE; ++i)
        {
            compute_norms2(forward_solution[i],uobserved_coeffs[i],diff2[i],fs2[i]);
        }
        logstream << "ue2 = [ue2;";
        for (unsigned int i=0; i< NUMBER_OF_TIMESTEPS_PLUS_ONE; ++i)
        {
            logstream << " " << diff2[i]; //l2_norm(forward_solution[i] - utrue_coeffs[i]);
        }
        logstream << "];\n";
        logstream << "u0 = [u0;";
        for (unsigned int i=0; i< NUMBER_OF_TIMESTEPS_PLUS_ONE; ++i)
        {
            logstream << " " << forward_solution[i].size();
        }
        logstream << "];\n";
        logstream << "u2 = [u2;";
        for (unsigned int i=0; i< NUMBER_OF_TIMESTEPS_PLUS_ONE; ++i)
        {
            logstream << " " << fs2[i];//l2_norm(forward_solution[i]);
        }
        logstream << "];\n";
    }
    
    
    
    template<unsigned int NUMBER_OF_TIMESTEPS, unsigned int NUMOFWAVELETSPERCOORDINATE>
    void plot_solutions(std::ofstream& plotstream, 
            const Array1D<InfiniteVector<double, int> >& w,
            const Array1D<InfiniteVector<double, int> >& w_true,
            const Array1D<InfiniteVector<double, int> >& w_limit)
    {
        // write plots of w to plotstream
#if _DIMENSION == 1
        FixedVector<double, (1<<(_HAAR_JMAX+1))> temp_haar_wav_coeffs, temp_haar_gen_coeffs;
#else
        FixedMatrix<double, (1<<(_HAAR_JMAX+1)), (1<<(_HAAR_JMAX+1))> temp_haar_wav_coeffs, temp_haar_gen_coeffs;
#endif
        for (unsigned int i=0; i<= NUMBER_OF_TIMESTEPS; ++i)
        {
            transform_IVTofixed(w[i], temp_haar_wav_coeffs);
            inverse_haar_wavelet_transform(temp_haar_wav_coeffs, temp_haar_gen_coeffs);
            plotstream << "w_plot{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE)*temp_haar_gen_coeffs[l];
            }
            plotstream << "];\n";
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*temp_haar_gen_coeffs.get_entry(l,k); // we do not need to reverse ordering in the 2nd component, since this is a matlab plot
                }
                plotstream << "; ";
            }
            plotstream << "];\n";
#endif

            transform_IVTofixed(w_true[i], temp_haar_wav_coeffs);
            inverse_haar_wavelet_transform(temp_haar_wav_coeffs, temp_haar_gen_coeffs);
            plotstream << "w_true_plot{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE)*temp_haar_gen_coeffs[l];
            }
            plotstream << "];\n";
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*temp_haar_gen_coeffs.get_entry(l,k); // we do not need to reverse ordering in the 2nd component, since this is a matlab plot
                }
                plotstream << "; ";
            }
            plotstream << "];\n";
#endif

            transform_IVTofixed(w_limit[i], temp_haar_wav_coeffs);
            inverse_haar_wavelet_transform(temp_haar_wav_coeffs, temp_haar_gen_coeffs);
            plotstream << "w_limit_plot{" << (i+1) << "} = [";
#if _DIMENSION == 1
            for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
            {
                plotstream << " " << sqrt(NUMOFWAVELETSPERCOORDINATE)*temp_haar_gen_coeffs[l];
            }
            plotstream << "];\n";
#else
            for (unsigned int k=0; k < NUMOFWAVELETSPERCOORDINATE; ++k)
            {
                for (unsigned int l=0; l < NUMOFWAVELETSPERCOORDINATE; ++l)
                {
                    plotstream << " " << NUMOFWAVELETSPERCOORDINATE*temp_haar_gen_coeffs.get_entry(l,k); // we do not need to reverse ordering in the 2nd component, since this is a matlab plot
                }
                plotstream << "; ";
            }
            plotstream << "];\n";
#endif
        } // end of write plots of w to logfile
    };

    template<unsigned int NUMBER_OF_TIMESTEPS_PLUS_ONE, unsigned int ONEDIMHAARCOUNT>
    void plot_solutions(std::ofstream& plotstream, 
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w_true,
            const FixedArray1D<Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> >, NUMBER_OF_TIMESTEPS_PLUS_ONE>& w_limit)
    {
        
        // write plots of w to plotstream
        FixedMatrix<double, ONEDIMHAARCOUNT> temp_haar_gen_coeffs;
        const unsigned int nop(w[0].size());
        //FixedMatrix<double, (1<<(_HAAR_JMAX+1)), (1<<(_HAAR_JMAX+1))> temp_haar_wav_coeffs, temp_haar_gen_coeffs;
        for (unsigned int t=0; t< NUMBER_OF_TIMESTEPS_PLUS_ONE; ++t)
        {
            for (unsigned int p=0; p<nop; ++p)
            {
                //transform_IVTofixed(w[i], temp_haar_wav_coeffs);
                inverse_haar_wavelet_transform(w[t][p], temp_haar_gen_coeffs);
                plotstream << "w_plot{" << (t+1)*(p+1) << "} = [";
                for (unsigned int k=0; k < ONEDIMHAARCOUNT; ++k)
                {
                    for (unsigned int l=0; l < ONEDIMHAARCOUNT; ++l)
                    {
                        plotstream << " " << ONEDIMHAARCOUNT*temp_haar_gen_coeffs.get_entry(k,l);
                    }
                    plotstream << "; ";
                }
                plotstream << "];\n";
                //transform_IVTofixed(w_true[i], temp_haar_wav_coeffs);
                inverse_haar_wavelet_transform(w_true[t][p], temp_haar_gen_coeffs);
                plotstream << "w_true_plot{" << (t+1)*(p+1) << "} = [";
                for (unsigned int k=0; k < ONEDIMHAARCOUNT; ++k)
                {
                    for (unsigned int l=0; l < ONEDIMHAARCOUNT; ++l)
                    {
                        plotstream << " " << ONEDIMHAARCOUNT*temp_haar_gen_coeffs.get_entry(k,l);
                    }
                    plotstream << "; ";
                }
                plotstream << "];\n";
                //transform_IVTofixed(w_limit[i], temp_haar_wav_coeffs);
                inverse_haar_wavelet_transform(w_limit[t][p], temp_haar_gen_coeffs);
                plotstream << "w_limit_plot{" << (t+1)*(p+1) << "} = [";
                for (unsigned int k=0; k < ONEDIMHAARCOUNT; ++k)
                {
                    for (unsigned int l=0; l < ONEDIMHAARCOUNT; ++l)
                    {
                        plotstream << " " << ONEDIMHAARCOUNT*temp_haar_gen_coeffs.get_entry(k,l); // we do not need to reverse ordering in the 2nd component, since this is a matlab plot
                    }
                    plotstream << "; ";
                }
                plotstream << "];\n";
            }
        } // end of write plots of w to logfile
    };
       
}

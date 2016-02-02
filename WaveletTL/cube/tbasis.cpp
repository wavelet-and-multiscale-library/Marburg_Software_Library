// implementation for tbasis.h

namespace WaveletTL
{
  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::TensorBasis()
    : bases_()
	{
    	// we only need one instance of IBASIS, without b.c.
    	IBASIS* b = new IBASIS();
    	bases_infact.push_back(b);
    	for (unsigned int i = 0; i < DIM; i++)
            bases_[i] = b;                    
        j0_[0] = b->j0();
        for (unsigned int i = 1; i < DIM; i++)
            j0_[i] = j0_[0];
    	delete_pointers = true;
  	}

/* Works only with DSBasis, but not with PBasis : 
  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::TensorBasis(const FixedArray1D<int,2*DIM>& s, const FixedArray1D<int,2*DIM>& sT)
	{
    	for (unsigned int i = 0; i < DIM; i++) {
      		// check whether the corresponding 1d basis already exists
      		IBASIS* b = 0;
      		for (typename list<IBASIS*>::const_iterator it(bases_infact.begin()); it != bases_infact.end(); ++it)
                {
                    if ((*it)->get_s0() == s[2*i]
                            && (*it)->get_s1() == s[2*i+1]
                            && (*it)->get_sT0() == sT[2*i]
                            && (*it)->get_sT1() == sT[2*i+1]) {
                        b = *it;
                        break;
                    }
      		}
      		if (b == 0) {
				//b = new IBASIS("",s[2*i], s[2*i+1], sT[2*i], sT[2*i+1]);
				b = new IBASIS(s[2*i], s[2*i+1], sT[2*i], sT[2*i+1]);
				bases_infact.push_back(b);
      		}
      		bases_[i] = b;
      		j0_[i] = b->j0();
    	}
    	delete_pointers = true;
  	}
*/
  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::TensorBasis(const FixedArray1D<int,2*DIM>& s) {
    	for (unsigned int i = 0; i < DIM; i++) {
      		// check whether the corresponding 1d basis already exists
      		IBASIS* b = 0;
      		for (typename list<IBASIS*>::const_iterator it(bases_infact.begin()); it != bases_infact.end(); ++it) {
      			if ((*it)->get_s0() == s[2*i]
	    			&& (*it)->get_s1() == s[2*i+1]) {
	  				b = *it;
	  				break;
				}
      		}
      		if (b == 0) {
				b = new IBASIS(s[2*i], s[2*i+1]);
				bases_infact.push_back(b);
      		}
      		bases_[i] = b;
      		j0_[i] = b->j0();
    	}
    	delete_pointers = true;
  	}

  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::TensorBasis(const FixedArray1D<bool,2*DIM>& bc) {
    	for (unsigned int i = 0; i < DIM; i++) {
      		// check whether the corresponding 1d basis already exists
      		IBASIS* b = 0;
      		for (typename list<IBASIS*>::const_iterator it(bases_infact.begin()); it != bases_infact.end(); ++it) {
				if (((*it)->get_s0()==1) == bc[2*i]
	    			&& ((*it)->get_s1()==1) == bc[2*i+1]) {
	  				b = *it;
	  				break;
				}
      		}
      		if (b == 0) {
				b = new IBASIS(bc[2*i], bc[2*i+1]);
				bases_infact.push_back(b);
      		}
      		bases_[i] = b;
      		j0_[i] = b->j0();
    	}
    	delete_pointers = true;
  	}

  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::TensorBasis(const FixedArray1D<IBASIS*,DIM> bases) {
    	for (unsigned int i = 0; i < DIM; i++) {
      		bases_[i] = bases[i];
      		j0_[i] = bases_[i]->j0();
    	}
    	delete_pointers = false;
  	}

  	template <class IBASIS, unsigned int DIM>
  	TensorBasis<IBASIS,DIM>::~TensorBasis()
  	{
	    if (delete_pointers) {
			for (typename list<IBASIS*>::const_iterator it(bases_infact.begin()); it != bases_infact.end(); ++it)
				delete *it;
            }
  	}


	template <class IBASIS, unsigned int DIM>
  	void
  	TensorBasis<IBASIS,DIM>::support(const Index& lambda, Support& supp) const
  	{
            // PERFORMANCE?
            typename Index::level_type temp_j = lambda.j();
            typename Index::type_type temp_e = lambda.e();
            //MultiIndex<int,DIM> temp_j=lambda.j(),temp_e=lambda.e();
            for (unsigned int i(0); i < DIM; i++) {
                supp.j[i] = temp_j[i]+temp_e[i];
                //supp.j[i] = lambda.j()[i] + lambda.e()[i];
                bases()[i]->support(temp_j[i],
                                    temp_e[i],
                                    lambda.k()[i],
                                    supp.a[i],
                                    supp.b[i]);
            }
  	}
  
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_generator() const
        {
            typename Index::type_type e;
            typename Index::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i]=0;
                k[i]=bases_[i]->DeltaLmin();
            }
            return Index(j0_, e, k, 0, this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_generator(const unsigned int j) const
        {
            typename Index::type_type e;
            typename Index::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i]=0;
                k[i]=bases_[i]->DeltaLmin();
            }
            return Index(j0_, e, k, 0,  this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_generator(const MultiIndex<int,DIM> j) const
        {
            typename Index::type_type e;
            typename Index::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i]=0;
                k[i]=bases_[i]->DeltaLmin();
            }
            return Index(j0_, e, k, 0, this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::last_generator() const
        {
            typename TensorIndex<IBASIS,DIM,TensorIndex<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
            for (unsigned int i = 0; i < DIM; i++)
            {
                e[i] = 0;
                k[i] = bases_[i]->DeltaRmax(j0_[i]);
            }
            int res=1;
            for (unsigned int i = 0; i < DIM; i++)
                res *= bases_[i]->Deltasize(j0_[i]);
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j0_, e, k, (res-1), this);
        }

#if _PRECOMPUTE_FIRSTLAST_WAVELETS
               
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_wavelet(const MultiIndex<int,DIM> j) const
        {
            MultiIndex<int,DIM> temp_mi(j);
            for (unsigned int i=0; i<DIM; ++i)
                temp_mi[i] -= j0_[i];
            return first_wavelets[temp_mi.number()];
        }
        
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_wavelet(const int levelsum) const
        {            
            MultiIndex<int,DIM> temp_mi;
            temp_mi[0] = levelsum - multi_degree(j0_);
            return first_wavelets[temp_mi.number()];
        }
        
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::last_wavelet(const MultiIndex<int,DIM> j) const
        {            
            MultiIndex<int,DIM> temp_mi(j);
            for (unsigned int i=0; i<DIM; ++i)
                temp_mi[i] -= j0_[i];
            return last_wavelets[temp_mi.number()];
        }
        
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::last_wavelet(const int levelsum) const
        {
            MultiIndex<int,DIM> temp_mi;
            temp_mi[0] = levelsum - multi_degree(j0_);
            return last_wavelets[temp_mi.number()];
        }
        
                
        template <class IBASIS, unsigned int DIM>
        void
        TensorBasis<IBASIS,DIM>::precompute_firstlast_wavelets()
        {
            MultiIndex<int,DIM> level_it;
            level_it[0] = jmax_ - multi_degree(j0_);
            const int numoflevels = level_it.number() + 1;
            level_it[0] = 0;
            first_wavelets.resize(numoflevels);
            last_wavelets.resize(numoflevels);
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::level_type temp_j(j0_);
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type temp_e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type temp_k;
            
            // first_level
            unsigned int l = 0;
            if (DIM == 1)
            {
                temp_e[0] = 1;
                temp_k[0] = bases_[0]->Nablamin();
            }
            else
            {
                for (unsigned int i = 0; i < DIM-1; i++) 
                {
                    temp_e[i] = 0;
                    temp_k[i] = bases_[i]->DeltaLmin();
                }
                temp_e[DIM-1] = 1;
                temp_k[DIM-1] = bases_[DIM-1]->Nablamin();
            }
            first_wavelets[0] = TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(temp_j, temp_e, temp_k, this);
            for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_e[i] = 1;
                    temp_k[i] = bases_[i]->Nablamax(temp_j[i]);
                }
            last_wavelets[0] = TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(temp_j, temp_e, temp_k, this);
            ++level_it;
            ++l;
            while (l < numoflevels)
            {
                //cout << "level_it == " << level_it << "; l == " << l <<  endl;
                for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_j[i] = j0_[i] + level_it[i];
                    if (level_it[i] == 0)
                    {
                        temp_e[i] = 0;
                        temp_k[i] = bases_[i]->DeltaLmin();
                    } else
                    {
                        temp_e[i] = 1;
                        temp_k[i] = bases_[i]->Nablamin();
                    }
                }
                //cout << "temp_j == " << temp_j << "; temp_e == " << temp_e << "; temp_k == " << temp_k << endl;
                first_wavelets[l] = TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(temp_j, temp_e, temp_k, this);
                for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_e[i] = 1;
                    temp_k[i] = bases_[i]->Nablamax(temp_j[i]);
                }
                //cout << "temp_j == " << temp_j << "; temp_e == " << temp_e << "; temp_k == " << temp_k << endl;
                last_wavelets[l] = TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(temp_j, temp_e, temp_k, this);
                ++level_it;
                ++l;
            }
        }
#else
        
        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_wavelet(const MultiIndex<int,DIM> j) const
        {
#if _TBASIS_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(j0_)
                    || (multi_degree(j) == multi_degree(j0_) && j >= j0_)
                  );
#endif
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
            bool first_level = true;
            for (unsigned int i = 0; i < DIM; i++) {
                if (j[i] == j0_[i])
                {
                    e[i] = 0;
                    k[i] = bases_[i]->DeltaLmin();
                } else
                {
                    e[i] = 1;
                    k[i] = bases_[i]->Nablamin();
                    first_level = false;
                }
            }
            if (first_level == true)
            {
                e[DIM-1] = 1;
                k[DIM-1] = bases_[DIM-1]->Nablamin();
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::first_wavelet(const int levelsum) const
        {            
#if _TBASIS_DEBUGLEVEL_ >= 1
            assert(levelsum >= multi_degree(j0_) );
#endif
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::level_type j;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;

            j[DIM-1] = j0_[DIM-1] + levelsum - multi_degree(j0_);
            e[DIM-1] = 1;
            k[DIM-1] = bases_[DIM-1]->Nablamin();

            for (unsigned int i = 0; i < DIM-1; i++)
            {
                j[i] = j0_[i];
                e[i] = 0;
                k[i] = bases_[i]->DeltaLmin();
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::last_wavelet(const MultiIndex<int,DIM> j) const
        {            
#if _TBASIS_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(j0_)
                    || ((multi_degree(j) == multi_degree(j0_)) && (j0_ <= j))
                  );
#endif
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i] = 1;
                k[i] = bases_[i]->Nablamax(j[i]);
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, this);
        }

        template <class IBASIS, unsigned int DIM>
        typename TensorBasis<IBASIS,DIM>::Index
        TensorBasis<IBASIS,DIM>::last_wavelet(const int levelsum) const
        {
#if _TBASIS_DEBUGLEVEL_ >= 1
            assert(levelsum >= multi_degree(j0_));
#endif
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::level_type j;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
            j[0]= j0_[0]+levelsum-multi_degree(j0_);
            e[0]=1;
            k[0]= bases_[0]->Nablamax(j[0]);
            for (unsigned int i = 1; i < DIM; i++) {
                j[i] = j0_[i];
                e[i] = 1;
                k[i] = bases_[i]->Nablamax(j0_[i]);
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, this);
        }
        
#endif // if _PRECOMPUTE_FIRSTLAST_WAVELETS
	template <class IBASIS, unsigned int DIM>
	void
	TensorBasis<IBASIS,DIM>::expand(const Function<DIM>* f,
                                        const bool primal,
                                        const MultiIndex<int,DIM> jmax,
                                        InfiniteVector<double,Index>& coeffs) const
  	{
            assert(primal == false); // only integrate against primal wavelets and generators
            for (unsigned int i=0; i<degrees_of_freedom(); i++)
            {
                const double coeff = integrate(f, full_collection[i]);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(full_collection[i], coeff);
            }
            /*
            for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
            {
                const double coeff = integrate(f, lambda);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(lambda, coeff);
                if (lambda == lambda_end)
                    break;
            }
            */
  	}

	template <class IBASIS, unsigned int DIM>
	void
	TensorBasis<IBASIS,DIM>::expand(const Function<DIM>* f,
                                        const bool primal,
                                        const unsigned int jmax,
                                        InfiniteVector<double,Index>& coeffs) const
  	{
            assert(primal == false); // only integrate against primal wavelets and generators
            /*
            for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
            {
                const double coeff = integrate(f, lambda);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(lambda, coeff);
                if (lambda == lambda_end)
                    break;
            }
            */
            for (unsigned int i=0; i<degrees_of_freedom(); i++)
            {
                const double coeff = integrate(f, full_collection[i]);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(full_collection[i], coeff);
            }
            
  	}
        
	template <class IBASIS, unsigned int DIM>
	double
	TensorBasis<IBASIS,DIM>::integrate(const Function<DIM>* f,
                                           const Index& lambda) const
  	{
            // f(v) = \int_0^1 g(t)v(t) dt
            double r = 0;
            // first compute supp(psi_lambda)
            Support supp;
            support(lambda, supp);
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = 5;
            FixedArray1D<double,DIM> h;
            for (unsigned int i=0; i<DIM;i++)
                h[i]=ldexp(1.0, -supp.j[i]); // granularity for the quadrature
                //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
            FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
            for (unsigned int i = 0; i < DIM; i++) {
                gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
                                = h[i]*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
                                = h[i]*GaussWeights[N_Gauss-1][n];
                    }
            }
            // compute the point values of the integrand (where we use that it is a tensor product)
            for (unsigned int i = 0; i < DIM; i++)
                bases()[i]->evaluate(0,
                                     lambda.j()[i],
                                     lambda.e()[i],
                                     lambda.k()[i],
                                     gauss_points[i],
                                     v_values[i]);
            // iterate over all points and sum up the integral shares
            int index[DIM]; // current multiindex for the point values
            for (unsigned int i = 0; i < DIM; i++)
                index[i] = 0;
            Point<DIM> x;
            while (true) {
                for (unsigned int i = 0; i < DIM; i++)
                    x[i] = gauss_points[i][index[i]];
                double share = f->value(x);
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
        }

        template <class IBASIS, unsigned int DIM>
  	double
  	TensorBasis<IBASIS,DIM>::evaluate(const unsigned int derivative, 
                                          const Index& lambda,
                                          const Point<DIM> x) const
  	{
            double value = 1.0;
            for (unsigned int i = 0; i < DIM; i++) // loop through components of the tensor product
                value *= bases_[i]->evaluate(derivative,
                                             lambda.j()[i], lambda.e()[i], lambda.k()[i],
                                             x[i]);
            return value;
        }


        template <class IBASIS, unsigned int DIM>
        void
        TensorBasis<IBASIS,DIM>::setup_full_collection()
        {
            if (jmax_ < multi_degree(j0_) ) {
                cout << "TensorBasis<IBASIS,DIM>::setup_full_collection(): the specified maximal level jmax is invalid. Specify a higher maximal level jmax_" << endl;
                cout << "jmax_ = " << jmax_ << "; j0_ = " << j0_ << endl;
                abort();
            }
            int degrees_of_freedom = last_wavelet_num<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(this, jmax_) +1; // +1 since numbering begins at 0
            cout << "total degrees of freedom between j0_ = " << j0_ << " and jmax_= " << jmax_ << " is " << degrees_of_freedom << endl;
            cout << "setting up collection of wavelet indices..." << endl;
            full_collection.resize(degrees_of_freedom);
            Index ind = first_generator();
            for (int k = 0; k < degrees_of_freedom; k++) {
                full_collection[k] = ind;
                ++ind;
            }
            /*
            int k=0;
            for (Index ind = first_generator(), itend = last_wavelet(jmax_); ind <= itend; ++ind) {
                full_collection[k] = ind;
                k++;
            }*/
            cout << "done setting up collection of wavelet indices..." << endl;
        }

}

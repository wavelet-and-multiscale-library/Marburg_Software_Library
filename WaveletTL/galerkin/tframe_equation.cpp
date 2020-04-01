// implementation for tframe_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>
#include <numerics/eigenvalues.h>
#include <cube/tframe_support.h>

namespace WaveletTL
{
//    //PERFORMANCE : wozu compute_rhs? welches jmax in den Konstruktoren? (wenn überhaupt)
//    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
//    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(EllipticBVP<DIM>* bvp,
//                                                           const FixedArray1D<bool,2*DIM>& bc,
//                                                           const bool precompute_rhs)
//    : bvp_(bvp), frame_(bc), normA(0.0), normAinv(0.0)
//    {
//        if (precompute_rhs == true)
//        {
//            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_.set_jpmax(multi_degree(frame_.j0()),0);
//            compute_rhs();
//        }
//    }

//    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
//    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const EllipticBVP<DIM>* bvp,
//                                                     const FixedArray1D<int,2*DIM>& bc,
//                                                     const bool precompute)
//    : bvp_(bvp), frame_(bc), normA(0.0), normAinv(0.0)
//    {
//        if (precompute == true)
//        {
//            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_.set_jpmax(multi_degree(frame_.j0()),0); // for a first quick hack
//            compute_rhs();
//        }
//    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
        TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const EllipticBVP<DIM>* bvp, const Frame* frame, const bool precompute_rhs)
    : nonzeroneumann_(false), bvp_(bvp), frame_(frame), normA(0.0), normAinv(0.0)
	{
#ifndef DYADIC
            compute_diagonal(); 
#endif
        if (precompute_rhs)
        {
//            cout << "Maximal level is set to "<<multi_degree(frame_->j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_->set_jpmax(multi_degree(frame_->j0()),0);
            compute_rhs();
        }
    }
    
    //constructor with given neumann-bc function
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
        TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const EllipticBVP<DIM>* bvp, const Frame* frame, const Function<DIM>* phi, const bool precompute_rhs)
    : nonzeroneumann_(true), bvp_(bvp), frame_(frame), phi_(phi), normA(0.0), normAinv(0.0)
	{
        cout<<"set non-zero neumann boundary conditions"<<endl;
#ifndef DYADIC
            compute_diagonal(); 
#endif
        if (precompute_rhs)
        {
//            cout << "Maximal level is set to "<<multi_degree(frame_->j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_->set_jpmax(multi_degree(frame_->j0()),0);
            compute_rhs();
        }
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::TensorFrameEquation(const TensorFrameEquation& eq)
    : bvp_(eq.bvp_), frame_(eq.frame_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
    {
#ifndef DYADIC
            compute_diagonal(); 
#endif
//            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_.set_jpmax(multi_degree(frame_.j0()),0);
    }

// TODO PERFORMANCE:: use setup_full_collection entries
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::compute_rhs()
    {
        cout << "TensorFrameEquation(): precompute right-hand side..." << endl;
        // precompute the right-hand side on a fine level
        InfiniteVector<double,Index> fhelp;
        InfiniteVector<double,int> fhelp_int;
        fnorm_sqr = 0;
        double coeff;
                
#if PARALLEL_RHS==1
        double fnorm_sqr_help = 0.;
        cout<<"parallel computing rhs"<<endl;
        
#pragma omp parallel 
{
//        cout<<"number of threads: "<<omp_get_num_threads()<<endl;
#pragma omp for private(coeff) schedule(static) reduction(+:fnorm_sqr_help)
        for (unsigned int i = 0; (int)i< frame_->degrees_of_freedom();i++)
        {
            coeff = f(frame_->get_quarklet(i)) / D(frame_->get_quarklet(i));
            if (fabs(coeff)>1e-15)
            {
#pragma omp critical
                {fhelp.set_coefficient(frame_->get_quarklet(i), coeff);
                fhelp_int.set_coefficient(i, coeff);}
                fnorm_sqr_help += coeff*coeff;
                
            }
        }
}
                fnorm_sqr=fnorm_sqr_help;
        
#else 
       for (int i = 0; i< frame_->degrees_of_freedom();i++)
        {
//           cout<<i<<": "<<*(frame_->get_quarklet(i))<<endl;
            coeff = f(*(frame_->get_quarklet(i))) / D(*(frame_->get_quarklet(i)));
            if (fabs(coeff)>1e-15)
            {
                fhelp.set_coefficient(*(frame_->get_quarklet(i)), coeff);
                fnorm_sqr += coeff*coeff;
                fhelp_int.set_coefficient(i, coeff);
            }
        } 
#endif
        cout << "... done, sort the entries in modulus..." << endl;
        // sort the coefficients into fcoeffs
        fcoeffs.resize(0); // clear eventual old values
        fcoeffs_int.resize(0);
        fcoeffs.resize(fhelp.size());
        fcoeffs_int.resize(fhelp_int.size());
        unsigned int id(0), id2(0);
        for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end());it != itend; ++it, ++id)
        {
            fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
//            cout << it.index() << ", " << *it << endl;
        }
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp_int.begin()), itend(fhelp_int.end()); it != itend; ++it, ++id2){
            fcoeffs_int[id2] = std::pair<int,double>(it.index(), *it);
//            cout << it.index() << ", " << *it << endl;
        }
        sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
        sort(fcoeffs_int.begin(), fcoeffs_int.end(), typename InfiniteVector<double,int>::decreasing_order());
        cout << "... done, all integrals for right-hand side computed!" << endl;
//        cout<<fcoeffs.size()<<endl;
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    inline
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::D(const Index& lambda) const
    {
#ifdef DYADIC
        double hspreconditioner(0), l2preconditioner(1);
        for (int i=0; i<space_dimension; i++){
            hspreconditioner+=pow(1+lambda.p()[i],8)*(1<<(2*lambda.j()[i]));
            l2preconditioner*=pow(1+lambda.p()[i],2);
        }
        double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
        return preconditioner;
#endif  
#ifdef TRIVIAL
        return 1;
#endif
#ifdef ENERGY
        return stiff_diagonal[lambda.number()];
#endif
#ifdef DYPLUSEN
        double hspreconditioner(0), l2preconditioner(1);
        for (int i=0; i<space_dimension; i++){
            hspreconditioner+=pow(1+lambda.p()[i],6);
            l2preconditioner*=pow(1+lambda.p()[i],2);
        }
        double preconditioner = sqrt(hspreconditioner)*l2preconditioner*stiff_diagonal[lambda.number()]/sqrt(space_dimension);
//        double preconditioner = sqrt(0.1)*stiff_diagonal[lambda.number()];
        return preconditioner;

#endif 
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    inline
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::a(const Index& lambda,
                                          const Index& nu) const
    {
        return a(lambda, nu, 2*IFRAME::primal_polynomial_degree());
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::a(const Index& la,
                                              const Index& nu,
                                              const unsigned int p) const
    {
        // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
        double r = 0;
//        Frame1D myframe;
        const Index* lambda = &la;
        const Index* mu     = &nu;
        // first decide whether the supports of psi_lambda and psi_mu intersect
        typedef typename Frame::Support Support;
        Support supp;
        if (intersect_supports(*frame_, *lambda, *mu, supp))
        {
//            Point<DIM> x;
//            const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0; //currently we do not support this
//            const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
            
            int N_Gauss [DIM];
            for (unsigned int i = 0; i < DIM; i++) {
               N_Gauss[i] = (p+1)/2+(lambda->p()[i]+mu->p()[i]+1)/2;
            }
            
            // loop over spatial direction
            for (unsigned int i = 0; i < DIM; i++) {
              double t = 1.;

              for (unsigned int j = 0; j < DIM; j++) {
                if (j == i)
                  continue;
//                IFRAME frame1d_la= frame_->frames(lambda->patch(),j);
//                IFRAME frame1d_mu= frame_->frames(mu->patch(),j);
                IndexQ1D<IFRAME> i1(IntervalQIndex<IFRAME> (
                                                          lambda->p()[j],lambda->j()[j],lambda->e()[j],lambda->k()[j],
                                                          frame_->frames()[j]
                                                          ),
                                   0,j,0 //patch not relevant
                                   );
                IndexQ1D<IFRAME> i2(IntervalQIndex<IFRAME> (mu->p()[j],mu->j()[j],mu->e()[j],mu->k()[j],
                                                          frame_->frames()[j]
                                                          ),
                                   0,j,0
                                   );


                t *= integrate(i1, i2, N_Gauss[j], j, supp);
//                cout << "Zwischenergebnis Richtung: " << j << "Wert: " << t << endl;
              }
              
//              IFRAME frame1d_la= frame_->frames(lambda->patch(),i);
//              IFRAME frame1d_mu= frame_->frames(mu->patch(),i);
              IndexQ1D<IFRAME> i1(IntervalQIndex<IFRAME> (
                                                        lambda->p()[i],lambda->j()[i],lambda->e()[i],lambda->k()[i],
                                                        frame_->frames()[i]
                                                        ),
                                 0,i,1 //patch not relevant
                                 );
              IndexQ1D<IFRAME> i2(IntervalQIndex<IFRAME> (mu->p()[i],mu->j()[i],mu->e()[i],mu->k()[i],
                                                        frame_->frames()[i]
                                                        ),
                                 0,i,1
                                 );

              t *= integrate(i1, i2, N_Gauss[i], i, supp);

              r += t;
            }
             
            
//            double integral[space_dimension], der_integral[space_dimension];
//            integral[0]=0, integral[1]=0, der_integral[0]=0, der_integral[1]=0;
//            FixedArray1D<Array1D<double>,DIM> psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
//                                              psi_mu_values,         // -"-, for psi_mu
//                                              psi_lambda_der_values, // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]
//                                              psi_mu_der_values;     // -"-, for psi_mu
//            
//
//            const int N_Gauss = (p+1)/2+(multi_degree(lambda->p())+multi_degree(mu->p())+1)/2;
//
//            FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights;
//
//            for (unsigned int i = 0; i < DIM; i++) {
//
//                
//                hi = ldexp(1.0, -supp.j[i]);
//                gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//                gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
//                for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
//                    for (int n = 0; n < N_Gauss; n++) {
//                        gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
//                                = hi*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
//                        gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
//                                = hi*GaussWeights[N_Gauss-1][n];
//                    }
//            
//            
//            
//                evaluate(*(frame_->frames()[i]), 
//                                                lambda->p()[i],
//                                                lambda->j()[i],
//                                                lambda->e()[i],
//                                                lambda->k()[i],
//                         gauss_points[i], psi_lambda_values[i], psi_lambda_der_values[i]);
//
//                evaluate(*(frame_->frames()[i]),
//                                                mu->p()[i],
//                                                mu->j()[i],
//                                                mu->e()[i],
//                                                mu->k()[i],
//                         gauss_points[i], psi_mu_values[i], psi_mu_der_values[i]);
//
//                    for (unsigned int ind = 0; ind < gauss_points[i].size(); ind++){
//
//                        integral[i] += psi_lambda_values[i][ind] * psi_mu_values[i][ind] * gauss_weights[i][ind];
//                        der_integral[i] += psi_lambda_der_values[i][ind] * psi_mu_der_values[i][ind] * gauss_weights[i][ind];
//                    }
//                
//            }
//                 
//                r = ax * (der_integral[0] * integral[1] + integral[0] * der_integral[1]) + qx * (integral[0] * integral[1]); 

        }
        return r;
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::f(const Index& lambda) const
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
        typedef typename Frame::Support Support;
        Support supp;
        support(*frame_, lambda, supp);
        // setup Gauss points and weights for a composite quadrature formula:
        const int N_Gauss = 6;
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
            evaluate(*(frame_->frames()[i]), 0,
                                            lambda.p()[i],
                                            lambda.j()[i],
                                            lambda.e()[i],
                                            lambda.k()[i],
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
//        cout<<"r="<<r<<endl;
//#ifdef NONZERONEUMANN
    if(nonzeroneumann_){
        //add  boundary integral
        double rrand=0;
//        cout<<"bin hier"<<endl;
        Array1D<double> phi_values;
        if(DIM==2){
            // loop over spatial direction
            for (int i = 0; i < 2; i++) {
              for (int j = 0; j < 2; j++) {
                if (j == i)
                  continue;
                
                phi_values.resize(N_Gauss*(supp.b[j]-supp.a[j]));
                for(int rand_punkt=0;rand_punkt<2;rand_punkt++){
                    double rand_value=frame_->frames()[i]->evaluate(0,lambda.p()[i],lambda.j()[i],lambda.e()[i],lambda.k()[i],rand_punkt);
//                    if((j==1 && rand_punkt==0) || (j==0 && rand_punkt==1)){
//                        
//                    }                   
                    for(int i2=0;i2<N_Gauss*(supp.b[j]-supp.a[j]);i2++){
                        x[i]=rand_punkt; x[j]=gauss_points[j][i2];
//                        cout<<x<<endl;
                        phi_values[i2]=phi(x);
                        rrand+=phi_values[i2]*rand_value*gauss_weights[j][i2]*v_values[j][i2]; 
                    }
//                    if(rand_value!=0){
//                        cout<<"weights:"<<gauss_weights[j]<<endl;
//                    cout<<"points:"<<gauss_points[j]<<endl;
//                    cout<<"v_values:"<<v_values[j]<<endl;
//                    cout<<"phi_values:"<<phi_values<<endl;
//                    }
                }
              }
//              cout<<"rrand="<<rrand<<endl;
            }
            
            
        } 
        r+=rrand;
    }
//#endif
//        cout<<"rrand="<<rrand<<endl;
#ifdef DELTADIS
      //r+= 4*frame_->frames()[0]->evaluate(0,lambda.p()[0],lambda.j()[0],lambda.e()[0],lambda.k()[0],0.5);  
      for(int i2=0;i2<N_Gauss*(supp.b[1]-supp.a[1]);i2++){
                        r+=4*frame_->frames()[0]->evaluate(0,lambda.p()[0],lambda.j()[0],lambda.e()[0],lambda.k()[0],0.5)*gauss_weights[1][i2]*v_values[1][i2]; 
                    }
#endif
        return r;
#endif
    }

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::RHS(const double eta,
                                                InfiniteVector<double, Index>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(fnorm_sqr - eta*eta);
//        typedef typename Frame::Index Index;
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
        typename Array1D<std::pair<int, double> >::const_iterator it(fcoeffs_int.begin());
        while (it != fcoeffs_int.end() && coarsenorm < bound)
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
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
    
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
        void
        TensorFrameEquation<IFRAME, DIM, TENSORFRAME>::compute_diagonal()
        {
          cout << " TensorFrameEquation(): precompute diagonal of stiffness matrix..." << endl;

          SparseMatrix<double> diag(1,frame_->degrees_of_freedom());

//          char filename[50];
//          char matrixname[50];
//      #ifdef ONE_D
//          int d = IFRAME::primal_polynomial_degree();
//          int dT = IFRAME::primal_vanishing_moments();
//      #else
//      #ifdef TWO_D
//          int d = IFRAME::primal_polynomial_degree();
//          int dT = IFRAME::primal_vanishing_moments();
//      #endif
//      #endif

          // prepare filenames for 1D and 2D case
      #ifdef ONE_D
          sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_interval_lap07_d", d, "_dT", dT);
          sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_1D_lap07_d", d, "_dT", dT);
      #endif
      #ifdef TWO_D
          //sprintf(filename, "%s%d%s%d", "stiff_diagonal_poisson_lshaped_lap1_d", d, "_dT", dT);
          //sprintf(matrixname, "%s%d%s%d", "stiff_diagonal_poisson_2D_lap1_d", d, "_dT", dT);
      #endif
      #ifndef PRECOMP_DIAG
          std::list<Vector<double>::size_type> indices;
          std::list<double> entries;
      #endif
      #ifdef PRECOMP_DIAG
          cout << "reading in diagonal of unpreconditioned stiffness matrix from file "
               << filename << "..." << endl;
          diag.matlab_input(filename);
          cout << "...ready" << endl;
      #endif
          stiff_diagonal.resize(frame_->degrees_of_freedom());
          for (int i = 0; i < frame_->degrees_of_freedom(); i++) {
      #ifdef PRECOMP_DIAG
            stiff_diagonal[i] = diag.get_entry(0,i);
      #endif 
      #ifndef PRECOMP_DIAG
      #ifdef FRAME
            stiff_diagonal[i] = sqrt(a(*(frame_->get_quarklet(i)),*(frame_->get_quarklet(i))));
      #endif
      #ifdef BASIS
            stiff_diagonal[i] = sqrt(a(*(frame_->get_wavelet(i)),*(frame_->get_wavelet(i))));
      #endif

            indices.push_back(i);
            entries.push_back(stiff_diagonal[i]);
      #endif
            //cout << stiff_diagonal[i] << " " << *(frame_->get_wavelet(i)) << endl;
          }
      #ifndef PRECOMP_DIAG
          diag.set_row(0,indices, entries);
          //diag.matlab_output(filename, matrixname, 1);
      #endif

          cout << "... done, diagonal of stiffness matrix computed" << endl;
        }

    // PERFORMANCE :: use setup_full_collection
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME,DIM,TENSORFRAME>::norm_A() const
    {
        if (normA == 0.0) {
            typedef typename Index::polynomial_type polynomial_type;
            int offsetj, offsetp;
            switch (space_dimension) {
                case 1:
                    offsetj = 2, offsetp = std::min((int)frame_->get_pmax(),2);
                    break;
                case 2:
//                    offsetj = 1, 
                    offsetj=std::min(1,(int)frame_->get_jmax()-(int)multi_degree(frame_->j0()));
                    offsetp = std::min((int)frame_->get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
            std::set<Index> Lambda;
            
            // cout << "tframe_equation.norm_A :: last quarklet = " << (frame_.last_quarklet(multi_degree(frame_.j0())+offset)) << endl;
            polynomial_type p;
            for (Index lambda = frame_->first_generator(frame_->j0(), p) ;;) {
                Lambda.insert(lambda);
                if (lambda == frame_->last_quarklet(multi_degree(frame_->j0())+offsetj, p)){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = frame_->first_generator(frame_->j0(), p);
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
            A_Lambda.compress(1e-10);
#if 1
//            double help;
//            unsigned int iterations;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//            normAinv = 1./help;
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
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
            typedef typename Frame::Index Index;
            typedef typename Index::polynomial_type polynomial_type;
            std::set<Index> Lambda;
            polynomial_type p;
            //const int j0 = frame().j0();
            //const int jmax = j0+3;
            int offsetj, offsetp;
            switch (space_dimension) {
                case 1:
                    offsetj = 2, offsetp = std::min((int)frame_->get_pmax(),2);
                    break;
                case 2:
//                    offsetj = 1, 
                      offsetj=std::min(1,(int)frame_->get_jmax()-(int)multi_degree(frame_->j0()));
                      offsetp = std::min((int)frame_->get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
            for (Index lambda = frame_->first_generator(frame_->j0()) ;;) {
                Lambda.insert(lambda);
                
                if (lambda == frame_->last_quarklet(multi_degree(frame_->j0())+offsetj, p) ){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = frame_->first_generator(frame_->j0(), p);
                }
                else
                    ++lambda;
                    
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
            A_Lambda.compress(1e-10);
#if 1
//            double help;
//            unsigned int iterations;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//            normAinv = 1./help;
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            cout << evals << endl;
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
#else
            Vector<double> xk(Lambda.size(), false);
            xk = 1;
            unsigned int iterations;
            normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
        }
        return normAinv;
    }
    
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    double
    TensorFrameEquation<IFRAME, DIM, TENSORFRAME>:: integrate(const IndexQ1D<IFRAME>& la,
                                                            const IndexQ1D<IFRAME>& nu,
                                                            const int N_Gauss,
                                                            const int dir,
                                                            typename Frame::Support supp) const
      {
        double res=0.;
//        cout<<"testing integrate"<<endl;
//        cout<<"direction="<<dir<<endl;
         // It makes sense to store the one dimensional
         // integrals arising when we make use of the tensor product structure. This costs quite
         // some memory, but really speeds up the algorithm!
        
        
        const IndexQ1D<IFRAME>* lambda= nu < la? &la : &nu;
        const IndexQ1D<IFRAME>* mu= nu < la ? &nu : &la;
//        if(lambda->number()==lambda->index().number())
//        cout << "truueee: " << lambda->index() << endl;
//        if(lambda->number()!=lambda->index().number())
//        cout << "faaaalsse: " << lambda->index() << ", " << lambda->number() << endl;
//        cout << mu->number() << endl;
        typename One_D_IntegralCache::iterator col_lb(one_d_integrals.lower_bound(*lambda));
        typename One_D_IntegralCache::iterator col_it(col_lb);
        if (col_lb == one_d_integrals.end() ||
            one_d_integrals.key_comp()(*lambda,col_lb->first))
          {
            // insert a new column
            typedef typename One_D_IntegralCache::value_type value_type;
            col_it = one_d_integrals.insert(col_lb, value_type(*lambda, Column1D()));
          }

        Column1D& col(col_it->second);

        typename Column1D::iterator lb(col.lower_bound(*mu));
        typename Column1D::iterator it(lb);
        if (lb == col.end() ||
            col.key_comp()(*mu, lb->first))
          {
    
            double h; // granularity for the quadrature
            Array1D<double> gauss_points, gauss_weights;
            
            //compute gauss points and weights
            double a =supp.a[dir];
            double b =supp.b[dir];
//            cout << "current support: ["<<a<<","<<b<<"]"<<endl;
            
            
//                int e = std::max(lambda->e()[i], mu->e()[i]); //correct granularity
            h = ldexp(1.0, -supp.j[dir]/*-e*/);
            gauss_points.resize(N_Gauss*(b-a));
            gauss_weights.resize(N_Gauss*(b-a));
            for (int interval = a; interval < b; interval++){
                for (int n = 0; n < N_Gauss; n++) {
                    gauss_points[(interval-a)*N_Gauss+n]= h*(2*interval+1+GaussPoints[N_Gauss-1][n])/2.;
                    gauss_weights[(interval-a)*N_Gauss+n]= h*GaussWeights[N_Gauss-1][n];
                }
            }
            
            
            Array1D<double> lambda_gauss_points(gauss_points), mu_gauss_points(gauss_points);
            Array1D<double> lambda_values,     // values of the components of psi_lambda at gauss_points[i]
                            mu_values;         // -"-, for psi_mu
            
            const IFRAME* frame1D_lambda = frame_->frames()[dir];
            const IFRAME* frame1D_mu     = frame_->frames()[dir];
            
            WaveletTL::evaluate(*frame1D_lambda, lambda->derivative(),  lambda->index(),         
                                                /*lambda->index().p(),
                                                lambda->index().j(),
                                                lambda->index().e(),
                                                lambda->index().k(),*/
                         lambda_gauss_points, lambda_values);
            
            WaveletTL::evaluate(*frame1D_mu, mu->derivative(), mu->index(),          
                                                /*mu->index().p(),
                                                mu->index().j(),
                                                mu->index().e(),
                                                mu->index().k(),*/
                         mu_gauss_points, mu_values);
            
            // - add all integral shares
            for (unsigned int ind = 0; ind < gauss_points.size(); ind++){
                   res += lambda_values[ind] * mu_values[ind] * gauss_weights[ind];
//                    der_integral += lambda_der_values[ind] * mu_der_values[ind] * gauss_weights[ind];
            }
            
            //store the calculated values
            typedef typename Column1D::value_type value_type;
//            cout << "Neuberechnung" << endl;
            it = col.insert(lb, value_type(*mu, res));
          }
        else {
//            cout << "aus dem Cache" << endl;
          res = it->second;
        }
    
//        (lambda->derivative()==0? return integral : return der_integral);
//        cout<<res<<endl;
        return res;
      }
    
}

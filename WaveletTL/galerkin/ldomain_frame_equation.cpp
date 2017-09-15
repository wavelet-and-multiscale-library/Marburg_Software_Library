// implementation for ldomain_frame_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>
#include <numerics/eigenvalues.h>
#include <Ldomain/ldomain_frame_support.h>


namespace WaveletTL
{
//    template <class IFRAME, class LDOMAINFRAME>
//    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp,
//					    const bool precompute_rhs)
//    : bvp_(bvp), frame_(), normA(0.0), normAinv(0.0)
//    {
//        cout << "hier2" << endl;
//        
//        if (precompute_rhs)
//        {
//            cout << "Maximal level is set to "<<multi_degree(frame_->j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_->set_jpmax(multi_degree(frame_->j0()),0);
//#ifndef DYADIC
//            compute_diagonal(); 
//#endif
//            compute_rhs();
//        }
//    }
    
    template <class IFRAME, class LDOMAINFRAME>
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp, const Frame* frame,
					    const bool precompute_rhs)
    : bvp_(bvp), frame_(frame), normA(0.0), normAinv(0.0)
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
  
    //template <class IFRAME, class LDOMAINFRAME>
    //LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp,
	//				   const FixedArray1D<bool,8>& bc,
	//				   const bool precompute_rhs)
    //: bvp_(bvp), frame_(bc), normA(0.0), normAinv(0.0)
    //{
    //    if (precompute_rhs)
    //    {
    //        cout << "Maximal level is set to "<<multi_degree(frame_->j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
    //        frame_->set_jpmax(multi_degree(frame_->j0()),0);
    //        compute_rhs();
    //   }
    //}
  
    template <class IFRAME, class LDOMAINFRAME>
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const LDomainFrameEquation& eq)
    : bvp_(eq.bvp_), frame_(eq.frame_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
    {
#ifndef DYADIC
            compute_diagonal(); 
#endif
//        cout << "Maximal level is set to "<<multi_degree(frame_->j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
//            frame_->set_jpmax(multi_degree(frame_->j0()),0);
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    void
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::compute_rhs()
    {
        cout << "LDomainFrameEquation(): precompute right-hand side..." << endl;
        // precompute the right-hand side on a fine level
        InfiniteVector<double,Index> fhelp;
        InfiniteVector<double,int> fhelp_int;
        fnorm_sqr = 0; 
        double coeff;
        
#if PARALLEL==1
        double fnorm_sqr_help = 0.;
        cout<<"parallel computing rhs"<<endl;
        
#pragma omp parallel 
{
//        cout<<"number of threads: "<<omp_get_num_threads()<<endl;
#pragma omp for  private(coeff) schedule(static) reduction(+:fnorm_sqr_help)
        for (int i = 0; i< frame_->degrees_of_freedom();i++)
        {
            coeff = f(*(frame_->get_quarklet(i))) / D(*(frame_->get_quarklet(i)));
            if (fabs(coeff)>1e-15)
            {
#pragma omp critical
                {fhelp.set_coefficient(*(frame_->get_quarklet(i)), coeff);}
                
                fnorm_sqr_help += coeff*coeff;
            }
        }
}
        fnorm_sqr=fnorm_sqr_help;
#else 
       for (int i = 0; i< frame_->degrees_of_freedom();i++)
        {
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
//        cout<<"number of threads: "<<omp_get_num_threads()<<endl;
        // sort the coefficients into fcoeffs
        fcoeffs.resize(0); // clear eventual old values
        fcoeffs_int.resize(0);
        fcoeffs.resize(fhelp.size());
        fcoeffs_int.resize(fhelp_int.size());
        unsigned int id(0), id2(0);
        for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end()); it != itend; ++it, ++id){
            fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
        }
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp_int.begin()), itend(fhelp_int.end()); it != itend; ++it, ++id2){
            fcoeffs_int[id2] = std::pair<int,double>(it.index(), *it);
        }
        sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
        sort(fcoeffs_int.begin(), fcoeffs_int.end(), typename InfiniteVector<double,int>::decreasing_order());
        cout << "... done, all integrals for right-hand side computed!" << endl;
    }

    template <class IFRAME, class LDOMAINFRAME>
    inline
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::D(const Index& lambda) const
    {
        //return ldexp(1.0, lambda.j());


        
#ifdef DYADIC
        double hspreconditioner(0), l2preconditioner(1);
        for (int i=0; i<space_dimension; i++){
            hspreconditioner+=pow(1+lambda.p()[i],8)*(1<<(2*lambda.j()[i]));
            l2preconditioner*=pow(1+lambda.p()[i],2);
        }
        double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
        return preconditioner;
#else  
#ifdef TRIVIAL
        return 1;
#else
        return stiff_diagonal[lambda.number()];
//        return sqrt(a(lambda, lambda));
#endif 
#endif
//        return 1;
    }

    template <class IFRAME, class LDOMAINFRAME>
    inline
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::a(const Index& lambda,
			     const Index& nu) const
    {
        return a(lambda, nu, IFRAME::primal_polynomial_degree()*IFRAME::primal_polynomial_degree());
    }

    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::a(const Index& la,
			     const Index& nu,
			     const unsigned int p) const
    {
        
        // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
        double r = 0.;
        const Index* lambda = &la;
        const Index* mu     = &nu;
        typedef typename Frame::Support Support;
        Support supp;
        if (intersect_supports(*frame_, *lambda, *mu, supp)) {
            //TODO: mit Fortsetzungen umgehen, integrate routine
             int N_Gauss [2];
             N_Gauss[0] = (p+1)/2+(lambda->p()[0]+mu->p()[0]+1)/2;
             N_Gauss[1] = (p+1)/2+(lambda->p()[1]+mu->p()[1]+1)/2;
             
        
        // loop over spatial direction
            for (int i = 0; i < 2; i++) {
              double t = 1.;

              for (int j = 0; j < 2; j++) {
                if (j == i)
                  continue;
//                IFRAME frame1d_la= frame_->frames(lambda->patch(),j);
//                IFRAME frame1d_mu= frame_->frames(mu->patch(),j);
                IndexQ1D<IFRAME> i1(IntervalQIndex<IFRAME> (
                                                          lambda->p()[j],lambda->j()[j],lambda->e()[j],lambda->k()[j],
                                                          frame_->frames(lambda->patch(),j)
                                                          ),
                                   lambda->patch(),j,0
                                   );
                IndexQ1D<IFRAME> i2(IntervalQIndex<IFRAME> (mu->p()[j],mu->j()[j],mu->e()[j],mu->k()[j],
                                                          frame_->frames(mu->patch(),j)
                                                          ),
                                   mu->patch(),j,0
                                   );


                t *= integrate(i1, i2, N_Gauss[j], j, supp);
//                cout << "Zwischenergebnis Richtung: " << j << "Wert: " << t << endl;
              }

//              IFRAME frame1d_la= frame_->frames(lambda->patch(),i);
//              IFRAME frame1d_mu= frame_->frames(mu->patch(),i);
              IndexQ1D<IFRAME> i1(IntervalQIndex<IFRAME> (
                                                        lambda->p()[i],lambda->j()[i],lambda->e()[i],lambda->k()[i],
                                                        frame_->frames(lambda->patch(),i)
                                                        ),
                                 lambda->patch(),i,1
                                 );
              IndexQ1D<IFRAME> i2(IntervalQIndex<IFRAME> (mu->p()[i],mu->j()[i],mu->e()[i],mu->k()[i],
                                                        frame_->frames(mu->patch(),i)
                                                        ),
                                 mu->patch(),i,1
                                 );

              t *= integrate(i1, i2, N_Gauss[i], i, supp);

              r += t;
            }
        }
        return r;
    }
//            
    
    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::f(const Index& lambda) const
    {
        // f(v) = \int_0^1 g(t)v(t) dt
//        cout << "LDomainEquation::f() called with lambda=" << lambda << endl;
        double r = 0;

        typedef typename Frame::Support Support;
        Support supp;
        support(*frame_,lambda,supp);
//        cout << "patch 0: [" << supp.xmin[0]  <<" , " <<supp.xmax[0] <<"]x["<<supp.ymin[0] <<" , "<<supp.ymax[0] <<"]"<<  endl;
//        cout << "patch 1: [" << supp.xmin[1]  <<" , " <<supp.xmax[1] <<"]x["<<supp.ymin[1] <<" , "<<supp.ymax[1] <<"]"<<  endl;
//        cout << "patch 2: [" << supp.xmin[2]  <<" , " <<supp.xmax[2] <<"]x["<<supp.ymin[2] <<" , "<<supp.ymax[2] <<"]"<<  endl;
        
        //determine on which patch lambda lies

//        const bool extended=(lambda->patch()>2) ? 1 : 0;

        const int mother_patch=(lambda.patch()==3) ? 0: (lambda.patch()==4 ? 2: lambda.patch());
        //cout << "extended: "<<extended<<endl;
        //cout << "mother_patch: "<<mother_patch<<endl;
    
        //setup gauss points
        //const int N_Gauss = (p+1)/2+(multi_degree(lambda.p()))/2;
        const int N_Gauss=6;
        double h; // granularity for the quadrature
        FixedArray1D<Array1D<double>,space_dimension> gauss_points, gauss_weights, v_values;
	
 	// per patch, collect all point values
 	for (int patch = 0; patch <= 2; patch++) {
            //cout << "current patch: "<<patch<<endl;
 	    if (supp.xmin[patch] != -1) { // psi_mu is nontrivial on patch p
                double a [2]; a[0]=supp.xmin[patch]; a[1]=supp.ymin[patch];
                double b [2]; b[0]=supp.xmax[patch]; b[1]=supp.ymax[patch];
                //cout << "current support: ["<<a[0]<<","<<b[0]<<"]x["<<a[1]<<","<<b[1]<<"]"<<endl;
                for (int i = 0; i < space_dimension; i++) {
                    // prepare Gauss points and weights 
                    h = ldexp(1.0,-supp.j[i]);
                    //cout << "current interval :" <<a[i]*h <<","<<b[i]*h<<endl;
                    gauss_points[i].resize(N_Gauss*(b[i]-a[i]));
                    gauss_weights[i].resize(N_Gauss*(b[i]-a[i]));
                    for (int interval = a[i]; interval < b[i]; interval++){
                        for (int n = 0; n < N_Gauss; n++) {
                            gauss_points[i][(interval-a[i])*N_Gauss+n] = h*(2*interval+1+GaussPoints[N_Gauss-1][n])/2.;
                            gauss_weights[i][(interval-a[i])*N_Gauss+n]= h*GaussWeights[N_Gauss-1][n];
                        }
                    }
                }
                
                
                
                //cout << "gauss points before reflection"<<endl;
                //cout << gauss_points[0]<<endl;
                //cout << gauss_points[1]<<endl;
                
                //special treat for patch 1 if the quarklet is extended:    
                FixedArray1D<Array1D<double>,space_dimension> lambda_gauss_points(gauss_points);
                if(patch==0){   // compute the point values of the integrand (where we use that it is a tensor product)
                    evaluate(*(frame_->frame1d_11()), 0,        //different 1d-frames for x- and y-direction due to boundary conditions
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                         lambda_gauss_points[0], v_values[0]);  
                    evaluate(*(frame_->frame1d_01()), 0,
                                            lambda.p()[1],
                                            lambda.j()[1],
                                            lambda.e()[1],
                                            lambda.k()[1],
                         lambda_gauss_points[1], v_values[1]); 
                }
                
                //reflect gauss points
                if(patch==1){
                    if(mother_patch==0){ //lambda has been extended from north so south, reflect gauss points in y-direction
                        for(int i=0;i<(int)lambda_gauss_points[1].size();i++){
                            lambda_gauss_points[1][i]=1-lambda_gauss_points[1][i];
                        }
                        evaluate(*(frame_->frame1d_11()), 0,        //different 1d-frames for x- and y-direction due to boundary conditions
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(*(frame_->frame1d_01()), 0,
                                            lambda.p()[1],
                                            lambda.j()[1],
                                            lambda.e()[1],
                                            lambda.k()[1],
                                        lambda_gauss_points[1], v_values[1]); 
                        for(int i=0;i<=(int)v_values[1].size()/2;i++){    //switch direction of values
                            //v_values[1].swap(i,v_values[1].size()-1-i);
                        }
                                
                    }
                    if(mother_patch==1){                //psi_lambda has not been extended, evaluate with zero bc
                        evaluate(*(frame_->frame1d_11()), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(*(frame_->frame1d_11()), 0,
                                            lambda.p()[1],
                                            lambda.j()[1],
                                            lambda.e()[1],
                                            lambda.k()[1],
                                        lambda_gauss_points[1], v_values[1]);
                    }
                    if(mother_patch==2){ //psi_lambda has been extended from east to west, reflect gauss points in x-direction
                        for(int i=0;i<(int)lambda_gauss_points[0].size();i++){
                            lambda_gauss_points[0][i]=1-lambda_gauss_points[0][i];
                        }
                        evaluate(*(frame_->frame1d_01()), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(*(frame_->frame1d_11()), 0,
                                            lambda.p()[1],
                                            lambda.j()[1],
                                            lambda.e()[1],
                                            lambda.k()[1],
                                        lambda_gauss_points[1], v_values[1]);
                        for(int i=0;i<=(int)v_values[0].size()/2;i++){    //switch direction of values
                            //v_values[0].swap(i,v_values[0].size()-1-i);
                        }
                    }
                } 
            
                if(patch==2){
                    evaluate(*(frame_->frame1d_01()), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                         lambda_gauss_points[0], v_values[0]);  
                    evaluate(*(frame_->frame1d_11()), 0,
                                            lambda.p()[1],
                                            lambda.j()[1],
                                            lambda.e()[1],
                                            lambda.k()[1],
                         lambda_gauss_points[1], v_values[1]); 
                }
                /*
                cout << "gauss points used for evaluation"<<endl;
                cout << lambda_gauss_points[0]<<endl;
                cout << lambda_gauss_points[1]<<endl;
                cout << "values:"<<endl;
                cout << v_values[0]<<endl;
                cout << v_values[1]<<endl;
                 */
                // iterate over all points and sum up the integral shares
                int index[space_dimension]; // current multiindex for the point values
                for (int i = 0; i < space_dimension; i++)
                    index[i] = 0;
                Point<space_dimension> x;
                while (true) {
                    //for (unsigned int i = 0; i < space_dimension; i++)
                    switch(patch) {
                        case 0:
                            x[0] = gauss_points[0][index[0]]-1;
                            x[1] = gauss_points[1][index[1]];
                            break;
                        case 1:
                            x[0] = gauss_points[0][index[0]]-1;
                            x[1] = gauss_points[1][index[1]]-1;
                            break;
                        case 2:
                            x[0] = gauss_points[0][index[0]];
                            x[1] = gauss_points[1][index[1]]-1;
                            break;
                        default:
                            x[0] = gauss_points[0][index[0]];
                            x[1] = gauss_points[1][index[1]];
                    }
                        
                    double share = bvp_->f(x);
                    for (int i = 0; i < space_dimension; i++)
                        share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    for (int i = 0; i < space_dimension; i++) {
                        if (index[i] == N_Gauss*(b[i]-a[i])-1) {
                            index[i] = 0;
                            exit = (i == space_dimension-1);
                        } else {
                            index[i]++;
                            break;
                        }
                    }
                    if (exit) break;
                }
        //return r;  
            }
	}    
        //cout << "LDomainEquation::f() result is r=" << r << endl;
	return r;  
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    void
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::RHS(const double eta,
                                    InfiniteVector<double,Index>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(fnorm_sqr - eta*eta);
        typename Array1D<std::pair<Index, double> >::const_iterator it(fcoeffs.begin());
        while (it != fcoeffs.end() && coarsenorm < bound){
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
            ++it;
        } 
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    void
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::RHS(const double eta,
                                    InfiniteVector<double,int>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(fnorm_sqr - eta*eta);
        typename Array1D<std::pair<int, double> >::const_iterator it(fcoeffs_int.begin());
        while (it != fcoeffs_int.end() && coarsenorm < bound){
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
            ++it;
        } 
    }
        

    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::norm_A() const
    {
        if (normA == 0.0) {
        //typedef typename Frame::Index Index;
        typedef typename Index::polynomial_type polynomial_type;
        int offsetj=std::min(1,frame_->get_jmax()-(int)multi_degree(frame_->j0()));
        int offsetp=std::min((int)frame_->get_pmax(),2);
        std::set<Index> Lambda;
      
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
//        cout << "Lambda size: "<<Lambda.size()<<endl;
//        typename  std::set<Index>::const_iterator itend(Lambda.end());
//        for(typename  std::set<Index>::const_iterator it=Lambda.begin();it!=itend;it++)
//            cout << *it << endl;
        
        SparseMatrix<double> A_Lambda;
        setup_stiffness_matrix(*this, Lambda, A_Lambda);
//        cout << A_Lambda << endl;
        A_Lambda.compress(1e-10);
#if 1
        //double help;
        //unsigned int iterations;
        //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
        //normAinv = 1./help;
//        cout << normA << ", " << normAinv << endl;
        Matrix<double> evecs;
        Vector<double> evals;
        SymmEigenvalues(A_Lambda, evals, evecs);
        int i = 0;
        while(abs(evals(i))<1e-2){
            ++i;
        }
        normA = evals(evals.size()-1);
        normAinv = 1./evals(i);
      
        cout << "Eigenwerte: " << evals << endl;
      //cout << "Eigenvektoren: " << endl << evecs << endl;
#else
        Vector<double> xk(Lambda.size(), false);
        xk = 1;
        unsigned int iterations;
        normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
        }
        return normA;
    }
   
    
    
    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::norm_Ainv() const
    {
        if (normAinv == 0.0) {
        typedef typename Frame::Index Index;
        typedef typename Index::polynomial_type polynomial_type;
        int offsetj=std::min(1,frame_->get_jmax()-(int)multi_degree(frame_->j0()));
        int offsetp=std::min((int)frame_->get_pmax(),2);
        std::set<Index> Lambda;
      
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
//        cout << "INVERS!!!" << endl;
//        cout << "Lambda size: "<<Lambda.size()<<endl;
//        
//        cout << "Lambda size: "<<Lambda.size()<<endl;
//        typename  std::set<Index>::const_iterator itend(Lambda.end());
//        for(typename  std::set<Index>::const_iterator it=Lambda.begin();it!=itend;it++)
//            cout << *it << endl;
        
        SparseMatrix<double> A_Lambda;
        setup_stiffness_matrix(*this, Lambda, A_Lambda);
        A_Lambda.compress(1e-10);
#if 1
        //double help;
        //unsigned int iterations;
        //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
        //normAinv = 1./help;
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
        normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
        }
        return normAinv;
    }
  

    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::s_star() const
    {
        // notation from [St04a]
        const double t = operator_order();
        const int n = 2;
        const int dT = Frame::primal_vanishing_moments();
        const double gamma = Frame::primal_regularity();
    
        return std::min((t+dT)/(double)n, (gamma-t)/(n-1.)); // [St04a, Th. 2.3]
    }


    template <class IFRAME, class LDOMAINFRAME>
        void
        LDomainFrameEquation<IFRAME, LDOMAINFRAME>::compute_diagonal()
        {
          cout << " LDomainFrameEquation(): precompute diagonal of stiffness matrix..." << endl;

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
    
#if 1
    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>:: integrate(const IndexQ1D<IFRAME>& la,
                                                            const IndexQ1D<IFRAME>& nu,
                                                            const int N_Gauss,
                                                            const int dir,
                                                            typename Frame::Support supp) const
      {

        double res=0.;

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
            int main_patch=1;           //the patch we will perform our computations on
            if(supp.xmin[0]!=-1){
                main_patch=0;
            }
            if(supp.xmin[2]!=-1){
                main_patch=2;
            }
            //cout << "Mainpatch: " << main_patch << endl;
                    
            bool mirrored_la, mirrored_mu;
            //are lambda, mu mirrored into the dircetion of the integration
            mirrored_la=((lambda->patch() == 3 && dir == 1) || (lambda->patch() == 4 && dir == 0)) ? 1 : 0;
            mirrored_mu=((mu->patch() == 3 && dir == 1) || (mu->patch() == 4 && dir == 0)) ? 1 : 0;
         
            
//            int mother_patch_la, mother_patch_mu;
//            mother_patch_la=(lambda->patch()==3) ? 0: (lambda->patch()==4 ? 2: lambda->patch());
//            mother_patch_mu=(mu->patch()==3) ? 0: (mu->patch()==4 ? 2: mu->patch());
            
            //compute gauss points and weights
            double a ; a=(dir==0? supp.xmin[main_patch] : supp.ymin[main_patch]);
            double b ; b=(dir==0? supp.xmax[main_patch] : supp.ymax[main_patch]);
            //cout << "current support: ["<<a[0]<<","<<b[0]<<"]x["<<a[1]<<","<<b[1]<<"]"<<endl;
            
            
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
                            
            
     
            const IFRAME* frame1D_lambda = frame_->frames(lambda->patch(),dir);
            const IFRAME* frame1D_mu     = frame_->frames(mu->patch(),dir);

            if(main_patch == 1){
                if(mirrored_la){
                    for(int i=0;i<(int)lambda_gauss_points.size();i++){
                             lambda_gauss_points[i]=1-lambda_gauss_points[i];
                    }
                }
                if(mirrored_mu){
                    for(int i=0;i<(int)mu_gauss_points.size();i++){
                             mu_gauss_points[i]=1-mu_gauss_points[i];
                    }
                }
            }
                
            WaveletTL::evaluate(*frame1D_lambda, lambda->derivative(),  lambda->index(),         //in x-direction evaluate with asymmetric boundary conditions
                                                /*lambda->index().p(),
                                                lambda->index().j(),
                                                lambda->index().e(),
                                                lambda->index().k(),*/
                         lambda_gauss_points, lambda_values);
            
            WaveletTL::evaluate(*frame1D_mu, mu->derivative(), mu->index(),          //in x-direction evaluate with asymmetric boundary conditions
                                                /*mu->index().p(),
                                                mu->index().j(),
                                                mu->index().e(),
                                                mu->index().k(),*/
                         mu_gauss_points, mu_values);
            
            if(main_patch == 1){
                if(mirrored_la && lambda->derivative() == 1){
                    for(int i=0;i<(int)lambda_values.size();i++){
                             lambda_values[i]=-lambda_values[i];
                    }
                }
                if(mirrored_mu && mu->derivative() == 1){
                    for(int i=0;i<(int)mu_values.size();i++){
                             mu_values[i]=-mu_values[i];
                    }
                }
            }

            // - add all integral shares
                
            for (unsigned int ind = 0; ind < gauss_points.size(); ind++){
                   res += lambda_values[ind] * mu_values[ind] * gauss_weights[ind];
//                    der_integral += lambda_der_values[ind] * mu_der_values[ind] * gauss_weights[ind];
            }
            
            
            //the result needs to be doubled if both quarklets are mirrored
            //into the direction of integration
            if(mirrored_la && mirrored_mu){
                res*=2;
//                der_integral*=2;
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
        return res;
      }
#endif
    }
    
    


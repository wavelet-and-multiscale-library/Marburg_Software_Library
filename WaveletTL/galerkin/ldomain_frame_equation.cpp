// implementation for ldomain_frame_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>
#include <numerics/eigenvalues.h>

namespace WaveletTL
{
    template <class IFRAME>
    LDomainFrameEquation<IFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp,
                                            const FixedArray1D<bool,8>& bc,
					    const bool precompute_rhs)
    : bvp_(bvp), frame(bc), normA(0.0), normAinv(0.0)
    {
        if (precompute_rhs)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
            compute_rhs();
        }
    }
  
    template <class IFRAME>
    LDomainFrameEquation<IFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp,
					   const FixedArray1D<int,8>& bc,
					   const bool precompute_rhs)
    : bvp_(bvp), basis_(basis), normA(0.0), normAinv(0.0)
    {
        if (precompute_rhs)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
            compute_rhs();
        }
    }
  
    template <class IFRAME>
    LDomainFrameEquation<IFRAME>::LDomainFrameEquation(const LDomainFrameEquation& eq)
    : bvp_(eq.bvp_), frame_(eq.frame_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
    {
        cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
    }
    
    template <class IFRAME>
    void
    LDomainFrameEquation<IFRAME>::compute_rhs()
    {
        cout << "LDomainFrameEquation(): precompute right-hand side..." << endl;
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
        cout << "... done, sort the entries in modulus..." << endl;

        // sort the coefficients into fcoeffs
        fcoeffs.resize(0); // clear eventual old values
        fcoeffs.resize(fhelp.size());
        unsigned int id(0);
        for (typename InfiniteVector<double,Index>::const_iterator it(fhelp.begin()), itend(fhelp.end()); it != itend; ++it, ++id){
            fcoeffs[id] = std::pair<Index,double>(it.index(), *it);
        }
        sort(fcoeffs.begin(), fcoeffs.end(), typename InfiniteVector<double,Index>::decreasing_order());
        cout << "... done, all integrals for right-hand side computed!" << endl;
    }

    template <class IFRAME>
    inline
    double
    LDomainFrameEquation<IFRAME>::D(const Index& lambda) const
    {
        //return ldexp(1.0, lambda.j());
        //return sqrt(a(lambda, lambda));
        double hspreconditioner(0), l2preconditioner(1);
        for (int i=0; i<space_dimension; i++){
            hspreconditioner+=pow(1+lambda.p()[i],8)*ldexp(1.0,2*lambda.j()[i]);
            l2preconditioner*=pow(1+lambda.p()[i],2);
        }
        double preconditioner = sqrt(hspreconditioner)*l2preconditioner;
        
        return preconditioner;
    }

    template <class IFRAME>
    inline
    double
    LDomainFrameEquation<IFRAME>::a(const Index& lambda,
			     const Index& nu) const
    {
        return a(lambda, nu, IFRAME::primal_polynomial_degree()*IFRAME::primal_polynomial_degree());
    }

    template <class IFRAME>
    double
    LDomainFrameEquation<IFRAME>::a(const Index& lambda,
			     const Index& mu,
			     const unsigned int p) const
    {
        // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
        double r = 0;
        // first compute the support intersection of psi_lambda and psi_mu
        typename Frame::Support supp;
        if (intersect_supports(frame_, lambda, mu, supp)) {
            //determine on which patches the support lies
            int main_patch=1;           //the patch we will perform our computations on
            if(supp.xmin[0]!=-1) main_patch=0;
            if(supp.xmin[2]!=-1) main_patch=2;
                    
            bool extended [2];
            extended[0]=(lambda.patch()>2) ? 1 : 0;
            extended[1]=(mu.patch()>2) ? 1 : 0;
            int mother_patch [2];
            mother_patch[0]=(lambda.patch()==3) ? 0: (lambda.patch()==4 ? 2: lambda.patch());
            mother_patch[1]=(mu.patch()==3) ? 0: (mu.patch()==4 ? 2: mu.patch());
            //due to symmetry we compute the bilinearform just on one patch
            
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = (p+1)/2+(multi_degree(lambda.p())+multi_degree(mu.p())+1)/2;
            double h; // granularity for the quadrature
            FixedArray1D<Array1D<double>,space_dimension> gauss_points, gauss_weights;
            //compute gauss points and weights for x- and y-direction (since the quarklets are anisotropic)
            double a [2]; a[0]=supp.xmin[main_patch]; a[1]=supp.ymin[main_patch];
            double b [2]; b[0]=supp.xmax[main_patch]; b[1]=supp.ymax[main_patch];
            
            //compute gauss points on main patch
            for(int i=0;i<space_dimension;i++){
                h = ldexp(1.0, -supp.j[i]);
                gauss_points[i].resize(N_Gauss*(b[i]-a[i]));
                gauss_weights[i].resize(N_Gauss*(b[i]-a[i]));
                for (int interval = a[i]; interval < b[i]; interval++){
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(interval-a[i])*N_Gauss+n]= h*(2*interval+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(interval-a[i])*N_Gauss+n]= h*GaussWeights[N_Gauss-1][n];
                    }
                }
            }
            
            //temp gauss points for reflection
            FixedArray1D<Array1D<double>,space_dimension> lambda_gauss_points(gauss_points), mu_gauss_points(gauss_points);
            //reflect gauss points
            if(main_patch==1){
                if(mother_patch[0]==0){ //lambda has been extended from north so south, reflect gauss points in y-direction
                    for(int i=0;i<lambda_gauss_points[1].size();i++){
                        lambda_gauss_points[1][i]=1-lambda_gauss_points[1][i];
                    }    
                }
                if(mother_patch[0]==2){ //lambda has been extended from east to west, reflect gauss points in x-direction
                    for(int i=0;i<lambda_gauss_points[0].size();i++){
                        lambda_gauss_points[0][i]=1-lambda_gauss_points[0][i];
                    }  
                }
                if(mother_patch[1]==0){ //mu has been extended from north so south, reflect gauss points in y-direction
                    for(int i=0;i<mu_gauss_points[1].size();i++){
                        mu_gauss_points[1][i]=1-mu_gauss_points[1][i];
                    }  
                }
                if(mother_patch[1]==2){ //mu has been extended from east to west, reflect gauss points in x-direction
                    for(int i=0;i<mu_gauss_points[0].size();i++){
                        mu_gauss_points[0][i]=1-mu_gauss_points[0][i];
                    } 
                }
            }
            
            
            // compute point values of the integrand (where we use that it is a tensor product)
            // evaluate method of the interval frame is called
            FixedArray1D<Array1D<double>,space_dimension> psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
                                              psi_mu_values,         // -"-, for psi_mu
                                              psi_lambda_der_values, // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]
                                              psi_mu_der_values;     // -"-, for psi_mu
            for (unsigned int i = 0; i < space_dimension; i++) {
                evaluate(*frame_.frame1d(), 0,
                                                lambda.p()[i],
                                                lambda.j()[i],
                                                lambda.e()[i],
                                                lambda.k()[i],
                         lambda_gauss_points[i], psi_lambda_values[i]);
                evaluate(*frame_.frame1d(), 1,
                                                lambda.p()[i],
                                                lambda.j()[i],
                                                lambda.e()[i],
                                                lambda.k()[i],
                         lambda_gauss_points[i], psi_lambda_der_values[i]);
                evaluate(*frame_.frame1d(), 0,
                                                mu.p()[i],
                                                mu.j()[i],
                                                mu.e()[i],
                                                mu.k()[i],
                         mu_gauss_points[i], psi_mu_values[i]);
                evaluate(*frame_.frame1d(), 1,
                                                mu.p()[i],
                                                mu.j()[i],
                                                mu.e()[i],
                                                mu.k()[i],
                         mu_gauss_points[i], psi_mu_der_values[i]);
            }
            
            // iterate over all points and sum up the integral shares
            int index[space_dimension]; // current multiindex for the point values
            for (unsigned int i = 0; i < space_dimension; i++)
                index[i] = 0;
            Point<space_dimension> x;
            const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
            const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
            double grad_psi_lambda[space_dimension], grad_psi_mu[space_dimension], weights;
            if (bvp_->constant_coefficients())
            {
                while (true)
                {
                    for (unsigned int i = 0; i < space_dimension; i++)
                        x[i] = gauss_points[i][index[i]];
                    // product of current Gauss weights
                    weights = 1.0;
                    for (unsigned int i = 0; i < space_dimension; i++)
                        weights *= gauss_weights[i][index[i]];
                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
                    for (unsigned int i = 0; i < space_dimension; i++)
                    {
                        grad_psi_lambda[i] = 1.0;
                        grad_psi_mu[i] = 1.0;
                        for (unsigned int s = 0; s < space_dimension; s++) 
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
                    for (unsigned int i = 0; i < space_dimension; i++)
                        share += grad_psi_lambda[i]*grad_psi_mu[i];
                    r += ax * weights * share;
                    // compute the share q(x)psi_lambda(x)psi_mu(x)
                    share = qx * weights;
                    for (unsigned int i = 0; i < space_dimension; i++)
                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    for(unsigned int i=0;i<space_dimension;i++){
                        if (index[i] == N_Gauss*(b[i]-a[i])-1)    
                        {
                            index[i] = 0;
                            exit = (i == space_dimension-1);
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
                    for (unsigned int i = 0; i < space_dimension; i++)
                        x[i] = gauss_points[i][index[i]];
                    // product of current Gauss weights
                    weights = 1.0;
                    for (unsigned int i = 0; i < space_dimension; i++)
                        weights *= gauss_weights[i][index[i]];
                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
                    for (unsigned int i = 0; i < space_dimension; i++) {
                        grad_psi_lambda[i] = 1.0;
                        grad_psi_mu[i] = 1.0;
                        for (unsigned int s = 0; s < space_dimension; s++) {
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
                    for (unsigned int i = 0; i < space_dimension; i++)
                        share += grad_psi_lambda[i]*grad_psi_mu[i];
                    r += bvp_->a(x) * weights * share;
                    // compute the share q(x)psi_lambda(x)psi_mu(x)
                    share = bvp_->q(x) * weights;
                    for (unsigned int i = 0; i < space_dimension; i++)
                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    
                    for(unsigned int i=0;i<space_dimension;i++){
                        if (index[i] == N_Gauss*(b[i]-a[i])-1)    
                        {
                            index[i] = 0;
                            exit = (i == space_dimension-1);
                        } else
                        {
                            index[i]++;
                            break;
                        }
                    }
                    
                    if (exit) break;
                }
            }
               
        }
        if(extended[0] && extended[1] && mother_patch[0]==mother_patch[1]) return 2*r;
        else return r;
    }
    
    template <class IFRAME>
    double
    LDomainFrameEquation<IFRAME>::f(const Index& lambda) const
    {
        // f(v) = \int_0^1 g(t)v(t) dt
        //cout << "LDomainEquation::f() called with lambda=" << lambda << endl;
        double r = 0;

        typedef typename Frame::Support Support;
        Support supp;
        support(frame_,lambda,supp);
        
        //determine on which patch lambda lies
        const bool extended=(lambda.patch()>2) ? 1 : 0;
        const int mother_patch=(lambda.patch()==3) ? 0: (lambda.patch()==4 ? 2: lambda.patch());
    
        //setup gauss points
        const int N_Gauss = (p+1)/2+(multi_degree(lambda.p())+multi_degree(mu.p())+1)/2;
        double h; // granularity for the quadrature
        FixedArray1D<Array1D<double>,space_dimension> gauss_points, gauss_weights, v_values;
	
 	// per patch, collect all point values
 	for (int patch = 0; patch <= 2; patch++) {
 	    if (supp.xmin[patch] != -1) { // psi_mu is nontrivial on patch p
                double a [2]; a[0]=supp.xmin[patch]; a[1]=supp.ymin[patch];
                double b [2]; b[0]=supp.xmax[patch]; b[1]=supp.ymax[patch];
                for (unsigned int i = 0; i < space_dimension; i++) {
                    // prepare Gauss points and weights 
                    h = ldexp(1.0,-supp.j[i]);
                    gauss_points[i].resize(N_Gauss*(b[i]-a[i]));
                    gauss_weights[i].resize(N_Gauss*(b[i]-a[i]));
                    for (int interval = a[i]; interval < b[i]; interval++){
                        for (int n = 0; n < N_Gauss; n++) {
                            gauss_points[i][(interval-a[i]])*N_Gauss = h*(2*interval+1+GaussPoints[N_Gauss-1][n])/2.;
                            gauss_weights[i][(interval-a[i])*N_Gauss+n]= h*GaussWeights[N_Gauss-1][n];
                        }
                    }
                }
                //special treat for patch 1 if the quarklet is extended:    
                FixedArray1D<Array1D<double>,space_dimension> lambda_gauss_points(gauss_points);
                //reflect gauss points
                if(patch==1){
                    if(mother_patch==0){ //lambda has been extended from north so south, reflect gauss points in y-direction
                        for(int i=0;i<lambda_gauss_points[1].size();i++){
                            lambda_gauss_points[1][i]=1-lambda_gauss_points[1][i];
                        }    
                    }
                    if(mother_patch==2){ //lambda has been extended from east to west, reflect gauss points in x-direction
                        for(int i=0;i<lambda_gauss_points[0].size();i++){
                            lambda_gauss_points[0][i]=1-lambda_gauss_points[0][i];
                        }  
                    }
                } 
            
                // compute the point values of the integrand (where we use that it is a tensor product)
                for (unsigned int i = 0; i < space_dimension; i++){
                    evaluate(*frame_.frame1d(), 0,
                                            lambda.p()[i],
                                            lambda.j()[i],
                                            lambda.e()[i],
                                            lambda.k()[i],
                         lambda_gauss_points[i], v_values[i]);
                }
                // iterate over all points and sum up the integral shares
                int index[space_dimension]; // current multiindex for the point values
                for (unsigned int i = 0; i < space_dimension; i++)
                    index[i] = 0;
                Point<space_dimension> x;
                while (true) {
                    for (unsigned int i = 0; i < space_dimension; i++)
                        x[i] = gauss_points[i][index[i]];
                    double share = bvp_->f(x);
                    for (unsigned int i = 0; i < space_dimension; i++)
                        share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    for (unsigned int i = 0; i < space_dimension; i++) {
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
    
    template <class IFRAME>
    void
    LDomainFrameEquation<IFRAME>::RHS(const double eta,
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
        

    template <class IFRAME>
    double
    LDomainFrameEquation<IFRAME>::norm_A() const
    {
        if (normA == 0.0) {
        //typedef typename Frame::Index Index;
        typedef typename Index::polynomial_type polynomial_type;
        int offsetj=1;
        int offsetp=std::min((int)frame_.get_pmax(),2);
        std::set<Index> Lambda;
      
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
   
    
    
    template <class IFRAME>
    double
    LDomainFrameEquation<IFRAMe>::norm_Ainv() const
    {
        if (normA == 0.0) {
        //typedef typename Frame::Index Index;
        typedef typename Index::polynomial_type polynomial_type;
        int offsetj=1;
        int offsetp=std::min((int)frame_.get_pmax(),2);
        std::set<Index> Lambda;
      
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
        return normAinv;
    }
  

    template <class IFRAME>
    double
    LDomainFrameEquation<IFRAME>::s_star() const
    {
        // notation from [St04a]
        const double t = operator_order();
        const int n = 2;
        const int dT = Frame::primal_vanishing_moments();
        const double gamma = Frame::primal_regularity();
    
        return std::min((t+dT)/(double)n, (gamma-t)/(n-1.)); // [St04a, Th. 2.3]
    }


}

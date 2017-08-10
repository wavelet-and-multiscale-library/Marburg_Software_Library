// implementation for ldomain_frame_equation.h

#include <cmath>
#include <time.h>
#include <utils/fixed_array1d.h>
#include <numerics/gauss_data.h>
#include <numerics/eigenvalues.h>
#include <Ldomain/ldomain_frame_support.h>


namespace WaveletTL
{
    template <class IFRAME, class LDOMAINFRAME>
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp,
					    const bool precompute_rhs)
    : bvp_(bvp), frame_(), normA(0.0), normAinv(0.0)
    {
        if (precompute_rhs)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
            compute_rhs();
        }
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const EllipticBVP<2>* bvp, const Frame& frame,
					    const bool precompute_rhs)
    : bvp_(bvp), frame_(frame), normA(0.0), normAinv(0.0)
    {
        if (precompute_rhs)
        {
            cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
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
    //        cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
    //        frame_.set_jpmax(multi_degree(frame_.j0()),0);
    //        compute_rhs();
    //   }
    //}
  
    template <class IFRAME, class LDOMAINFRAME>
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::LDomainFrameEquation(const LDomainFrameEquation& eq)
    : bvp_(eq.bvp_), frame_(eq.frame_),
      fcoeffs(eq.fcoeffs), fnorm_sqr(eq.fnorm_sqr),
      normA(eq.normA), normAinv(eq.normAinv)
    {
        cout << "Maximal level is set to "<<multi_degree(frame_.j0())<< ". Maximal polynomial to " << 0 << ". You may want to increase that." << endl;
            frame_.set_jpmax(multi_degree(frame_.j0()),0);
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    void
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::compute_rhs()
    {
        cout << "LDomainFrameEquation(): precompute right-hand side..." << endl;
        // precompute the right-hand side on a fine level
        InfiniteVector<double,Index> fhelp;
        fnorm_sqr = 0;
        for (int i = 0; (int)i< frame_.degrees_of_freedom();i++)
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

    template <class IFRAME, class LDOMAINFRAME>
    inline
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::D(const Index& lambda) const
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
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::a(const Index& lambda,
			     const Index& mu,
			     const unsigned int p) const
    {
        
        // a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
        double r = 0;
        double integral[space_dimension], der_integral[space_dimension];
        integral[0]=0, integral[1]=0, der_integral[0]=0, der_integral[1]=0;
        // first compute the support intersection of psi_lambda and psi_mu
        typedef typename Frame::Support Support;
        Support supp;
        
        if (intersect_supports(frame_, lambda, mu, supp)) {
            
//            cout << "support in a-routine for index " << lambda << ", " << mu << endl;
//    cout << "j[0]: " << supp.j[0] << ", j[1]: " << supp.j[1] << endl;
//    cout << "patch 0: [" << supp.xmin[0] <<" , " <<supp.xmax[0]<<"]x["<<supp.ymin[0]<<" , "<<supp.ymax[0]<<"]"<<  endl;
//    cout << "patch 1: [" << supp.xmin[1] <<" , " <<supp.xmax[1]<<"]x["<<supp.ymin[1]<<" , "<<supp.ymax[1]<<"]"<<  endl;
//    cout << "patch 2: [" << supp.xmin[2] <<" , " <<supp.xmax[2]<<"]x["<<supp.ymin[2]<<" , "<<supp.ymax[2]<<"]"<<  endl;
//            cout << "Supports intersect!" << endl;
//            cout << "Vergleichswerte: " << supp.xmin[0] << ", " << supp.xmin[1] << ", " << supp.xmin[2] << endl;
            //determine on which patches the support lies
            int main_patch=1;           //the patch we will perform our computations on
            if(supp.xmin[0]!=-1){
                main_patch=0;
            }
            if(supp.xmin[2]!=-1){
                main_patch=2;
            }
            //cout << "Mainpatch: " << main_patch << endl;
                    
            bool extended[2];
            extended[0]=(lambda.patch()>2) ? 1 : 0;
            extended[1]=(mu.patch()>2) ? 1 : 0;
            int mother_patch [2];
            mother_patch[0]=(lambda.patch()==3) ? 0: (lambda.patch()==4 ? 2: lambda.patch());
            mother_patch[1]=(mu.patch()==3) ? 0: (mu.patch()==4 ? 2: mu.patch());
            //due to symmetry we compute the bilinearform just on one patch
            //if(main_patch==1 && mother_patch[0]==mother_patch[1]) cout << mother_patch[0] << endl;
            
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = (p+1)/2+(multi_degree(lambda.p())+multi_degree(mu.p())+1)/2;
            //const int N_Gauss=10;
            double h; // granularity for the quadrature
            FixedArray1D<Array1D<double>,space_dimension> gauss_points, gauss_weights;
            //compute gauss points and weights for x- and y-direction (since the quarklets are anisotropic)
            double a [2]; a[0]=supp.xmin[main_patch]; a[1]=supp.ymin[main_patch];
            double b [2]; b[0]=supp.xmax[main_patch]; b[1]=supp.ymax[main_patch];
            //cout << "current support: ["<<a[0]<<","<<b[0]<<"]x["<<a[1]<<","<<b[1]<<"]"<<endl;
            
            //compute gauss points on main patch
            for(int i=0;i<space_dimension;i++){
//                int e = std::max(lambda.e()[i], mu.e()[i]); //correct granularity
                h = ldexp(1.0, -supp.j[i]/*-e*/);
                gauss_points[i].resize(N_Gauss*(b[i]-a[i]));
                gauss_weights[i].resize(N_Gauss*(b[i]-a[i]));
                for (int interval = a[i]; interval < b[i]; interval++){
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(interval-a[i])*N_Gauss+n]= h*(2*interval+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(interval-a[i])*N_Gauss+n]= h*GaussWeights[N_Gauss-1][n];
                    }
                }
            }
            //cout << "Gauss Points LDomain vor Änderung xrichtung: " << gauss_points[0] << endl;
            //cout << "Gauss Points LDomain vor Änderung yrichtung: " << gauss_points[1] << endl;
            
            // compute point values of the integrand (where we use that it is a tensor product)
            //temp gauss points for reflection
            FixedArray1D<Array1D<double>,space_dimension> lambda_gauss_points(gauss_points), mu_gauss_points(gauss_points);
            // evaluate method of the interval frame is called
            FixedArray1D<Array1D<double>,space_dimension> psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
                                              psi_mu_values,         // -"-, for psi_mu
                                              psi_lambda_der_values, // values of the 1st deriv. of the components of psi_lambda at gauss_points[i]
                                              psi_mu_der_values;     // -"-, for psi_mu
            //cout << main_patch <<endl;
            switch(main_patch){
                case 0:            //both functions lie on patch 0, evaluate with 0,1 boundary conditions
                    evaluate(frame_.frame1d_11(),            //in x-direction evaluate with homogeneous boundary conditions
                                                lambda.p()[0],
                                                lambda.j()[0],
                                                lambda.e()[0],
                                                lambda.k()[0],
                         gauss_points[0], psi_lambda_values[0], psi_lambda_der_values[0]);

                    evaluate(frame_.frame1d_11(),
                                                mu.p()[0],
                                                mu.j()[0],
                                                mu.e()[0],
                                                mu.k()[0],
                         gauss_points[0], psi_mu_values[0], psi_mu_der_values[0]);

                
                    evaluate(frame_.frame1d_01(),                //in y-direction evaluate with asymmetric boundary conditions
                                                lambda.p()[1],
                                                lambda.j()[1],
                                                lambda.e()[1],
                                                lambda.k()[1],
                         gauss_points[1], psi_lambda_values[1], psi_lambda_der_values[1]);

                    evaluate(frame_.frame1d_01(),
                                                mu.p()[1],
                                                mu.j()[1],
                                                mu.e()[1],
                                                mu.k()[1],
                         gauss_points[1], psi_mu_values[1], psi_mu_der_values[1]);

                    break;
            
            
            //reflect gauss points
                case 1:
                //first assume that both functions stem from patch 1, evaluate with zero bc
                //cout << "bin hier" << endl;
                //cout << "points and values before reflection"<<endl;
                //cout << "mu_gauss_points[0]: "<<mu_gauss_points[0]<<endl;
                //cout << "mu_gauss_points[1]: "<<mu_gauss_points[1]<<endl;
                    evaluate(frame_.frame1d_11(),            
                                                lambda.p()[0],
                                                lambda.j()[0],
                                                lambda.e()[0],
                                                lambda.k()[0],
                                        lambda_gauss_points[0], psi_lambda_values[0], psi_lambda_der_values[0]);
                //cout << "Gauss Points LDomain: " << lambda_gauss_points[0] << endl;
               // cout << "psi_lambda_values[0]: "<<psi_lambda_values[0] << endl;
//                    
                //cout <<"psi_lambda_der_values[0]: "<< psi_lambda_der_values[0] << endl;
                    evaluate(frame_.frame1d_11(),
                                                mu.p()[0],
                                                mu.j()[0],
                                                mu.e()[0],
                                                mu.k()[0],
                                        mu_gauss_points[0], psi_mu_values[0], psi_mu_der_values[0]);
                //cout << "psi_mu_values[0]: "<<psi_mu_values[0] << endl;
//                    
                //cout << "psi_mu_der_values[0]: "<<psi_mu_der_values[0] << endl;
                
                    evaluate(frame_.frame1d_11(),                
                                                lambda.p()[1],
                                                lambda.j()[1],
                                                lambda.e()[1],
                                                lambda.k()[1],
                                        lambda_gauss_points[1], psi_lambda_values[1], psi_lambda_der_values[1]);
                //cout << "psi_lambda_values[1]: "<<psi_lambda_values[1] << endl;
//                    
                //cout <<"psi_lambda_der_values[1]"<< psi_lambda_der_values[1] << endl;
                    evaluate(frame_.frame1d_11(),
                                                mu.p()[1],
                                                mu.j()[1],
                                                mu.e()[1],
                                                mu.k()[1],
                                        mu_gauss_points[1], psi_mu_values[1], psi_mu_der_values[1]);
                //cout << "psi_mu_values[1]: "<<psi_mu_values[1] << endl;
//                    
                //cout << "psi_mu_der_values[1]: "<<psi_mu_der_values[1] << endl;
                    if(mother_patch[0]==0 && mother_patch[1]!=0){ //lambda has been extended from north so south, reflect gauss points in y-direction
                        for(int i=0;i<(int)lambda_gauss_points[1].size();i++){
                            lambda_gauss_points[1][i]=1-lambda_gauss_points[1][i];
                        }    
                        //evaluate lambda with asymmetric bc in y-direction
                        evaluate(frame_.frame1d_01(),                
                                                lambda.p()[1],
                                                lambda.j()[1],
                                                lambda.e()[1],
                                                lambda.k()[1],
                                        lambda_gauss_points[1], psi_lambda_values[1], psi_lambda_der_values[1]);
//                        
                        for(int i=0;i<(int)psi_lambda_values[1].size();i++){
                            psi_lambda_der_values[1][i]=-psi_lambda_der_values[1][i];
                        }
                          
                    }
                    if(mother_patch[0]==2 && mother_patch[1]!=2){ //lambda has been extended from east to west, reflect gauss points in x-direction
                        for(int i=0;i<(int)lambda_gauss_points[0].size();i++){
                            lambda_gauss_points[0][i]=1-lambda_gauss_points[0][i];
                        }    
                        //evaluate lambda with asymmetric bc in x-direction, mu with zero boundary conditions
                        evaluate(frame_.frame1d_01(),            
                                                lambda.p()[0],
                                                lambda.j()[0],
                                                lambda.e()[0],
                                                lambda.k()[0],
                                        lambda_gauss_points[0], psi_lambda_values[0], psi_lambda_der_values[0]);
//                        
                        for(int i=0;i<(int)psi_lambda_values[0].size();i++){
                            psi_lambda_der_values[0][i]=-psi_lambda_der_values[0][i];
                        }
                               
                    }
                    if(mother_patch[1]==0 && mother_patch[0]!=0){ //mu has been extended from north so south, reflect gauss points in y-direction
                        for(int i=0;i<(int)mu_gauss_points[1].size();i++){
                            mu_gauss_points[1][i]=1-mu_gauss_points[1][i];
                        }  
                        //cout << "bin hier"<<endl;
                        //evaluate mu with asymmetric bc in y-direction, lambda with zero boundary conditions
                        evaluate(frame_.frame1d_01(),
                                                mu.p()[1],
                                                mu.j()[1],
                                                mu.e()[1],
                                                mu.k()[1],
                                        mu_gauss_points[1], psi_mu_values[1], psi_mu_der_values[1]);
                        //cout << "psi_mu_der_values[1]: "<<psi_mu_der_values[1] << endl;
//                        
                        for(int i=0;i<(int)psi_mu_values[1].size();i++){
                            psi_mu_der_values[1][i]=-psi_mu_der_values[1][i];
                        }
                        //cout << "psi_mu_der_values[1]: "<<psi_mu_der_values[1] << endl;
                        
                      
                    }
                    if(mother_patch[1]==2 && mother_patch[0]!=2){ //mu has been extended from east to west, reflect gauss points in x-direction
                        for(int i=0;i<(int)mu_gauss_points[0].size();i++){
                            mu_gauss_points[0][i]=1-mu_gauss_points[0][i];
                        }    
                        //evaluate mu with asymmetric bc in x-direction, lambda with zero boundary conditions
                        evaluate(frame_.frame1d_01(),
                                                mu.p()[0],
                                                mu.j()[0],
                                                mu.e()[0],
                                                mu.k()[0],
                                        mu_gauss_points[0], psi_mu_values[0], psi_mu_der_values[0]);
//                        
                        for(int i=0;i<(int)psi_mu_values[0].size();i++){
                            psi_mu_der_values[0][i]=-psi_mu_der_values[0][i];
                        }
                         
                    }
                    else{}
                    break;
                    
            
                case 2:  //both functions lie on patch 2, evaluate with 0,1 boundary conditions
                    evaluate(frame_.frame1d_01(),            //in x-direction evaluate with asymmetric boundary conditions
                                                lambda.p()[0],
                                                lambda.j()[0],
                                                lambda.e()[0],
                                                lambda.k()[0],
                         gauss_points[0], psi_lambda_values[0], psi_lambda_der_values[0]);
//                cout << psi_lambda_values[0] << endl;
//                    
//                cout << psi_lambda_der_values[0] << endl;
                    evaluate(frame_.frame1d_01(),
                                                mu.p()[0],
                                                mu.j()[0],
                                                mu.e()[0],
                                                mu.k()[0],
                         gauss_points[0], psi_mu_values[0], psi_mu_der_values[0]);
//                cout << psi_mu_values[0] << endl;
//                    
//                cout << psi_mu_der_values[0] << endl;
                
                    evaluate(frame_.frame1d_11(),                //in y-direction evaluate with homogeneous boundary conditions
                                                lambda.p()[1],
                                                lambda.j()[1],
                                                lambda.e()[1],
                                                lambda.k()[1],
                         gauss_points[1], psi_lambda_values[1], psi_lambda_der_values[1]);
//                cout << psi_lambda_values[1] << endl;
//                    
//                cout << psi_lambda_der_values[1] << endl;
                    evaluate(frame_.frame1d_11(),
                                                mu.p()[1],
                                                mu.j()[1],
                                                mu.e()[1],
                                                mu.k()[1],
                         gauss_points[1], psi_mu_values[1], psi_mu_der_values[1]);
//                cout << psi_mu_values[1] << endl;
//                    
//                cout << psi_mu_der_values[1] << endl;
                    break;
            }
            //cout << "gauss points used for evaluation for a"<<endl;
            //    cout << lambda_gauss_points[0]<<endl;
            //    cout << lambda_gauss_points[1]<<endl;
            //cout << "points and values after reflection and evaluation"<<endl;
            //cout << "mu_gauss_points[0]: "<<mu_gauss_points[0]<<endl;
            //cout << "mu_gauss_points[1]: "<<mu_gauss_points[1]<<endl;
            //cout << "psi_lambda_values[0]: "<<psi_lambda_values[0] << endl;
            //cout << "psi_lambda_der_values[0]: "<<psi_lambda_der_values[0] << endl;
            //cout << "psi_mu_values[0]: "<<psi_mu_values[0] << endl;
            //cout << "psi_mu_der_values[0]: "<<psi_mu_der_values[0] << endl;
            //cout << "psi_lambda_values[1]: "<<psi_lambda_values[1] << endl;
            //cout << "psi_lambda_der_values[1]: "<<psi_lambda_der_values[1] << endl;
            //cout << "psi_mu_values[1]: "<<psi_mu_values[1] << endl;
            //cout << "psi_mu_der_values[1]: "<<psi_mu_der_values[1] << endl;
            
            
            
            
            // iterate over all points and sum up the integral shares
            int index[space_dimension]; // current multiindex for the point values
            for (int i = 0; i < space_dimension; i++)
                index[i] = 0;
            Point<space_dimension> x;
            const double ax = bvp_->constant_coefficients() ? bvp_->a(x) : 0.0;
            const double qx = bvp_->constant_coefficients() ? bvp_->q(x) : 0.0;
            double grad_psi_lambda[space_dimension], grad_psi_mu[space_dimension], weights;
            if (bvp_->constant_coefficients())
            {
                // - add all integral shares
                
                for (int i = 0; i < space_dimension; i++){
//                    cout << endl << "gauss_weights[" << i << "]: " << gauss_weights[i] << endl;
//                    cout << "psi_lambda_values[" << i << "]: " << psi_lambda_values[i] << endl;
//                    cout << "psi_mu_values[" << i << "]: " << psi_mu_values[i] << endl;
                    for (unsigned int ind = 0; ind < gauss_points[i].size(); ind++){
//                        if(i==1)
//                        cout << "Zwischenwert integral: " << integral[i] << endl;
                        integral[i] += psi_lambda_values[i][ind] * psi_mu_values[i][ind] * gauss_weights[i][ind];
                        der_integral[i] += psi_lambda_der_values[i][ind] * psi_mu_der_values[i][ind] * gauss_weights[i][ind];
                    }
                }
                
//                cout << "der_integral[0]: " << der_integral[0] << endl;
//                cout << "der_integral[1]: " << der_integral[1] << endl;
//                cout << "integral[0]: " << integral[0] << endl;
//                cout << "integral[1]: " << integral[1] << endl;
                
                r = ax * (der_integral[0] * integral[1] + integral[0] * der_integral[1]) + qx * (integral[0] * integral[1]);
                
                
//                while (true)
//                {
////                    switch(main_patch) {
////                        case 0:
////                            x[0] = gauss_points[0][index[0]]-1;
////                            x[1] = gauss_points[1][index[1]];
////                            break;
////                        case 1:
////                            x[0] = gauss_points[0][index[0]]-1;
////                            x[1] = gauss_points[1][index[1]]-1;
////                            break;
////                        case 2:
////                            x[0] = gauss_points[0][index[0]];
////                            x[1] = gauss_points[1][index[1]]-1;
////                            break;
////                        default:
////                            x[0] = gauss_points[0][index[0]];
////                            x[1] = gauss_points[1][index[1]];
////                    }
//                    // product of current Gauss weights
//                    weights = 1.0;
//                    for (int i = 0; i < space_dimension; i++)
//                        weights *= gauss_weights[i][index[i]];
//                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
//                    for (int i = 0; i < space_dimension; i++)
//                    {
//                        grad_psi_lambda[i] = 1.0;
//                        grad_psi_mu[i] = 1.0;
//                        for (int s = 0; s < space_dimension; s++) 
//                        {
//                            if (i == s)
//                            {
//                                grad_psi_lambda[i] *= psi_lambda_der_values[i][index[i]];
//                                grad_psi_mu[i]     *= psi_mu_der_values[i][index[i]];
//                            } else
//                            {
//                                grad_psi_lambda[i] *= psi_lambda_values[s][index[s]];
//                                grad_psi_mu[i] *= psi_mu_values[s][index[s]];
//                            }
//                        }
//                    }
//                    double share = 0;
//                    if(ax != 0){
//                        for (int i = 0; i < space_dimension; i++)
//                            share += grad_psi_lambda[i]*grad_psi_mu[i];
//                        r += ax * weights * share;
//                    }
//                    if(qx != 0){
//                    // compute the share q(x)psi_lambda(x)psi_mu(x)
//                    share = qx * weights;
//                    for (int i = 0; i < space_dimension; i++)
//                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
//                    r += share;
//                    }
//                    // "++index"
//                    bool exit = false;
//                    for(int i=0;i<space_dimension;i++){
//                        if (index[i] == N_Gauss*(b[i]-a[i])-1)    
//                        {
//                            index[i] = 0;
//                            exit = (i == space_dimension-1);
//                        } else
//                        {
//                            index[i]++;
//                            break;
//                        }
//                    }
//                    
//                    if (exit) break;
//                }
            } else // coefficients are not constant:
            {
                while (true) {
                    switch(main_patch) {
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
                    // product of current Gauss weights
                    weights = 1.0;
                    for (int i = 0; i < space_dimension; i++)
                        weights *= gauss_weights[i][index[i]];
                    // compute the share a(x)(grad psi_lambda)(x)(grad psi_mu)(x)
                    for (int i = 0; i < space_dimension; i++) {
                        grad_psi_lambda[i] = 1.0;
                        grad_psi_mu[i] = 1.0;
                        for (int s = 0; s < space_dimension; s++) {
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
                    for (int i = 0; i < space_dimension; i++)
                        share += grad_psi_lambda[i]*grad_psi_mu[i];
                    r += bvp_->a(x) * weights * share;
                    // compute the share q(x)psi_lambda(x)psi_mu(x)
                    share = bvp_->q(x) * weights;
                    for (int i = 0; i < space_dimension; i++)
                        share *= psi_lambda_values[i][index[i]] * psi_mu_values[i][index[i]];
                    r += share;
                    // "++index"
                    bool exit = false;
                    
                    for(int i=0;i<space_dimension;i++){
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
            if(extended[0] && extended[1] && mother_patch[0]==mother_patch[1]){
                r= 2*r; 
            }
        }
        
        return r;
    }
    
    template <class IFRAME, class LDOMAINFRAME>
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::f(const Index& lambda) const
    {
        // f(v) = \int_0^1 g(t)v(t) dt
//        cout << "LDomainEquation::f() called with lambda=" << lambda << endl;
        double r = 0;

        typedef typename Frame::Support Support;
        Support supp;
        support(frame_,lambda,supp);
//        cout << "patch 0: [" << supp.xmin[0]  <<" , " <<supp.xmax[0] <<"]x["<<supp.ymin[0] <<" , "<<supp.ymax[0] <<"]"<<  endl;
//        cout << "patch 1: [" << supp.xmin[1]  <<" , " <<supp.xmax[1] <<"]x["<<supp.ymin[1] <<" , "<<supp.ymax[1] <<"]"<<  endl;
//        cout << "patch 2: [" << supp.xmin[2]  <<" , " <<supp.xmax[2] <<"]x["<<supp.ymin[2] <<" , "<<supp.ymax[2] <<"]"<<  endl;
        
        //determine on which patch lambda lies

//        const bool extended=(lambda.patch()>2) ? 1 : 0;

        const int mother_patch=(lambda.patch()==3) ? 0: (lambda.patch()==4 ? 2: lambda.patch());
        //cout << "extended: "<<extended<<endl;
        //cout << "mother_patch: "<<mother_patch<<endl;
    
        //setup gauss points
        //const int N_Gauss = (p+1)/2+(multi_degree(lambda.p()))/2;
        const int N_Gauss=10;
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
                    evaluate(frame_.frame1d_11(), 0,        //different 1d-frames for x- and y-direction due to boundary conditions
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                         lambda_gauss_points[0], v_values[0]);  
                    evaluate(frame_.frame1d_01(), 0,
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
                        evaluate(frame_.frame1d_11(), 0,        //different 1d-frames for x- and y-direction due to boundary conditions
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(frame_.frame1d_01(), 0,
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
                        evaluate(frame_.frame1d_11(), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(frame_.frame1d_11(), 0,
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
                        evaluate(frame_.frame1d_01(), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                                        lambda_gauss_points[0], v_values[0]);  
                        evaluate(frame_.frame1d_11(), 0,
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
                    evaluate(frame_.frame1d_01(), 0,        
                                            lambda.p()[0],
                                            lambda.j()[0],
                                            lambda.e()[0],
                                            lambda.k()[0],
                         lambda_gauss_points[0], v_values[0]);  
                    evaluate(frame_.frame1d_11(), 0,
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
    double
    LDomainFrameEquation<IFRAME, LDOMAINFRAME>::norm_A() const
    {
        if (normA == 0.0) {
        //typedef typename Frame::Index Index;
        typedef typename Index::polynomial_type polynomial_type;
        int offsetj=std::min(1,frame_.get_jmax()-(int)multi_degree(frame_.j0()));
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
      
//        cout << "Eigenwerte: " << evals << endl;
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
        int offsetj=std::min(1,frame_.get_jmax()-(int)multi_degree(frame_.j0()));
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


}

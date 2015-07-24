// implementation for cdd2.h

#include <cmath>
#include <set>

#if _WAVELETTL_USE_TBASIS == 1
#include <adaptive/apply_tensor.h>
#else
#include <adaptive/apply.h>
#endif

using std::set;

namespace WaveletTL
{
    template <class PROBLEM>
    void CDD2_SOLVE(PROBLEM& P, const double nu, const double epsilon,
            InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
            const unsigned int maxlevel)
    {
        typedef typename PROBLEM::WaveletBasis::Index Index;

        // compute optimal relaxation parameter omega
        const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
        cout << "CDD2_SOLVE: omega=" << omega << endl;

        // compute spectral norm rho
        const double cond_A = P.norm_A() * P.norm_Ainv();
        const double rho = (cond_A - 1.0) / (cond_A + 1.0);
        cout << "CDD2_SOLVE: rho=" << rho << endl;

        // desired error reduction factor theta < 1/3
    //     const double theta = 2.0/7.0;
        const double theta = 0.333;
        cout << "CDD2_SOLVE: theta=" << theta << endl;

        // compute minimal K such that 3*rho^K < theta
        const int K = (int) ceil(log10(theta/3.0) / log10(rho));
        cout << "CDD2_SOLVE: K=" << K << endl;

        u_epsilon.clear();

        double epsilon_k = nu, eta;
        InfiniteVector<double,Index> f, v, Av;
#if _WAVELETTL_USE_TBASIS == 1
        Array1D<int> jp_guess(0);
#endif
        while (epsilon_k > epsilon) {
            epsilon_k *= 3*pow(rho, K) / theta;
            cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
            eta = theta * epsilon_k / (6*omega*K)*10;
            cout << "eta = " << eta << endl;
            P.RHS(eta, f);
            v = u_epsilon;
            for (int j = 1; j <= 1/*K*/; j++) {
#if _WAVELETTL_USE_TBASIS == 1
                APPLY_TENSOR(P, v, eta, Av, maxlevel, tensor_simple, true);
                //APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
#else
                APPLY(P, v, eta, Av, maxlevel, CDD1);
#endif
                //Av.COARSE(eta, Av);
                cout << "Number of degrees of freedom " << Av.size() << endl;
                cout << "current residual error ||f-Av||=" << l2_norm(f - Av) << endl;
                v += 0.5 * omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
            }
            cout << "coarse tol = " << (1-theta)*epsilon_k << endl;
            //v.COARSE((1-theta)*epsilon_k, u_epsilon);
            //v.COARSE(1.0e-6, u_epsilon);
            v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
            cout << "CDD2:: v.size() = " << v.size() << endl;
        } 
    }
    
    template <class PROBLEM>
    void CDD2_SOLVE(PROBLEM& P, const double nu, const double epsilon,
            InfiniteVector<double, int>& u_epsilon,
            const unsigned int maxlevel)
    {
//typedef typename PROBLEM::WaveletBasis::Index Index;

        // compute optimal relaxation parameter omega
        const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
        cout << "CDD2_SOLVE: omega=" << omega << endl;

        // compute spectral norm rho
        const double cond_A = P.norm_A() * P.norm_Ainv();
        const double rho = (cond_A - 1.0) / (cond_A + 1.0);
        cout << "CDD2_SOLVE: rho=" << rho << endl;

        // desired error reduction factor theta < 1/3
        const double theta = 0.333;
        cout << "CDD2_SOLVE: theta=" << theta << endl;

        // compute minimal K such that 3*rho^K < theta
        const int K = (int) ceil(log10(theta/3.0) / log10(rho));
        cout << "CDD2_SOLVE: K=" << K << endl;

        u_epsilon.clear();

        double epsilon_k = nu, eta;
        InfiniteVector<double,int> f, v, Av;
#if _WAVELETTL_USE_TBASIS == 1
        Array1D<int> jp_guess(0);
#endif
        while (epsilon_k > epsilon) 
        {
            epsilon_k *= 3*pow(rho, K) / theta;
            cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
            eta = theta * epsilon_k / (6*omega*K)*10;
            cout << "eta = " << eta << endl;
            P.RHS(eta, f);
            //cout << "l2norm(f) = " << l2_norm(f) << endl;
            v = u_epsilon;
            for (int j = 1; j <= 1/*K*/; j++) 
            {
#if _WAVELETTL_USE_TBASIS == 1
                APPLY_TENSOR(P, v, eta, Av, maxlevel, tensor_simple, true);
                //APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
#else
                APPLY(P, v, eta, Av, maxlevel, CDD1);
#endif
                //Av.COARSE(eta, Av);
                cout << "Number of degrees of freedom " << Av.size() << endl;
                cout << "current residual error ||f-Av||=" << l2_norm(f - Av) << endl;
                v += 0.5 * omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
            }
            cout << "coarse tol = " << (1-theta)*epsilon_k << endl;
            v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
            cout << "CDD2:: v.size() = " << v.size() << endl;
        } 
    }
  
  
}

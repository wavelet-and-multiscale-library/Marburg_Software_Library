// implementation for cdd2.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>

#if _WAVELETTL_USE_TBASIS == 1
#include <adaptive/apply_tensor.h>
#else
#include <adaptive/apply.h>
#endif

#ifndef PRIMALORDER
#define PRIMALORDER 3
#endif
#ifndef DUALORDER
#define DUALORDER 3
#endif
#ifndef JMAX
#define JMAX 8
#endif
#ifndef PMAX
#define PMAX 0
#endif


using std::set;
using std::map;

namespace WaveletTL
{
  template <class PROBLEM>
  void CDD2_SOLVE(const PROBLEM& P, const double nu, const double epsilon,
                  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
                  const unsigned int maxlevel, CompressionStrategy strategy)
  {
    typedef typename PROBLEM::WaveletBasis::Index Index;

    // compute optimal relaxation parameter omega
    const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 0.2;
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
    cout << "CDD2_SOLVE: K=" << K << endl << endl;

    u_epsilon.clear();
    
//    cout << "TEST_epsilon_k: " << 
    
    double epsilon_k = nu, eta;
        InfiniteVector<double,Index> f, v, Av, tempAv;
#if _WAVELETTL_USE_TBASIS == 1
        Array1D<int> jp_guess(0);
#endif
        while (epsilon_k > epsilon) {
        cout << "CDD2:: u.size() = " << u_epsilon.size() << endl;
        epsilon_k *= 3*pow(rho, K) / theta;
        cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
//        eta = theta * epsilon_k / (6*omega*K) ;//WAS SOLL HIER DIE 10@PHK
        eta = theta * epsilon_k / (6*omega*K)*10;
        cout << "eta = " << eta << endl;
        P.RHS(eta, f);
        cout << "CDD2:: f.size() = " << f.size() << endl;
        //cout << f << endl;
      
        v = u_epsilon;
        
        //cout << "2.CDD2:: v.size() = " << v.size() << endl;
      
        for (int j = 1; j <= 1 /*K*/; j++) {
#if _WAVELETTL_USE_TBASIS == 1
          APPLY(P, v, eta, Av, maxlevel, tensor_simple);
          //APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
#else
          APPLY_COARSE(P, v, eta, Av, 0.5, maxlevel, strategy);
//          APPLY(P, v, eta, Av, maxlevel, strategy);
          //APPLY with successive COARSE @PHK
//         APPLY(P, v, eta, tempAv, maxlevel, strategy, pmax, a, b);
//          tempAv.COARSE(1e-4, Av);
         

          
          
          
#endif
	//Av.COARSE(eta, tempAv);
        //Av = tempAv;
        /////cout << tempAv << endl;
        //cout << Av << endl;
//	cout << "Number of degrees of freedom (before coarsening) " << tempAv.size() << endl;
        cout << "Number of degrees of freedom " << Av.size() << endl;
        cout << "current residual error ||f-Av||=" << l2_norm(f - Av) << endl;
        cout << "coarse tol = " << (1-theta)*epsilon_k << endl;
	v += 0.5 *  omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
      }
      
      
      
      
      cout << "CDD2:: v.size() = " << v.size() << endl << endl;
//      v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
      v.COARSE((1-theta)*epsilon_k, u_epsilon);
//      v.COARSE(1.0e-6, u_epsilon);
      
      
      
      
//      cout << "f:" << endl<< f << endl;
//      cout << "Av:" << endl << Av << endl;
//      cout << "v:" << endl << v << endl;
//      cout << "u:" << endl <<u_epsilon << endl;
      
      
      
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
  
  template <class PROBLEM>
  void CDD2_QUARKLET_SOLVE(const PROBLEM& P, const double nu, const double epsilon,
                  InfiniteVector<double, typename PROBLEM::QuarkletFrame::Index>& u_epsilon,
                  const unsigned int maxlevel, CompressionStrategy strategy,
                  const int pmax, const double a, const double b)
  {
    typedef typename PROBLEM::QuarkletFrame::Index Index;

    // compute optimal relaxation parameter omega
    const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 0.2;
    cout << "CDD2_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    const double cond_A = P.norm_A() * P.norm_Ainv();
    cout << "CDD2_SOLVE: cond_A=" << cond_A << endl;
    const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    cout << "CDD2_SOLVE: rho=" << rho << endl;

    // desired error reduction factor theta < 1/3
    //     const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "CDD2_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log10(theta/3.0) / log10(rho));
    cout << "CDD2_SOLVE: K=" << K << endl << endl;

    u_epsilon.clear();
    
    double epsilon_k = nu, eta;
        InfiniteVector<double,Index> f, v, Av, tempAv;
#if _WAVELETTL_USE_TFRAME == 1
        Array1D<int> jp_guess(0);
#endif
        while (epsilon_k > epsilon) {
        cout << "CDD2:: u.size() = " << u_epsilon.size() << endl;
        epsilon_k *= 3*pow(rho, K) / theta;
        cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
        eta = theta * epsilon_k / (6*omega*K) ;//WAS SOLL HIER DIE 10@PHK
        //eta = theta * epsilon_k / (6*omega*K)*10;
        cout << "eta = " << eta << endl;
        P.RHS(eta, f);
        cout << "CDD2:: f.size() = " << f.size() << endl;
//        cout << "f: " << endl << f << endl;
      
        v = u_epsilon;
        
        //cout << "2.CDD2:: v.size() = " << v.size() << endl;
      
        for (int j = 1; j <= 1 /*K*/; j++) {
//#if _WAVELETTL_USE_TFRAME == 1
//          APPLY_QUARKLET(P, v, eta, Av, maxlevel, tensor_simple, pmax, a, b);
//          //APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
//#else
          //APPLY_COARSE(P, v, eta, Av, 0.5, maxlevel, CDD1);
//            cout << "v: " << endl << v << endl;
          APPLY_QUARKLET(P, v, eta, Av, maxlevel, strategy, pmax, a, b);
//          cout << "Av: " << endl << Av << endl;
          //APPLY with successive COARSE @PHK
//         APPLY(P, v, eta, tempAv, maxlevel, strategy, pmax, a, b);
//          tempAv.COARSE(1e-4, Av);
         

          
          
          
//#endif
	//Av.COARSE(eta, tempAv);
        //Av = tempAv;
        /////cout << tempAv << endl;
        //cout << Av << endl;
//	cout << "Number of degrees of freedom (before coarsening) " << tempAv.size() << endl;
        cout << "Number of degrees of freedom " << Av.size() << endl;
        cout << "current residual error ||f-Av||=" << l2_norm(f - Av) << endl;
        cout << "coarse tol = " << (1-theta)*epsilon_k << endl;
	v += 0.5 *  omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
      }
      
      
      
      
      cout << "CDD2:: v.size() = " << v.size() << endl << endl;
      v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
//      v.COARSE((1-theta)*epsilon_k, u_epsilon);
//      v.COARSE(1.0e-6, u_epsilon);
      
      
      
      
//      cout << "f:" << endl<< f << endl;
//      cout << "Av:" << endl << Av << endl;
//      cout << "v:" << endl << v << endl;
//      cout << "u:" << endl <<u_epsilon << endl;
      
      
      
    }
        
        
  }
  
  
  template <class PROBLEM>
  void CDD2_QUARKLET_SOLVE(const PROBLEM& P, const double nu, const double epsilon,
                  InfiniteVector<double, int>& u_epsilon,
                  const unsigned int maxlevel, CompressionStrategy strategy,
                  const int pmax, const double a, const double b)
  {
//    typedef typename PROBLEM::QuarkletFrame::Index Index;

    // compute optimal relaxation parameter omega
    const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    //const double omega = 0.2;
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
    cout << "CDD2_SOLVE: K=" << K << endl << endl;

    u_epsilon.clear();
    
    double epsilon_k = nu, eta;
    unsigned int loops(0), iter(0);
        InfiniteVector<double,int> f, v, Av, tempAv;
#if _WAVELETTL_USE_TFRAME == 1
        Array1D<int> jp_guess(0);
#endif
        while (epsilon_k > epsilon) {
            ++loops;
            cout << "CDD2:: u.size() = " << u_epsilon.size() << endl;
            epsilon_k *= 3*pow(rho, K) / theta;
            cout << "CDD2_SOLVE: epsilon_k=" << epsilon_k << endl;
            eta = theta * epsilon_k / (6*omega*K) ;//WAS SOLL HIER DIE 10@PHK
            //eta = theta * epsilon_k / (6*omega*K)*10;
            cout << "eta = " << eta << endl;
            P.RHS(eta, f);
            cout << "CDD2:: f.size() = " << f.size() << endl;
    //        cout << "f: " << endl << f << endl;

            v = u_epsilon;
        
        //cout << "2.CDD2:: v.size() = " << v.size() << endl;
      
        for (int j = 1; j <= 1 /*K*/; j++) {
//#if _WAVELETTL_USE_TFRAME == 1
//          APPLY_QUARKLET(P, v, eta, Av, maxlevel, tensor_simple, pmax, a, b);
//          //APPLY(P, v, eta, jp_guess, Av, maxlevel, tensor_simple);
//#else
          //APPLY_COARSE(P, v, eta, Av, 0.5, maxlevel, CDD1);
//            cout << "v: " << endl << v << endl;
        
	++iter;

	
//	
          APPLY_QUARKLET(P, v, eta, Av, maxlevel, strategy, pmax, a, b);
//          cout << "Av: " << endl << Av << endl;
          //APPLY with successive COARSE @PHK
//         APPLY(P, v, eta, tempAv, maxlevel, strategy, pmax, a, b);
//          tempAv.COARSE(1e-4, Av);
         

          
          
          
//#endif
	//Av.COARSE(eta, tempAv);
        //Av = tempAv;
        /////cout << tempAv << endl;
        //cout << Av << endl;
//	cout << "Number of degrees of freedom (before coarsening) " << tempAv.size() << endl;
        cout << "Number of degrees of freedom " << Av.size() << endl;
        cout << "current residual error ||f-Av||=" << l2_norm(f - Av) << endl;
        cout << "coarse tol = " << (1-theta)*epsilon_k << endl;
	v += 0.5 *  omega * (f - Av); // the factor 0.5 is needed in case the computed value of normA or normAinv isn't accurate enough
      }
      
      
      
      
      cout << "CDD2:: v.size() = " << v.size() << endl;
      cout << "loop: " << loops << endl;
      cout << "Calls of APPLY_QUARKLET: " << iter << endl << endl;
      v.COARSE(std::min((1-theta)*epsilon_k,1.0e-6), u_epsilon);
//      v.COARSE((1-theta)*epsilon_k, u_epsilon);
//      v.COARSE(1.0e-6, u_epsilon);
      
      
      
      
//      cout << "f:" << endl<< f << endl;
//      cout << "Av:" << endl << Av << endl;
//      cout << "v:" << endl << v << endl;
//      cout << "u:" << endl <<u_epsilon << endl;
      
      
      
    }
        
        
  }
  
  
  template <class PROBLEM>
  void richardson_QUARKLET_SOLVE(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon, CompressionStrategy strategy,
                        const double a, const double b)
  {


    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
//    typedef PQFrame<d,dT> FRAME1D;

//     Singularity1D_RHS_2<double> sing1D;
//     Singularity1D_2<double> exactSolution1D;
     
//      Point<2> origin;
//      origin[0] = 0.0;
//      origin[1] = 0.0;

//      CornerSingularity sing2D(origin, 0.5, 1.5);
//      CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    const unsigned int jmax = JMAX;
    const unsigned int pmax = PMAX;

    const double nu = P.norm_Ainv()*P.F_norm();
    typedef typename PROBLEM::Index Index;

    cout << "Rich_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    
    // ####### 1D #######
//    const double omega = 1;
//    const double omega = 0.2; //for p=4
//    const double omega = 0.4;
//    const double omega = 0.15;

    // ####### 2D #######
    // bad one
//    const double omega = 0.05;
    // good one
    const double omega = 0.25;

    //const double omega = 0.3;


    //const double omega = 2.0/2.47-0.5;
    cout << "Rich_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    //const double cond_A = P.norm_A() * P.norm_Ainv();
    //const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.8;
    cout << "Rich_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "Rich_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "Rich_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;
    map<double,double> weak_ell_tau_norms;
    
    bool exit = 0;
    unsigned int loops = 0;

    double epsilon_k = nu;
    InfiniteVector<double,Index> f, v, Av;

    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

//     EvaluateFrame<Basis1D,2,2> evalObj;
    
    
     while (epsilon_k > epsilon) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "Rich_SOLVE: epsilon_k=" << epsilon_k << endl;
      double eta = theta * epsilon_k / (6*omega*K);
//      cout << "eta= " << eta << endl;
      cout << "Rich_SOLVE: eta=" << eta << endl;
      P.RHS(eta, f);
      for (int j = 1; j <= 1/*K*/; j++) {
	APPLY_QUARKLET_COARSE(P, v, eta, Av, jmax, strategy, pmax, a, b, 1.0e-6);

	v += omega * (f - Av);

	++loops;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;
	

        // ############ output #############
	P.RHS(1.0e-6, f);
 	APPLY_QUARKLET(P, v, 1.0e-6, Av, jmax, strategy, pmax, a, b);
  	double residual_norm = l2_norm(f - Av);
 	double tmp1 = log10(residual_norm);
	cout << "current residual error ||f-Av||=" << residual_norm << endl;
#if 1	
	asymptotic[log10( (double)v.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;


#ifdef ONE_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./2.5); // d=3, dT=3
#endif
#ifdef TWO_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./1.5); // d=3, dT=3
#endif


	cout << "active indices: " << v.size() << endl;
 	cout << "loop: " << loops << endl;


	
	char name1[128];
	char name2[128];
	char name3[128];

#ifdef ONE_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results/rich1D_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results/rich1D_time_asymptotic_P_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_results/rich1D_weak_ell_tau_norms_P_jmax18_d", d, "_dT", dT, ".m");
#endif	

#ifdef TWO_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results_2D/rich2D_asymptotic_P_jmax", jmax,",_pmax", pmax,"_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results_2D/rich2D_time_asymptotic_P_jmax18_d", jmax,",_pmax", pmax,"_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_results_2D/rich2D_weak_ell_tau_norms_P_jmax18_d", jmax,",_pmax", pmax,"_d", d, "_dT", dT, ".m");
#endif

	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();

	std::ofstream os3(name3);
	matlab_output(weak_ell_tau_norms,os3);
	os3.close();

	// ############ end output #############	
#endif
	tstart = clock();

#ifdef ONE_D
	if (residual_norm < 3.1623e-04 || loops == 5000) 
#endif
#ifdef TWO_D
	if (residual_norm < 0.01 || loops == 5000) 
#endif
        {
	  u_epsilon = v;
	  exit = true;
          cout << "bin hier"  << endl;
	  break;
	}
      }
      //v.COARSE((1-theta)*epsilon_k, u_epsilon);
      
      if (exit)         
        break;
      
      
    } 
    
    

  }
    
    
  template <class PROBLEM>
  void richardson_QUARKLET_SOLVE(const PROBLEM& P, const double epsilon,
			InfiniteVector<double,int>& u_epsilon, const unsigned int maxiter,
                        CompressionStrategy strategy,
                        const double a, const double b, const double shrink)
  {


    const int d  = PRIMALORDER;
    const int dT = DUALORDER;
//    typedef PQFrame<d,dT> FRAME1D;

//     Singularity1D_RHS_2<double> sing1D;
//     Singularity1D_2<double> exactSolution1D;
     
//      Point<2> origin;
//      origin[0] = 0.0;
//      origin[1] = 0.0;

//      CornerSingularity sing2D(origin, 0.5, 1.5);
//      CornerSingularityRHS singRhs(origin, 0.5, 1.5);

    const unsigned int jmax = JMAX;
    const unsigned int pmax = PMAX;

    const double nu = P.norm_Ainv()*P.F_norm();
    

    cout << "Rich_SOLVE: nu=" << nu << endl;

    // compute optimal relaxation parameter omega
    //const double omega = 2.0 / (P.norm_A() + 1.0/P.norm_Ainv());
    
    // ####### 1D #######
//    const double omega = 1.;
//    const double omega = 0.1; //for p=4
//    const double omega = 0.4;
//    const double omega = 0.15;

    // ####### 2D #######
    // bad one
//    const double omega = 0.05;
    // good one
    const double omega = 0.25;

//    const double omega = 0.05;


    //const double omega = 2.0/2.47-0.5;
    cout << "Rich_SOLVE: omega=" << omega << endl;

    // compute spectral norm rho
    //const double cond_A = P.norm_A() * P.norm_Ainv();
    //const double rho = (cond_A - 1.0) / (cond_A + 1.0);
    const double rho = 0.8;
    cout << "Rich_SOLVE: rho=" << rho << endl;
    
    // desired error reduction factor theta < 1/3
    //const double theta = 2.0/7.0;
    const double theta = 0.333;
    cout << "Rich_SOLVE: theta=" << theta << endl;

    // compute minimal K such that 3*rho^K < theta
    const int K = (int) ceil(log(theta/3.0) / log(rho));
    cout << "Rich_SOLVE: K=" << K << endl;
    
    u_epsilon.clear();

    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;
    map<double,double> weak_ell_tau_norms;
    map<double,double> loop_asymptotic;
    
    bool exit = 0;
    unsigned int loops = 0;

    double epsilon_k = nu;
    InfiniteVector<double,int> f, v, Av;

    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

//     EvaluateFrame<Basis1D,2,2> evalObj;
    
    
     while (epsilon_k > epsilon) {
//    while (true) {
      epsilon_k *= 3*pow(rho, K) / theta;
      cout << "Rich_SOLVE: epsilon_k=" << epsilon_k << endl;
      double eta = theta * epsilon_k / (6*omega*K);
      cout << "eta= " << eta << endl;
      cout << "Rich_SOLVE: eta=" << eta << endl;
      P.RHS(eta, f);
      for (int j = 1; j <= 1/*K*/; j++) {
	APPLY_QUARKLET_COARSE(P, v, eta, Av, jmax, strategy, pmax, a, b, 1.0e-6);

	v += 0.5 * omega * (f - Av);
//        v.shrinkage(shrink);

	++loops;
	tend = clock();
	time += (double)(tend-tstart)/CLOCKS_PER_SEC;
	

        // ############ output #############
	P.RHS(1.0e-6, f);
 	APPLY_QUARKLET(P, v, 1.0e-6, Av, jmax, strategy, pmax, a, b);
  	double residual_norm = l2_norm(f - Av);
 	double tmp1 = log10(residual_norm);
	cout << "current residual error ||f-Av||=" << residual_norm << endl;
#if 1	
	asymptotic[log10( (double)v.size() )] = tmp1;
	time_asymptotic[log10(time)] = tmp1;
        loop_asymptotic[log10((double) loops)] = tmp1;


#ifdef ONE_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./2.5); // d=3, dT=3
#endif
#ifdef TWO_D
	weak_ell_tau_norms[loops] = v.weak_norm(1./1.5); // d=3, dT=3
#endif


	cout << "active indices: " << v.size() << endl;
 	cout << "loop: " << loops << endl;


	
	char name1[128];
	char name2[128];
	char name3[128];
        char name4[128];

#ifdef ONE_D
	sprintf(name1, "%s%d%s%d%s", "./Richardson_results/rich1D_asymptotic_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s", "./Richardson_results/rich1D_time_asymptotic_jmax18_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s", "./Richardson_results/rich1D_weak_ell_tau_norms_jmax18_d", d, "_dT", dT, ".m");
        sprintf(name4, "%s%d%s%d%s", "./Richardson_results/rich1D_loop_asymptotic_jmax18_d", d, "_dT", dT, ".m");
#endif	

#ifdef TWO_D
	sprintf(name1, "%s%d%s%d%s%d%s%d%s", "./Richardson_results_2D/rich2D_asymptotic_jmax", jmax,"_pmax", pmax,"_d", d, "_dT", dT, ".m");
	sprintf(name2, "%s%d%s%d%s%d%s%d%s", "./Richardson_results_2D/rich2D_time_asymptotic_jmax", jmax,"_pmax", pmax,"_d", d, "_dT", dT, ".m");
	sprintf(name3, "%s%d%s%d%s%d%s%d%s", "./Richardson_results_2D/rich2D_weak_ell_tau_norms_jmax", jmax,"_pmax", pmax,"_d", d, "_dT", dT, ".m");
        sprintf(name4, "%s%d%s%d%s%d%s%d%s", "./Richardson_results_2D/rich2D_loop_asymptotic_jmax", jmax,"_pmax", pmax,"_d", d, "_dT", dT, ".m");
#endif

	std::ofstream os1(name1);
	matlab_output(asymptotic,os1);
	os1.close();
	
	std::ofstream os2(name2);
	matlab_output(time_asymptotic,os2);
	os2.close();

	std::ofstream os3(name3);
	matlab_output(weak_ell_tau_norms,os3);
	os3.close();
        
        std::ofstream os4(name4);
	matlab_output(loop_asymptotic,os4);
	os4.close();

	// ############ end output #############	
#endif
	tstart = clock();

#ifdef ONE_D
	if (residual_norm < 1e-4 || loops == 5000 || (shrink!=0 && residual_norm < 0.29)) 
#endif
#ifdef TWO_D
	if (residual_norm < 0.01 || loops == maxiter) 
#endif
        {
	  u_epsilon = v;
	  exit = true;
	  break;
	}
      }
      //v.COARSE((1-theta)*epsilon_k, u_epsilon);
      
      if (exit)
	break;
      
    } 
    
    

  }
  
  
  
  
  
  
}
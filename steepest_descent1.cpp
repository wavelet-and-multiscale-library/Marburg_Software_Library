// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>
// #include <frame_index.h>
// #include <biharmonic_equation.h>
// #include <cube/cube_basis.h>
// #include <cube/cube_index.h>

using std::set;

namespace FrameTL{



  template <class PROBLEM>
  void steepest_descent1_SOLVE(const PROBLEM& P,  const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			       int jmax, 
			       InfiniteVector<double, typename PROBLEM::Index> rhs)
  {
    //u_epsilon wird mit der "Lösung" beschrieben

    /*

    Ax + b = 0 lösen

    ersetze dies durch Minimierung der Funktion

    F(v) := 1/2 <v,Av> + <b,v> 

    Setze 

    p = -(Ax + b) = -r (Richtungsvektor)

    b=-f

    t = - <p,r>/<p,Ap> (Minimalität in Richtung p)

    v' = v + t*p; (F(v') ist min) 

    F(v') = 1/2 * t^2 <p,Ap> + t* (p,r) + F(v)
    
    */
    typedef typename PROBLEM::Index Index;

    InfiniteVector<double, Index> w, r, b, Ap, p, f, Av, help;
    double t, eps, delta;
    unsigned int loops = 0;
    delta = 1.0;
    eps = 0.01;
    map<double,double> log_10_residual_norms;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;

    //double a_inv     = P.norm_Ainv();

    // double omega_i   = 90*1.0/10000;
    double omega_i = 1000;
    //cout << "||A^-1||=" << a_inv << endl;
    double time = 0.;
    clock_t tstart, tend;
    tstart = clock();

    APPLY_COARSE(P, w,omega_i, r, 0.00001, jmax, CDD1);
    r*=-1;
    cout << "r:" << r << endl;
     P.RHS(omega_i,b);
    r+=b;
    cout << "b:" << b << endl;
    //if(f.size()==0) 

    while(loops<40000){
      cout << "omega_i:" << omega_i << endl;
      loops++;
      cout << "Scheife:" << loops << endl;
      APPLY_COARSE(P, r ,omega_i*eps*eps, Ap , 0.00001, jmax, CDD1);
      t=(r*r)/(r*Ap);
      r*=t;
      w+=r;
      
      P.RHS(0.000000001,f);
      APPLY(P, w, .0, Av, jmax, CDD1);
      help = f-Av;
  
      omega_i*=0.99; 
      APPLY_COARSE(P, w, omega_i, r, 0.00001, jmax, CDD1);
      r*=-1;
      P.RHS(omega_i,b);
      r+=b;
//     cout << "loop: " << loops << " nu = " 
// 	 << nu_i << " epsilon = " << omega_i/((1+3.*mu)*a_inv) << endl;
//     cout << "xi: " << xi_i << endl; 
    
      double tmp = l2_norm(help);
      double tmp1 = log10(tmp);
      cout << "residual norm = " << tmp << endl;
      asymptotic[log10( (double)w.size() )] = tmp1;
      log_10_residual_norms[loops] = tmp1;
      
      tend = clock();
    
      time += ((double) (tend-tstart))/((double) CLOCKS_PER_SEC);
    
      time_asymptotic[log10(time)] = tmp1;
    cout << "active indices: " << w.size() << endl;
    
    
    tstart = clock();
    
    if (tmp < 0.00001 || loops ==10000) {
      u_epsilon = w;
      break;
    }
	  std::ofstream os3("steep1D_asymptotic_13_2000_b3.m");
	  matlab_output(asymptotic,os3);
	  os3.close();
	  
	  std::ofstream os4("steep1D_time_asymptotic_13_2000_b3.m");
	  matlab_output(time_asymptotic,os4);
	  os4.close();
    } //end while 
  }
}


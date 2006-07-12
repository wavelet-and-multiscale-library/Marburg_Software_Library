// implementation for steepest_descent.h

#include <cmath>
#include <set>
#include <utils/plot_tools.h>
#include <adaptive/apply.h>
#include <numerics/corner_singularity.h>
#include <frame_evaluate.h>

using std::set;

namespace FrameTL
{


/*!
*/
template<class VALUE = double>
class Singularity1D_RHS_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_RHS_2() {};
  virtual ~Singularity1D_RHS_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    return -sin(3.*M_PI*p[0])*9.*M_PI*M_PI - 4.;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};

/*!
  special function with steep gradients
  near the right end of the interval
*/
template<class VALUE = double>
class Singularity1D_2
  : public Function<1, VALUE>
{
public:
  Singularity1D_2() {};
  virtual ~Singularity1D_2() {};
  VALUE value(const Point<1>& p,
	      const unsigned int component = 0) const
  {
    if (0. <= p[0] && p[0] < 0.5)
      return -sin(3.*M_PI*p[0]) + 2.*p[0]*p[0];

    if (0.5 <= p[0] && p[0] <= 1.0)
      return -sin(3.*M_PI*p[0]) + 2.*(1-p[0])*(1-p[0]);

    return 0.;

  }
  
  void vector_value(const Point<1> &p,
		    Vector<VALUE>& values) const { ; }
  
};


  template <class PROBLEM>
  void  multiplicative_Schwarz_SOLVE(const PROBLEM& P,  const double epsilon,
				     InfiniteVector<double, typename PROBLEM::Index>& u_epsilon)
  {
    typedef DSBasis<4,6> Basis1D;
    //typedef PBasis<3,3> Basis1D;	

    Point<2> origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    
    CornerSingularity sing2D(origin, 0.5, 1.5);
    CornerSingularityRHS singRhs(origin, 0.5, 1.5);



    unsigned int loops = 0;
    const int jmax = 6;
    typedef typename PROBLEM::Index Index;

    double a_inv     = P.norm_Ainv();
    double omega_i   = a_inv*P.F_norm();

    InfiniteVector<double, Index> u_k, u_k_1_2, r, help, f, Av, r_exact;
    
    
    
    map<double,double> log_10_residual_norms;
    map<double,double> degrees_of_freedom;
    map<double,double> asymptotic;
    map<double,double> time_asymptotic;
    map<double,double> log_10_L2_error;

    bool exit = 0;
    double time = 0;
    clock_t tstart, tend;
    tstart = clock();

    EvaluateFrame<Basis1D,1,1> evalObj;

    double eta = 1.;


    const double alpha = 0.9;

    unsigned int global_iterations = 0;
    double tmp = 5.;
    while (tmp > 1.0e-3) {

      //approximate residual
      P.RHS(eta, f);
      APPLY_COARSE(P, u_k, eta, help, 0.00000001, jmax, CDD1);
      r = f - help;

      //approximate residual
      P.RHS(0., f);
      APPLY(P, u_k, 0., help, jmax, CDD1);
      r_exact = f - help;
      
      
      
      tmp = l2_norm(r_exact);
      cout << "residual norm = " << tmp  << endl;
      double tmp1 = log10(tmp);
      if (u_k.size() != 0)
	asymptotic[log10( (double)u_k.size() )] = tmp1;
      std::ofstream os3("steep_asymptotic_46_2D_ds_0507.m");
      matlab_output(asymptotic,os3);
      os3.close();

      // setup local index set
      set<Index> Lambda;
      r.support(Lambda);    
      
      set<Index> Lambda_1;
      typename set<Index>::const_iterator it = Lambda.begin();
      for (; it != Lambda.end(); ++it) {
	if ((*it).p() == 0) {
	  Lambda_1.insert(*it);
	}
      }

      // setup local stiffness matrix
      cout << "setting up full stiffness matrix..." << endl;
      SparseMatrix<double> A_Lambda_1;
      WaveletTL::setup_stiffness_matrix(P, Lambda_1, A_Lambda_1);


      // setup local right hand side
      InfiniteVector<double, Index> r_1;
      typename InfiniteVector<double, Index>::const_iterator it2 = r.begin();
      for (; it2 != r.end(); ++it2) {
	if ((it2.index()).p() == 0) {
	  r_1.set_coefficient(it2.index(),*it2);
	}
      }



      cout << "setting up full right hand side..." << endl;
      Vector<double> F_Lambda_1(r_1.size());
      unsigned int id = 0;
      typename InfiniteVector<double, Index>::const_iterator it3 = r_1.begin();
      for (; it3 != r_1.end(); ++it3, ++id) {
	F_Lambda_1[id] = *it3;
      }

      // compute approximation to local problem
      Vector<double> xk(Lambda_1.size());
                 cout << "r_1 size = " << r_1.size() << endl;
      unsigned int iterations = 0;
      CG(A_Lambda_1, F_Lambda_1, xk, 0.0001, 300, iterations);
      cout << "CG done!!!!" << endl;           

      help.clear();

      id = 0;
      for (typename set<Index>::const_iterator it = Lambda_1.begin(), itend = Lambda_1.end();
	   it != itend; ++it, ++id)
	help.set_coefficient(*it, xk[id]);

      u_k_1_2 = u_k + help;
      //u_k = u_k + alpha*help;



      // ############## first half step completed #####################
      

      //approximate residual
      P.RHS(eta, f);
      APPLY_COARSE(P, u_k_1_2, eta, help, 0.00000001, jmax, CDD1);
      r = f - help;


      // setup local index set
      Lambda.clear();
      r.support(Lambda);    
      
//       cout << "r " << r.size() << endl; 
//       cout << "Lambda " << Lambda.size() << endl; 

      set<Index> Lambda_2;
      it = Lambda.begin();
      for (; it != Lambda.end(); ++it) {
	//	cout << "patch = " << (*it).p() << endl;
	if ((*it).p() == 1) {
	  Lambda_2.insert(*it);
	}
      }

      cout << "lambda_2 " << Lambda_2.size() << endl; 

      // setup local stiffness matrix
      cout << "setting up full stiffness matrix..." << endl;
      SparseMatrix<double> A_Lambda_2;

      WaveletTL::setup_stiffness_matrix(P, Lambda_2, A_Lambda_2);

  
      // setup local right hand side
      InfiniteVector<double, Index> r_2;
      it2 = r.begin();
      for (; it2 != r.end(); ++it2) {
	if ((it2.index()).p() == 1) {
	  r_2.set_coefficient(it2.index(),*it2);
	}
      }

      
      
      cout << "setting up full right hand side..." << endl;
      Vector<double> F_Lambda_2(r_2.size());
      id = 0;
      it3 = r_2.begin();
      for (; it3 != r_2.end(); ++it3, ++id)
	F_Lambda_2[id] = *it3;

      // compute approximation to local problem
      iterations = 0;
     
      Vector<double> yk(Lambda_2.size());
      CG(A_Lambda_2, F_Lambda_2, yk, 0.0001, 300, iterations);
      cout << "CG done!!!!" << endl;
      

      help.clear();

      id = 0;
      for (typename set<Index>::const_iterator it = Lambda_2.begin(), itend = Lambda_2.end();
	   it != itend; ++it, ++id)
	help.set_coefficient(*it, yk[id]);

      u_k = u_k_1_2 + help;
      //u_k = u_k + alpha*help;

      global_iterations++;

      eta *= 0.85;



      u_epsilon = u_k;

      cout << ++loops << " loops completed" << endl;

//       if (global_iterations % 1 == 0) {
// 	u_epsilon = u_k;
// 	u_epsilon.scale(&P,-1);
// 	char filename1[50];
// 	char filename2[50];
	
// 	sprintf(filename1, "%s%d%s%d%s", "approx_sol_steep22_2D_out_", loops, "_nactive_", u_k.size(),".m");
// 	sprintf(filename2, "%s%d%s%d%s", "error_steep_2D_out_", loops, "_nactive_", u_k.size(),".m");
// 	cout << "...plotting approximate solution" << endl;
// 	Array1D<SampledMapping<2> > U = evalObj.evaluate(P.basis(), u_epsilon, true, 5);//expand in primal basis
	
// 	std::ofstream ofs5(filename1);
// 	matlab_output(ofs5,U);
// 	ofs5.close();



// // 	Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(P.basis(), u_epsilon, exactSolution1D, 6);
// // 	cout << "...plotting error" << endl;
//       }
  }

 




  }
}

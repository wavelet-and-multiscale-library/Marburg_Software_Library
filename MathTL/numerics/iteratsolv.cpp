// implementation of iterative solvers

#include <numerics/preconditioner.h>
#include <cmath>
#include <set>
#include <map>
#include <utils/plot_tools.h>

namespace MathTL
{
  template <class VECTOR, class MATRIX>
  void Richardson(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
		  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // first residual
    VECTOR rk(A.row_dimension(), false);
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	xk.add(-omega, rk);

	A.apply(xk, rk);
	rk -= b;
	error = l2_norm(rk);
      }
  }
  
  template <class VECTOR, class MATRIX>
  void Richardson_sparse(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega, const double mu,
		  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // first residual
    VECTOR rk(A.row_dimension(), false);
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);
    double shrink = mu;
    int dof;
    map<double,double> residual_norm;
    map<double,double> degrees_of_freedom;

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	xk.add(-omega, rk);
//        if(iterations==5000)
//            cout << "x vor shrinkage: " << endl << xk << endl;
        xk.shrinkage(shrink);
//        if(shrink > 0.01*mu && iterations >50)
//            shrink*=9./10.;
//        if(iterations==5000)
//            cout << "x nach shrinkage: " << endl << xk << endl;

	A.apply(xk, rk);
	rk -= b;	
        //////////////////////Output////////////////////////
        error = l2_norm(rk);        
        double tmp1 = log10(error);
        residual_norm[log10( (double)iterations )]=tmp1;
        
//        cout << "iteration: " << iterations <<  ", residual norm = " << error << endl;
        
        
            for(unsigned int i = 0;i<xk.size();i++){
                if(xk[i]!=0)
                    dof++;            
            }
//        cout << "Degrees of freedom: " << dof << endl;
        degrees_of_freedom[log10( (double)iterations )]=log10( (double)dof );
        dof=0;
      }
    // perform matlab output
    std::ofstream os1("residual_plot.m");
    matlab_output(residual_norm,os1);
    os1.close();

    std::ofstream os2("degrees_of_freedom_plot.m");
    matlab_output(degrees_of_freedom,os2);
    os2.close();
  }
  
  template <class VECTOR, class MATRIX>
  void Landweber(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
		  const double tol, const unsigned int maxiter, unsigned int& iterations, const double& mu)
  {
    // first residual
    VECTOR rk(A.row_dimension(), false), Ark(A.row_dimension(), false);
    
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	
        A.apply(rk,Ark);        
        xk.add(-omega, Ark);
        xk.shrinkage(mu);


	A.apply(xk, rk);
	rk -= b;
	error = l2_norm(rk);
        cout << "error: " << error << endl;
      }
  }
  

  template <class VECTOR, class MATRIX>
  void Jacobi(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	      const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());
    
    // first residual
    VECTOR rk(A.row_dimension(), false);
    A.apply(xk, rk);
    rk -= b; // rk=A*xk-b
    double error = l2_norm(rk);

    for (iterations = 1; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename VECTOR::size_type i(0); i < n; i++)
//	  xk[i] -= rk[i] / A(i, i);
          xk[i] -= rk[i] / A.get_entry(i, i);

	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
        cout <<"Iteration: " << iterations << ", Error: " << error << endl;
      }
  }

  template <class VECTOR, class MATRIX>
  void GaussSeidel(const MATRIX &A, const VECTOR &b, VECTOR &xk,
		   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());

    VECTOR rk(A.row_dimension(), false);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    xk[i] = b[i];
	    for (typename MATRIX::size_type j(0); j < n; j++)
	      if (j != i)
		xk[i] -= A.get_entry(i, j) * xk[j];
	    xk[i] /= A.get_entry(i, i);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());

    VECTOR rk(A.row_dimension(), false);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    double sigma(b[i]);
	    for (typename MATRIX::size_type j(0); j < n; j++)
	      if (j != i)
		sigma -= A(i, j) * xk[j];
	    sigma /= A(i, i);
	    xk[i] += omega * (sigma - xk[i]);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  void SSOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    const typename MATRIX::size_type n(A.column_dimension());
    
    VECTOR rk(A.row_dimension(), false), xkpuf(xk);
    double error = 2*tol;
    
    for (iterations = 0; iterations <= maxiter && error > tol; iterations++)
      {
	for (typename MATRIX::size_type i(0); i < n; i++)
	  {
	    double sigma(0);
	    for (typename MATRIX::size_type j(0); j+1 <= i; j++)
	      sigma += A(i, j) * xkpuf[j];
	    for (typename MATRIX::size_type j(i+1); j < n; j++)
	      sigma += A(i, j) * xk[j];
	    sigma = (b[i]-sigma) / A(i, i);
	    xkpuf[i] = xk[i] + omega * (sigma - xk[i]);
	  }

	for (typename MATRIX::size_type i(n); i >= 1; i--)
	  {
	    double sigma(0);
	    for (typename MATRIX::size_type j(0); j+1 <= i-1; j++)
	      sigma += A(i-1, j) * xkpuf[j];
	    for (typename MATRIX::size_type j(i); j < n; j++)
	      sigma += A(i-1, j) * xk[j];
	    sigma = (b[i-1]-sigma) / A(i-1, i-1);
	    xk[i-1] = xkpuf[i-1] + omega * (sigma - xkpuf[i-1]);
	  }
	
	A.apply(xk, rk);
	rk -= b;

	error = l2_norm(rk);
      }
  }

  template <class VECTOR, class MATRIX>
  inline
  bool CG(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
//    cout << b << endl;
    IdentityPreconditioner<MATRIX,VECTOR> I(A);
    return PCG(A, b, I, xk, tol, maxiter, iterations);

//     iterations = 0;
    
//     VECTOR gk(A.row_dimension(), false),
//       dk(A.row_dimension(), false),
//       Adk(A.row_dimension(), false);
    
//     // first residual
//     A.apply(xk, gk);
//     gk -= b; // gk = A*xk-b
    
//     // first direction
//     dk = gk;
//     dk.scale(-1.0); // dk = -gk
    
//     double error = l2_norm(gk);
//     for (iterations = 1; error > tol && iterations <= maxiter; iterations++)
//       {
// 	const double denom = gk * gk;
// 	A.apply(dk, Adk);
// 	const double alpha = denom / (dk * Adk);
// 	xk.add(alpha, dk);
// 	gk.add(alpha, Adk);
// 	error = l2_norm(gk);
// 	const double beta = (gk * gk) / denom;
// 	dk *= beta; 
// 	dk -= gk;
//       }
    
//     return (iterations <= maxiter);
  }

  template <class VECTOR, class MATRIX, class PREC>
  bool PCG(const MATRIX &A, const VECTOR &b, const PREC& P, VECTOR &xk,
	   const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
//    cout << b << endl;
    // see: "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods"

    VECTOR rk(A.row_dimension(), false),
      zk(A.row_dimension(), false),
      pk(A.row_dimension(), false),
      Apk(A.row_dimension(), false);

//    cout << A << endl;
    //A.matlab_output("AMAT", "A", 1);
    // first (negative) residual
    A.apply(xk, rk);
    rk.subtract(b);
    const double normr0 = l2_norm_sqr(rk);
//    cout << "normr0" << endl;
    double normrk = normr0, rhok = 0, oldrhok = 0;
    for (iterations = 1; normrk/normr0 > tol*tol && iterations <= maxiter; iterations++)
      //for (iterations = 1; normrk > tol*tol && iterations <= maxiter; iterations++)
      {
	P.apply_preconditioner(rk, zk);
	rhok = rk * zk;

	if (iterations == 1) // TODO: shift this case upwards!
	  pk = zk;
	else
	  pk.sadd(rhok/oldrhok, zk);

	A.apply(pk, Apk);
	const double alpha = rhok/(pk*Apk);
	xk.add(-alpha,  pk);
	rk.add(-alpha, Apk);
	normrk = l2_norm_sqr(rk);

//	cout << "normrk = " << sqrt(normrk) << endl;
	oldrhok = rhok;
      }

    return (iterations <= maxiter);
  }
}

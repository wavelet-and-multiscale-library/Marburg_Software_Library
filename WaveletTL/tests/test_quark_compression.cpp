//Gives an Estimation for the norm of the Laplacian of the compresses matrices A_J


#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <algebra/vector.h>
#include <Rd/r_q_index.h>
#include <Rd/r_index.h>
#include <Rd/quarklet_frame.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/i_q_index.h>
#include <interval/periodic_frame.h>
#include <interval/periodic.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/matrix.h>
#include <algebra/matrix_norms.h>
#include <algebra/sparse_matrix.h>
#include <time.h>
#include <geometry/grid.h>


using namespace std;
using namespace WaveletTL;
using namespace MathTL;


int main(){
    
    
    
  //Testreihe evaluate
    /*
  const int d  = 3;
  const int dt = 3;
  typedef PeriodicFrame<QuarkletFrame<d,dt> > Frame;
  typedef PeriodicBasis<CDFBasis<d,dt> > Basis;
  typedef Frame::Index Index;
  typedef Basis::Index BIndex;
  Frame frame;
  Basis basis;
  Index lambda3(0,3,1,0), mu3(0,3,1,14);
  BIndex lambda4(3,1,2), mu4(6,1,2);
  cout << frame.integrate(1,lambda3, mu3) << endl;
  cout << basis.integrate(lambda4, mu4) << endl;
   */ 
   /* cout << "Test CDFBasis class..." << endl;

  const int d = 3;
  const int dt = 3;

  CDFBasis<d,dt> basis;

  cout << "A CDF basis with d=" << d << " and dt=" << dt
       << " has the following data:" << endl;
  cout << "  + primal mask:" << endl;
  cout << basis.a() << endl;
  cout << "  + dual mask:" << endl;
  cout << basis.aT() << endl;
  double b[8];
  Matrix<double> phiwerte(3,3, "0.5 -1 0.5 1 1 -2 0.5 0.5 2");
  cout << phiwerte << endl;

  b[0]=0;
  b[1]=1;
  cout << "  + wavelet coefficients:" << endl;
  int i = 2;
  for (int k = 1-CDFBasis<d,dt>::dual_mask::end();
       k <= 1-CDFBasis<d,dt>::dual_mask::begin();
       k++) {
      b[i] = basis.b(k);
      i++;
    cout << "b(" << k << ")=" << basis.b(k) << endl;
  }
  b[i]=0;
  b[i+1]=0;
  Vector<double> temp(3);
  Matrix<double> solution(3,10);
  for(int k= 0; k<10; k++){
      solution(0,k)=1;
      solution(1,k)=2;
      solution(2,k)=3;
  }*/
  
    /*################################################################
####################################################################
NUMERISCHE EXPERIMENTE: MATRIXNORM*/
    
#if 1
 ofstream ausgabe("../../Desktop/numerische_exp_13_08.txt", ios::out|ios::app);
 
 
 //Setup the frame
  const int d  = 3;
  const int dt = 3;
  
  typedef PeriodicFrame<QuarkletFrame<d,dt> > Frame;
  typedef Frame::Index Index;
  Frame frame, basis2;
  
  
  
  
  // Parameter
  const int j0=frame.j0();
  const int jmax = 11;
  const int pmax = 6;
  const int delta = 3;
  const double a = 2;
  const double b = 2;
  const int maxJ = 25;
  //const int tau = delta + 3 - d;
  
  ausgabe << endl << "Periodic Quarklet frame of order " << d <<
    " with " << dt << " vanishing moments, a=" << a << ", b=" << b << ", tau=" << delta+3-d<< endl;
#endif
#if 0  
  //Testreihe evaluate
  ofstream ausgabe4("../../Arbeitsfläche/numerische_testreihe_evaluate_slow.txt");
  QuarkletFrame<3,3> myquark;
  CDFBasis<3,3> mycdf;
  RIndex lambda4(0,1,0);
  RQIndex lambda3(1,0,1,0);
  double x = -2;
  for(int i=0; i<=50; i++){
    ausgabe4 << "(" << x <<"|" << myquark.evaluate(1,lambda3, x) << ")" << endl;
    x+=0.1;
  }
#endif  
  
#if 0
  //Testreihe Integrate
  ofstream ausgabe2("../../Arbeitsfläche/testreihe_integrate_fast.txt");
  int p2 = 0;
  Index lambda2, mu2;
  lambda2 = frame.first_generator(frame.j0());
  //lambda2 = RQIndex(1,4,1,0);
  mu2 = frame.first_generator(frame.j0());
  for(int i=0; i<100 ; i++){
  double temp3;
  temp3 = frame.integrate(1,lambda2, mu2);
  ausgabe2 << lambda2 << ", " << mu2 << ", " << temp3 << endl;
  
  if(mu2.k() == (1<<jmax)-1){
          ++p2;
          mu2 = Index(p2, j0, 0, 0);
          }
   else
      ++mu2;
  
  }
  ausgabe2.close();
#endif
#if 1
  //Setup the matrix
  int matdim = (pmax+1) * pow(2,jmax+1);
  ausgabe << "jmax=" << jmax
	  << ", pmax="<< pmax << ", Matrixdimension: " << matdim << endl;
  SparseMatrix<double> ARes, tempmat;
  tempmat.resize(matdim, matdim);
  Index lambda, mu;
  Vector<double> xval(maxJ+1), yval(maxJ+1);
  
  double time1 = 0.0, tstart;
  tstart = clock();
  
  //helper values
  double ln2rez = 1./log(2);
  Vector<double> log_table(pmax+1);
  for(int i(0); i<=pmax; ++i)
      log_table[i] = log(i+1);
  
 
#endif 
#if 1
  for(int J(0); J<=maxJ;J++){
  //int J = 0;
  ausgabe << "J=" << J << ": ";
  ARes.resize(matdim, matdim);
  lambda = frame.first_generator(frame.j0());
  int besetzung = 0;
  int p1 = 0;
  int p2 = 0;

  //Initialization of the marix;
  for (int row(0); row<matdim; row++){
    cout << row << endl;
    cout << lambda << endl;
    p2 = p1;
    mu = lambda;
    for(int column(row); column<matdim; column++){
      //If Condition is fullfilled, the entry is 0
      if(a * log_table[abs(lambda.p()-mu.p())] * ln2rez + b * abs(lambda.j()-mu.j()) <= J) ;
      else{
         //Preconditioning factors
	 double factor;
         factor = ldexp(1.0,-lambda.j()-mu.j()) * pow(lambda.p()+1, -2-delta) * pow(mu.p()+1, -2-delta);
	 //Evaluate entry
         double tempone = frame.integrate(1,lambda, mu);
	 double temp = tempone * factor;
	 if(abs(temp)>10e-12){
         //if(temp!=0){
	  //Set entry
            if(row==column){
	    ARes.set_entry(row, column, temp);
            //tempmat.set_entry(row, column, temp);
	    besetzung++;
             }
	    else{ 
             ARes.set_entry(row, column, temp);
	     ARes.set_entry(column, row, temp);
             //tempmat.set_entry(row, column, temp);
	     //tempmat.set_entry(column, row, temp);
	     besetzung+=2;
             }
            
          }
        }
      //Determine the next mu
      if(mu.k() == (1<<jmax)-1 && mu.e() == 1){
          ++p2;
          mu = Index(p2, j0, 0, 0);
          }
          else
          ++mu;
      
      
    }
    //Determine the next lambda
    if(lambda.k() == (1<<jmax)-1 && lambda.e() == 1){
        ++p1;
	lambda = Index(p1, j0, 0, 0);
          }
          else
          ++lambda;   
  
  }
  //cout << ARes << endl;
  //Manipulation of the matrix
  //ARes(matdim,matdim) = -0.03;

  //Compute Matrixnorm
  cout << "- calculating the maximal absolute eigenvalue of ARes ..." << endl;
  unsigned int iterations(0);
  double lambdamin(0), lambdamax(0);
  LanczosIteration(ARes, 1e-6, lambdamin, lambdamax, matdim, iterations);
  double norm = max(abs(lambdamin), abs(lambdamax));
  ausgabe << "LanczosIteration yields (after " << iterations << " iterations) ||A-A_J||_l_2 = " << norm  << endl;
  xval(J)=J;
  yval(J)=norm;

  double anzeintr = pow(matdim,2);
  ausgabe << "Besetzungsgrad der Matrix=" << ARes.size()/anzeintr << endl;
  
  time1 += clock() - tstart;     // end
 
  time1 = time1/CLOCKS_PER_SEC;  // rescale to seconds

  ausgabe << "  time(after step J=" << J << ") = " << time1 << " sec." << endl;
  }
  
#if 0  
for(int J(1); J<=maxJ;J++){
  ausgabe << "J=" << J << ": ";
  ARes.resize(matdim, matdim);
  lambda = frame.first_generator(frame.j0());
  int besetzung = 0;
  int p1 = 0;
  int p2 = 0;

  //Initialization of the marix;
  for (int row(0); row<matdim; row++){
    cout << row << endl;
    cout << lambda << endl;
    p2 = p1;
    mu = lambda;
    for(int column(row); column<matdim; column++){
      //If Condition is fullfilled, the entry is 0
      if(a * log_table[abs(lambda.p()-mu.p())] * ln2rez + b * abs(lambda.j()-mu.j()) <= J) ;
      else{
          //Set entry
          double temp2 = tempmat.get_entry(row, column);
          if(temp2 != 0){
            if(row==column){
	    ARes.set_entry(row, column, temp2);
            besetzung++;
             }
	    else{ 
             ARes.set_entry(row, column, temp2);
	     ARes.set_entry(column, row, temp2);
             besetzung+=2;
             }
          }
        }
      //Determine the next mu
      if(mu.k() == (1<<jmax)-1){
          ++p2;
          mu = Index(p2, j0, 0, 0);
          }
          else
          ++mu;
      
      
    }
    //Determine the next lambda
    if(lambda.k() == (1<<jmax)-1){
        ++p1;
	lambda = Index(p1, j0, 0, 0);
          }
          else
          ++lambda;   
  
  }
  //cout << ARes << endl;
  //Manipulation of the matrix
  //ARes(matdim,matdim) = -0.03;

  //Compute Matrixnorm
  cout << "- calculating the maximal absolute eigenvalue of ARes ..." << endl;
  unsigned int iterations(0);
  double lambdamin(0), lambdamax(0);
  LanczosIteration(ARes, 1e-6, lambdamin, lambdamax, matdim, iterations);
  double norm = max(abs(lambdamin), abs(lambdamax));
  ausgabe << "LanczosIteration yields (after " << iterations << " iterations) ||A-A_J||_l_2 = " << norm  << endl;
  xval(J)=J;
  yval(J)=norm;

  double anzeintr = pow(matdim,2);
  ausgabe << "Besetzungsgrad der Matrix=" << ARes.size()/anzeintr << endl;
  
  time1 += clock() - tstart;     // end
 
  time1 = time1/CLOCKS_PER_SEC;  // rescale to seconds

  ausgabe << "  time(after step J=" << J << ") = " << time1 << " sec." << endl;
}
#endif
  ausgabe.close();

  ofstream ausgabe2("../../Desktop/plotter13_08.m");
  ausgabe2 << "x=" << xval << ";" << endl;
  ausgabe2 << "y=" << yval << ";" << endl;
  ausgabe2 << "figure;\nsemilogy(x,y);"
              << "title('myplot');" << endl;
  ausgabe2.close();
  
  #endif
  #if 0
  int k1, k2;
  basis2.support(RQIndex(0,3,1,0), k1, k2);
  Index myindex(0,1,0,0);
  cout << basis2.first_generator(basis2.j0()) << endl;
  cout << pow(2,-4) * k1 << "," << pow(2, -4) * k2 << endl;

  cout << "- j0=" << basis2.j0() << endl;
  cout << "- the default quarklet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_q_generator(&basis2, basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_q_generator(&basis2, basis2.j0()) << endl;
  cout << "- leftmost quarklet on the coarsest level: " << first_quarklet(&basis2, basis2.j0()) << endl;
  cout << "- rightmost quarklet on the coarsest level: " << last_quarklet(&basis2, basis2.j0()) << endl;


  QuarkletFrame<3, 3>  qfr2;
  PeriodicFrame<QuarkletFrame<3, 3> > qfr;
  //cout << qfr2.evaluate(0, RQIndex(0,2,0,0), 0.0) << endl;

  PeriodicBasis<CDFBasis<3, 3> > cdfb;
  //cout << cdfb.evaluate(0, RIndex(2,0,0), 0.0) << endl;
  const int p = 1;
  const int j = 2;
  const int e = 0;
  const int k = 1;
  //cdfb.support(RIndex(j,e,k),k1,k2);
  //cout << cdfb.j0() << endl;
  //cout << pow(2,-j-e) * k1 << ", " << pow(2,-j-e) * k2 << endl;
  ofstream raus("../../Arbeitsfläche/perfkt.m");
  qfr.evaluate(Index(p,j,e,k), 10).matlab_output(raus);
  raus << "figure;\nplot(x,y);"
              << "title('Superplot');" << endl;
  raus.close();
  //cout << "Ein paar Werte: " << qfr2.evaluate(0,Index(p,j,e,k), 0.15) << endl;
  cout << qfr.evaluate(1, Index(p,j,e,k), 0.5) << endl;
  cout << cdfb.evaluate(1, RIndex(j,e,k), 0.8) << endl;
  //cout << cdfb.integrate(RIndex(j,e,k), RIndex(j,e,k+1)) << endl;
  //cout << qfr.integrate(0, Index(p,j,e,k), Index(p,j,e,k)) << endl;
    #endif
   

}


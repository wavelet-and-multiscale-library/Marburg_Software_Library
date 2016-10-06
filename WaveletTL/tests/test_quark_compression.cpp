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
#include <interval/pq_frame.h>
#include <interval/pq_expansion.h>
#include <interval/periodic.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/matrix.h>
#include <algebra/matrix_norms.h>
#include <algebra/sparse_matrix.h>
#include <time.h>
#include <geometry/grid.h>

#define PERIODIC
#undef PRIMBSQUARK


using namespace std;
using namespace WaveletTL;
using namespace MathTL;


int main(){
  //NUMERISCHE EXPERIMENTE: MATRIXNORM*/
    
#if 1
 ofstream ausgabe("numerische_exp.txt", ios::out|ios::app);
 
 
 //Setup the frame
  const int d  = 3;
  const int dT = 3;
  
#ifdef PERIODIC
  typedef PeriodicFrame<QuarkletFrame<d,dT> > Frame;
  Frame frame;
#endif
#ifdef PRIMBSQUARK
  typedef PQFrame<d,dT> Frame;
  Frame frame(true,true,false);
#endif
  typedef Frame::Index Index;
  
  
  const int j0=frame.j0();
  const int jmax = 12;
  const int pmax = 6;
  const int delta = 3;
  const double a = 2;
  const double b = 2;
  const int maxJ = 25;
  //const int tau = delta + 3 - d;
  
  ausgabe << endl << "Periodic Quarklet frame of order " << d <<
    " with " << dT << " vanishing moments, a=" << a << ", b=" << b << ", tau=" << delta+3-d<< endl;
#endif

#if 1
  //Determine the index set
  int p = 0;
  set<Index> Lambda;
  
  for (Index lambda = frame.first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == frame.last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == frame.last_wavelet(jmax,p)){
            ++p;
            lambda = frame.first_generator(j0,p);
        }
        else
            ++lambda;
      }
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it){
       cout << *it << endl;
  }
  //Setup the matrix
//  int matdim = (pmax+1) * pow(2,jmax+1);
  int matdim = Lambda.size();
  ausgabe << "jmax=" << jmax
	  << ", pmax="<< pmax << ", Matrixdimension: " << matdim << endl;
  SparseMatrix<double> ARes, tempmat;
  tempmat.resize(matdim, matdim);
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
  for(int J(0); J<=2;J++){
  ausgabe << "J=" << J << ": ";
  ARes.resize(matdim, matdim);
  int besetzung = 0;
  int row = 0;
  int column = 0;
  
  //Initialization of the marix;
  //noch einbauen @PHK
  for (typename std::set<Index>::const_iterator it1(Lambda.begin()), it1end(Lambda.end());
	 it1 != it1end; ++it1, ++row){
    cout << "row: " << *it1 << endl;
    column=0;
    for (typename std::set<Index>::const_iterator it2(Lambda.begin()), it2end(Lambda.end());
	 it2 != it2end; ++it2, ++column){
            cout << "column: " << *it2 << endl;
//    //If Condition is fullfilled, the entry is 0
      if(a * log_table[abs((*it1).p()-(*it2).p())] * ln2rez + b * abs((*it1).j()-(*it2).j()) <= J) ;
      else{
         //Preconditioning factors
	 double factor;
         factor = ldexp(1.0,-(*it1).j()-(*it2).j()) * pow((*it1).p()+1, -2-delta) * pow((*it2).p()+1, -2-delta);
	 //Evaluate entry
#ifdef PERIODIC
         double tempone = frame.integrate(1,*it1, *it2);
#endif
#ifdef PRIMBSQUARK
         double tempone = integrate(frame,*it1, *it2, 1);
#endif
         double temp = tempone * factor;
	 if(abs(temp)>10e-12){
            ARes.set_entry(row, column, temp);
            besetzung+=1;
         }
            
      }
    }
      
  }
  
  

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
  

  ausgabe.close();

  ofstream ausgabe2("plotter_matrixnorm.m");
  ausgabe2 << "x=" << xval << ";" << endl;
  ausgabe2 << "y=" << yval << ";" << endl;
  ausgabe2 << "figure;\nsemilogy(x,y);"
              << "title('myplot');" << endl;
  ausgabe2.close();
  
  #endif
  


}


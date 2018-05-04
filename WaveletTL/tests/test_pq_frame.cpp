#undef BASIS
#define FRAME

#define DYADIC
#undef TRIVIAL

#define _DIM 1

#ifdef FRAME
#include <interval/pq_frame.h>
#include <interval/pq_support.h>
#include <interval/pq_expansion.h>
#include <interval/pq_evaluate.h>
#include <interval/periodic_frame.h>
#include <Rd/quarklet_frame.h>
#include <galerkin/cached_quarklet_problem.h>
#endif

#include <iostream>
#include <interval/p_basis.h>
#include <interval/p_evaluate.h>
#include <interval/p_expansion.h>
#include <interval/ds_basis.h>
#include <numerics/bvp.h>
#include <cube/cube_basis.h>
#include <galerkin/cube_equation.h>
#include <geometry/grid.h>
#include <numerics/sturm_bvp.h>
#include <galerkin/cached_problem.h>
#include <adaptive/cdd2.h>
#include <utils/function.h>
#include <algebra/vector.h>
//#include <Example_CDD2/TestProblem.h>
#include <galerkin/sturm_equation.h>

using WaveletTL::CubeEquation;
using WaveletTL::CubeBasis;
using MathTL::IdentityBVP;


using namespace std;
using namespace WaveletTL;
using namespace MathTL;

template <unsigned int N>
class TestProblem : public SimpleSturmBVP
{
    public:
        double p(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 0;
                    break;
                case 2:
                    return 1;
                    break;
                case 3:
                    return 1;
                    break;
                case 4:
                    return 0;
                    break;
                case 5:
                    return 1;
                    break;   
                default:
                    return 0;
                    break;
            }
        }

        double p_prime(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 0;
                    break;
                case 2:
                    return 0;
                    break;
                case 3:
                    return 0;
                    break;
                case 4:
                    return 0;
                    break;
                case 5:
                    return 0;
                    break;
                default:
                    return 0;
                    break;
            }
        }

        double q(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 1;
                    break;
                case 2:
                    return 0;
                    break;
                case 3:
                    return 0;
                    break;
                case 4:
                    return 3.1;
                    break;
                case 5:
                    return 0;
                    break;    
                default:
                    return 0;
                    break;
            }
        }

        double g(const double t) const
 
        {
            switch(N)
            {
                case 1:
                    return t*(1-t);
                    break;
                case 2:
                    //return exp(-100*(t-0.5)*(t-0.5));
                    return (200-(200*t-100)*(200*t-100))*exp(-100*(t-0.5)*(t-0.5));
                    break;
                case 3:
                    return ( -100*exp(5*t)*(1-(exp(5*t)-1)/(exp(5.)-1)) /
                             (exp(5.)-1)+200*exp(10*t)/((exp(5.)-1) *
                             (exp(5.)-1))+100*(exp(5*t)-1)*exp(5*t) /
                             ((exp(5.)-1)*(exp(5.)-1)) );
                    break;
                case 4:
                    return 4;
                    break;
                case 5:
                    return 2;
                    break;
                default:
                    return 0;
                    break;
            }
        }

        string toString() const
        {
            switch(N)
            {
                case 1:
                    return "Gramian Test Problem";
                    break;
                case 2:
                    return "y(t) = exp(-100*(x-0.5)^2), -y''(t) = (200-(200x-100)^2)*exp(-100*(x-0.5)^2), y(0) = y(1) = 0";
                    break;
                case 3:
                    return "1D example from [BBCCDDU]";
                    break;
                case 4:
                    return "konstante Funktion";
                    break;
                case 5:
                    return "y(t) = t*(1-t), -y''(t)=2, y(0)=y(1)=0";
                    break;default:
                    return "TestProblem: N not defined.";
                    break;
            }
        }

  bool bc_left() const { return true; }   /* boundary conditions left: true */
  bool bc_right() const { return true; }  /* boundary conditions right: true */
};


//class NewFunction : public Function<1>
//{
//    public:
//        inline double value(const Point<1>& p,
//                            const unsigned int component = 0) const
//        {
//            double t = 0;
//            if(0<=p[0] && p[0]<=1){
//                //t=p[0]*(1-p[0]);
//                t=1;
//            }
//            return t;
//            //return p[0]*p[0]*p[0]-1.5*p[0]*p[0]+0.5*p[0];
//        }
//
//        void vector_value(const Point<1> &p,
//                          Vector<double>& values) const
//        {
//            values.resize(1, false);
//            values[0] = value(p);
//        }
//};



int main()
{
  cout << "Testing wavelet frame from [DFKRS] ..." << endl;
  

  const int d  = 3;
  const int dT = 3;
  
//  typedef PeriodicFrame<QuarkletFrame<d,dT> >PerFrame;
//  typedef PerFrame::Index PerFrIndex;
//  
//  
//  
//  
//  
//  PerFrame perframe;main.cpp:409:5: error: ‘sm3’ was not declared in this scope
//  
//  
//  FrIndex lambda (1,2,0,2, &frame);
//  PerFrIndex mu(1,2,0,2);
//  cout << lambda << endl;
//  cout<< (1<<2)-ell1<d>()-d << endl;
  
  
//  double epsilon1 = 1e-6;
  const int jmax = 3;
  
  TestProblem<2> testproblem;
  
  
#ifdef FRAME
  typedef PQFrame<d,dT> Frame;
  typedef Frame::Index Index;
  Frame frame(false,false, true);
  const int pmax = 1;
  cout << "hier0" << endl;
  frame.set_jpmax(jmax,pmax);
  cout << "hier1" << endl;
  SturmEquation<Frame> problem(testproblem, frame);
  cout << "hier2" << endl;
  CachedQuarkletProblem<SturmEquation<Frame> > cproblem(&problem);
  
  
  cout << "NormA: "<< cproblem.norm_A()<< endl;
  cout << "NormAinv: "<< cproblem.norm_Ainv()<< endl;
  
  
  
  
  cout << "Deltas" << endl;
  cout << frame.DeltaLmin(0) << endl;
  cout << frame.DeltaRmax(3,0) << endl;
  cout << frame.Nablamin(0) << endl;
  cout << frame.Nablamax(3,0) << endl;
  Index eta(1,5,1,9, &frame);
  Index teta(1,5,1,9, &frame);
  cout << eta.number() << endl;
//  cout << "Hau raus: " << cproblem.number(eta,jmax) << endl;
  cout << "Deltasize: " << frame.Deltasize(4,1) << endl;
  cout << "DeltaNablasize: " << frame.DeltaNablasize(4,1) << endl;
  cout << cproblem.a(eta,teta) << endl;
  cout << problem.a(eta,teta) << endl;
  std::list<typename PQFrame<d,dT>::Index> intersecting;
  
  
  typedef std::list<Index> IntersectingList;
  IntersectingList nus;
  intersecting_quarklets(frame,eta,5, 0,nus, 1);
  
  for (typename IntersectingList::const_iterator it2(nus.begin()), itend2(nus.end());
		 it2 != itend2; ++it2) {
                
	      cout << *it2 <<  endl;
  }
  
#endif
#ifdef BASIS
  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index; 
  Basis basis(0,0); // 1st order complementary b.c.'s at x=0 and x=1
  Basis basis10(1,0);
  Basis basis01(0,1);
  Basis basis11(true,true);
  basis11.set_jmax(jmax);
  basis.set_jmax(jmax);
  SturmEquation<Basis> problem(testproblem, basis11);
  CachedProblem<SturmEquation<Basis> > cproblem(&problem);
  Index eta(3,0,0, &basis11);
  Index teta(5,1,0, &basis11);
  cout << "nummer1: " << eta.number() << endl;
  cout << "nummer2: " << teta.number() << endl;
  cout << "Hier: " << problem.a(eta,teta) << endl;
  cout << "Hier: " << cproblem.a(eta,teta) << endl;
#endif
#ifdef basis
 Index phi(6,1,0,&basis); 
 InfiniteVector<double, Index> c,c2;
 basis.reconstruct_1(phi, 7, c);
 basis.decompose_1(phi,4,c2);
  cout << "two-scale: " << endl << c << endl; 
  cout << "decompose" << endl << c2 << endl;
#endif
  
  
  
  
#ifdef FRAME
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
#if 0 //für Randanpassung
  Matrix<double> A;
  Index lambda(0,6,0,-1, &frame);
  system_matrix(frame, A, lambda, 0);
  cout << "LGS linke Seite Ansatz 2" << endl;
  cout << A << endl;
  A.matlab_output("../../Desktop/plots/Ansatz2matrix", "A", 0);
  const int r = 0;
  cout << "rechte Seite: " <<rightside(frame, lambda, r, 0) << endl;
  Vector<double> coeffs2;
  rightsidevector(frame, lambda, coeffs2, 0);
  cout << "RHS: " << coeffs2 << endl;
  coeffs2.matlab_output("../../Desktop/plots/Ansatz2vector", "v");
#endif  
  SparseMatrix<double> mj1;
 frame.assemble_Mj0(3, mj1);
 cout << "primale Zweiskalenrelation generator/generator" << endl;
  cout << mj1 << endl;
  cout << "primale Zweiskalenrelation generator/wavelet" << endl;
  cout << frame.get_Mj1() << endl;
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  //cproblem.norm_Ainv();
  //cproblem.norm_A();
  MonomeFunction newfunction(2);
  Index mu(0,4,0,-1, &frame);
  Index phi(5,3,1,0,&frame);
  InfiniteVector<double, Index> c;
  InfiniteVector<double, int> c2;
  frame.reconstruct_1(phi, 4, c);
  frame.reconstruct_1(phi.p(), phi.j(), phi.e(), phi.k(), 4, c2);
  cout << "two-scale: " << endl << c << endl;
  cout << "alt. two-scale: " << endl << c2[5] << endl;
  double mysumme=0;
  for(int k=5;k<=12;k++){
      mysumme+=pow(k+37,2)*c2[k];
  cout << "Summe: " << mysumme << endl;
  }
  double summe=0;
  //cout << "Size: " << frame.Deltasize(4) << endl;
  cout << "Test auf verschwindende Momente: " << integrate(&newfunction, frame, phi) <<endl;
  
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
 
#if 1
  SparseMatrix<double> A_Lambda;
//  setup_factor_matrix(frame, A_Lambda, 3, 0);
//  A_Lambda.matlab_output("../../Desktop/plots/my2matrix", "A", 0);
//  cout << A_Lambda << endl;
  //cout << A_Lambda.row_dimension() << endl;
  Vector<double> x(A_Lambda.column_dimension()), Ax(A_Lambda.row_dimension());
//  const double wert = 1;
  x = 1;
  A_Lambda.apply(x, Ax);
  //for(int k=0; k<=x.size();k++){
  //x(4) += wert;
  //}
  cout << x << endl;
  cout << Ax << endl;
  Ax.compress();
  //cout << Ax << endl;
  for(int l = -1; l<=6; l++){
//   cout << mu << ", " << integrate(&newfunction, frame, mu) <<endl;
//   summe+=c[mu]*integrate(&newfunction, frame, mu);
   double Faktor=factor(frame, l, 0, 2, 0, 3);
   //cout << Faktor << endl;
   summe+= Faktor;
   

   
   //cout << "Zwischensumme: " << summe << endl;
  }
#endif
  //cout << "Summe: " << summe << endl;
//  cout << "Test auf verschwindende Momente: " << integrate(&newfunction, frame, mu) <<endl;
//  Index lambda(0,3, 0,0, &basis11);
//  cout << "Skalarprodukt: " << cproblem.f(lambda) << endl;
 // Index lambda(0,8,1,255,&frame);
 // Index lambda(8,1,255,&basis11);
 // cout << "Integral: " << integrate(&newfunction, frame, mu) << endl;
  //cout << evaluate(frame, 0, lambda, 0.999908) << endl;
  //cout << evaluate(basis11, 0, lambda, 0.999908) << endl;
//  cout << evaluate(basis11, 0, 9,0,basis11.DeltaLmin()+511, 0.999908) << endl;
  //cout << evaluate(frame, 0, 0, 9,0,frame.DeltaLmin()+511, 0.999908) << endl;
  //cout << "Integral: " << integrate(&newfunction, basis11, lambda) << endl; 
  
  
#endif
  
  InfiniteVector<double,Index> f, v, Av;
  InfiniteVector<double, Index> F1_eta, coeffs, rhs, precrhs;
#ifdef FRAME
  expand(&newfunction, frame, 0, jmax, coeffs, pmax);
  expand(&newfunction, frame, 1, jmax, rhs, pmax);
  //coeffs.scale(&cproblem,1);
  //cproblem.RHS(1e-6, precrhs);
  //cout << frame.last_generator(3) << endl;
  //cout << (frame.last_wavelet(jmax)).number() << endl;
  //cout << first_q_generator(&frame, 3) << endl;
  //Index lambda(2,3,0,1,&frame);
  //cout << lambda.number() + lambda.p() * ((frame.last_wavelet(jmax)).number()+1) << endl;
  //Index lambda(1,4,1,0,&frame);
  //cout << "Test: " << cproblem.f(lambda)/cproblem.D(lambda) << endl;
  //APPLY(cproblem, coeffs, 1e-6, Av, jmax, DKOR, pmax, 2,2);
#endif
#ifdef BASIS
//  expand(&newfunction, basis11, 0, jmax, coeffs);
//  expand(&newfunction, basis11, 1, jmax, rhs);
//  coeffs.scale(&cproblem,1);
//  cproblem.RHS(1e-6, precrhs);
////  APPLY(cproblem, coeffs, 1e-6, Av, jmax, CDD1);
#endif
#if 1  
//  Index lambda(0,3,0,1, &frame);
//  Index nu(0,5,1,31, &frame);
//  cproblem.a(lambda, nu);
//  cproblem.a(*frame.get_wavelet(111), nu);
//  
//  cout << *(frame.get_wavelet(62)) << endl;
  coeffs.scale(&cproblem, 1); 
//  cout << "linke Seite: " << endl;
//  cout << coeffs << endl;
//  
//  cout << "rechte Seite: " << endl;
//  cout << rhs << endl;factor(const int l, const int k, const int r, const int q, const int j)
//  cout << "Applyergebnis: " << endl;
//  cout << Av << endl;
  cproblem.RHS(1e-6, F1_eta);
//  cout << F1_eta;
//  const double nu1 = cproblem.norm_Ainv() * l2_norm(F1_eta);
//  cout << "nu = " << nu1 << endl;
  InfiniteVector<double, Index> solution1_epsilon;
#endif
#ifdef BASIS
  //CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon, jmax);
#endif
#ifdef FRAME  
  //CDD2_SOLVE(cproblem, nu1, epsilon1, solution1_epsilon, jmax, DKOR, pmax, 2, 2);
#endif
//  cout << solution1_epsilon << endl;
//  cout << "exakte Lösung: " << endl;
//  cout << coeffs << endl;
//  
  /* plot solution graph */
//    coeffs.scale(&cproblem, -1); /* scaling because ... */
#if 1
  coeffs.clear();
    Index chi(1,3,0,0,&frame);
    coeffs[chi]=1;
   //coeffs.scale(&cproblem, -1);
    SampledMapping<1> sm3(evaluate(cproblem.frame(), coeffs, true, 8,0));
    //solution1_epsilon.scale(&cproblem, -1); /* scaling because ... */
    //SampledMapping<1> sm3(evaluate(cproblem.basis(), solution1_epsilon, true, 2*jmax)); //" Increase last parameter, if Assertion `resolution >= 0' failed."
    std::ofstream u_stream3("plotthis3.m");
    sm3.matlab_output(u_stream3);
    u_stream3 << "figure;\nplot(x,y);"
              << "title('quarkgraph');" << endl;
    u_stream3.close();

#endif 
//  InfiniteVector<double,int> gcoeffs;
//  frame.reconstruct_1(0,2,1,2, 3, gcoeffs); 
//  cout << gcoeffs;
  
#if 0
//  SampledMapping<1> sm1(evaluate(frame, mu, 1,8));
//  SampledMapping<1> sm2(perframe.evaluate(mu, 8, 0));
//  cout << "Periodisches Integral: " << perframe.integrate(0, mu, mu, 0) << endl;
//  
//  cout << "Primbs-Integral: " << integrate(frame, lambda, lambda) << endl;
//  cout << "Integral aus Sturmklasse: " << problem.a(lambda, lambda) << endl;
  
  std::ofstream u_stream1("../../Desktop/plotthis.m");
  sm1.matlab_output(u_stream1);
  u_stream1 << "figure;\nplot(x,y);"
            << "title('primbsquark')" << endl;
  u_stream1.close();
  
//  std::ofstream u_stream2("../../Desktop/plotthis2.m");
//  sm2.matlab_output(u_stream2);
//  u_stream2 << "figure;\nplot(x,y);"
//            << "title('periodicquark')" << endl;
//  u_stream2.close();
//  //SampledMapping<1> sm1(frame.evaluate;
//  
//  cout << frame.first_generator(frame.j0(), 0) << endl;
#endif
  
  //  Basis basis; // no b.c.'s
  
  
//  cout << basis.DeltaLmin() << endl;
//  cout << basis.DeltaLmax() << endl;
//  cout << basis.DeltaRmin(basis.j0()) << endl;
//  cout << basis.DeltaRmax(basis.j0()) << endl;
  
  //cout << basis.first_generator(3) << endl;
  //   Basis basis(1, 0); // complementary b.c. at x=0
  //   Basis basis(0, 1); // complementary b.c. at x=1
  //Array1D<int>() Ar;
  //evaluate(Basis,2,Index,Ar,Ar);
//  cout << "- d=" << d << ", dT=" << dT << endl;
//  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
//  cout << "- the default wavelet index: " << Index() << endl;
//  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
//  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
//  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
//  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;

//   basis.set_jmax(basis.j0()+4);
//   for (int i = 0; i < basis.Deltasize(basis.j0()+5); i++) {
//     const Index* ind = basis.get_wavelet(i);
//     cout << (*ind) << endl;
//   }


  //  abort();  
  
//  frame.set_jpmax(8,1);
//
//    FixedArray1D<Basis*,4> BasisArray;
//    BasisArray[0] = &basis;
//    BasisArray[1] = &basis10;
//    BasisArray[2] = &basis01;
//    BasisArray[3] = &basis11;
//    
//    unsigned int levelrange(8);
//    for (unsigned int b=0; b < 4; ++b)
//    {
//        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
//    }
        
#if 0
    
    cout << "improving the log2 Function" << endl;
    /*
    for (unsigned int i = 0; i < (-1); i++)
    {
        log2(i);
        //cout << "i = " << i << "; log2(i) = " << log2(i) << "; floor(log2(i)) = " << floor(log2(i)) << endl;
    }
    cout << "check 1 done" << endl;
    */
    
    
    unsigned int v;  // 32-bit value to find the log2 of 
    const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}; //, 0xFFFFFFFF00000000};
    const unsigned int S[] = {1, 2, 4, 8, 16}; //, 32};
    int i;
undefine
    
#if 1
    for (i = 4; i >= 0; i--)
    {
        cout << "b[" << i <<"] = " << b[i] << endl;
    }
        
    for (unsigned int num = 1; num < -1; ++num)
    {
        v = num;
        register unsigned int r = 0; // result of log2(v) will go here
        for (i = 4; i >= 0; i--) // unroll for speed...
        {
          if (v & b[i])
          {
            v >>= S[i];
            r |= S[i];
          } 
        }
        assert (r == log2(num));
        //assert (r == floor(log2((double)num)));
        //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
    }
    
    for (int num = 1; num > 0; ++num)
    {
        v = num;
        register unsigned int r = 0; // result of log2(v) will go here
        for (i = 4; i >= 0; i--) // unroll for speed...
        {
          if (v & b[i])
          {
            v >>= S[i];
            r |= S[i];
          } 
        }
        assert (r == log2((unsigned int)num));
        //assert (r == floor(log2((double)num)));
        //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
    }
    
  cout << "sizeof(int) = " << sizeof(int) << endl;
  cout << "sizeof(unsigned int) = " << sizeof(unsigned int) << endl;
    cout << "done" << endl;
    abort();
#endif
    
    cout << "Compare speed of floor(log2(double)) and log2(int)" << endl;
  
  int repetitions(1);
  clock_t tstart, tend;
  double time1(0), time2(0);
 
  
  tstart = clock();
  for (unsigned int rep = 0; rep < repetitions; ++rep)
  {
      for (unsigned int num = 1; num< -1; ++num)
      {
          
          v = num;
          register unsigned int r = 0; // result of log2(v) will go here
          for (i = 4; i >= 0; i--) // unroll for speed...
          {
              if (v & b[i])
              {
                  v >>= S[i];
                  r |= S[i];
              } 
          }
          //v = log2(num);
            //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
        }
  }
  tend = clock();
  time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  
  tstart = clock();
  for (unsigned int rep = 0; rep < repetitions; ++rep)
  {
      for (unsigned int num = 1; num< -1; ++num)
      {
          //v = floor(log2(num));
          v = log2(num);
      }
  }
  tend = clock();
  time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "Repetitions = " << repetitions << endl
          << "test1 " << time1 << "sec; test2 " << time2 << "sec" << endl;

  
/*
    static const char LogTable256[256] = 
    {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
        -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
        LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
        LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
    };
    unsigned int v; // 32-bit word to find the log of
    unsigned r;     // r will be lg(v)
    register unsigned int t, tt; // temporaries

    for (unsigned int i = 1; i<-1; ++i)
    {
        if (tt = v >> 16)
        {
          r = (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
        }
        else 
        {
          r = (t = v >> 8) ? 8 + LogTable256[t] : LogTable256[v];
        }
        cout << "i = " << i << "; r = " << r << "; log2(i) = " << log2(i) << "; floor(log2(i)) = " << floor(log2(i)) << endl;
        //assert (r == floor(log2(i)));
    }
    
  */  
    
#endif
    
#if 0
    const int levelrange(12);
    
    cout << "improve speed of IntervalIndex(int)" << endl;
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis No " << b << endl;
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
        for (unsigned int num=0; num < BasisArray[b]->degrees_of_freedom(); num++)
        {
            ///*
            int j_, e_, k_;
            int j0 = BasisArray[b]->j0();
            int nabla = BasisArray[b]->Nablasize(j0);
            double tmp;
            int delta = BasisArray[b]->Deltasize(j0);
            // using that delta(j0) + sum_{l=0}^{j} 2^l *nabla = Deltasize(j+j0) 
            if (num < delta) 
            {
              j_ = j0;
              e_ = 0;
            }
            else 
            {
              tmp = num - delta;
              tmp = (double)tmp/nabla +1.;
              j_ = floor(log2(tmp)) + j0;
              e_ = 1;
            }

            if (e_ == 0)
              k_ = BasisArray[b]->DeltaLmin() + num;
            else
              k_ = BasisArray[b]->Nablamin()  + num -delta +nabla -(1<<(j_-j0))*nabla ;
             /*
            cout << "N = " << num << " = (" << j_ << "," << e_ << "," << k_ << ")" 
                    << "; num -delta = " << (num-delta) 
                    << "; (num - delta)/nabla = " << ((num - delta)/nabla) 
                    << "; (double)(num - delta)/nabla = " << ((double)(num - delta)/nabla) 
                    << "; floor(log2(tmp)) = " << floor(log2(((double)(num - delta)/nabla))) << endl;
              * */
            //*/
            // ----------
            int j2_, e2_, k2_;
            int j02 = BasisArray[b]->j0();
            int nabla2 = BasisArray[b]->Nablasize(j02);
            int tmp2;
            int delta2 = BasisArray[b]->Deltasize(j02);
            // using that delta(j0) + sum_{l=0}^{j} 2^l *nabla = Deltasize(j+j0) 
            if (num < delta2) 
            {
              j2_ = j02;
              e2_ = 0;
            }
            else 
            {
              tmp2 = num - delta2;
              tmp2 = tmp2/nabla2 +1;
              j2_ = log2((unsigned int)tmp2) + j02;
              e2_ = 1;
            }

            if (e2_ == 0)
              k2_ = BasisArray[b]->DeltaLmin() + num;
            else
              k2_ = BasisArray[b]->Nablamin()  + num -delta2 +nabla2 -(1<<j2_) ;
            //-----
            assert (((j_ == j2_) && (e_ == e2_) ) && (k_ == k2_) );
        }
    }
    
    cout << "Compare speed of old and new code" << endl;
    int repetitions(30);
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    cout << "test 1 - old code" << endl;
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; ++rep)
    {
        for (unsigned int b=0; b < 4; ++b)
        {
            cout << "Basis No " << b << endl;
            //BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
            for (unsigned int num=0; num < BasisArray[b]->degrees_of_freedom(); num++)
            {
                int j_, e_, k_;
                int j0 = BasisArray[b]->j0();
                int nabla = BasisArray[b]->Nablasize(j0);
                double tmp;
                int delta = BasisArray[b]->Deltasize(j0);
                // using that delta(j0) + sum_{l=0}^{j} 2^l *nabla = Deltasize(j+j0) 
                if (num < delta) 
                {
                  j_ = j0;
                  e_ = 0;
                }
                else 
                {
                  tmp = num - delta;
                  tmp = (double)tmp/nabla +1.;
                  j_ = floor(log2(tmp)) + j0;
                  e_ = 1;
                }

                if (e_ == 0)
                  k_ = BasisArray[b]->DeltaLmin() + num;
                else
                  k_ = BasisArray[b]->Nablamin()  + num -delta +nabla -(1<<(j_-j0))*nabla ;
            }
        }
    }
    tend = clock();
  
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  
  
    tstart = clock();
  
    cout << "test 2 - new code" << endl;
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; ++rep)
    {
        for (unsigned int b=0; b < 4; ++b)
        {
            cout << "Basis No " << b << endl;
            //BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
            for (unsigned int num=0; num < BasisArray[b]->degrees_of_freedom(); num++)
            {
                int j2_, e2_, k2_;
                int j02 = BasisArray[b]->j0();
                int nabla2 = BasisArray[b]->Nablasize(j02);
                int tmp2;
                int delta2 = BasisArray[b]->Deltasize(j02);
                // using that delta(j0) + sum_{l=0}^{j} 2^l *nabla = Deltasize(j+j0) 
                if (num < delta2) 
                {
                    j2_ = j02;
                    e2_ = 0;
                }
                else 
                {
                    tmp2 = num - delta2;
                    tmp2 = tmp2/nabla2 +1;
                    j2_ = log2((unsigned int)tmp2) + j02;
                    e2_ = 1;
                }
                if (e2_ == 0)
                    k2_ = BasisArray[b]->DeltaLmin() + num;
                else
                    k2_ = BasisArray[b]->Nablamin()  + num -delta2 +nabla2 -(1<<j2_) ;
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "Repetitions = " << repetitions << endl
            << "test1 " << time1 << "sec; test2 " << time2 << "sec" << "; time1/time2 = " << (time1/time2) << endl;
#endif
    
#if 0
    cout << "Check implementation of log2(int) in IntervalIndex(num, basis)" << endl;
    int repetitions(40);
    const unsigned int levelrange(12);
    clock_t tstart, tend;
    double time1(0), time2(0);
    cout << "test 1 - old code" << endl;
    tstart = clock();
    for (unsigned int b=0; b < 4; ++b)
    {
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
        cout << "Basis b = " << b << "; dof = " << BasisArray[b]->degrees_of_freedom() << endl;
        for (unsigned int rep = 0; rep < repetitions; ++rep)
        {
            for (unsigned int num=0; num < BasisArray[b]->degrees_of_freedom(); num++)
            {
                Index temp_ind(num, BasisArray[b]);
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "time needed " << time1 << " sec" << endl;
#endif
    
#if 0
  
    //const int levelrange(8);
    //basis.set_jmax(basis.j0()+levelrange);
    cout << "Testing reconstruct_1" << endl;
    cout << "Compare output of reconstruct_1(lambda, InfiniteVector) and reconstruct_1(j,e,k, int)" << endl;
    Index temp_ind;
    InfiniteVector<double, Index> temp_iv1;
    InfiniteVector<double, int> temp_iv2;
    
    int minlambdanum(0), maxlambdanum;
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis b = " << b << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int j=BasisArray[b]->j0(); j < BasisArray[b]->j0()+levelrange; ++j)
        {
            for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
            {
                temp_ind = BasisArray[b]->get_wavelet(i);
                //cout << "j = " << j << "; temp_ind = " << temp_ind << endl;
                //cout << "basis.DeltaLmin() = " << basis.DeltaLmin() << ";";
                //cout << "basis.DeltaSize(" << temp_ind.j() << ") = " << basis.Deltasize(temp_ind.j()) << ";";
                //cout << "basis.Nablamin() = " << basis.Nablamin() << endl;
                temp_iv1.clear();
                BasisArray[b]->reconstruct_1(temp_ind,j, temp_iv1);
                //cout << "temp_iv1 = " << temp_iv1 << endl;
                //cout << "j = " << j << "; temp_ind = " << temp_ind << "; temp_ind.j() = " << temp_ind.j() << "; temp_ind.e() = " << temp_ind.e() << "; temp_ind.k() = " << temp_ind.k() << endl;
                temp_iv2.clear();
                BasisArray[b]->reconstruct_1(temp_ind.j(), temp_ind.e(), temp_ind.k(), j, temp_iv2);
                //cout << "temp_iv2 = " << temp_iv2 << endl;
                //cout.flush();
                assert (temp_iv1.size() == temp_iv2.size());

                InfiniteVector<double, int>::const_iterator it2(temp_iv2.begin());

                for (InfiniteVector<double, Index>::const_iterator it(temp_iv1.begin()), itend(temp_iv1.end()); it != itend; ++it, ++it2 )
                {
                    //assert (it.index() == it2.index());
                    if (!(it.index().number() == it2.index()))
                    {
                        cout << "Problem detected!" << endl;
                        cout << "i = " << i << "; j = " << j << "; temp_ind = " << temp_ind << "; temp_ind.j() = " << temp_ind.j() << "; temp_ind.e() = " << temp_ind.e() << "; temp_ind.k() = " << temp_ind.k() << endl;
                        cout << "temp_iv1 = " << temp_iv1 << endl;
                        cout << "temp_iv2 = " << temp_iv2 << endl;
                        cout << it.index() << "; .number() = " << it.index().number() << "; *it = " << *it << "; ";
                        cout << it2.index() << "; *it2 = " << *it << endl;
                    }
                    assert (it.index().number() == it2.index());
                    assert (*it == *it2);
                    //cout << it.index() << "; .number() = " << it.index().number() << "; *it = " << *it << "; ";
                    //cout << it2.index() << "; *it2 = " << *it << endl;
                    //cout << "it = (" << it << ", " << *it << "); it2 = (" << it2 << ", " << *it2 << "); .number() = " << endl;
                }
            }
        }
    }
    cout << "Compare speed of reconstruct_1(lambda,Index) and reconstruct_1(j,e,k,int)" << endl;
  
    int repetitions(10);
    Index lambda_run;
    clock_t tstart, tend;
    double time1(0), time2(0);
 
    tstart = clock();
    
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis b = " << b << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int rep = 0; rep < repetitions; ++rep)
        {
            for (unsigned int j=BasisArray[b]->j0(); j < BasisArray[b]->j0()+levelrange; ++j)
            {
                for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
                {
                    temp_ind = BasisArray[b]->get_wavelet(i);
                    temp_iv1.clear();
                    BasisArray[b]->reconstruct_1(temp_ind,j, temp_iv1);
                }
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  
    tstart = clock();
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis b = " << b << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int rep = 0; rep < repetitions; ++rep)
        {
            for (unsigned int j=BasisArray[b]->j0(); j < BasisArray[b]->j0()+levelrange; ++j)
            {
                for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
                {
                    temp_ind = BasisArray[b]->get_wavelet(i);
                    temp_iv2.clear();
                    BasisArray[b]->reconstruct_1(temp_ind.j(),temp_ind.e(),temp_ind.k(),j, temp_iv2);
                }
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "Repetitions = " << repetitions << "; levelrange = " << levelrange << "; minlambdanum = " << minlambdanum << "; maxlambdanum = " << maxlambdanum << endl
            << "reconstruct_1(lambda,Index) " << time1 << "sec; reconstruct_1(j,e,k,int) " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
#endif
  
  
#if 0
  //const int levelrange(8);
  //basis.set_jmax(basis.j0()+levelrange);
  cout << "Testing evaluate" << endl;
  cout << "Compare output of evaluate(lambda, x) and evaluate(j,e,k, x)" << endl;

  // check all values x\in [a, a+h,..., a+steps*h]
  double a(-1), b(2); 
  int steps (100);
  double h = (b-a) / steps;
  
  int minlambdanum(0), maxlambdanum;
  double x1,x2,x3,x4;
  
  Index temp_ind;
  InfiniteVector<double, Index> temp_iv1, temp_iv2;
  
  
  for (unsigned int b=0; b < 4; ++b)
  {
      cout << "Basis b = " << b << endl;
      maxlambdanum = BasisArray[b]->degrees_of_freedom();
      for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
      {
          //cout << "i = " << i << endl;
          for (unsigned int j=0; j<steps; ++j)
          {
              //cout << "j = " << j << endl;
              temp_ind = BasisArray[b]->get_wavelet(i);
              //cout << "j = " << j << "; temp_ind = " << temp_ind << endl;
              x1 = BasisArray[b]->evaluate(0,temp_ind,a+h*j);
              x2 = BasisArray[b]->evaluate(1,temp_ind,a+h*j);
              x3 = BasisArray[b]->evaluate(0,temp_ind.j(), temp_ind.e(), temp_ind.k(),a+h*j);
              x4 = BasisArray[b]->evaluate(1,temp_ind.j(), temp_ind.e(), temp_ind.k(),a+h*j);
              if (x1 != x3)
              {
                  cout << "x1 = " << x1 << "; x3 = " << x3 << "; diff = " << (x1-x3) << endl;
              }
              if (x2 != x4)
              {
                  cout << "x2 = " << x2 << "; x4 = " << x4 << "; diff = " << (x2-x4) << endl;
              }
              //assert (x1 == x3);
              //assert (x2 == x4);
          }
      }
  }
  
  
  cout << "Compare speed of evaluate(lambda) and evaluate(j,e,k)" << endl;
  int repetitions(2);
  clock_t tstart, tend;
  double time1(0), time2(0);
 
  tstart = clock();
  for (unsigned int b=0; b < 4; ++b)
  {
      cout << "Basis b = " << b << endl;
      maxlambdanum = BasisArray[b]->degrees_of_freedom();
      for (unsigned int rep = 0; rep < repetitions; ++rep)
      {
          for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
          {
              //cout << "i = " << i << endl;
              for (unsigned int j=0; j<steps; ++j)
              {
                  temp_ind = BasisArray[b]->get_wavelet(i);
                  //cout << "j = " << j << "; temp_ind = " << temp_ind << endl;
                  x1 = BasisArray[b]->evaluate(0,temp_ind,a+h*j);
                  x2 = BasisArray[b]->evaluate(1,temp_ind,a+h*j);
              }
          }
      }
  }
  tend = clock();
  time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  
  tstart = clock();
  for (unsigned int b=0; b < 4; ++b)
  {
      cout << "Basis b = " << b << endl;
      maxlambdanum = BasisArray[b]->degrees_of_freedom();
      for (unsigned int rep = 0; rep < repetitions; ++rep)
      {
          for (unsigned int i=minlambdanum; i< maxlambdanum; ++i)
          {
              //cout << "i = " << i << endl;
              for (unsigned int j=0; j<steps; ++j)
              {
                  temp_ind = basis.get_wavelet(i);
                  x3 = basis.evaluate(0,temp_ind.j(), temp_ind.e(), temp_ind.k(),a+h*j);
                  x4 = basis.evaluate(1,temp_ind.j(), temp_ind.e(), temp_ind.k(),a+h*j);
              }
          }
      }
  }
  tend = clock();
  time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "Repetitions = " << repetitions << "; levelrange = " << levelrange << "; minlambdanum = " << minlambdanum << "; maxlambdanum = " << maxlambdanum << endl
          << "evaluate(lambda) " << time1 << "sec; evaluate(j,e,k) " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
#endif  
  
#if 0
    cout << "Check code modification in p_support::evaluate(basis, derivative, Index, points, values)" << endl;
    //const unsigned int levelrange(12);
    
    double aa(-0.1), bb(1.1); 
    int steps (1000);
    double h = (bb-aa) / steps;
    Array1D<double> points, values, values2;
    points.resize(steps);
    for (unsigned int i=0; i<steps; ++i)
        points[i] = aa+i*h;
    Index temp_ind;
    
#if 0 //test equality of evaluate(Index,points) and evaluate(j,e,k,points)
     // test a suspicious wavelet
    {
        BasisArray[0]->set_jmax(BasisArray[0]->j0()+levelrange);
        temp_ind = BasisArray[0]->get_wavelet(514);
        evaluate(*BasisArray[0], 0, temp_ind, points, values);
        evaluate(*BasisArray[0], 1, temp_ind, points, values2);
        for (unsigned int i = 0; i< steps; ++i)
        {
            double temp_d = evaluate(*BasisArray[0],0,temp_ind,points[i]);
            double temp_d2 = evaluate(*BasisArray[0],1,temp_ind,points[i]);
            if (((abs(temp_d) != 0) || (abs(temp_d2) != 0) )
                || ((abs(values[i]) != 0) || (abs(values2[i]) != 0) ))
                cout << "i=" << i << "; points["<<i<<"] = " << points[i] <<"; evaluate(0,x) = " << temp_d << "; evaluate(0,points)["<<i<<"] = " << values[i] <<"; evaluate(1,x) = " << temp_d2 << "; evaluate(1,points)["<<i<<"] = " << values2[i] << endl;
        }
    }
    cout.flush();
    //return 0;
    
    for (unsigned int b=0; b < 4; ++b)
    {
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
        cout << "Basis b = " << b << "; dof = " << BasisArray[b]->degrees_of_freedom() << endl;
        cout << "get_evaluate_with_pre_computation() = " << BasisArray[b]->get_evaluate_with_pre_computation() << endl;
        for (unsigned int num=0; num < BasisArray[b]->degrees_of_freedom(); num++)
        {
            cout << "num = " << num << endl;
            temp_ind = BasisArray[b]->get_wavelet(num);
            evaluate(*BasisArray[b], 0, temp_ind, points, values);
            evaluate(*BasisArray[b], 0, temp_ind.j(), temp_ind.e(), temp_ind.k(), points, values2);
            assert (values.size() == steps);
            assert (values2.size() == steps);
            
            
            double temp_d(0), temp_d2(0);
            for (unsigned int i=0; i<steps; ++i)
            {
                temp_d += abs(values[i]);
                temp_d2 += abs(values2[i]);
                if (num == 514) cout << "i = " << i << "; v[" << i << "] = " << values[i] << "; v2[" << i << "] = " << values2[i] << endl;
                assert (values[i] == values2[i]);
            }
            if (num == 514)
            {
                cout << "val[" << 675 << "] = " << values[675] << "; val2[" << 675 << "] = " << values2[675] << "; temp_d = " << temp_d << "; temp_d2 = " << temp_d2 << endl;
                //return 0;
            }
            if ((temp_d == 0) || (temp_d2 == 0))
            {
                cout << "Problem: evaluation of temp_ind = " << temp_ind << "results in 0 on the whole interval!" << endl;
                cout << "basis b = " << b << "; num = " << num << "; temp_ind.number() = " << temp_ind.number() << endl;
                steps *= 3;
                cout << "resolution too coarse?! Increase steps and retry!. New stepsize = " << steps << endl;
                cout << "values = " << endl << values << endl << endl << "values2 = " << endl << values2 << endl;
                h = (bb-aa) / steps;
                points.resize(steps);
                for (unsigned int i=0; i<steps; ++i)
                    points[i] = aa+i*h;
                num = num - 1;
            }
            
            evaluate(*BasisArray[b], 1, temp_ind, points, values);
            evaluate(*BasisArray[b], 1, temp_ind.j(), temp_ind.e(), temp_ind.k(), points, values2);
            assert (values.size() == steps);
            assert (values2.size() == steps);
            temp_d=0;
            temp_d2=0;
            for (unsigned int i=0; i<steps; ++i)
            {
                temp_d += abs(values[i]);
                temp_d2 += abs(values2[i]);
                //cout << "i = " << i << "; v[" << i << "] = " << values[i] << "; v2[" << i << "] = " << values2[i] << endl;
                assert (values[i] == values2[i]);
            }
            if (num == 514)
            {
                cout << "val[" << 675 << "] = " << values[675] << "; val2[" << 675 << "] = " << values2[675] << "; temp_d = " << temp_d << "; temp_d2 = " << temp_d2 << endl;
                return 0;
            }
            if ((temp_d == 0) || (temp_d2 == 0))
            {
                cout << "Problem: evaluation of temp_ind = " << temp_ind << "results in 0 on the whole interval!" << endl;
                cout << "basis b = " << b << "; num = " << num << "; temp_ind.number() = " << temp_ind.number() << endl;
                steps *= 3;
                cout << "resolution too coarse?! Increase steps and retry!. New stepsize = " << steps << endl;
                cout << "values = " << endl << values << endl << endl << "values2 = " << endl << values2 << endl;
                h = (bb-aa) / steps;
                points.resize(steps);
                for (unsigned int i=0; i<steps; ++i)
                    points[i] = aa+i*h;
                num = num - 1;
            }
        }
    }
    cout << "done" << endl;
    return 0;
#endif
#if 1
    cout << "Compare speed of evaluate(lambda, points) and evaluate(j,e,k, points)" << endl;
    int repetitions(2);
    unsigned int minlambdanum(0), maxlambdanum;
    clock_t tstart, tend;
    double time1(0), time2(0);

    double temp_d(0), temp_d2(0);

    steps = 900;
    h = (bb-aa) / steps;
    points.resize(steps);
    for (unsigned int i=0; i<steps; ++i)
        points[i] = aa+i*h;
            
    for (unsigned int b=0; b < 4; ++b)
    {
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
    }
    
    tstart = clock();
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis b = " << b << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int rep = 0; rep < repetitions; ++rep)
        {
            /*
            steps = 100;
            h = (bb-aa) / steps;
            points.resize(steps);
            for (unsigned int i=0; i<steps; ++i)
                points[i] = aa+i*h;
             * */
            for (unsigned int num=minlambdanum; num< maxlambdanum; ++num)
            {
                //cout << "num = " << num << endl;
                temp_ind = BasisArray[b]->get_wavelet(num);
                evaluate(*BasisArray[b], 0, temp_ind, points, values);
                evaluate(*BasisArray[b], 1, temp_ind, points, values2);

                /*
                temp_d = 0;
                temp_d2 = 0;
                for (unsigned int i=0; i<steps; ++i)
                {
                    temp_d += abs(values[i]);
                    temp_d2 += abs(values2[i]);

                }
                if ((temp_d == 0) || (temp_d2 == 0))
                {
                    //cout << "Problem: evaluation of temp_ind = " << temp_ind << "results in 0 on the whole interval!" << endl;
                    //cout << "basis b = " << b << "; num = " << num << "; temp_ind.number() = " << temp_ind.number() << endl;
                    steps *= 3;
                    //cout << "resolution too coarse?! Increase steps and retry!. New stepsize = " << steps << endl;
                    //cout << "values = " << endl << values << endl << endl << "values2 = " << endl << values2 << endl;
                    h = (bb-aa) / steps;
                    points.resize(steps);
                    for (unsigned int i=0; i<steps; ++i)
                        points[i] = aa+i*h;
                    num = num - 1;
                }
                */
                
            }
        }
        //cout << "maximal stepsize for basis b = " << b << " was " << steps << endl;
    }
    tend = clock();
    
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    for (unsigned int b=0; b < 4; ++b)
    {
        cout << "Basis b = " << b << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int rep = 0; rep < repetitions; ++rep)
        {
            /*
            steps = 100;
            h = (bb-aa) / steps;
            points.resize(steps);
            for (unsigned int i=0; i<steps; ++i)
                points[i] = aa+i*h;
             */
            for (unsigned int num=minlambdanum; num< maxlambdanum; ++num)
            {
               
                temp_ind = BasisArray[b]->get_wavelet(num);
                evaluate(*BasisArray[b], 0, temp_ind.j(), temp_ind.e(), temp_ind.k(), points, values);
                evaluate(*BasisArray[b], 1, temp_ind.j(), temp_ind.e(), temp_ind.k(), points, values2);

                /*
                temp_d = 0;
                temp_d2 = 0;
                for (unsigned int i=0; i<steps; ++i)
                {
                    temp_d += abs(values[i]);
                    temp_d2 += abs(values2[i]);

                }
                if ((temp_d == 0) || (temp_d2 == 0))
                {
                    //cout << "Problem: evaluation of temp_ind = " << temp_ind << "results in 0 on the whole interval!" << endl;
                    //cout << "basis b = " << b << "; num = " << num << "; temp_ind.number() = " << temp_ind.number() << endl;
                    steps *= 3;
                    //cout << "resolution too coarse?! Increase steps and retry!. New stepsize = " << steps << endl;
                    //cout << "values = " << endl << values << endl << endl << "values2 = " << endl << values2 << endl;
                    h = (bb-aa) / steps;
                    points.resize(steps);
                    for (unsigned int i=0; i<steps; ++i)
                        points[i] = aa+i*h;
                    num = num - 1;
                }
                 */
            }
        }
        //cout << "maximal stepsize for basis b = " << b << " was " << steps << endl;
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; levelrange = " << levelrange << "; minlambdanum = " << minlambdanum << "; maxlambdanum = " << maxlambdanum << endl
            << "evaluate(lambda,points) " << time1 << "sec; evaluate(j,e,k,points) " << time2 << "sec; time1/time2 = " << (time1/time2) << endl;
#endif 
#endif //end of evaluate(points)
  
#if 0
    cout << "Check code modification in p_support::evaluate(basis, Index, points, fktvalues, dervalues)" << endl;
    //const unsigned int levelrange(12);
    
    double aa(-0.1), bb(1.1); 
    int steps (900);
    double h = (bb-aa) / steps;
    Array1D<double> points, funcvalues1, funcvalues2, dervalues1, dervalues2;;
    points.resize(steps);
    for (unsigned int i=0; i<steps; ++i)
        points[i] = aa+i*h;
    Index temp_ind;
    

    //test equality of evaluate(Index,points) and evaluate(j,e,k,points)
    // no speed check. should be the same speedup as above
    /*
     // test a suspicious wavelet
    {
        BasisArray[0]->set_jmax(BasisArray[0]->j0()+levelrange);
        temp_ind = BasisArray[0]->get_wavelet(10);
        evaluate(*BasisArray[0], 0, temp_ind, points, values);
        evaluate(*BasisArray[0], 1, temp_ind, points, values2);
        for (unsigned int i = 0; i< steps; ++i)
        {
            double temp_d = evaluate(*BasisArray[0],0,temp_ind,points[i]);
            double temp_d2 = evaluate(*BasisArray[0],1,temp_ind,points[i]);
            if (((abs(temp_d) != 0) || (abs(temp_d2) != 0) )
                || ((abs(values[i]) != 0) || (abs(values2[i]) != 0) ))
                cout << "i=" << i << "; points["<<i<<"] = " << points[i] <<"; evaluate(0,x) = " << temp_d << "; evaluate(0,points)["<<i<<"] = " << values[i] <<"; evaluate(1,x) = " << temp_d2 << "; evaluate(1,points)["<<i<<"] = " << values2[i] << endl;
        }
    }
    cout.flush();
    //return 0;
    */
    
    unsigned int minlambdanum(0), maxlambdanum;
    
    for (unsigned int b=0; b < 4; ++b)
    {
        BasisArray[b]->set_jmax(BasisArray[b]->j0()+levelrange);
        cout << "Basis b = " << b << "; dof = " << BasisArray[b]->degrees_of_freedom() << endl;
        cout << "get_evaluate_with_pre_computation() = " << BasisArray[b]->get_evaluate_with_pre_computation() << endl;
        maxlambdanum = BasisArray[b]->degrees_of_freedom();
        for (unsigned int num=minlambdanum; num< maxlambdanum; ++num)
        {
             /*
            if ( (num == 11))
            {
                return 0;
            }
              * **/
            cout << "num = " << num << endl;
            temp_ind = BasisArray[b]->get_wavelet(num);
            evaluate(*BasisArray[b], temp_ind, points, funcvalues1, dervalues1);
            //evaluate2(*BasisArray[b], temp_ind, points, funcvalues2, dervalues2);
            evaluate(*BasisArray[b], temp_ind.j(), temp_ind.e(), temp_ind.k(), points, funcvalues2, dervalues2);
            
            //cout << "funcvalues1 = " << endl << funcvalues1 << endl << endl << "funcvalues2 = " << endl << funcvalues2 << endl;
            //cout << "dervalues1 = " << endl << dervalues1 << endl << endl << "dervalues2 = " << endl << dervalues2 << endl;
                
            assert (funcvalues1.size() == steps);
            assert (funcvalues2.size() == steps);
            assert (dervalues1.size() == steps);
            assert (dervalues2.size() == steps);
            
            
            double temp_d1(0), temp_d2(0), temp_d3(0), temp_d4(0);
            for (unsigned int i=0; i<steps; ++i)
            {
                temp_d1 += abs(funcvalues1[i]);
                temp_d2 += abs(funcvalues2[i]);
                temp_d3 += abs(dervalues1[i]);
                temp_d4 += abs(dervalues2[i]);
                //if (num == 514) cout << "i = " << i << "; v[" << i << "] = " << values[i] << "; v2[" << i << "] = " << values2[i] << endl;
                if (funcvalues1[i] != funcvalues2[i])
                {
                    cout << "i = " << i << "; funcvalues1[" << i << "] = " << funcvalues1[i] << "; funcvalues2[" << i << "] = " << funcvalues2[i] << endl;
                }
                //assert (funcvalues1[i] == funcvalues2[i]);
                //assert (dervalues1[i] == dervalues2[i]);
            }
            /*
            if (num == 514)
            {
                cout << "val[" << 675 << "] = " << values[675] << "; val2[" << 675 << "] = " << values2[675] << "; temp_d = " << temp_d << "; temp_d2 = " << temp_d2 << endl;
                //return 0;
            }
             * */
            if ( ((temp_d1 == 0) || (temp_d2 == 0)) || ((temp_d3 == 0) || (temp_d4 == 0)) )
            {
                cout << "Problem: evaluation of temp_ind = " << temp_ind << "results in 0 on the whole interval!" << endl;
                cout << "basis b = " << b << "; num = " << num << "; temp_ind.number() = " << temp_ind.number() << endl;
                steps *= 3;
                cout << "resolution too coarse?! Increase steps and retry!. New stepsize = " << steps << endl;
                cout << "funcvalues1 = " << endl << funcvalues1 << endl << endl << "funcvalues2 = " << endl << funcvalues2 << endl;
                cout << "dervalues1 = " << endl << dervalues1 << endl << endl << "dervalues2 = " << endl << dervalues2 << endl;
                h = (bb-aa) / steps;
                points.resize(steps);
                for (unsigned int i=0; i<steps; ++i)
                    points[i] = aa+i*h;
                num = num - 1;
            }
            /*
            if ( (num == 11))
            {
                //cout << "val[" << 675 << "] = " << values[675] << "; val2[" << 675 << "] = " << values2[675] << "; temp_d = " << temp_d << "; temp_d2 = " << temp_d2 << endl;
                return 0;
            }
             */
            
        }
    }
    cout << "done" << endl;
    return 0;


#endif //end of evaluate(points, fktval, derval)
    
    
  
#if 0
  {
    SparseMatrix<double> M;
    const int j = basis.j0();

    basis.assemble_Mj0(j, M);
    M.scale(M_SQRT2);
    cout << "* Mj0 (without factor 1/sqrt(2)) for j=j0=" << j << ":" << endl;
    cout << M;
    std::ofstream result("Mj0.m");
    result << "M=";
    print_matrix(M, result);
    result << endl;
    result.close();

    basis.assemble_Mj1(j, M);
    M.scale(M_SQRT2);
    cout << "* Mj1 (without factor 1/sqrt(2)) for j=j0=" << j << ":" << endl;
    cout << M;
    result.open("Mj1.m");
    result << "M=";
    print_matrix(M, result);
    result << endl;
    result.close();
    
    basis.assemble_Mj0T(j, M);
    M.scale(M_SQRT2);
    cout << "* Mj0T (without factor 1/sqrt(2)) for j=j0=" << j << ":" << endl;
    cout << M;
    result.open("Mj0T.m");
    result << "M=";
    print_matrix(M, result);
    result << endl;
    result.close();

    basis.assemble_Mj1T(j, M);
    M.scale(M_SQRT2);
    cout << "* Mj1T (without factor 1/sqrt(2)) for j=j0=" << j << ":" << endl;
    cout << M;
    result.open("Mj1T.m");
    result << "M=";
    print_matrix(M, result);
    result << endl;
    result.close();
  }
#endif


#if 0
  
  typedef CubeBasis<Basis,1> CBasis;
  typedef CBasis::Index CIndex;

  //CBasis.set_jmax()

  FixedArray1D<int,2> bc;
  bc[0] = 0;
  bc[1] = 0;
  
  CBasis cbasis(bc);
  cbasis.set_jmax(cbasis.j0());

  set<CIndex> Lambda;
  for (CIndex lambda(first_generator<Basis,1,CBasis>(&cbasis, cbasis.j0()));; ++lambda) {
    cout << lambda << endl;
    Lambda.insert(lambda);
    if (lambda == last_generator<Basis,1,CBasis>(&cbasis, cbasis.j0())) break;
  }


   Vector<double> value(1);
   value[0] = 1;
   ConstantFunction<1> const_fun(value);
   IdentityBVP<1> trivial_bvp(&const_fun);
   trivial_bvp.set_f(&const_fun);

   CubeEquation<Basis,1,CBasis> eq(&trivial_bvp, bc);
   

   cout << "setting up full right hand side..." << endl;
   Vector<double> rh;
   WaveletTL::setup_righthand_side(eq, Lambda, rh);
   
   cout << rh << endl;
   
   cout << "...done setting up full right hand side" << endl;

   InfiniteVector<double, Index> coeff;
   Index index(first_generator(&basis, basis.j0()));
   for (int i = 0;; ++index, i++) {
     cout << index << endl;
     coeff.set_coefficient(index, rh[i]);
     if (index == last_generator(&basis, basis.j0())) break;
   }
   
   //coeff.scale(&eq,-1);

   //cout << index << endl;
//    ++index;
//    ++index;
//    ++index;
//    ++index;
//    ++index;
//    coeff.set_coefficient(index, 1.);



   SampledMapping<1> res = evaluate(basis, coeff, false, 6);

   cout << "...done evaluating expansion in dual basis" << endl;


   std::ofstream ofs5("reproduced_function.m");
   res.matlab_output(ofs5);
   ofs5.close();

      abort();

#endif


#if 0
  cout << "- checking biorthogonality of Mj0, Mj0T for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++)
    {
      SparseMatrix<double> mj0_t, mj0T;
      basis.assemble_Mj0_t(level, mj0_t);
      basis.assemble_Mj0T(level, mj0T);

      SparseMatrix<double> T = mj0_t * mj0T;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0^T*Mj0T-I||_infty: " << row_sum_norm(T) << endl;

      SparseMatrix<double> mj0, mj0T_t;
      basis.assemble_Mj0(level, mj0);
      basis.assemble_Mj0T_t(level, mj0T_t);

      T = mj0T_t * mj0;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0T^T*Mj0-I||_infty: " << row_sum_norm(T) << endl;
    }
#endif

#if 0
  cout << "- checking biorthogonality of Mj<->Gj and MjT<->GjT for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+1; level++)
    {
      SparseMatrix<double> mj0, mj1;
      basis.assemble_Mj0(level, mj0);
      basis.assemble_Mj1(level, mj1);
      SparseMatrix<double> mj(mj0.row_dimension(), mj0.row_dimension());
      mj.set_block(0, 0, mj0);
      mj.set_block(0, mj0.column_dimension(), mj1);

      SparseMatrix<double> mj0T_t, mj1T_t;
      basis.assemble_Mj0T_t(level, mj0T_t);
      basis.assemble_Mj1T_t(level, mj1T_t);
      SparseMatrix<double> gj(mj0T_t.column_dimension(), mj0T_t.column_dimension());
      gj.set_block(0, 0, mj0T_t);
      gj.set_block(mj0T_t.row_dimension(), 0, mj1T_t);

      SparseMatrix<double> T = mj * gj;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj*Gj-I||_infty: " << row_sum_norm(T) << endl;

      T = gj * mj;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Gj*Mj-I||_infty: " << row_sum_norm(T) << endl;
      
      SparseMatrix<double> mj0T, mj1T;
      basis.assemble_Mj0T(level, mj0T);
      basis.assemble_Mj1T(level, mj1T);
      SparseMatrix<double> mjt(mj.row_dimension(), mj.row_dimension());
      mjt.set_block(0, 0, mj0T);
      mjt.set_block(0, mj0T.column_dimension(), mj1T);

      SparseMatrix<double> mj0_t, mj1_t;
      basis.assemble_Mj0_t(level, mj0_t);
      basis.assemble_Mj1_t(level, mj1_t);
      SparseMatrix<double> gjt(mj.row_dimension(), mj.row_dimension());
      gjt.set_block(0, 0, mj0_t);
      gjt.set_block(mj0_t.row_dimension(), 0, mj1_t);

      T = mjt * gjt;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||MjT*GjT-I||_infty: " << row_sum_norm(T) << endl;

      T = gjt * mjt;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||GjT*MjT-I||_infty: " << row_sum_norm(T) << endl;
    }
#endif

#if 0
  cout << "- checking access to single rows of the M_{j,i} matrices:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++)
    {
      InfiniteVector<double, Vector<double>::size_type> v, w;
      double maxerr = 0.0;
      SparseMatrix<double> mj0, mj0_t, mj1, mj1_t, mj0T, mj0T_t, mj1T, mj1T_t;
      basis.assemble_Mj0(level, mj0); mj0_t = transpose(mj0);
      basis.assemble_Mj1(level, mj1); mj1_t = transpose(mj1);

      //cout << "#########################" << endl;
      //cout << mj1 << endl;

      basis.assemble_Mj0T(level, mj0T); mj0T_t = transpose(mj0T);
      basis.assemble_Mj1T(level, mj1T); mj1T_t = transpose(mj1T);
      for (size_t row = 0; row < mj0.row_dimension(); row++)
	{
	  mj0.get_row(row, v);
	  basis.Mj0_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0_t.row_dimension(); row++)
	{
	  mj0_t.get_row(row, v);
	  basis.Mj0_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1.row_dimension(); row++)
	{
	  mj1.get_row(row, v);
	  
	  basis.Mj1_get_row(level, row, w);
	  
	  maxerr = max(maxerr, linfty_norm(v-w));
	  if ( linfty_norm(v-w)> 1.0e-10) {
	    cout << "########" << endl;
	    cout << "row = " << row << endl;
	    cout << "maxerr= " << maxerr << endl;
	    cout << "v= " << v << endl;
	    cout << "w= " << w << endl;
	  }
	}
      cout << "* j=" << level << ", max. error in Mj1: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1_t.row_dimension(); row++)
	{
	  mj1_t.get_row(row, v);
	  basis.Mj1_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0T.row_dimension(); row++)
	{
	  mj0T.get_row(row, v);
	  basis.Mj0T_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0T: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0T_t.row_dimension(); row++)
	{
	  mj0T_t.get_row(row, v);
	  basis.Mj0T_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0T_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1.row_dimension(); row++)
	{
	  mj1T.get_row(row, v);
	  basis.Mj1T_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1T: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1T_t.row_dimension(); row++)
	{
	  mj1T_t.get_row(row, v);
	  basis.Mj1T_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1T_t: " << maxerr << endl;
    }
#endif

#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(basis.first_generator(level));
      for (;; ++index)
	{
	  InfiniteVector<double, Index> origcoeff;
	  origcoeff[index] = 1.0;
	  
  	  cout << "* original coeffs:" << endl << origcoeff;

	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose(origcoeff, basis.j0(), wcoeff);

  	  cout << "* after decompose():" << endl << wcoeff;
	  
	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct(wcoeff, level, transformcoeff);

	  cout << "* after reconstruct():" << endl << transformcoeff;
	  
	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == basis.last_generator(level)) break;
	}
    }
#endif

#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose_t() and reconstruct_t() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(first_generator(&basis, level));
      for (;; ++index)
	{
	  InfiniteVector<double, Index> origcoeff;
	  origcoeff[index] = 1.0;
	  
	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose_t(origcoeff, basis.j0(), wcoeff);
	  
	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct_t(wcoeff, level, transformcoeff);
	  
	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis, level)) break;
	}
    }
#endif

#if 0
  cout << "- evaluating some primal generators:" << endl;
  for (Index lambda(first_generator(&basis, basis.j0()));; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, true, 7).matlab_output(cout);
    if (lambda == last_generator(&basis, basis.j0())) break;
  }
  
  cout << "- evaluating some primal wavelets:" << endl;
  for (Index lambda = first_wavelet(&basis, basis.j0());; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, true, 6).matlab_output(cout);
    if (lambda == last_wavelet(&basis, basis.j0())) break;
  }
#endif

#if 0
  cout << "- evaluating some dual generators:" << endl;
  for (Index lambda(first_generator(&basis, basis.j0()));; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, false, 7).matlab_output(cout);
    if (lambda == last_generator(&basis, basis.j0())) break;
  }
#endif

#if 0
  InfiniteVector<double, Index> coeff;
  Index index(last_wavelet(&basis, basis.j0()));
  coeff.set_coefficient(index,1.0);
 //  for (int i = 0;; ++index, i++) {
//     //cout << index << endl;
//     coeff.set_coefficient(index, rh[i]);
//     if (index == last_generator(&basis, basis.j0())) break;
  //}
  int dil=8;
  Array1D<double> grid;
  Array1D<double> funval;
  grid.resize((1<<dil)+1);
  for(unsigned int k=0;k<grid.size();k++)
  	grid[k]=k*(1.0/(1<<dil));
  Grid<1>gitter(grid);
  evaluate(basis,0,index,grid,funval);  
  SampledMapping<1> res (gitter, funval);
  std::ofstream ofs5("reproduced_function.m");
  res.matlab_output(ofs5);
  ofs5.close();
   

#endif

  return 0;
}

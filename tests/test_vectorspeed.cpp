#include <iostream>
#include <cmath>
#include <time.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing speed of vector operations..." << endl;

  clock_t tstart, tend;
  double time;

  const int d = 10;
  const int N = 10000000;

  Vector<double> v2(d), w2(v2);
  v2(0) = v2(4) = 1e-6;
  cout << "- adding a sparsely populated " << d << "-d RawVector (dense representation) " << N << " times ..." << endl;
  tstart = clock();
  for (int n(0); n < N; n++)
    w2 += v2;
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time: " << time << "s" << endl;

//   RawVector<double,SparseArray1D<double> > v3(d), w3(v3);
//   v3(0) = v3(4) = 1e-6;
//   cout << "- adding a sparsely populated " << d << "-d RawVector (sparse representation) " << N << " times ..." << endl;
//   tstart = clock();
//   for (int n(0); n < N; n++)
//     w3 -= v3;
//   tend = clock();
//   time = (double)(tend-tstart)/CLOCKS_PER_SEC;
//   cout << "  ... done, time: " << time << "s" << endl;

  return 0;
}

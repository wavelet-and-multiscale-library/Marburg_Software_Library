#include <iostream>
#include <fstream>
#include <cassert>

#include <utils/random.h>
#include <ctime>

#include <interval/s_basis.h>
#include <interval/s_support.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

typedef std::list<SBasis::Index> SupportList;
typedef std::list<std::pair<SBasis::Index, SBasis::Support> > SupportPairList;

// prototype declaration of helper functions
void print_list(SupportList nu);
void print_pair_list(SupportPairList nu);


int main()
{
  cout << "Testing the S support routines ..." << endl << endl;

  SBasis basis;
  SBasis::Index lambda;
  int j, e;
  int k1, k2;

  // initialize random generator
  #ifdef HAVE_CLOCK_GETTIME
  struct timespec d;
  clock_gettime(CLOCK_REALTIME, &d);
  srand(d.tv_nsec);
  #else // !HAVE_CLOCK_GETTIME
  time_t rawtime;
  tm* gt;
  rawtime = time(NULL);
  gt = gmtime(&rawtime);
  srand(gt->tm_sec);
  #endif // !HAVE_CLOCK_GETTIME

  j = basis.j0();// + random_integer(0, 5);
  e = random_integer(0, 1);
  if (e == E_GENERATOR) {
    k1 = basis.DeltaLmin();
    k2 = basis.DeltaRmax(j);
  }
  else { // e == E_WAVELET
    k1 = basis.Nablamin();
    k2 = basis.Nablamax(j);
  }
  lambda = SBasis::Index(j, e, random_integer(k1,k2), random_integer(0, basis.number_of_components-1), &basis);
  assert(lambda.is_valid());
  j = basis.j0() + random_integer(0, 5);
  SupportList intersecting;
  cout << "Intersecting level-"<<j<<"-generators of \\Psi_" << lambda << ":" << endl;
  intersecting_wavelets(basis, lambda, j, true, intersecting);
  print_list(intersecting);
  cout << "Intersecting level-"<<j<<"-wavelets of \\Psi_" << lambda << ":" << endl;
  intersecting_wavelets(basis, lambda, j, false, intersecting);
  print_list(intersecting);

  cout << "Also give the intersections:" << endl;
  cout << "lambda = " << lambda << endl;
  SupportPairList nu;
  intersecting_wavelets(basis, lambda, j, true, nu);
  print_pair_list(nu);
  intersecting_wavelets(basis, lambda, j, false, nu);
  print_pair_list(nu);

  return 0;
}


void print_list(SupportList nu)
{
  for (SupportList::const_iterator it(nu.begin()); it != nu.end(); ++it) {
    cout << *it << " ";
  }
  cout << endl;
}

void print_pair_list(SupportPairList nu)
{
  for (SupportPairList::const_iterator it(nu.begin()); it != nu.end(); ++it) {
    cout << "    nu=" << it->first
         << " with support intersection "
         << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
  }
}

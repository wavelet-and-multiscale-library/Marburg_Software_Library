// tests for the class AdaptedBasis

#include <interval/s_basis.h>
#include <interval/s_support.h>
#include <interval/interval_evaluate.h>

#include <interval/adapted_basis.h>

#include <utils/random.h>
#include <ctime>


typedef WaveletTL::SBasis MultiBasis; // multiwavelet basis to adapt
typedef WaveletTL::AdaptedBasis<MultiBasis> Basis; // adapted basis
typedef WaveletTL::AdaptedBasis<MultiBasis>::Index Index; // index of the adapted basis

int main(void)
{
  MultiBasis multi_basis;
  Basis basis(&multi_basis);
  Basis::Index lambda;
  int k1, k2;
  int j, number;

  cout << "Testing AdaptedBasis ..." << endl;

  cout << "j0 = " << basis.j0() << endl;
  cout << lambda << endl;
  cout << basis.first_generator(basis.j0()) << endl;

  for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_generator(basis.j0()); ++lambda) {
    basis.primal_support(lambda, k1, k2);
    cout << "\\phi_" << lambda << "(multi-index " << *lambda.multi_index() << ") has support [" << k1 << ", " << k2 << "]" << endl;
  }
  for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_generator(basis.j0()); ++lambda) {
    basis.dual_support(lambda, k1, k2);
    cout << "\\tilde\\phi_" << lambda << "(multi-index " << *lambda.multi_index() << ") has support [" << k1 << ", " << k2 << "]" << endl;
  }

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

  for (k1 = 0; k1 < 5; k1++) { // do this test several times
    j = basis.j0() + random_integer(10);
    lambda = Index(j, E_WAVELET, random_integer(basis.DeltaLmin(),basis.DeltaRmax(j)), &basis);
    number = lambda.number();
    cout << lambda << " has number " << number << ", has index " << Index(number, &basis) << endl;
    cout << "Multi-index: " << *lambda.multi_index() << " re-constructed from number " << number << ": " << MultiBasis::Index(number, &multi_basis) << endl;
  }
  // test it the other way
  for (k1 = 0; k1 < 5; k1++) { // do this test several times
    number = random_integer(10000);
    lambda = Index(number, &basis);
    cout << "Random number " << number << " corresponds to index " << lambda << ", which is number " << lambda.number() << "." << endl;
  }

  return 0;
}

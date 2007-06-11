// tests for the class AdaptedBasis

#include <interval/s_basis.h>
#include <interval/s_support.h>
#include <interval/interval_evaluate.h>

#include <interval/adapted_basis.h>


typedef WaveletTL::SBasis MultiBasis; // multiwavelet basis to adapt
typedef WaveletTL::AdaptedBasis<MultiBasis> Basis; // adapted basis

int main(void)
{
  MultiBasis multi_basis;
  Basis adapted_basis(&multi_basis);
  Basis::Index lambda;
  int k1, k2;

  cout << "Testing AdaptedBasis ..." << endl;

  cout << "j0 = " << adapted_basis.j0() << endl;
  cout << lambda << endl;
  cout << adapted_basis.first_generator(adapted_basis.j0()) << endl;

  for (lambda = adapted_basis.first_generator(adapted_basis.j0()); lambda <= adapted_basis.last_generator(adapted_basis.j0()); ++lambda) {
    adapted_basis.primal_support(lambda, k1, k2);
    cout << "\\phi_" << lambda << "(multi-index " << *lambda.multi_index() << ") has support [" << k1 << ", " << k2 << "]" << endl;
  }
  for (lambda = adapted_basis.first_generator(adapted_basis.j0()); lambda <= adapted_basis.last_generator(adapted_basis.j0()); ++lambda) {
    adapted_basis.dual_support(lambda, k1, k2);
    cout << "\\tilde\\phi_" << lambda << "(multi-index " << *lambda.multi_index() << ") has support [" << k1 << ", " << k2 << "]" << endl;
  }

  return 0;
}

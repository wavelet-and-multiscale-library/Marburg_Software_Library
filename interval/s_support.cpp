// implementation for s_support.h

namespace WaveletTL
{
  inline
  void support(const SBasis& basis,
               const SBasis::Index& lambda,
               int& k1, int& k2)
  {
    basis.primal_support(lambda, k1, k2);
  }

  bool intersect_supports(const SBasis& basis,
                          const SBasis::Index& lambda,
                          const int m, const int a, const int b,
                          int& j, int& k1, int& k2)
  {
    // universal routine
    const int j_lambda = lambda.j() + lambda.e();
    j = std::max(j_lambda, m);
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    const int jmb = (1<<(j-m)) * b;
    const int jjk1 = (1<<(j-j_lambda)) * k1_lambda;
    if (jjk1 >= jmb)
      return false;
    else {
      const int jma = (1<<(j-m)) * a;
      const int jjk2 = (1<<(j-j_lambda)) * k2_lambda;
      if (jma >= jjk2)
        return false;
      else {
        k1 = std::max(jjk1, jma);
        k2 = std::min(jjk2, jmb);
      }
    }
    return true;
  }

  bool intersect_supports(const SBasis& basis,
                          const SBasis::Index& lambda,
                          const SBasis::Index& nu,
                          int& j, int& k1, int& k2)
  {
    // compute support of \psi_\lambda
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);

    // call intersect_supports routine for the interval 2^{-j} [k1_lambda, k2_lambda]
    return intersect_supports(basis, nu, lambda.j()+lambda.e(), k1_lambda, k2_lambda, j, k1, k2);
  }

  inline
  bool intersect_supports(const SBasis& basis,
                          const SBasis::Index& lambda,
                          const SBasis::Index& nu,
                          SBasis::Support& supp)
  {
    // universal routine
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
                              supp.j, supp.k1, supp.k2);
  }

  void intersecting_wavelets(const SBasis& basis,
                             const SBasis::Index& lambda,
                             const int j, const bool generators,
                             int& k1, int& k2)
  {
    // compute support of \psi_\lambda
    int k1_lambda, k2_lambda;
    support(basis, lambda, k1_lambda, k2_lambda);
    //const int j_lambda = lambda.j() + lambda.e(); // wavelets have granularity lambda.j()+1
    if (lambda.e() == E_WAVELET) { // re-norm with respect to lambda.j()
      k1_lambda /= 2;
      k2_lambda /= 2;
    }

    // compute leftmost and rightmost generator/wavelet that intersects with \psi_\lambda
    // k1 and k2 respective to lambda.j()
    int delta_k;
    if (j >= lambda.j()) {
      delta_k = (1<<(j-lambda.j()));
//      cout << "j-lambda.j = "<<j-lambda.j()<<", delta_k = "<<delta_k<<" -> ["<<k1_lambda*delta_k<<", "<<k2_lambda*delta_k<<"]"<<endl;
      // on level j, \psi_\lambda has support [k1_lambda*delta_k, k2_lambda*delta_l]
      // look for the first index which function support intersects this interval for at least one j-unit
      if (generators) {
        k1 = max(basis.DeltaLmin(), k1_lambda*delta_k); // a generator reaches 1 unit to the right
        k2 = min(basis.DeltaRmax(j), k2_lambda*delta_k); // a generator reaches 1 unit to the left
      }
      else { // wavelets
        k1 = max(basis.Nablamin(), k1_lambda*delta_k - 1); // a wavelets reaches 2 lambda.j()-units to the right
        k2 = min(basis.Nablamax(j), k2_lambda*delta_k); // a wavelets reaches 1 lambda.j()-unit to the left
      }
    }
    else { // (j < lambda.j())
      // a maximum of 2 k-indices can intersect
      delta_k = (1<<(lambda.j()-j));
      if (generators) {
        k1 = max(basis.DeltaLmin(), (k1_lambda-1)/delta_k); // ceil
        k2 = min(basis.DeltaRmax(j), (k2_lambda+1)/delta_k + (k2_lambda % delta_k == 0) ? 0 : 1); // models floor
      }
      else { // wavelets
        k1 = max(basis.Nablamin(), (k1_lambda-1)/delta_k); // ceil
        k2 = min(basis.Nablamax(j), (k2_lambda+2)/delta_k + (k2_lambda % delta_k == 0) ? 0 : 1);  // models floor
      }
    }
//    cout << "lambda = "<<lambda<<", j = "<<j<<": ["<<k1_lambda<<","<<k2_lambda<<"] -> "<<k1<<"..."<<k2 << endl;
  }

  void intersecting_wavelets(const SBasis& basis,
                             const SBasis::Index& lambda,
                             const int j, const bool generators,
                             std::list<SBasis::Index>& intersecting)
  {
    typedef SBasis::Index Index;
    typedef SBasis::Support Support;

    intersecting.clear();

    int k1, k2;
    intersecting_wavelets(basis, lambda, j, generators, k1, k2);

    // create list (all indices between k1 and k2)
    for (int k = k1; k <= k2; k++)
      for (int c = 0; c < (int)basis.number_of_components; c++)
        intersecting.push_back(SBasis::Index(j, !generators, k, c, &basis));
  }

  void intersecting_wavelets
    (const SBasis& basis,
     const SBasis::Index& lambda,
     const int j, const bool generators,
     std::list<std::pair<SBasis::Index, SBasis::Support> >& intersecting)
  {
    typedef SBasis::Index Index;
    typedef SBasis::Support Support;

    intersecting.clear();

    int k1, k2;
    intersecting_wavelets(basis, lambda, j, generators, k1, k2);

    // create list (all indices between k1 and k2)
    Support supp;
    Index nu;
    for (int k = k1; k <= k2; k++) {
      for (int c = 0; c < (int)basis.number_of_components; c++) {
        nu = Index(j, !generators, k, c, &basis);
        #ifdef NDEBUG // no debuging, final version
        intersect_supports(basis, nu, lambda, supp.j, supp.k1, supp.k2);
        #else // !NDEBUG (debuging)
        if (!intersect_supports(basis, nu, lambda, supp.j, supp.k1, supp.k2))
          cerr << "intersecting_wavelets: Error: intersect_supports returned false. lambda = "<<lambda<<", nu = "<<nu << endl;
        #endif // !NDEBUG
        intersecting.push_back(std::make_pair(nu, supp));
      }
    }
  }
}

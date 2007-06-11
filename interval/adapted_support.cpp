// implementation for adapted_support.h

namespace WaveletTL
{
  template <class IBASIS>
  inline
  void
  support(const AdaptedBasis<IBASIS>& basis,
          const typename AdaptedBasis<IBASIS>::Index& lambda,
          int& k1, int& k2)
  {
    basis.primal_support(lambda, k1, k2);
  }

  template <class IBASIS>
  bool intersect_supports(const AdaptedBasis<IBASIS>& basis,
                          const typename AdaptedBasis<IBASIS>::Index& lambda,
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

  template <class IBASIS>
  inline
  bool intersect_supports(const AdaptedBasis<IBASIS>& basis,
                          const typename AdaptedBasis<IBASIS>::Index& lambda,
                          const typename AdaptedBasis<IBASIS>::Index& nu,
                          typename AdaptedBasis<IBASIS>::Support& supp)
  {
    // universal routine
    int k1_nu, k2_nu;
    support(basis, nu, k1_nu, k2_nu);
    return intersect_supports(basis, lambda, nu.j()+nu.e(), k1_nu, k2_nu,
                              supp.j, supp.k1, supp.k2);
  }

  template <class IBASIS>
  void intersecting_wavelets(const AdaptedBasis<IBASIS>& basis,
                             const typename AdaptedBasis<IBASIS>::Index& lambda,
                             const int j, const bool generators,
                             std::list<typename AdaptedBasis<IBASIS>::Index>& intersecting)
  {
    std::list<typename IBASIS::Index>& intersecting_m;
    intersecting_wavelets(*basis.multi_basis(), *lambda.multi_index(), j, generators, intersecting_m);
    // convert list of multi-indices to list of adapted indices
    intersecting.clear();
    for (typename std::list<typename IBASIS::Index>::const_iterator iter = intersecting_m.begin(), it_end = intersecting_m.end(); iter != it_end; ++iter) {
      intersecting.push_back(AdaptedIndex(*iter, &basis));
    }
  }

  template <class IBASIS>
  void intersecting_wavelets(const AdaptedBasis<IBASIS>& basis,
                             const typename AdaptedBasis<IBASIS>::Index& lambda,
                             const int j, const bool generators,
                             std::list<std::pair<typename AdaptedBasis<IBASIS>::Index, typename AdaptedBasis<IBASIS>::Support> >& intersecting)
  {
    std::list<std::pair<typename IBASIS::Index, typename IBASIS::Support> >& intersecting_m;
    intersecting_wavelets(*basis.multi_basis(), *lambda.multi_index(), j, generators, intersecting_m);
    // convert list of multi-indices to list of adapted indices
    intersecting.clear();
    for (typename std::list<typename IBASIS::Index>::const_iterator iter = intersecting_m.begin(), it_end = intersecting_m.end(); iter != it_end; ++iter) {
//      IBASIS::Supp supp_m = *iter.second();
      intersecting.push_back(std::make_pair(AdaptedIndex(*iter.first(), &basis), *iter.second()));
    }
  }

}

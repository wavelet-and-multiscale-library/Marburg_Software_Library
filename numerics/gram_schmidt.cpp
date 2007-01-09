namespace MathTL
{
  template <class V>
  void gramSchmidtProcess(Array1D<V>& vectors)
  {
    V* vtmp;
    unsigned int i, n;

    for (n = 0; n < vectors.size(); n++) {
      // orthogonalize the n-th vector respective to the previous ones
      for (i = 0; i < n; i++) {
        vtmp = new V(vectors[i]);
        vtmp->scale(-vectors[i].inner_product(vectors[n]));
        vectors[n].add(*vtmp);
        delete vtmp;
        vtmp = NULL;
      }
      // norm the n-th vector
      vectors[n].scale(1/sqrt(vectors[n].inner_product(vectors[n])));
    }
  }
}

// implementations of Function<DIM>:: inline functions

#include <cassert>

namespace MathTL
{
  template <unsigned int DIM>
  inline
  Function<DIM>::Function(const unsigned int n_components,
			  const double initial_time)
    : FunctionTime(initial_time), n_components(n_components)
  {
    assert(n_components > 0);
  }
  
  template <unsigned int DIM>
  inline
  Function<DIM>::~Function()
  {
  } 
  
  template <unsigned int DIM>
  inline
  const typename Function<DIM>::size_type
  Function<DIM>::memory_consumption() const
  {
    return sizeof(*this);
  }
  
  template <unsigned int DIM>
  inline
  ZeroFunction<DIM>::ZeroFunction(const unsigned int n_components)
    : Function<DIM>(n_components)
  {
  }
  
  template <unsigned int DIM>
  inline
  ZeroFunction<DIM>::~ZeroFunction()
  {
  } 
  
  template <unsigned int DIM>
  inline
  double ZeroFunction<DIM>::value(const Point<DIM>& p,
				  const unsigned int component) const
  {
    return 0.0;
  }

  template <unsigned int DIM>
  inline
  void ZeroFunction<DIM>::vector_value(const Point<DIM> &p,
				       Array1D<double>& values) const
  {
    for (unsigned int i(0); i < n_components; i++)
      values[i] = value(p, i);
  }
}

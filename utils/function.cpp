// implementations of Function<DIM>:: inline functions

#include <cassert>

namespace MathTL
{
  template <unsigned int DIM, class VALUE>
  inline
  Function<DIM, VALUE>::Function(const unsigned int n_components,
				 const double initial_time)
    : FunctionTime(initial_time), n_components(n_components)
  {
    assert(n_components > 0);
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  Function<DIM, VALUE>::~Function()
  {
  } 
  
  template <unsigned int DIM, class VALUE>
  inline
  const typename Function<DIM, VALUE>::size_type
  Function<DIM, VALUE>::memory_consumption() const
  {
    return sizeof(*this);
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  ZeroFunction<DIM, VALUE>::ZeroFunction(const unsigned int n_components)
    : Function<DIM, VALUE>(n_components)
  {
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  ZeroFunction<DIM, VALUE>::~ZeroFunction()
  {
  } 
  
  template <unsigned int DIM, class VALUE>
  inline
  VALUE ZeroFunction<DIM, VALUE>::value(const Point<DIM,VALUE>& p,
					const unsigned int component) const
  {
    return 0.0;
  }

  template <unsigned int DIM, class VALUE>
  inline
  void ZeroFunction<DIM, VALUE>::vector_value(const Point<DIM,VALUE> &p,
					      Vector<VALUE>& values) const
  {
    values = 0;
  }

  template <unsigned int DIM, class VALUE>
  inline
  ConstantFunction<DIM, VALUE>::ConstantFunction(const Vector<VALUE>& value)
    : Function<DIM, VALUE>(value.size()), c(value)
  {
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  ConstantFunction<DIM, VALUE>::ConstantFunction(const VALUE& value)
    : Function<DIM, VALUE>((unsigned int)1), c(value)
  {
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  ConstantFunction<DIM, VALUE>::~ConstantFunction()
  {
  } 
  
  template <unsigned int DIM, class VALUE>
  inline
  VALUE ConstantFunction<DIM, VALUE>::value(const Point<DIM,VALUE>& p,
					    const unsigned int component) const
  {
    return c[component];
  }

  template <unsigned int DIM, class VALUE>
  inline
  void ConstantFunction<DIM, VALUE>::vector_value(const Point<DIM,VALUE> &p,
						  Vector<VALUE>& values) const
  {
    values = c;
  }

  template <unsigned int DIM, class VALUE>
  inline
  ProductFunction<DIM, VALUE>::ProductFunction(const Function<DIM,VALUE>* f1,
					       const Function<DIM,VALUE>* f2)
    : Function<DIM, VALUE>(f1->n_components), f1_(f1), f2_(f2)
  {
  }
  
  template <unsigned int DIM, class VALUE>
  inline
  ProductFunction<DIM, VALUE>::~ProductFunction()
  {
  } 
  
  template <unsigned int DIM, class VALUE>
  inline
  VALUE ProductFunction<DIM, VALUE>::value(const Point<DIM,VALUE>& p,
					   const unsigned int component) const
  {
    return f1_->value(p) * f2_->value(p);
  }

  template <unsigned int DIM, class VALUE>
  inline
  void ProductFunction<DIM, VALUE>::vector_value(const Point<DIM,VALUE> &p,
						 Vector<VALUE>& values) const
  {
    Vector<VALUE> v2;
    f1_->vector_value(p, values);
    f2_->vector_value(p, v2);
    for (unsigned int i = 0; i < DIM; i++)
      values[i] *= v2[i];
  }
  
}

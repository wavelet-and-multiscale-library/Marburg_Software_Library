// implementation of MathTL::Vector inline functions

#include <io/vector_io.h>
#include <cassert>
#include <algorithm>

namespace MathTL
{
  template <class C>
  inline
  Vector<C>::Vector()
    : Array1D<C>()
  {
  }

  template <class C>
  inline
  Vector<C>::Vector(const size_type s, const bool initialize)
    : Array1D<C>(s)
  {
    if (initialize)
      *this = 0;
  }

  template <class C>
  inline
  Vector<C>& Vector<C>::operator = (const C& c)
  {
    assert(size() > 0);
    std::fill(begin(), end(), c);
    return *this;
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const Vector<C>& v)
  {
    print_vector(v, os);
  }
}

// implementation for row_method.h

namespace MathTL
{
  template <class VECTOR, class IVP>
  ROWMethod<VECTOR, IVP>::ROWMethod(const typename WMethod<VECTOR, IVP>::Method method)
    : WMethod<VECTOR, IVP>(method, this)
  {
  }
  
}

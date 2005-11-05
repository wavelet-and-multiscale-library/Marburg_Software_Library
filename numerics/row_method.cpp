// implementation for row_method.h

namespace MathTL
{
  template <class VECTOR>
  ROWMethod<VECTOR>::ROWMethod(const typename WMethod<VECTOR>::Method method)
    : WMethod<VECTOR>(method, this)
  {
  }
  
}

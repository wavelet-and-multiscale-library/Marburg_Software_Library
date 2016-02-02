// implementation for ivp.h

namespace MathTL
{
  template <unsigned int DIM>
  IVP<DIM>::~IVP()
  {
  }

  template <class VECTOR>
  AbstractIVP<VECTOR>::~AbstractIVP()
  {
  }
  
  template <class VECTOR>
  AbstractCachedIVP<VECTOR>::~AbstractCachedIVP()
  {
  }
}

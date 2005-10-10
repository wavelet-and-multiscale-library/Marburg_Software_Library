// implementation for sturm_bvp.h

namespace MathTL
{
  inline
  SturmBVP::~SturmBVP()
  {
  }

  inline
  void
  SturmBVP::apply_A(const Point<2>& v, Point<2>& result) const
  {
    result[0] = alpha0()*v[0] + alpha1()*p(0)*v[1];
    result[1] = 0;
  }
   
  inline
  void
  SturmBVP::apply_B(const Point<2>& v, Point<2>& result) const
  {
    result[0] = 0;
    result[1] = beta0()*v[0] + beta1()*p(1)*v[1];
  }

  inline
  void
  SturmBVP::apply_f(const double t, const Point<2>& v, Point<2>& result) const
  {
    result[0] = v[1];
    result[1] = (-p_prime(t)*v[1]+q(t)*v[0]-g(t))/p(t);
  }

  inline
  simpleSturmBVP::~simpleSturmBVP()
  {
  }

  inline
  periodicSturmBVP::~periodicSturmBVP()
  {
  }
}

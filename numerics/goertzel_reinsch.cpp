// implementation for goertzel_reinsch.h

namespace MathTL
{
  template <class C>
  C
  Goertzel<C>::value(const Point<1,C>&p, const unsigned int component) const
  {
    // adjoint summation for the three-term recursion
    //   p_k = 2*cos(x) * p_{k-1} - p_{k-2}

    assert(coeffs_.size() >= 1);
    
    C uk(0), ukplus1(0), ukplus2(0);
    
    const C gamma = cos(p[0]);
    for (int k(coeffs_.size()-1); k >= 1; --k) {
      uk = 2*gamma*ukplus1 - ukplus2 + coeffs_[k];
      if (k > 0) {
	ukplus2 = ukplus1;
	ukplus1 = uk;
      }
    }
    
    return (sine_ ? uk*sin(p[0]) : coeffs_[0]+uk*cos(p[0])-ukplus1);
  }
  
  template <class C>
  C
  GoertzelReinsch<C>::value(const Point<1,C>&p, const unsigned int component) const
  {
    assert(coeffs_.size() >= 1);
    
    C uk(0), ukplus1(0), wk(0), wkplus1(0), lambda(0);

    if (cos(p[0])>0) {
      lambda = -4*sin(p[0]/2.)*sin(p[0]/2.);
      for (int k(coeffs_.size()-1); k >= 0; --k) {
	uk = wk + ukplus1;
	wkplus1 = wk;
	wk = (lambda*uk + wkplus1) + coeffs_[k];
	if (k > 0) {
	  ukplus1 = uk;
	}
      }
    } else {
      lambda = 4*cos(p[0]/2.)*cos(p[0]/2.);
      for (int k(coeffs_.size()-1); k >= 0; --k) {
	uk = wk - ukplus1;
	wkplus1 = wk;
	wk = (lambda*uk - wkplus1) + coeffs_[k];
	if (k > 0) {
	  ukplus1 = uk;
	}
      }
    }
    
    return sine_ ? uk*sin(p[0]) : wk - wkplus1*(lambda/2);
  }
  
}

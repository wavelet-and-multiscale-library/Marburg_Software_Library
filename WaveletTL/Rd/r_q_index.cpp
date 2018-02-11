// implementation for r_q_index.h

namespace WaveletTL
{
  inline
  RQIndex::RQIndex(const RQIndex& lambda)
    : p_(lambda.p()), j_(lambda.j()), e_(lambda.e()), k_(lambda.k())
  {
  }

  inline
  RQIndex::RQIndex(const int p, const int j, const int e, const int k)
    : p_(p), j_(j), e_(e), k_(k)
  {
  }

  inline
  RQIndex& RQIndex::operator = (const RQIndex& lambda)
  {
    p_ = lambda.p();
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();

    return *this;
  }

  inline
  bool RQIndex::operator == (const RQIndex& lambda) const
  {
    return (p_ == lambda.p() &&
	    j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  inline
  bool RQIndex::operator < (const RQIndex& lambda) const
  {
    return  (p_ < lambda.p() ||
	    (p_ == lambda.p() && j_ < lambda.j()) ||
	    (p_ == lambda.p() && j_ == lambda.j() && e_ < lambda.e()) ||
	    (p_ == lambda.p() && j_ == lambda.j() && 
	     e_ == lambda.e() && k_ < lambda.k()));
  }
}

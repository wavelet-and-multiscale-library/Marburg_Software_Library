// implementation for r_index.h

namespace WaveletTL
{
  inline
  RIndex::RIndex(const RIndex& lambda)
    : j_(lambda.j()), e_(lambda.e()), k_(lambda.k())
  {
  }

  inline
  RIndex::RIndex(const int j, const int e, const int k)
    : j_(j), e_(e), k_(k)
  {
  }

  inline
  RIndex& RIndex::operator = (const RIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();

    return *this;
  }

  inline
  bool RIndex::operator == (const RIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  inline
  bool RIndex::operator < (const RIndex& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && e_ < lambda.e()) ||
	    (j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k()));
  }
}

// implementation for r_index.h

namespace WaveletTL
{
  RIndex::RIndex(const RIndex& lambda)
    : j_(lambda.j()), e_(lambda.e()), k_(lambda.k())
  {
  }

  RIndex::RIndex(const int j, const int e, const int k)
    : j_(j), e_(e), k_(k)
  {
  }

  RIndex& RIndex::operator = (const RIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();

    return *this;
  }

  bool RIndex::operator == (const RIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  bool RIndex::operator < (const RIndex& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && e_ < lambda.e()) ||
	    (j_ == lambda.j() && e_ == lambda.e() && k_ < lambda.k()));
  }
}

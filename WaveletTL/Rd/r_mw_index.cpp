// implementation for r_mw_index.h

namespace WaveletTL
{
  RMWIndex::RMWIndex(const RMWIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    c_ = lambda.c();
    k_ = lambda.k();
  }

  RMWIndex::RMWIndex(const int j, const int e, const int c, const int k)
  {
    j_ = j;
    e_ = e;
    c_ = c;
    k_ = k;
  }

  RMWIndex& RMWIndex::operator = (const RMWIndex& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    c_ = lambda.c();
    k_ = lambda.k();

    return *this;
  }

  bool RMWIndex::operator == (const RMWIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    c_ == lambda.c() &&
	    k_ == lambda.k());
  }
  
  bool RMWIndex::operator < (const RMWIndex& lambda) const
  {
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && (c_ < lambda.c() ||
							(c_ == lambda.c() && k_< lambda.k()))))));
  }
}

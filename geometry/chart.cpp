// implementation for chart.h

#include <cassert>

namespace MathTL
{
  template <unsigned int DIM>
  inline
  void
  IdentityMapping<DIM>::map_point(const Point<DIM>& x, Point<DIM>& y) const {
    y = x;
  }
    
  template <unsigned int DIM>
  inline
  void
  IdentityMapping<DIM>::map_point_inv(const Point<DIM>& x, Point<DIM>& y) const {
    y = x;
  }

  template <unsigned int DIM>
  inline
  const double
  IdentityMapping<DIM>::Gram_factor(const Point<DIM>&) const {
    return 1.0;
  }

  template <unsigned int DIM>
  inline
  const double
  IdentityMapping<DIM>::Gram_D_factor(const unsigned int i, const Point<DIM>& x) const {
    return 0.0;
  }

  template <unsigned int DIM>
  inline
  const double
  IdentityMapping<DIM>::Dkappa_inv(const unsigned int i, const unsigned int j,
				   const Point<DIM>& x) const {
    return (i == j ? 1.0 : 0.0);
  }
  
  template <unsigned int DIM>
  const bool
  IdentityMapping<DIM>::in_patch(const Point<DIM>& x) const {
    for (unsigned int i = 0; i < DIM; i++)
      if (x[i] < 0.0 || x[i] > 1.0) return false;
    return true;
  }

  template <unsigned int DIM>
  AffineLinearMapping<DIM>::AffineLinearMapping(const Matrix<double>& A, const Point<DIM>& b)
    : A_(A), A_inv(DIM, DIM), b_(b)
  {
    assert(A.row_dimension() == DIM && A.column_dimension() == DIM);
    assert(DIM == 1 || DIM == 2);
    switch (DIM) {
    case 1:
      det_A = A_(0,0);
      A_inv(0,0) = 1./det_A;
      break;
    case 2:
      det_A = A_(0,0)*A_(1,1)-A_(1,0)*A_(0,1);
      A_inv(0,0) = A_(1,1)/det_A;
      A_inv(0,1) = -A_(0,1)/det_A;
      A_inv(1,0) = -A_(1,0)/det_A;
      A_inv(1,1) = A_(0,0)/det_A;
      break;
    }
  }

  template <unsigned int DIM>
  inline
  void
  AffineLinearMapping<DIM>::map_point(const Point<DIM>& x, Point<DIM>& y) const {
    A_.apply(x, y);
    y += b_;
  }
    
  template <unsigned int DIM>
  inline
  void
  AffineLinearMapping<DIM>::map_point_inv(const Point<DIM>& x, Point<DIM>& y) const {
    A_inv.apply(Point<DIM>(x-b_), y);
  }

  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Gram_factor(const Point<DIM>& x) const {
    return fabs(det_A);
  }

  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Gram_D_factor(const unsigned int i, const Point<DIM>& x) const {
    return 0.0;
  }
  
  template <unsigned int DIM>
  inline
  const double
  AffineLinearMapping<DIM>::Dkappa_inv(const unsigned int i, const unsigned int j,
				       const Point<DIM>& x) const
  {
    return A_inv(i, j);
  }

  template <unsigned int DIM>
  const bool
  AffineLinearMapping<DIM>::in_patch(const Point<DIM>& x) const {
    Point<DIM> y;
    map_point_inv(x, y);
    for (unsigned int i = 0; i < DIM; i++)
      if (y[i] < 0 || y[i] > 1) return false;
    return true;
  }
}

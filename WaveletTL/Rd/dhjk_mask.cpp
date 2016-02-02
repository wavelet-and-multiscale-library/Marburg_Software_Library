// implementation for dhjk_mask.h

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>
#include <utils/multiindex.h>
#include <algebra/fixed_matrix.h>

namespace WaveletTL
{
  DHJKMask_primal::DHJKMask_primal()
  {
    set_coefficient(-1, FixedMatrix<double, 2>("0.5 0.75 -0.125 -0.125"));
    set_coefficient(0, FixedMatrix<double, 2>("1.0 0 0 0.5"));
    set_coefficient(1, FixedMatrix<double, 2>("0.5 -0.75 0.125 -0.125"));
  }

  DHJKMask_dual::DHJKMask_dual()
  {
    set_coefficient(-2, FixedMatrix<double, 2>("-14 -10 87 62", 128, true));
    set_coefficient(-1, FixedMatrix<double, 2>("16 6 -99 -37", 32, true));
    set_coefficient(0, FixedMatrix<double, 2>("39 0 0 60", 32, true));
    set_coefficient(1, FixedMatrix<double, 2>("16 -6 99 -37", 32, true));
    set_coefficient(2, FixedMatrix<double, 2>("-14 10 -87 62", 128, true));
  }
}

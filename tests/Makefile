# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

MATHTL_DIR = $(HOME)/c++/MathTL

CXX = g++
CXXFLAGS += -I$(MATHTL_DIR) -Wall -pipe

EXEOBJF = \
  test_random.o\
  test_tensor.o test_point.o test_array1d.o test_fixed_array1d.o\
  test_vector.o test_infinite_vector.o test_vectorspeed.o test_matrix.o\
  test_function.o test_polynomial.o test_laurent_polynomial.o\
  test_iteratsolv.o test_eigenvalues.o test_decomp.o test_decomposable_matrix.o\
  test_ortho_poly.o test_quadrature.o test_gauss_quadrature.o\
  test_recursion.o\
  test_grid.o test_sampled_mapping.o\
  test_splines.o\
  test_rosenbrock.o\
  test_differences.o\
  test_multiindex.o test_multi_lp.o\
  test_sturm_bvp.o

EXES = $(EXEOBJF:.o=)

all:: tests

tests:: $(EXES)

clean::
	rm -f $(EXEOBJF) $(EXES)

veryclean:: clean
	rm -f *~

$(EXEOBJF): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<

$(EXES): %: %.o
	$(CXX) $< -g -o $@

# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

MATHTL_DIR = $(HOME)/c++/MathTL

CXX = g++-3.4
CXXFLAGS += -I$(MATHTL_DIR) -Wall

# mature test programs
EXEOBJF = \
  test_random.o\
  test_tensor.o test_point.o test_array1d.o\
  test_vector.o test_infinite_vector.o test_vectorspeed.o test_matrix.o\
  test_function.o test_polynomial.o test_laurent_polynomial.o\
  test_iteratsolv.o test_eigenvalues.o test_decomp.o\
  test_ortho_poly.o test_quadrature.o test_gauss_quadrature.o\
  test_recursion.o\
  test_grid.o test_sampled_mapping.o

EXES = $(EXEOBJF:.o=)

# test programs under construction
EXEOBJF2 = \
  test_splines.o\
  test_rosenbrock.o

EXES2 = $(EXEOBJF2:.o=)

all:: tests1 tests

tests1:: $(EXES)

tests:: $(EXES2)

clean::
	rm -f $(EXEOBJF) $(EXES) $(EXEOBJF2) $(EXES2)

veryclean:: clean
	rm -f *~

$(EXEOBJF): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<

$(EXEOBJF2): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<

$(EXES): %: %.o
	$(CXX) $< -g -o $@

$(EXES2): %: %.o
	$(CXX) $< -g -o $@

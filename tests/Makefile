# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2004                                             |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

MATHTL_DIR = $(HOME)/c++/MathTL

CXX = g++
CXXFLAGS += -I$(MATHTL_DIR) -Wall

EXEOBJF = \
  test_tensor.o test_point.o test_array1d.o\
  test_function.o test_polynomial.o\
  test_vector.o test_vectorspeed.o test_matrix.o\
  test_iteratsolv.o test_eigenvalues.o test_decomp.o

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

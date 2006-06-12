# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2002-2006                                        |
# | Thorsten Raasch, Manuel Werner                                 |
# +----------------------------------------------------------------+

# choose between GNU/Intel/Portland Group compiler
COMPILER = gnu
#COMPILER = intel
#COMPILER = pgi

ifeq ($(COMPILER),gnu)
CXX = g++
#CXX = g++-2.95
#CXX = $(HOME)/gcc-4.1.0/bin/g++
#CXX = $(HOME)/gcc-4.1.1/bin/g++

ifeq ($(CXX),$(HOME)/gcc-4.1.0/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.0/lib
else
ifeq ($(CXX),$(HOME)/gcc-4.1.1/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.1/lib
else
LDFLAGS =
endif
endif

CXXFLAGS += -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(MATHTL_DIR) -Wall -pipe
endif

ifeq ($(COMPILER),intel)
CXX = icpc
CXXFLAGS += -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(MATHTL_DIR) -Wall -pipe

# disable some annoying warnings of icpc
CXXFLAGS += -wd383 -wd981
endif

ifeq ($(COMPILER),pgi)
CXX = pgCC
CXXFLAGS += -I$(MATHTL_DIR)
#CXXFLAGS += -O2 -I$(MATHTL_DIR)
#CXXFLAGS += -fastsse -I$(MATHTL_DIR)
endif

EXEOBJF = \
  test_random.o\
  test_tensor.o test_point.o test_array1d.o test_fixed_array1d.o\
  test_vector.o test_infinite_vector.o test_vectorspeed.o test_matrix.o\
  test_block_matrix.o\
  test_preconditioner.o\
  test_function.o test_polynomial.o test_laurent_polynomial.o\
  test_iteratsolv.o test_eigenvalues.o test_decomp.o test_decomposable_matrix.o\
  test_ortho_poly.o test_quadrature.o test_gauss_quadrature.o test_extrapolation.o\
  test_recursion.o\
  test_grid.o test_sampled_mapping.o test_colormap.o\
  test_splines.o\
  test_one_step_scheme.o test_rosenbrock.o\
  test_differences.o\
  test_multiindex.o test_multi_lp.o\
  test_sturm_bvp.o test_bvp.o\
  test_chart.o test_atlas.o\
  test_corner_singularity.o

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
	$(CXX) $(LDFLAGS) $< -g -o $@

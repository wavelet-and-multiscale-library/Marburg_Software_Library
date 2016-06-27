# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2002-2009                                        |
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
#CXX = $(HOME)/gcc-4.2-20060610/bin/g++
#CXX = $(HOME)/gcc-4.2.0/bin/g++

ifeq ($(CXX),$(HOME)/gcc-4.1.0/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.0/lib
else
ifeq ($(CXX),$(HOME)/gcc-4.1.1/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.1/lib
else
ifeq ($(CXX),$(HOME)/gcc-4.2-20060610/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.2-20060610/lib
else
ifeq ($(CXX),$(HOME)/gcc-4.2.0/bin/g++)
#LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.2.0/lib
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.2.0/lib64
else
LDFLAGS =
endif
endif
endif
endif

TARGETDIR:=bin

#CXXFLAGS += -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O2 -I$(MATHTL_DIR) -Wall -pipe
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
 test_multiindex.o\
 test_atlas.o\
 test_multi_lp.o\
 test_random.o test_tools.o\
 test_tensor.o test_point.o test_array1d.o test_fixed_array1d.o\
 test_vector.o test_infinite_vector.o test_vectorspeed.o test_matrix.o\
 test_block_matrix.o test_qs_matrix.o test_qs_matrixspeed.o\
 test_preconditioner.o\
 test_function.o test_polynomial.o test_laurent_polynomial.o\
 test_iteratsolv.o test_eigenvalues.o test_decomp.o test_decomposable_matrix.o\
 test_ortho_poly.o test_goertzel_reinsch.o\
 test_quadrature.o\
 test_gauss_quadrature.o\
 test_extrapolation.o\
 test_recursion.o\
 test_grid.o test_sampled_mapping.o test_colormap.o\
 test_splines.o test_bezier.o test_up_function.o\
 test_rosenbrock.o\
 test_differences.o\
 test_sturm_bvp.o test_bvp.o\
 test_chart.o\
 test_corner_singularity.o\
 test_map_tools.o\
 test_corner_singularity_biharmonic.o\
 test_gram_schmidt.o test_piecewise.o\
 test_fixed_vector.o test_fixed_matrix.o\
 test_cardinalsplines.o
 

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

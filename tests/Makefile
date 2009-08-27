# +----------------------------------------------------------------+
# | Makefile for some FrameTL test programs                        |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Manuel Werner                                                  |
# +----------------------------------------------------------------+

#DEBUGMODE = on
DEBUGMODE = off

# choose between GNU/Intel/Portland Group compiler
COMPILER = gnu
#COMPILER = intel
#COMPILER = pgi
#COMPILER = mpi

ifeq ($(COMPILER),mpi)
CXX = mpiCC.mpich
CXXFLAGS +=  -O2 -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS +=  -O2 -I$(FRAMETL_DIR_PAR) -I$(WAVELETTL_DIR_PAR) -I$(MATHTL_DIR_PAR) -Wall -pipe
endif

ifeq ($(COMPILER),gnu)
CXX = g++
#CXX = g++-2.95
#CXX = $(HOME)/gcc-4.1.0/bin/g++
#CXX = $(HOME)/gcc-4.1.1/bin/g++
#CXX = $(HOME)/gcc-4.2-20060610/bin/g++
ifeq ($(CXX),g++-2.95)
CXXFLAGS += -ftemplate-depth-40
endif

ifeq ($(CXX),$(HOME)/gcc-4.1.0/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.0/lib
#LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.0/lib64
else
ifeq ($(CXX),$(HOME)/gcc-4.1.1/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.1/lib
#LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.1.1/lib64
else
ifeq ($(CXX),$(HOME)/gcc-4.2-20060610/bin/g++)
LDFLAGS += -Xlinker -rpath -Xlinker $(HOME)/gcc-4.2-20060610/lib
else
LDFLAGS =
endif
endif
endif

#CXX = mpiCC

#CXX = icpc
CXXFLAGS +=  -O2 -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -mcpu=pentium4 -msse3 -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -fast -mcpu=7450 -faltivec -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
endif


EXEOBJF = test_aggregated_frame.o\
test_biharmonic_1d_testcase.o\
test_elliptic_equation_1D.o\
test_elliptic_equation.o\
test_simple_elliptic_equation.o\
test_simple_biharmonic_equation_1D.o\
test_simple_biharmonic_equation_2D.o\
test_biharmonic_equation_1D.o\
test_biharmonic_equation.o\
test_additive_Schwarz33_2D.o\
test_multiplicative_Schwarz33_biharmonic_2D_v2.o\
test_multiplicative_Schwarz46_biharmonic_2D.o\
test_multiplicative_Schwarz33_2D_ring.o\
test_additive_Schwarz33_2D_ring.o\
test_additive_Schwarz33_1D.o\
test_multiplicative_Schwarz22_2D.o\
test_multiplicative_Schwarz33_2D_Full.o\
test_multiplicative_Schwarz33_1D.o\
test_multiplicative_Schwarz22_1D.o\
test_multiplicative_Schwarz46_1D.o\
test_multiplicative_Schwarz33_2D.o\
test_multiplicative_Schwarz46_2D.o\
test_multiplicative_Schwarz33_2D_j03.o\
test_steepest_descent_biharmonic_1D.o\
test_steepest_descent33.o\
test_steepest_descent46.o\
test_steepest_descent33_basis.o\
test_richardson_2D.o\
test_steepest_descent22_2D.o\
test_steepest_descent46_2D.o\
test_richardson.o\
test_steepest_descent22.o\
test_additive_Schwarz.o\
test_steepest_descent_2D.o\
test_multiplicative_Schwarz.o\
test_apply_1D.o \
test_apply.o \
test_compute_dual_frame_elements.o


EXES = $(EXEOBJF:.o=)

all:: tests

tests:: $(EXES)

clean::
	rm -f $(EXEOBJF) $(EXES)

veryclean:: clean
	rm -f *~

$(EXEOBJF): %.o: %.cpp
#	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(EXES): %: %.o
#	$(CXX) $< -g -o $@
	$(CXX) $< -o $@

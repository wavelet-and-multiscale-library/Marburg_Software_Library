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
#CXXFLAGS += -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -I /home/werner/PROMOTION/Cuba-1.2 -Wall -pipe
#CXXFLAGS += -O3 -mcpu=pentium4 -msse3 -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -fast -mcpu=7450 -faltivec -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
endif


EXEOBJF = test_steepest_descent.o\
test_steepest_descent_2D.o\
test_simple_elliptic_equation.o\
test_multiplicative_Schwarz.o\
test_multiplicative_Schwarz_1D.o\
test_elliptic_equation_1D.o\
test_elliptic_equation.o\
test_biharmonic_equation_1D.o



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

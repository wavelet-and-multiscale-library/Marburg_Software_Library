# +----------------------------------------------------------------+
# | Makefile for some FrameTL test programs                        |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Manuel Werner                                                  |
# +----------------------------------------------------------------+

CXX = g++
#CXX = mpiCC

#CXX = icpc
#CXXFLAGS += -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -I /home/werner/PROMOTION/Cuba-1.2 -Wall -pipe
#CXXFLAGS += -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -I /home/werner/PROMOTION/Cuba-1.2 -Wall -pipe
#CXXFLAGS += -O3 -mcpu=pentium4 -msse3 -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -fast -mcpu=7450 -faltivec -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

//CXXFLAGS += -L /ia32/usr/lib

# ++++++++ additionals for Cuba quadrature library +++++++
#CXXFLAGS += -I /home/werner/PROMOTION/Cuba-1.2
#CXXFLAGS += -L /home/werner/PROMOTION/Cuba-1.2/libcuba.a
#CXXFLAGS += -L -lm

EXEOBJF = test_biharmonic_equation_1D.o\
test_elliptic_equation_1D.o\
test_elliptic_equation.o

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
	$(CXX) $< -o $@ #/home/werner/PROMOTION/Cuba-1.2/libcuba.a

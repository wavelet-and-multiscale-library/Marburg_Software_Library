# +----------------------------------------------------------------+
# | Makefile for some WaveletTL test programs                      |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

CXX = g++
#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -funroll-loops -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -fast -mcpu=7450 -faltivec -funroll-loops -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

EXEOBJF = \
  test_refinable.o\
  test_dm.o\
  test_r_basis.o test_cdf_expansion.o\
  test_ds_basis.o test_ds_support.o test_ds_expansion.o\
  test_periodic.o\
  test_sturm_bvp.o test_cdd1.o test_cdd2.o test_duv.o

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

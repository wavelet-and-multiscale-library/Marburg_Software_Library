# +----------------------------------------------------------------+
# | Makefile for some WaveletTL test programs                      |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

include ../config

CXX = g++
CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall

EXEOBJF = \
  test_refinable.o\
  test_r_basis.o\
  test_periodic.o\
  test_dku.o

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

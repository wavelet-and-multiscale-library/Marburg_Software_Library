# +----------------------------------------------------------------+
# | Makefile for some MathTL test programs                         |
# |                                                                |
# | Copyright (c) 2004                                             |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

MATHTL_DIR = $(HOME)/c++/MathTL

CXX = g++
CXXFLAGS += -I$(MATHTL_DIR)

EXEOBJF = \
  test_tensor.o test_point.o test_array1d.o

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

# +----------------------------------------------------------------+
# | Makefile for some FrameTL test programs                        |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Manuel Werner                                                  |
# +----------------------------------------------------------------+

CXX = g++
#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O2 -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -mcpu=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -fast -mcpu=7450 -faltivec -funroll-loops -I$(FRAMETL_DIR) -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

EXEOBJF = \
  test_parametrization.o test_atlas.o test_frame_index.o \
  test_aggregated_frame.o


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

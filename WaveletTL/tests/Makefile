# +----------------------------------------------------------------+
# | Makefile for some WaveletTL test programs                      |
# |                                                                |
# | Copyright (c) 2002-2009                                        |
# | Thorsten Raasch, Manuel Werner                                 |
# +----------------------------------------------------------------+

DEBUGMODE = on
#DEBUGMODE = off

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

#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -msse2 -ffast-math -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -ffast-math -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse2 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse2 -ffast-math -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -mtune=opteron -msse2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -mtune=opteron -msse2 -ffast-math -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -march=pentium4 -msse2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -march=pentium4 -msse2 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -march=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -mtune=7450 -faltivec -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -mtune=7450 -faltivec -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
endif

ifeq ($(COMPILER),intel)
CXX = icpc
#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -mp1 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

# disable some annoying warnings of icpc
CXXFLAGS += -wd383 -wd981
endif

ifeq ($(COMPILER),pgi)
CXX = pgCC
#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR)
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR)
#CXXFLAGS += -O3 -tp amd64e -fastsse -I$(WAVELETTL_DIR) -I$(MATHTL_DIR)
CXXFLAGS += -O3 -Mvect=assoc,cachesize:1048576 -Mcache_align -fast -I$(WAVELETTL_DIR) -I$(MATHTL_DIR)
endif

# set 1 of test programs: stuff on R and R^d
EXEOBJF1 = \

  

# set 2 of test programs: wavelet bases on the interval ([DS],[P],[JL],[A],[S])
EXEOBJF2 = \
  
  

# set 2a of test programs: wavelet bases on the interval with improved performance (periodic)
EXEOBJF2a = \
  test_periodic_gramian.o\
  test_periodic_laplacian.o
  

# set 3 of test programs: wavelet bases on general higher-dim. domains ((mapped) cube, tensor prod.)
EXEOBJF3 = \
  test_tbasis.o\
  test_tbasis_support.o\
  test_tbasis_index.o\
  test_tbasis_indexplot.o\
  test_tbasis_index2.o\
  test_tbasis_evaluate.o\
  solve1.o\
  test_tframe_index.o\

# set 3a of test programs: wavelet bases on the ring domain
EXEOBJF3a = \
   

# set 3b of test programs: generic tensor product wavelet bases
EXEOBJF3b = \
  

# set 4 of test programs: wavelet bases on the L-domain
EXEOBJF4 = \
  

# set 5 of test programs: adaptive wavelet schemes for elliptic equations
EXEOBJF5 = \
  
  
  


# set 6 of test programs: adaptive wavelet schemes for parabolic equations
EXEOBJF6 = \
  

# set 7 of test programs: nonadaptive wavelet solvers for elliptic problems	
EXEOBJF7 = \
  test_sturm_bvp.o\
  

# set 8 of test programs: ExpandASPP	
EXEOBJF8 = \
  
  

EXES1 = $(EXEOBJF1:.o=)
EXES2 = $(EXEOBJF2:.o=)
EXES2a = $(EXEOBJF2a:.o=)
EXES3 = $(EXEOBJF3:.o=)
EXES3a = $(EXEOBJF3a:.o=)
EXES3b = $(EXEOBJF3b:.o=)
EXES4 = $(EXEOBJF4:.o=)
EXES5 = $(EXEOBJF5:.o=) $(EXEOBJF5b:.o=)
EXES6 = $(EXEOBJF6:.o=)
EXES7 = $(EXEOBJF7:.o=)
EXES8 = $(EXEOBJF8:.o=) $(EXEOBJF5b:.o=)

all:: tests

tests:: tests1 tests2 tests2a tests3 tests3a tests3b tests4 tests5 tests6 tests7 tests8

tests1:: $(EXES1)
tests2:: $(EXES2)
tests2a:: $(EXES2a)
tests3:: $(EXES3)
tests3a:: $(EXES3a)
tests3b:: $(EXES3b)
tests4:: $(EXES4)
tests5:: $(EXES5)
tests6:: $(EXES6)
tests7:: $(EXES7)
tests8:: $(EXES8)

clean::
	rm -f $(EXEOBJF1) $(EXES1)
	rm -f $(EXEOBJF2) $(EXES2)
	rm -f $(EXEOBJF2a) $(EXES2a)
	rm -f $(EXEOBJF3) $(EXES3)
	rm -f $(EXEOBJF3a) $(EXES3a)
	rm -f $(EXEOBJF3b) $(EXES3b)
	rm -f $(EXEOBJF4) $(EXES4)
	rm -f $(EXEOBJF5) $(EXES5)
	rm -f $(EXEOBJF6) $(EXES6)
	rm -f $(EXEOBJF7) $(EXES7)
	rm -f $(EXEOBJF8) $(EXES8)

veryclean:: clean
	rm -f *~

$(EXEOBJF1): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES1): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF2): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES2): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF2a): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES2a): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF3): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES3): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF3a): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES3a): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF3b): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES3b): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF4): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES4): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF5): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES5): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF6): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES6): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF7): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES7): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF8): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -DNDEBUG -o $@ $<
endif

$(EXES8): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif


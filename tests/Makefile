# +----------------------------------------------------------------+
# | Makefile for some WaveletTL test programs                      |
# |                                                                |
# | Copyright (c) 2002-2006                                        |
# | Thorsten Raasch, Manuel Werner                                 |
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

ifeq ($(CXX),g++-2.95)
CXXFLAGS += -ftemplate-depth-40
endif

#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse2 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
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
CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

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
  test_refinable.o test_regularity.o test_dm.o\
  test_r_basis.o test_cdf_expansion.o

# set 2 of test programs: wavelet bases on the interval ([DS],[P],periodic)
EXEOBJF2 = \
  test_ds_basis.o test_ds_support.o test_ds_expansion.o\
  test_p_basis.o test_p_support.o test_p_expansion.o\
  test_periodic.o\

# set 3 of test programs: wavelet bases on higher-dim. domains ((mapped) cube, tensor prod., L-domain)
EXEOBJF3 = \
  test_cube_index.o test_cube_basis.o test_mapped_cube_basis.o\
  test_tp_basis.o\
  test_lshaped_basis.o

# set 4 of test programs: adaptive wavelet schemes for elliptic equations
EXEOBJF4 = \
  test_sturm_bvp.o test_cdd1.o test_cdd1_cube.o test_cdd2.o test_duv.o\
  test_cube_equation.o

# set 5 of test programs: adaptive wavelet schemes for parabolic equations
EXEOBJF5 = \
  test_lin_par_equation.o

EXES1 = $(EXEOBJF1:.o=)
EXES2 = $(EXEOBJF2:.o=)
EXES3 = $(EXEOBJF3:.o=)
EXES4 = $(EXEOBJF4:.o=)
EXES5 = $(EXEOBJF5:.o=)

all:: tests

tests:: tests1 tests2 tests3 tests4 tests5

tests1:: $(EXES1)
tests2:: $(EXES2)
tests3:: $(EXES3)
tests4:: $(EXES4)
tests5:: $(EXES5)

clean::
	rm -f $(EXEOBJF1) $(EXES1)
	rm -f $(EXEOBJF2) $(EXES2)
	rm -f $(EXEOBJF3) $(EXES3)
	rm -f $(EXEOBJF4) $(EXES4)
	rm -f $(EXEOBJF5) $(EXES5)

veryclean:: clean
	rm -f *~

$(EXEOBJF1): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES1): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

$(EXEOBJF2): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES2): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

$(EXEOBJF3): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES3): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

$(EXEOBJF4): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES4): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

$(EXEOBJF5): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES5): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

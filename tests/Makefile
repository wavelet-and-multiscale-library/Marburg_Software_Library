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

CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
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
  test_refinable.o test_regularity.o test_dm.o\
  test_r_basis.o test_cdf_expansion.o

# set 2 of test programs: wavelet bases on the interval ([DS],[P],[JL],periodic)
EXEOBJF2 = \
  test_ds_basis.o test_ds_support.o test_ds_expansion.o\
  test_p_basis.o test_p_support.o test_p_expansion.o\
  test_jl_basis.o test_jl_support.o\
  test_periodic.o\
  test_spline_basis.o

# set 3 of test programs: wavelet bases on general higher-dim. domains ((mapped) cube, tensor prod.)
EXEOBJF3 = \
  test_cube_index.o test_cube_support.o test_cube_basis.o test_mapped_cube_basis.o\
  test_tp_basis.o

# set 4 of test programs: wavelet bases on the L-domain
EXEOBJF4 = \
  test_ldomain_index.o test_ldomain_support.o test_ldomain_basis.o\
  test_ldomain_equation.o

# set 5 of test programs: adaptive wavelet schemes for elliptic equations
EXEOBJF5 = \
  test_sturm_bvp.o test_cdd1.o test_cdd1_cube.o test_cdd2.o test_duv.o\
  test_cube_equation.o test_cdd1_ldomain.o\
  test_cdd1_gramian.o test_cdd1_poisson1d.o test_cdd1_helm1d.o\
  test_helmholtz_equation.o

# set 6 of test programs: adaptive wavelet schemes for parabolic equations
EXEOBJF6 = \
  test_lin_par_equation.o test_row_stage_equation.o test_row_stage_equation_2d.o

# set 7 of test programs: nonadaptive wavelet solvers for elliptic problems	
EXEOBJF7 = \
  test_full_laplacian.o test_full_laplacian_speed.o\
  test_full_gramian.o test_full_helmholtz.o		


EXES1 = $(EXEOBJF1:.o=)
EXES2 = $(EXEOBJF2:.o=)
EXES3 = $(EXEOBJF3:.o=)
EXES4 = $(EXEOBJF4:.o=)
EXES5 = $(EXEOBJF5:.o=)
EXES6 = $(EXEOBJF6:.o=)
EXES7 = $(EXEOBJF7:.o=)

all:: tests

tests:: tests1 tests2 tests3 tests4 tests5 tests6 tests7

tests1:: $(EXES1)
tests2:: $(EXES2)
tests3:: $(EXES3)
tests4:: $(EXES4)
tests5:: $(EXES5)
tests6:: $(EXES6)
tests7:: $(EXES7)

clean::
	rm -f $(EXEOBJF1) $(EXES1)
	rm -f $(EXEOBJF2) $(EXES2)
	rm -f $(EXEOBJF3) $(EXES3)
	rm -f $(EXEOBJF4) $(EXES4)
	rm -f $(EXEOBJF5) $(EXES5)
	rm -f $(EXEOBJF6) $(EXES6)
	rm -f $(EXEOBJF7) $(EXES7)

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
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF2): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES2): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF3): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES3): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

$(EXEOBJF4): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
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
	$(CXX) $(CXXFLAGS) -c -o $@ $<
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
	$(CXX) $(CXXFLAGS) -c -o $@ $<
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
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES7): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $(LDFLAGS) $< -g -o $@
else
	$(CXX) $(LDFLAGS) $< -o $@
endif

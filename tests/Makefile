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

ifeq ($(CXX),g++-2.95)
CXXFLAGS += -ftemplate-depth-40
endif

#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -march=pentium4 -msse2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -march=pentium4 -msse3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -msse3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O4 -march=pentium4 -msse3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -mtune=7450 -faltivec -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -mtune=7450 -faltivec -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
endif

ifeq ($(COMPILER),intel)
CXX = icpc
#CXXFLAGS += -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O2 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
#CXXFLAGS += -O3 -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe
CXXFLAGS += -O3 -funroll-loops -fomit-frame-pointer -I$(WAVELETTL_DIR) -I$(MATHTL_DIR) -Wall -pipe

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


EXEOBJF = \
  test_refinable.o test_regularity.o\
  test_dm.o\
  test_r_basis.o test_cdf_expansion.o\
  test_ds_basis.o test_ds_support.o test_ds_expansion.o\
  test_periodic.o\
  test_p_basis.o test_p_support.o test_p_expansion.o\
  test_sturm_bvp.o test_cdd1.o test_cdd1_cube.o test_cdd2.o test_duv.o\
  test_tp_basis.o test_cube_basis.o test_mapped_cube_basis.o\
  test_cube_index.o test_cube_equation.o\
  test_lin_par_equation.o

EXES = $(EXEOBJF:.o=)

all:: tests

tests:: $(EXES)

clean::
	rm -f $(EXEOBJF) $(EXES)

veryclean:: clean
	rm -f *~

$(EXEOBJF): %.o: %.cpp
ifeq ($(DEBUGMODE),on)
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<
else
	$(CXX) $(CXXFLAGS) -c -o $@ $<
endif

$(EXES): %: %.o
ifeq ($(DEBUGMODE),on)
	$(CXX) $< -g -o $@
else
	$(CXX) $< -o $@
endif

# +----------------------------------------------------------------+
# | Makefile for WaveletTL - Wavelet Template Library for C++      |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Thorsten Raasch                                                |
# +----------------------------------------------------------------+

all:: tests

tests::
	cd tests; $(MAKE)

doc::
	cd doc; $(MAKE)

clean::
	cd tests; $(MAKE) clean

veryclean:: clean
	cd tests; $(MAKE) veryclean

countlines::
	wc -l Rd/*.{h,cpp} interval/*.{h,cpp} galerkin/*.{h,cpp} adaptive/*.{h,cpp} generic/*.{h,cpp} cube/*.{h,cpp} tests/*.cpp

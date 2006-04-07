# +----------------------------------------------------------------+
# | Makefile for WaveletTL - Wavelet Template Library for C++      |
# |                                                                |
# | Copyright (c) 2002-2006                                        |
# | Thorsten Raasch, Manuel Werner                                 |
# +----------------------------------------------------------------+

all:: tests

tests::
	cd tests; $(MAKE)

tests1::
	cd tests; $(MAKE) tests1

tests2::
	cd tests; $(MAKE) tests2

tests3::
	cd tests; $(MAKE) tests3

tests4::
	cd tests; $(MAKE) tests4

tests5::
	cd tests; $(MAKE) tests5

doc::
	cd doc; $(MAKE)

clean::
	cd tests; $(MAKE) clean

veryclean:: clean
	cd tests; $(MAKE) veryclean

countlines::
	wc -l Rd/*.{h,cpp} interval/*.{h,cpp} galerkin/*.{h,cpp} adaptive/*.{h,cpp} generic/*.{h,cpp} cube/*.{h,cpp} tests/*.cpp

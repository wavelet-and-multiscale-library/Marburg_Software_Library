# +----------------------------------------------------------------+
# | Makefile for FrameTL - Mathematical Template Library for C++   |
# |                                                                |
# | Copyright (c) 2002-2005                                        |
# | Manuel Werner                                                  |
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
	wc -l  tests/*.cpp ./*.{h,cpp}

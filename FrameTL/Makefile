# +----------------------------------------------------------------+
# | Makefile for FrameTL - Mathematical Template Library for C++   |
# |                                                                |
# | Copyright (c) 2002-2010                                        |
# | Manuel Werner                                                  |
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


doc::
	cd doc; $(MAKE) doc

docclean::
	cd doc; $(MAKE) docclean

clean::
	cd tests; $(MAKE) clean

clean1::
	cd tests; $(MAKE) clean1
clean2::
	cd tests; $(MAKE) clean2
clean3::
	cd tests; $(MAKE) clean3



veryclean:: clean
	cd tests; $(MAKE) veryclean

countlines::
	wc -l  tests/*.cpp ./*.{h,cpp}

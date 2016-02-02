# +----------------------------------------------------------------+
# | Makefile for MathTL - Mathematical Template Library for C++    |
# |                                                                |
# | Copyright (c) 2002-2009                                        |
# | Thorsten Raasch, Manuel Werner                                 |
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
	wc -l algebra/*.{h,cpp} geometry/*.{h,cpp} io/*.h numerics/*.{h,cpp} tests/*.cpp utils/*.{h,cpp}

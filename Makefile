# +----------------------------------------------------------------+
# | Makefile for MathTL - Mathematical Template Library for C++    |
# |                                                                |
# | Copyright (c) 2004                                             |
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

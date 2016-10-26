all:: repository

repository::
	cd MathTL; $(MAKE)
	cd WaveletTL; $(MAKE)
	cd Examples; $(MAKE)
	
	

clean::
	cd MathTL; $(MAKE) clean
	cd WaveletTL; $(MAKE) clean
	cd Examples; $(MAKE) clean

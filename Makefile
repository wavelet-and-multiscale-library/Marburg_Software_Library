all:: repository

repository::
	cd MathTL; $(MAKE)
	cd WaveletTL; $(MAKE)
	cd FrameTL; $(MAKE)
	cd Examples; $(MAKE)
	
	

clean::
	cd MathTL; $(MAKE) clean
	cd WaveletTL; $(MAKE) clean
	cd FrameTL; $(MAKE) clean
	cd Examples; $(MAKE) clean

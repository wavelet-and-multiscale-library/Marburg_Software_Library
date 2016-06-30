all:: repository

repository::
	cd MathTL; $(MAKE)
	cd Examples; $(MAKE)
	
	

clean::
	cd MathTL; $(MAKE) clean
	cd Examples; $(MAKE) clean

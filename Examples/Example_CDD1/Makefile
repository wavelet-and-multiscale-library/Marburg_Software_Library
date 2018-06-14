CXX = g++

CXXFLAGS = -O2 -Wall -pipe -I./ -I./../../WaveletTL -I./../../MathTL 
LDFLAGS =


EXEOBJF = \
 Example_CDD1.o


EXES = $(EXEOBJF:.o=)

all:: tests

tests:: $(EXES)

clean::
	rm -f $(EXEOBJF) $(EXES)
	rm -f *~

$(EXEOBJF): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -g -fopenmp -o $@ $<

$(EXES): %: %.o
	$(CXX) $(LDFLAGS) $< -g -fopenmp -o $@


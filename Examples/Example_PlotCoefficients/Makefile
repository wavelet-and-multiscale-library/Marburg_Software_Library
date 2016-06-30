CXX = g++

CXXFLAGS = -O2 -Wall -pipe -I./ -I./../../WaveletTL -I./../../MathTL 
LDFLAGS =


EXEOBJF = \
 Example_PlotCoefficients.o


EXES = $(EXEOBJF:.o=)

all:: tests

tests:: $(EXES)

clean::
	rm -f $(EXEOBJF) $(EXES)
	rm -f *~

$(EXEOBJF): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -g -o $@ $<

$(EXES): %: %.o
	$(CXX) $(LDFLAGS) $< -g -o $@


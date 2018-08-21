SDIR = source
ODIR = object

CXX = g++
CXXFLAGS = -g -O2 -Wcomments -Wall -std=c++11 -fopenmp -Wextra -lm -lquadmath -fext-numeric-literals
LDFLAGS  = -g -O2 -Wcomments -Wall -std=c++11 -fopenmp -Wextra -lm -lquadmath -fext-numeric-literals

_OBJ = gtopw.o auxfun1.o auxfun2.o caps.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

xgtopw: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(ODIR)/gtopw.o: $(SDIR)/gtopw.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/auxfun1.o: $(SDIR)/auxfun1.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/auxfun2.o: $(SDIR)/auxfun2.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/caps.o: $(SDIR)/caps.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)


clean:
	rm -f $(ODIR)/*.o
	rm -f xgtopw

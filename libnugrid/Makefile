CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11

all: main


main: grid.cpp grid.h
	$(CXX) $(CXXFLAGS) -c grid.cpp -o libnugrid.a
	
.PHONY: clean 
clean:
	-rm -f *.o main *.a

.PHONY: install
install:
	cp libnugrid.a ../tests
	cp grid.h ../tests/libnugrid.h


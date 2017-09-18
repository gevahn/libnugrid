CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11 -O3
INCLUDE = -Iinclude
LIBS = -lnugrid

all: test

test: test.cpp
	$(CXX) $(CXXFLAGS)  test.cpp -o test $(INCLUDE) $(LIBS)
	
.PHONY: clean 
clean:
	-rm -f *.o test *.a


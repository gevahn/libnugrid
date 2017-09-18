CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11 -O3
INCLUDE = -Iinclude
LIBS = -Llib -lnugrid

all: test

test: test.cpp
	$(CXX) $(CXXFLAGS)  test.cpp -o test.exe $(INCLUDE) $(LIBS)
	
.PHONY: clean 
clean:
	-rm -f *.o test.exe *.a


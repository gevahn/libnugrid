CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11

all: main

main: main.cpp libnugrid.h
	$(CXX) $(CXXFLAGS)  main.cpp -o main -L. -lnugrid
	
.PHONY: clean 
clean:
	-rm -f *.o main *.a


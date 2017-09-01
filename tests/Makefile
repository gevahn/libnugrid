CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11

all: main

main: main.cpp
	$(CXX) $(CXXFLAGS)  main.cpp -o main -lnugrid
	
.PHONY: clean 
clean:
	-rm -f *.o main *.a


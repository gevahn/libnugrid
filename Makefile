CXX = g++ 
CXXFLAGS = -g -Wall -std=c++11

all: main

main: main.o grid.o
	$(CXX) $(CXXFLAGS) -o main main.o grid.o

main.o: main.cpp grid.h
	$(CXX) $(CXXFLAGS) -c main.cpp

grid.o: grid.h
	$(CXX) $(CXXFLAGS) -c grid.cpp
	
clean:
	rm *.o main


#include "nugrid.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

enum {X, Y, Z, W};

double gaussian(std::vector<double> meanPoint, double variance, std::vector<double> evalPoint) {	

	double distance = 0;
	for (int i = 0; i < 3; i++) {
		distance += (meanPoint[i] - evalPoint[i]) * (meanPoint[i] - evalPoint[i]);
	}
	
	return 1 / (std::pow(2 * M_PI * variance, 3.0/2.0)) * std::exp(-distance/(2 * variance));
}

double integrate(std::vector< std::vector<double> > grid, std::vector< std::vector<double> > atomPos) {

	double sum = 0;
	double variance = 0.1;
	for (int i = 0; i < grid[0].size(); i++) {
		for (int j = 0; j < atomPos.size(); j++) {
			sum += gaussian(atomPos[j], variance, {grid[X][i], grid[Y][i], grid[Z][i]}) * grid[W][i];
		}
	}
	return sum;
}

int main(int argc, char * argv[]) {
	std::vector<int> atomNames;
	std::vector<std::vector<double> > atomPos;

	atomNames.push_back(8);
	atomNames.push_back(1);
	atomNames.push_back(1);

	atomPos.push_back({0,0,0});
	atomPos.push_back({0,1,0});
	atomPos.push_back({0,0,1});

	std::vector<std::vector<double> > grid;
	try{
	std::vector<double> paramVector =  std::vector<double> (atomPos.size()* 2);
	for (int i = 0; i < paramVector.size(); i++) {
		paramVector[i] = std::stoi(argv[1]);
		paramVector[i + 1] = std::stoi(argv[2]);
	}
	grid = libnugrid::genGrid(atomNames, atomPos, libnugrid::gridTypes::LEB_EM, paramVector);
	}
	catch (char const* param) { std::cout << param << std::endl;}

//	libnugrid::printGrid(grid);
	double gaussIntgrated = integrate(grid, atomPos);
	std::cout << "Integrating" << std::endl;
	std::cout << gaussIntgrated << std::endl;
	std::cout << gaussIntgrated - atomPos.size() << std::endl;
	assert(abs(gaussIntgrated - atomPos.size()) < 0.00001);
	return 0;

}

#include "grid.h"
#include <vector>
#include <iostream>

int main() {
	std::vector<int> atomNames;
	std::vector<std::vector<double> > atomPos;

	atomNames.push_back(8);
	atomNames.push_back(1);
	atomNames.push_back(1);

	atomPos.push_back({0,0,0});
	atomPos.push_back({0,1,0});
	atomPos.push_back({0,0,1});


	std::vector<std::vector<double> > grid;
	grid = libnugrid::genGrid(atomNames, atomPos, "Rect", std::vector<double> (1,0.5));

	libnugrid::printGrid(grid);
}

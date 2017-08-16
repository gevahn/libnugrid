#include "grid.h"
#include <limits>
#include <cmath>
#include <iostream>

using namespace libnugrid;

std::vector< std::vector< double >> libnugrid::genGrid(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, 
             std::string gridType, std::vector<double> gridParams) {

	if (gridType == "Rect")  return genGridRect(atomName, atomPos, gridParams);
	throw "Unknown grid type";	
}


std::vector< std::vector< double >> libnugrid::genGridRect(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams) {


	//get Axis Aligned Bounding Box (AABB)
	double inf = std::numeric_limits<double>::infinity();
	double maxX = -inf, minX = inf, maxY = -inf, minY = inf, maxZ = -inf, minZ = inf;

	for (auto& pos: atomPos) {
		if (pos[0] > maxX) maxX = pos[0];
		if (pos[0] < minX) minX = pos[0];
		if (pos[1] > maxY) maxY = pos[1];
		if (pos[1] < minY) minY = pos[1];
		if (pos[2] > maxZ) maxZ = pos[2];
		if (pos[2] < minZ) minZ = pos[2];
	}

	// add padding to bounding box of 2 angstrom
	maxX = maxX + 2;	
	minX = minX - 2;	
	maxY = maxY + 2;	
	minY = minY - 2;	
	maxZ = maxZ + 2;	
	minZ = minZ - 2;	

	// set grid spacing
	int gridDimX, gridDimY, gridDimZ;
	if (gridParams.size() == 1) {
		gridParams.push_back(gridParams[0]);	
		gridParams.push_back(gridParams[0]);	
	}
	gridDimX = std::round((maxX - minX) / gridParams[0]);	
	gridDimY = std::round((maxY - minY) / gridParams[1]);	
	gridDimZ = std::round((maxZ - minZ) / gridParams[2]);	

	int gridSize = gridDimX * gridDimY * gridDimZ;

	//generate grid
	std::vector<std::vector<double> > gridVector(4, std::vector<double> (gridSize, 0));

	for (int i = 0; i < gridDimX; i++) {
		for (int j = 0; j < gridDimY; j++) {
			for (int k = 0; k < gridDimZ; k++) {
				int currentIdx = k + j * gridDimY + i * (gridDimZ * gridDimY);
				gridVector[0][currentIdx] = i * gridParams[0];
				gridVector[1][currentIdx] = j * gridParams[1];
				gridVector[2][currentIdx] = k * gridParams[2];
				gridVector[3][currentIdx] = gridParams[0] * gridParams[1] * gridParams[2];
			}	
		}	
	}	
	return gridVector;
}


//grid printer
void libnugrid::printGrid(const std::vector<std::vector<double> >& grid) {
        for (int i = 0; i < grid[0].size(); i++) {
                std::cout << grid[0][i] << '\t';
                std::cout << grid[1][i] << '\t';
                std::cout << grid[2][i] << '\t';
                std::cout << grid[3][i] << std::endl;
        }

}

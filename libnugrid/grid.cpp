#include "grid.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace libnugrid;

std::vector< std::vector< double >> libnugrid::genGrid(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, 
             gridTypes gridType, std::vector<double> gridParams) {

	if (gridType == RECT)  return genGridRect(atomName, atomPos, gridParams);
	throw "Unknown grid type";	
}


std::vector< std::vector< double >> libnugrid::genGridRect(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams) {


	//get Axis Aligned Bounding Box (AABB)
	double inf = std::numeric_limits<double>::infinity();
	double maxX = -inf, minX = inf, maxY = -inf, minY = inf, maxZ = -inf, minZ = inf;

	for (auto& pos: atomPos) {
		if (pos[X] > maxX) maxX = pos[X];
		if (pos[X] < minX) minX = pos[X];
		if (pos[Y] > maxY) maxY = pos[Y];
		if (pos[Y] < minY) minY = pos[Y];
		if (pos[Z] > maxZ) maxZ = pos[Z];
		if (pos[Z] < minZ) minZ = pos[Z];
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
		gridParams.push_back(gridParams[X]);	
		gridParams.push_back(gridParams[X]);	
	}
	gridDimX = std::round((maxX - minX) / gridParams[X]);	
	gridDimY = std::round((maxY - minY) / gridParams[Y]);	
	gridDimZ = std::round((maxZ - minZ) / gridParams[Z]);	

	int gridSize = gridDimX * gridDimY * gridDimZ;

	//generate grid
	std::vector<std::vector<double> > gridVector(4, std::vector<double> (gridSize, 0));

	for (int i = 0; i < gridDimX; i++) {
		for (int j = 0; j < gridDimY; j++) {
			for (int k = 0; k < gridDimZ; k++) {
				int currentIdx = k + j * gridDimY + i * (gridDimZ * gridDimY);
				gridVector[X][currentIdx] = minX + i * gridParams[X];
				gridVector[Y][currentIdx] = minY + j * gridParams[Y];
				gridVector[Z][currentIdx] = minZ + k * gridParams[Z];
				gridVector[W][currentIdx] = gridParams[X] * gridParams[Y] * gridParams[Z];
			}	
		}	
	}	
	return gridVector;
}


std::vector< std::vector< double >> libnugrid::genGridLebedev(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams) {

	int gridSize = 0;
	for (auto& param : gridParams) {
		gridSize += lebedevGetGridSize(std::round(param));	
	}
 
	//generate grid
	std::vector<std::vector<double> > gridVector(4, std::vector<double> (gridSize, 0));

	fstream gridFile;
	gridFile.open("lebedev_grids/lebedev_%3d.txt", std::fstream::in);

	int currentPoint = 0;
	for (auto& param : gridParams) {
		int size = lebedevGetGridSize(std::round(param));	
		for (int i = 0; i < size; i++) {
			double theta, phi, w;
			gridFile >> theta;
			gridFile >> phi;
			gridFile >> w;
			gridVector[X][currentPoint] = cos(theta) * sin(phi);
			gridVector[Y][currentPoint] = sin(theta) * sin(phi);
			gridVector[Z][currentPoint] = cos(phi);
			gridVector[W][currentPoint] = 4 * M_PI * w;
		}

	}

	
}

// Utils
// Grid printer
void libnugrid::printGrid(const std::vector<std::vector<double> >& grid) {
        for (int i = 0; i < grid[0].size(); i++) {
                std::cout << grid[X][i] << '\t';
                std::cout << grid[Y][i] << '\t';
                std::cout << grid[Z][i] << '\t';
                std::cout << grid[W][i] << std::endl;
        }

}

int libnugrid::lebedevGetGridSize(int gridPrecision) {
	
	vector<int> lebedevGridSizes = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810}

	vector<int> lebedevGridPrecisions = {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131}

	auto gridPrecision = std::find(lebedevGridPrecisions.begin(), lebedevGridPrecisions.end(), gridPrecision);
	if (gridPrecision == lebedevGridPrecisions.end()) {
		throw "Requested Grid Precision is not present";
	}
	int index = std::distance(lebedevGridPrecisions.begin(), gridPrecision);

	return lebedevGridSizes[index];
}

#include "nugrid.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include "lebedev_grids/lebedevData.h"

using namespace libnugrid;

std::vector< std::vector< double >> libnugrid::genGrid(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, 
             gridTypes gridType, std::vector<double> gridParams) {


	if (gridType == RECT)  return genGridRect(atomName, atomPos, gridParams);
	if (gridType == LEB_EM) {
		std::vector< std::vector< double > > gridToReturn = genGridLebedev(atomName, atomPos, gridParams);
		addBeckeWeightToGrid(atomName, atomPos, gridToReturn);
		return gridToReturn;

	}

	throw "Unknown grid type1";	
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
	for (int i = 0; i < gridParams.size(); i += 2) {
		gridSize += std::round(gridParams[i]) * std::round(gridParams[i + 1]);	
	}

 
	//generate grid
	std::vector<std::vector<double> > gridVector(5, std::vector<double> (gridSize, 0));


	int currentPoint = 0;
	for (int atomIdx = 0; atomIdx < atomName.size(); atomIdx++) {
		int sphericalNumPoints = gridParams[2*atomIdx];
		int radialNumPoints = gridParams[2*atomIdx + 1];
		std::vector< std::vector<double> > lebedevGrid = getLebedevGridByNumPoints(sphericalNumPoints);
		double alpha = getBraggSlaterRad(atomName[atomIdx]);
		for (int r = 1; r < radialNumPoints; r++) {
			double radialWeight = 1.0 / radialNumPoints;
			for (int i = 0; i < sphericalNumPoints; i++) {
				double r0 = EMradialFunc((double)r / radialNumPoints);
				double rprime = EMradialDervFunc((double)r / radialNumPoints);
				gridVector[X][currentPoint] = lebedevGrid[X][i] * r0 * alpha + atomPos[atomIdx][X];
				gridVector[Y][currentPoint] = lebedevGrid[Y][i] * r0 * alpha + atomPos[atomIdx][Y];
				gridVector[Z][currentPoint] = lebedevGrid[Z][i] * r0 * alpha + atomPos[atomIdx][Z];
				gridVector[W][currentPoint] = lebedevGrid[W][i] * r0 * r0 * alpha * alpha * alpha * rprime * radialWeight;
				gridVector[parentAtom][currentPoint] = atomIdx;
				currentPoint++;
			}
		}
	}
	return gridVector;
	
}

void libnugrid::addBeckeWeightToGrid(const std::vector<int>& atomName, 
             const std::vector<std::vector<double> >& atomPos, std::vector<std::vector<double> >& grid) {

	int natom = atomName.size();
	std::vector<std::vector<double> > inv_dist(natom, std::vector<double>(natom));
	std::vector<std::vector<double> > sizeRatioMatrix(natom, std::vector<double>(natom));
 	for (int i = 0; i < natom; i++) {
		Point thisAtomPos = {atomPos[i][X], atomPos[i][Y], atomPos[i][Z]};
		for (int j = 0; j < i; j++) {
			Point thatAtomPos = {atomPos[j][X], atomPos[j][Y], atomPos[j][Z]};
			inv_dist[i][j] = 1 / thatAtomPos.distanceTo(thisAtomPos);
			inv_dist[j][i] = inv_dist[i][j];
			double radiusRatio = getBraggSlaterRad(atomName[i]) / getBraggSlaterRad(atomName[j]);
			double chi = radiusRatio;
			sizeRatioMatrix[i][j] = getAfromChi(chi);
			sizeRatioMatrix[j][i] = -sizeRatioMatrix[i][j];
		}
	}

	
	
	for (int idx = 0; idx < grid[0].size(); idx++) {
		Point currentGridPoint = {grid[X][idx], grid[Y][idx], grid[Z][idx]};		
		std::vector<double> dist(natom);
	 	for (int i = 0; i < natom; i++) {
			Point thisAtomPos = {atomPos[i][X], atomPos[i][Y], atomPos[i][Z]};
			dist[i] = currentGridPoint.distanceTo(thisAtomPos);				
		}
		

		double numerator = 0;
		double denominator = 0;
 		for (int i = 0; i < natom; i++) {
			double prod = 1;
			for (int j = 0; j < natom; j++) {
				if (i == j)
					continue;
				double mu = (dist[i] - dist[j]) * inv_dist[i][j];
				double nu = mu + sizeRatioMatrix[i][j]*(1-mu*mu);
				double s = BeckeStepFunction(nu);
				prod *= s;
			}
		        if (i == grid[parentAtom][idx]) numerator = prod; 
		        denominator += prod;
		 }
	 	 grid[W][idx] *= numerator/denominator;
	}
}

// Utils

// Eular Maclauren transform function
double libnugrid::EMradialFunc(double q) {
	double temp = q / (1-q);
	return temp * temp;
}
double libnugrid::EMradialDervFunc(double q) {
	double temp = (1 - q);
	return 2 * q / (temp * temp * temp);
}


double libnugrid::getAfromChi(double chi) {
	double a = (1-chi*chi)/(4*chi);
	return (a < -0.5) ? -0.5 : (a > 0.5) ? 0.5 : a;
}

double libnugrid::BeckeStepFunction(double x)
{
    double   px =   x*(3-  x*x  )/2;
    double  ppx =  px*(3- px*px )/2;
    double pppx = ppx*(3-ppx*ppx)/2;
    return (1 - pppx)/2;
}

// Bragg Slater radius
double libnugrid::getBraggSlaterRad(int Z) {
	const static std::vector<double> radii = {1.000, // ghost atom, and moves the table to 1-index
        1.000,
        1.001,                                                                                                                 1.012,
        0.825, 1.408,                                                                       1.485, 1.452, 1.397, 1.342, 1.287, 1.243,
        1.144, 1.364,                                                                       1.639, 1.716, 1.705, 1.683, 1.639, 1.595,
        1.485, 1.474, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.650, 1.727, 1.760, 1.771, 1.749, 1.727,
        1.628, 1.606, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.672, 1.804, 1.881, 1.892, 1.892, 1.881};
	if (Z < radii.size())	return radii[Z];
	else return radii[radii.size() - 1];
}


// Grid printer
void libnugrid::printGrid(const std::vector<std::vector<double> >& grid) {
        for (int i = 0; i < grid[0].size(); i++) {
                std::cout << grid[X][i] << '\t';
                std::cout << grid[Y][i] << '\t';
                std::cout << grid[Z][i] << '\t';
                std::cout << grid[W][i] << std::endl;
        }

}
/*
int libnugrid::lebedevGetGridSize(int reqGridPrecision) {
	
	static const std::vector<int> lebedevGridSizes = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};

	static const std::vector<int> lebedevGridPrecisions = {3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};

	auto gridPrecision = std::find(lebedevGridPrecisions.begin(), lebedevGridPrecisions.end(), reqGridPrecision);
	if (gridPrecision == lebedevGridPrecisions.end()) {
		throw "Requested Grid Precision is not present";
	}
	int index = std::distance(lebedevGridPrecisions.begin(), gridPrecision);

	return lebedevGridSizes[index];
}*/

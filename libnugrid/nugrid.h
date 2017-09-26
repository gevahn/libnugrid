#pragma once

#include <vector>
#include <string>
#include <cmath>


namespace libnugrid {

enum {X, Y, Z, W, parentAtom};

enum gridTypes {RECT, LEBEDEV, LEB_RECT, LEB_EM};

// Master function to call all other types of grid generators
std::vector< std::vector< double > > genGrid(const std::vector<int>& atomName,
					    const std::vector<std::vector<double> >& atomPos,
					    gridTypes gridType, std::vector<double> gridParams);

// Rectalinear equal spacing grid
std::vector< std::vector< double > > genGridRect(const std::vector<int>& atomName,
					    const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams);

// Lebedev on a unit sphere
std::vector< std::vector< double > > genGridLebedev(const std::vector<int>& atomName,
					    const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams);

void addBeckeWeightToGrid(const std::vector<int>& atomName,
             const std::vector<std::vector<double> >& atomPos, std::vector<std::vector<double> >& grid);

// Utils
double EMradialFunc(double q); 
double EMradialDervFunc(double q); 
double getAfromChi(double chi);  
double BeckeStepFunction(double x);
double getBraggSlaterRad(int Z);
struct Point {
	double x;
	double y;
	double z;
	inline double distanceTo(const Point& otherPoint) { return sqrt((x-otherPoint.x) * (x-otherPoint.x)
							       + (y-otherPoint.y) * (y-otherPoint.y)
							       + (z-otherPoint.z) * (z-otherPoint.z));};
};

// Print a grid to cout
void printGrid(const std::vector<std::vector<double> >& grid);

// Get the number of points in a levedev grid of the given precision
int lebedevGetGridSize(int gridPrecision);
}

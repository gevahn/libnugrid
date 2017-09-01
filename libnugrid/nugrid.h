#include <vector>
#include <string>

#pragma once

namespace libnugrid {

enum {X, Y, Z, W};

enum gridTypes {RECT, LEBEDEV, LEB_RECT, LEB_GAUSS};

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

// Utils
// Print a grid to cout
void printGrid(const std::vector<std::vector<double> >& grid);

// Get the number of points in a levedev grid of the given precision
int lebedevGetGridSize(int gridPrecision);
}
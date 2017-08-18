#include <vector>
#include <string>

#pragma once

namespace libnugrid {

//Master function to call all other types of grid generators
std::vector< std::vector< double > > genGrid(const std::vector<int>& atomName,
					    const std::vector<std::vector<double> >& atomPos,
					    std::string gridType, std::vector<double> gridParams);

//Rectalinear equal spacing grid
std::vector< std::vector< double > > genGridRect(const std::vector<int>& atomName,
					    const std::vector<std::vector<double> >& atomPos, std::vector<double> gridParams);



//Print a grid to cout
void printGrid(const std::vector<std::vector<double> >& grid);

}

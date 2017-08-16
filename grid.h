#include <vector>
#include <string>

#pragma once

namespace libnugrid {

//Master function to call all other types of grid generators
std::vector< std::vector< double >> genGrid(const std::vector<int>& atomName,
					    const std::vector<std::vector<double>>& atomPos,
					    std::sting gridType, std::vector<double> gridParams);

//Rectalinear equal spacing grid
std::vector< std::vector< double >> genGridRect(const std::vector<int>& atomName,
					    const std::vector<std::vector<double>>& atomPos, std::vector<double> gridParams);

}

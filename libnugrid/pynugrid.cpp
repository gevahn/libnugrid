#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "nugrid.h"

namespace py = pybind11;


// how to grid type...?
py::array_t<double> gen_grid(py::array_t<double> atomPos, py::array_t<int> atomName, py::array_t<double> gridParams) {
    py::buffer_info param_info = gridParams.request();
    py::buffer_info pos_info = atomPos.request();
    py::buffer_info name_info = atomNames.request();
    
    if (param_info.ndim != 1)
        throw std::runtime_error("parameters is not a 1D vector)");
    if (pos_info.ndim != 2)
        throw std::runtime_error("atom positions is not a matrix");
    if (name_info.shape[0] != pos_info.shape[0])
        throw std::runtime_error("atom positions length does not match atom names length");

    std::vector<double> passParams(param_info.shape[0]);
    std::vector<std::vector<double>> passPos(pos_info.shape[0], std::vector<double>(pos_info.shape[1]);
    std::vector<int> passNames(name_info.shape[0]);

    passParams.assign(param_info.data(), param_info.data() + param_info.shape[0]);
    passNames.assign(name_info.data(), name_info.data() + name_info.shape[0]);

    const double * pos_data = pos_info.data();

    for (size_t i = 0; i < pos_info.size(); i+=3) {
        passPos[i][0] = pos_data[i + 0];
        passPos[i][1] = pos_data[i + 1];
        passPos[i][2] = pos_data[i + 2];
    }

    // generate a grid from libnugrid
    std::vector<std::vector<double>> gridArray = libnugrid::genGrid(passNames, passPos, libnugrid::LEB_EM, passParams);


    //build pybind numpy data object
    py::buffer_info grid_info = {
        gridArray.data(),
        sizeof(double),
        py::format_descriptor<double>::format(),
        2,
        {gridArray[0].size(), 5},
        {&gridArray[1] - &gridArray[0], sizeof(double)}};

   return py::array(grid_info);
}



PYBIND11_PLUGIN(nugrid) {
    py::module m("pynugrid", "The non-uniform grid module");

    m.def("genGrid", &gen_grid, "Generats an atom centered grid");

    return m.ptr();
}

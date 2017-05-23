#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "result.hpp"

result::result
(
    latticeBase &lb,
    fluidField &field
)
: lb_ (lb),
  field_ (field)
{
    auto results = result::initializeCleanFolder();
    if (results != 0) throw std::runtime_error("Error in folder initialization");
}

int result::initializeCleanFolder()
{
    // creates output folders if they don't exist
    auto vtk_folder = system("mkdir -p vtk_fluid");
    auto old_vtk_fluid_files = system("rm -f vtk_fluid/*");
    return vtk_folder + old_vtk_fluid_files;
}

void result::writeResultVTK(int time)
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    std::ofstream vtk_file;
    vtk_file.open("vtk_fluid/fluid_t" + std::to_string(time) + ".vtk");
    vtk_file << "# vtk DataFile Version 3.0" << std::endl;
    vtk_file << "fluid_state" << std::endl;
    vtk_file << "ASCII" << std::endl;
    vtk_file << "DATASET RECTILINEAR_GRID" << std::endl;
    vtk_file << "DIMENSIONS " << nx << " " << ny << " 1" << std::endl;

    // Write x, y, z coordinates. z set to be 1 since it's 2D
    vtk_file << "X_COORDINATES " << nx << " float" << std::endl;
    for (auto x = 0u; x < nx; ++x) vtk_file << x << " ";
    vtk_file << std::endl;
    vtk_file << "Y_COORDINATES " << ny << " float" << std::endl;
    for (auto y = 0u; y < ny; ++y) vtk_file << y << " ";
    vtk_file << std::endl;
    vtk_file << "Z_COORDINATES " << 1 << " float" << std::endl;
    vtk_file << 0 << std::endl;
    vtk_file << "POINT_DATA " << nx * ny << std::endl;

    // Write relative pressure
    vtk_file << "SCALARS relative_pressure float" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (auto pressure : field_.p) vtk_file << pressure << std::endl;

    // Write velocity as vectors
    vtk_file << "VECTORS velocity_vector float" << std::endl;
    for (auto v : field_.u) vtk_file << v[0] << " " << v[1] << " 0" << std::endl;
    vtk_file.close();
}

#include <iostream>
#include <stdexcept>
#include <vector>

#include "algorithm.hh"

#include "latticeBase.hpp"
#include "collisionBase.hpp"

collisionBase::collisionBase
(
    latticeBase &lb,
    double initial_density
)
: lb_ (lb),
  tau_ {0},
  c_ {lb.getLatticeSpeed()}
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto nz = lb_.getNumberOfNz();
    const auto nc = lb_.getNumberOfDirections();
    const auto lat_size = nx * ny * nz;
    eqdf.assign(lat_size, std::vector<double>(nc, 0.0));
    rho.assign(lat_size, initial_density);
    skip.assign(lat_size, false);
    computeEq();
}

collisionBase::collisionBase
(
    latticeBase &lb,
    const std::vector<double> &initial_density
)
: eqdf {},
  rho {initial_density},
  skip {},
  lb_ (lb),
  tau_ {0},
  c_ {lb.getLatticeSpeed()}
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto nz = lb_.getNumberOfNz();
    const auto nc = lb_.getNumberOfDirections();
    const auto lat_size = nx * ny * nz;
    eqdf.assign(lat_size, std::vector<double>(nc, 0.0));
    skip.assign(lat_size, false);
    computeEq();
}

std::vector<double> collisionModel::computeRho
(
    const std::vector<std::vector<double>> &df
)
{
    std::vector<double> rho(df.size(), 0.0);
    auto it_rho = begin(rho);
    for (auto i_df : df) (*it_rho++) = getZerothMoment(i_df);
    return rho;
}



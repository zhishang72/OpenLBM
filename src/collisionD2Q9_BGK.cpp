#include <vector>

#include "momentComputing.h"

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "collisionD2Q9_BGK.hpp"

collisionD2Q9_BGK::collisionD2Q9_BGK
(
    latticeBase &lb,
    double kinematic_viscosity,
    double initial_density,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: collisionBase(lb, initial_density),
  field_ (field),
  D2Q9_ (D2Q9),
  tau_ {0},
  skip {}
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto lat_size = nx * ny;
    skip.assign(lat_size, false);
    const auto dt = lb_.getTimeStep();
    // BGK tau_ formula from "Discrete lattice effects on the forcing term in
    // the lattice Boltzmann method" Guo2002
    tau_ = 0.5 + kinematic_viscosity / (cs_sqr_ * dt);  //BGK
}

collisionD2Q9_BGK::collisionD2Q9_BGK
(
    latticeBase &lb,
    double kinematic_viscosity,
    const std::vector<double> &initial_density,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: collisionBase(lb, initial_density),
  tau_ {0},
  field_ (field),
  D2Q9_ (D2Q9),
  skip {}
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto lat_size = nx * ny;
    skip.assign(lat_size, false);
    const auto dt = lb_.getTimeStep();
    // BGK tau_ formula from "Discrete lattice effects on the forcing term in
    // the lattice Boltzmann method" Guo2002
    tau_ = 0.5 + kinematic_viscosity / (cs_sqr_ * dt);  //BGK
}

void collisionD2Q9_BGK::computeEq()
{
    auto nx = lb_.getNumberOfNx();
    auto ny = lb_.getNumberOfNy();
    auto nc = lb_.getNumberOfDirections();
    for (auto n = 0u; n < nx * ny; ++n)
    {
        double u_sqr = innerProduct(field_.u[n], field_.u[n]);
        u_sqr /= 2.0 * cs_sqr_;
        for (auto i = 0u; i < nc; ++i)
        {
            double c_dot_u = innerProduct(D2Q9_.e.at(i), field_.u[n]);
            c_dot_u /= cs_sqr_;
            eqdf[n][i] = D2Q9_.weight[i] * rho_[n] *
                         (1.0 + c_dot_u * (1.0 + c_dot_u / 2.0) - u_sqr);
        }  // i
    }  // n
}

std::vector<double> collisionD2Q9_BGK::computeRho
(
    const std::vector<std::vector<double>> &df
)
{
    std::vector<double> rho_update(df.size(), 0.0);
    auto it_rho = begin(rho_update);
    for (auto it_df : df) (*it_rho++) = getZerothMoment(it_df);
    return rho_update;
}

std::vector<std::vector<double>> collisionD2Q9_BGK::computeU
(
    const std::vector<std::vector<double>> &df
)
{
    std::vector<std::vector<double>> rhou;
    auto index = 0u;
    for (auto it_df : df) rhou.push_back(getFirstMoment(it_df, D2Q9_.e)); //now is rho*u
    for (auto &it_rhou : rhou)
    {
        for (auto &it_u : it_rhou) it_u /= rho_[index]; //now is rho*u/rho
        ++index;
    }  // node
    return rhou; //now is only u
}

void collisionD2Q9_BGK::computeMacroscopicProperties
(
    const std::vector<std::vector<double>> &df
)
{
    rho_ = computeRho(df);
    auto rho = rho_;
    for (auto &it_rho : rho) it_rho *= cs_sqr_;  //now rho is pressure
    field_.p = rho;
    field_.u = computeU(df);
}

void collisionD2Q9_BGK::addNodeToSkip(std::size_t n)
{
  skip[n] = true;
}

void collisionD2Q9_BGK::collide
(
    std::vector<std::vector<double>> &df_lattice
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto nc = lb_.getNumberOfDirections();
    for (auto n = 0u; n < nx * ny; ++n)
    {
        if (!skip[n])
        {
            for (auto i = 0u; i < nc; ++i)
            {
                df_lattice[n][i] += (eqdf[n][i] - df_lattice[n][i]) / tau_;
            }  // i
        }
    }  // n
}

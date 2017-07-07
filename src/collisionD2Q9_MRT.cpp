#include <vector>

#include "momentComputing.h"

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "collisionD2Q9_MRT.hpp"

collisionD2Q9_MRT::collisionD2Q9_MRT
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
  s_ {},
  m_ {},
  mEq_ {},
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
    auto nc = lb_.getNumberOfDirections();
    s_.assign(nc, 0.0);
    s_[7] = 1.0 / tau_;
    s_[1] = 1.6;
    s_[2] = 1.8;
    s_[4] = 8.0*(2.0-s_[7])/(8.0-s_[7]);
    s_[6] = s_[4];
    s_[8] = s_[7];

}

collisionD2Q9_MRT::collisionD2Q9_MRT
(
    latticeBase &lb,
    double kinematic_viscosity,
    const std::vector<double> &initial_density,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: collisionBase(lb, initial_density),
  tau_ {0},
  s_{},
  m_ {},
  mEq_ {},
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
    auto nc = lb_.getNumberOfDirections();
    s_.assign(nc, 0.0);
    s_[7] = 1.0 / tau_;
    s_[1] = 1.6;
    s_[2] = 1.8;
    s_[4] = 8.0*(2.0-s_[7])/(8.0-s_[7]);
    s_[6] = s_[4];
    s_[8] = s_[7];
}

void collisionD2Q9_MRT::computefEq()
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
            eqdf[n][i] = D2Q9_.weight[i] * rho_[n] * (1.0 + c_dot_u * (1.0 + c_dot_u / 2.0) - u_sqr);
        }  // i
    }  // n
}

void collisionD2Q9_MRT::computemEq()
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    auto nc = lb_.getNumberOfDirections();
    const auto lat_size = nx * ny;
    mEq_.assign(lat_size, std::vector<double>(nc, 0.0));
    for (auto n = 0u; n < nx * ny; ++n)
    {
        double jx = rho_[n] * field_.u[n][0];
        double jy = rho_[n] * field_.u[n][1];

        mEq_[n][0] = rho_[n];
        mEq_[n][1] = -2.0 * rho_[n] + 3.0 * (jx * jx + jy * jy);
        mEq_[n][2] = rho_[n] - 3.0 * (jx * jx + jy * jy);
        mEq_[n][3] = jx;
        mEq_[n][4] = -jx;
        mEq_[n][5] = jy;
        mEq_[n][6] = -jy;
        mEq_[n][7] = (jx * jx - jy * jy);
        mEq_[n][8] = jx * jy;
    }  // n
}

void collisionD2Q9_MRT::computeM
(
    const std::vector<std::vector<double>> &df
)
{
    auto nx = lb_.getNumberOfNx();
    auto ny = lb_.getNumberOfNy();
    auto nc = lb_.getNumberOfDirections();
    const auto lat_size = nx * ny;
    m_.assign(lat_size, std::vector<double>(nc, 0.0));
    computemEq();
    for (auto n = 0u; n < nx * ny; ++n)
    {
        m_.at(n) = MVmultiplyProduct(D2Q9_.M, df.at(n));
        for (auto i = 0u; i < nc; ++i)
        {
            m_[n][i] += s_[i] * (mEq_[n][i] - m_[n][i]);
        }
    }  // n
}

std::vector<double> collisionD2Q9_MRT::computeRho
(
    const std::vector<std::vector<double>> &df
)
{
    std::vector<double> rho_update(df.size(), 0.0);
    auto it_rho = begin(rho_update);
    for (auto it_df : df) (*it_rho++) = getZerothMoment(it_df);
    return rho_update;
}

std::vector<std::vector<double>> collisionD2Q9_MRT::computeU
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

void collisionD2Q9_MRT::computeMacroscopicProperties
(
    const std::vector<std::vector<double>> &df
)
{
    rho_ = computeRho(df);
    auto rho = rho_;
    for (auto &it_rho : rho) it_rho = cs_sqr_ * (it_rho - 1.0);  //now rho is pressure
    field_.p = rho;
    field_.u = computeU(df);
}

void collisionD2Q9_MRT::addNodeToSkip(std::size_t n)
{
  skip[n] = true;
}

void collisionD2Q9_MRT::collide
(
    std::vector<std::vector<double>> &df_lattice
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    auto nc = lb_.getNumberOfDirections();

    computeM(df_lattice);

    for (auto n = 0u; n < nx * ny; ++n)
    {
        if (!skip[n])
        {
            df_lattice.at(n) = MVmultiplyProduct(D2Q9_.Minv, m_.at(n));
        }
    }  // n
}

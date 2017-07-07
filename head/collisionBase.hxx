#ifndef COLLISIONBASE_HXX_INCLUDED
#define COLLISIONBASE_HXX_INCLUDED

#include <vector>

#include "latticeBase.hpp"

class collisionBase
{
    public:
        // Constructor: Creates collision base with the same density at each node
        // param lm lattice model used for simulation
        // param initial_density initial density of the lattice
        collisionBase
        (
            latticeBase &lb,
            double initial_density
        )
        : eqdf {},
          lb_ (lb),
          rho_ {},
          c_ {lb.getLatticeSpeed()}
        {
            const auto nx = lb_.getNumberOfNx();
            const auto ny = lb_.getNumberOfNy();
            const auto nc = lb_.getNumberOfDirections();
            const auto lat_size = nx * ny;
            eqdf.assign(lat_size, std::vector<double>(nc, 0.0));
            rho_.assign(lat_size, initial_density);
        };
        // Constructor: Creates collision base with the same density at each node
        // param lm lattice model used for simulation
        // param initial_density initial density of the lattice
        collisionBase
        (
            latticeBase &lb,
            const std::vector<double> &initial_density
        )
        : eqdf {},
          lb_ (lb),
          rho_ {initial_density},
          c_ {lb.getLatticeSpeed()}
        {
            const auto nx = lb_.getNumberOfNx();
            const auto ny = lb_.getNumberOfNy();
            const auto nc = lb_.getNumberOfDirections();
            const auto lat_size = nx * ny;
            eqdf.assign(lat_size, std::vector<double>(nc, 0.0));
        };
        // https://stackoverflow.com/questions/353817/should-every-class-have-a-
        // virtual-destructor
        // Virtual destructor since we are deriving from this class
        virtual ~collisionBase()= default;
        // Calculates equilibrium distribution function according to LBIntro
        virtual void computefEq() = 0;
        // Compute density at each node by summing up its distribution functions
        // param lattice 2D vector containing distribution functions
        // return density of lattice stored row-wise in a 1D vector
        virtual std::vector<double> computeRho
        (
            const std::vector<std::vector<double>> &df
        ) = 0;
        // Pure virtual function to compute the macroscopic properties of the lattice
        // depending on the equation, density and velocity for Navier-Stokes, only
        // density for Convection-diffusion equation
        // This is used to unify function calling in LatticeBoltzmann takeStep()
        // method
        // param df Particle distribution functions of the lattice stored row-wise
        // in a 2D vector
        virtual void computeMacroscopicProperties
        (
            const std::vector<std::vector<double>> &df
        ) = 0;
        // Adds a node to exclude it from the collision step
        // param n index of the node in the lattice
        virtual void addNodeToSkip
        (
            std::size_t n
        ) = 0;
        // Pure virtual function to compute collision step and apply force step
        // according to "A new scheme for source term in LBGK model for
        // convection  diffusion equation" and Guo2002
        // param lattice 2D vector containing distribution functions
        virtual void collide
        (
            std::vector<std::vector<double>> &df
        ) = 0;
        // Density stored row-wise in a 1D vector
        std::vector<double> rho_;
        // Equilibrium distribution function stored row-wise in a 2D vector
        std::vector<std::vector<double>> eqdf;
    protected:
        // Lattice model to handle number of rows, columns, dimensions, directions,
        // velocity/
        latticeBase &lb_;
        // Speed of sound in lattice/
        double c_;
        // Square of speed of sound in lattice, used to simplify computations in the
        // collision step
        double cs_sqr_ = c_ * c_ / 3.0;
};

#endif // COLLISIONBASE_HXX_INCLUDED

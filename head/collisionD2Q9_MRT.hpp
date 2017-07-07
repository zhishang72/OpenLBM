#ifndef COLLISIOND2Q9_MRT_HPP_INCLUDED
#define COLLISIOND2Q9_MRT_HPP_INCLUDED

#include "latticeModel.hxx"
#include "collisionBase.hxx"

#include "latticeBase.hpp"

class collisionD2Q9_MRT: public collisionBase
{
    public:
        // Constructor: Creates collision model for NS equation with the same density
        // at each node
        // param lm lattice model used for simulation
        // param kinematic viscosity
        // param initial_density_f initial density of NS lattice
        collisionD2Q9_MRT
        (
            latticeBase &lb,
            double kinematic_viscosity,
            double initial_density_f,
            latticeModelD2Q9 &D2Q9,
            fluidField &field
        );
        // Constructor: Creates collision model for NS equation with the same density
        // at each node
        // param lm lattice model used for simulation
        // param kinematic viscosity
        // param initial_density_f initial density of NS lattice
        collisionD2Q9_MRT
        (
            latticeBase &lb,
            double kinematic_viscosity,
            const std::vector<double> &initial_density_f,
            latticeModelD2Q9 &D2Q9,
            fluidField &field
        );
        // Virtual destructor since we may be deriving from this class
        virtual ~collisionD2Q9_MRT() = default;
        // Calculates equilibrium distribution function according to LBIntro
        void computefEq();
        // Calculates momentem equilibrium distribution function according to LBIntro
        void computemEq();
        // Calculates momentem equilibrium distribution function according to LBIntro
        void computeM
        (
            const std::vector<std::vector<double>> &df
        );
        // Compute density at each node by summing up its distribution functions
        // param lattice 2D vector containing distribution functions
        // return density of lattice stored row-wise in a 1D vector
        std::vector<double> computeRho
        (
            const std::vector<std::vector<double>> &df
        );
        // Calculated velocity for NS equation without body force based on formula in Guo2002
        // param df 2D vector containing distribution functions of the NS equation
        // return 2D vector containing velocity at each node of lattice
        std::vector<std::vector<double>> computeU
        (
            const std::vector<std::vector<double>> &df
        );
        // Computes the macroscopic properties based on the collision model used, both
        // velocity and density in this case. Based on "Discrete lattice effects on
        // the forcing term in the lattice Boltzmann method"
        // param df lattice distribution functions stored row-wise in a 2D vector
        void computeMacroscopicProperties
        (
            const std::vector<std::vector<double>> &df
        );
        // Adds a node to exclude it from the collision step
        // param n index of the node in the lattice
        void addNodeToSkip
        (
            std::size_t n
        );
        // Collides according to Guo2002
        // param lattice 2D vector containing distribution functions
        void collide
        (
            std::vector<std::vector<double>> &df_lattice
        );
    private:
        // define fluid field;
        fluidField &field_;
        // define lattice model
        latticeModelD2Q9 &D2Q9_;
        // Relaxation time for BGK/
        double tau_;
        // Relaxation rate for MRT/
        std::vector<double> s_;
        // Relaxation rate for MRT/
        std::vector<std::vector<double>> m_;
        // Relaxation rate for MRT/
        std::vector<std::vector<double>> mEq_;
        // Skips the collision step for the node if it is a full-way bounceback node
        std::vector<bool> skip;
};

#endif // COLLISIOND2Q9_MRT_HPP_INCLUDED

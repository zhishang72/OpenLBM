#ifndef LATTICEBOLTZMANN_HPP_INCLUDED
#define LATTICEBOLTZMANN_HPP_INCLUDED

#include <vector>

#include "latticeBase.hpp"
#include "collisionBase.hxx"
#include "streamBase.hxx"
#include "boundaryNode.hxx"

class latticeBoltzmann
{
    public:
        // Constructor: Creates a LatticeBoltzmann object
        // param lm lattice model which contains information on the number of rows,
        //       columns, dimensions, discrete directions and lattice velocity
        // param cm collision model used by the lattice: Convection-diffusion,
        //       Navier-Stokes and Navier-Stokes with force
        // param sm stream mode used by the lattice: Periodic stream, non-periodic streaming
        latticeBoltzmann
        (
            latticeBase &lb,
            collisionBase &cb,
            streamBase &sb
        );
        latticeBoltzmann(const latticeBoltzmann&) = default;
        ~latticeBoltzmann() = default;
        // Adds a boundary condition to the lattice
        // param bn pointer to the boundary condition to be added
        void addBoundaryNode
        (
            boundaryNode *bn
        );
        // Performs one cycle of evolution equation, computes the relevant macroscopic
        // properties such as velocity and density
        void takeStep();
        // by reference, similar to by pointer
        // https://stackoverflow.com/questions/9285627/is-it-possible-to-pass-derived-
        // classes-by-reference-to-a-function-taking-base-cl
    private:
        // Lattice distribution function stored row-wise in a 2D vector
        std::vector<std::vector<double>> df;
        // Lattice model which contains information on the number of rows, columns,
        // dimensions, discrete directions and lattice velocity
        latticeBase &lb_;
        // Collision models to perform the collision step based on the model chosen
        collisionBase &cb_;
        // Stream model to perform the streaming step based on the model chosen
        streamBase &sb_;
        // Pointers to boundary conditions in the lattice stored in a vector.
        // References cannot be used as it is not possible to store a vector of references
        std::vector<boundaryNode*> bn_;
};

#endif // LATTICEBOLTZMANN_HPP_INCLUDED

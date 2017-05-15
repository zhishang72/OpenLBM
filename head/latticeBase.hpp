#ifndef LATTICEBASE_HPP_INCLUDED
#define LATTICEBASE_HPP_INCLUDED

#include <iostream>
//#include <vector>

class latticeBase
{
    public:
        // Constructor: creates lattice model with the same velocity at each node
        // param num_dims number of dimensions
        // param num_dirs number of discrete directions
        // param dx space step
        // param dt time step
        latticeBase
        (
            std::size_t nx,
            std::size_t ny,
            std::size_t num_dims,
            std::size_t num_dirs,
            double dl,
            double dt
        );
        // Constructor: creates lattice model with the same velocity at each node
        // param num_dims number of dimensions
        // param num_dirs number of discrete directions
        // param dx space step
        // param dt time step
        latticeBase
        (
            std::size_t nx,
            std::size_t ny,
            std::size_t nz,
            std::size_t num_dims,
            std::size_t num_dirs,
            double dl,
            double dt
        );
        // Virtual destructor since the deriving from this class, see collision.hpp
        virtual ~latticeBase() = default;
        // Get the number of grid along x coordinate
        // return number of grid along x coordinate
        std::size_t getNumberOfNx() const;
        // Get the number of grid along y coordinate
        // return number of grid along y coordinate
        std::size_t getNumberOfNy() const;
        // Get the number of grid along z coordinate
        // return number of grid along z coordinate
        std::size_t getNumberOfNz() const;
        // Get the number of dimensions of the lattice. 2 for 2D and 3 for 3D.
        // return number of dimensions of the lattice
        std::size_t getNumberOfDimensions() const;
        // Get the number of discrete velocities of the lattice, specified by the
        // model used. 9 for Q9.
        // return number of discrete velocities of the lattice
        std::size_t getNumberOfDirections() const;
        // Get the space step (dl) of the model
        // return space step of the model
        double getSpaceStep() const;
        // Get the time step (dt) of the model
        // return time step of the model
        double getTimeStep() const;
        // Get the lattice speed (c) of the model
        // return lattice speed of the model
        double getLatticeSpeed() const;
        // Checks if input parameters for lattice base is valid, prevents creation of
        // invalid lattice base, such as a size 0 x 0 lattice
        // return validity of lattice base
        //        TRUE: input parameters are valid
        //        FALSE: input parameters are invalid
        bool checkInput();
    private:
        // Number of grid along x coordinate
        std::size_t number_of_nx_;
        // Number of grid along y coordinate
        std::size_t number_of_ny_;
        // Number of grid along z coordinate
        std::size_t number_of_nz_;
        // Number of dimensions of the lattice model, 2 for D2Q9 model
        std::size_t number_of_dimensions_;
        // Number of discrete directions of the lattice model, 9 for D2Q9 model
        std::size_t number_of_directions_;
        // Space step of the lattice model, dl
        double space_step_;
        // Time step of the lattice model, dt
        double time_step_;
        // Propagation speed on the lattice. Based on "Introduction to Lattice Boltzmann Methods"
        double c_ = space_step_ / time_step_;
};

#endif // LATTICEBASE_HPP_INCLUDED

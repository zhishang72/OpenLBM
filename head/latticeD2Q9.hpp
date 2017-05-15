#ifndef LATTICED2Q9_HPP_INCLUDED
#define LATTICED2Q9_HPP_INCLUDED

#include <vector>

#include "latticeModel.hxx"

#include "latticeBase.hpp"

class latticeD2Q9: public latticeBase
{
    public:
        // Constructor: Create lattice model for D2Q9 with the same velocity at each node
        // param num_rows number of nx
        // param num_cols number of ny
        // param dl space step
        // param dt time step
        // param initial model of the lattice
        latticeD2Q9
        (
            std::size_t num_nx,
            std::size_t num_ny,
            double dl,
            double dt,
            latticeModelD2Q9 &D2Q9
        );
        // Destructor
        virtual ~latticeD2Q9() = default;
    private:
        // define lattice model
        latticeModelD2Q9 &D2Q9_;
};

#endif // LATTICED2Q9_HPP_INCLUDED

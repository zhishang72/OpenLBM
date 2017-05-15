#ifndef STREAMD2Q9_HPP_INCLUDED
#define STREAMD2Q9_HPP_INCLUDED

#include <vector>

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "streamBase.hxx"

class streamD2Q9: public streamBase
{
    public:
        // Constructor: Creates a non-periodic streaming model for D2Q9 lattice model
        // param lm Lattice model which contains information on the number of rows,
        // columns, dimensions, discrete directions and lattice velocity
        streamD2Q9
        (
            latticeBase &lb,
            latticeModelD2Q9 &D2Q9
        );
        // Destructor
        ~streamD2Q9() = default;
        // Performs the streaming function based on "Introduction to Lattice Boltzmann
        // Methods". Distribution functions which require off-lattice streaming are
        // unchanged
        // param df lattice distribution functions stored row-wise in a 2D vector
        std::vector<std::vector<double>> stream
        (
            const std::vector<std::vector<double>> &df
        );
    private:
        // define lattice model
        latticeModelD2Q9 &D2Q9_;
};

#endif // STREAMD2Q9_HPP_INCLUDED

#ifndef STREAMBASE_HXX_INCLUDED
#define STREAMBASE_HXX_INCLUDED

#include <vector>

#include "latticeBase.hpp"

class streamBase
{
    public:
        // Constructor: Base class for stream models
        // param lm lattice model which contains information on the number of rows,
        // columns, dimensions, discrete directions and lattice velocity
        streamBase
        (
            latticeBase &lb
        )
        : lb_ (lb)
        {};
        //Virtual destruction since we are deriving from this class
        virtual ~streamBase() = default;
        // Pure virtual function for the streaming function
        // param df lattice distribution function stored row-wise in a 2D vector
        // df[n][i]: n is the index of grid; i is index of lattice velocity at grid
        virtual std::vector<std::vector<double>> stream
        (
            const std::vector<std::vector<double>> &df
        ) = 0;
    protected:
        // Lattice model to handle number of rows, columns, dimensions, directions,
        // velocity/
        latticeBase &lb_;
};

#endif // STREAMBASE_HXX_INCLUDED

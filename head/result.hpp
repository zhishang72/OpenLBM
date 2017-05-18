#ifndef RESULT_HPP_INCLUDED
#define RESULT_HPP_INCLUDED

#include <string>
#include <vector>

#include "latticeBase.hpp"

class result
{
    public:
        // Constructor: Creates results class with reference to LatticeModel for
        // information on number of rows, columns, space step, time step and lattice
        // velocity. Creates and cleans the folders for output results as well.
        // Throws exception if folder initialization fails
        // param lm reference to LatticeModel
        result
        (
            latticeBase &lb,
            fluidField &field

        );
        result(const result&) = default;
        result& operator= (const result&) = default;
        // Destructor
        ~result() = default;
        // Creates output folders if they don't already exist and make sure old
        // results are deleted
        // return sum of status codes return by the called commands, 0 if successful
        int initializeCleanFolder();
        // Writes results at a particular time point to .vtk files for post-processing
        // with ParaView. Currently writes: coordinates, density difference
        // velocity in x- and y- direction, for NS only. Will throw exception if NS
        // collision model is not registered
        // param time time point
        void writeResultVTK
        (
            int time
        );
    private:
        // Reference to LatticeModel
        latticeBase &lb_;
        // define fluid field;
        fluidField &field_;
};

#endif // RESULT_HPP_INCLUDED

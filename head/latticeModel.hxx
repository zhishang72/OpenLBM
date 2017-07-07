#ifndef LATTICEMODEL_HXX_INCLUDED
#define LATTICEMODEL_HXX_INCLUDED

#include <iostream>
#include <vector>

struct fluidField
{
    // pressure 1D with n grids * pressure scalar in the field
    std::vector<double> p;
    // velocity 2D with n grids * velocity vector in the field
    std::vector<std::vector<double>> u;
    // Constructor: Create lattice model for D2Q9 with variable velocity at each node
    // param num_rows number of nx
    // param num_cols number of ny
    // param initial_velocity initial velocity of the lattice
    fluidField
    (
        std::size_t num_nx,
        std::size_t num_ny,
        const std::vector<double> &initial_velocity
    )
    : p {}
    {
       u.assign(num_nx * num_ny, initial_velocity);
    };
    // Constructor: Create lattice model for D2Q9 with variable velocity at each node
    // param num_rows number of nx
    // param num_cols number of ny
    // param initial_velocity initial velocity of the lattice
    fluidField
    (
        const std::vector<std::vector<double>> &initial_velocity
    )
    : p {}
    {
        u = initial_velocity;
    };
    // Destructor
    virtual ~fluidField() = default;
 };

struct latticeModelD2Q9
{
    // 6  2  5  ^ y
    //  \ | /   |
    // 3--0--1  |
    //  / | \   |       x
    // 7  4  8  +------->
    // velocities D2Q9 due to the velocities will be updated by lattice velocity
    std::vector<std::vector<double>> e =
    {
        {0.0, 0.0},                             //at 0
        {1.0, 0.0}, { 0.0, 1.0}, {-1.0,  0.0}, {0.0, -1.0},   //at 1,2,3,4
        {1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}  //at 5,6,7,8
    };
    // Enumeration for discrete directions D2Q9 to be used with distribution functions
    enum e_directions
    {
        E = 1,
        N,
        W,
        S,
        NE,
        NW,
        SW,
        SE
    };
    // weight for D2Q9
    const std::vector<double> weight =
    {
        16.0 / 36.0,                                      //at 0
        4.0 / 36.0,  4.0 / 36.0, 4.0 / 36.0, 4.0 / 36.0,  //at 1,2,3,4
        1.0 / 36.0,  1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0   //at 5,6,7,8
    };
    // convertion matrix for D2Q9 MRT
    std::vector<std::vector<double>> M =
    {
        { 1.0,  1.0,  1.0,  1.0,  1.0, 1.0,  1.0,  1.0,  1.0},
        {-4.0, -1.0, -1.0, -1.0, -1.0, 2.0,  2.0,  2.0,  2.0},
        { 4.0, -2.0, -2.0, -2.0, -2.0, 1.0,  1.0,  1.0,  1.0},
        { 0.0,  1.0,  0.0, -1.0,  0.0, 1.0, -1.0, -1.0,  1.0},
        { 0.0, -2.0,  0.0,  2.0,  0.0, 1.0, -1.0, -1.0,  1.0},
        { 0.0,  0.0,  1.0,  0.0, -1.0, 1.0,  1.0, -1.0, -1.0},
        { 0.0,  0.0, -2.0,  0.0,  2.0, 1.0,  1.0, -1.0, -1.0},
        { 0.0,  1.0, -1.0,  1.0, -1.0, 0.0,  0.0,  0.0,  0.0},
        { 0.0,  0.0,  0.0,  0.0,  0.0, 1.0, -1.0,  1.0, -1.0}
    };
    // convertion inverse matrix for D2Q9 MRT
    std::vector<std::vector<double>> Minv =
    {
        {1.0/9.0, -4.0/36.0,  4.0/36.0,  0.0/6.0,  0.0/12.0,  0.0/6.0,  0.0/12.0,  0.0/4.0,  0.0/4.0},
        {1.0/9.0, -1.0/36.0, -2.0/36.0,  1.0/6.0, -2.0/12.0,  0.0/6.0,  0.0/12.0,  1.0/4.0,  0.0/4.0},
        {1.0/9.0, -1.0/36.0, -2.0/36.0,  0.0/6.0,  0.0/12.0,  1.0/6.0, -2.0/12.0, -1.0/4.0,  0.0/4.0},
        {1.0/9.0, -1.0/36.0, -2.0/36.0, -1.0/6.0,  2.0/12.0,  0.0/6.0,  0.0/12.0,  1.0/4.0,  0.0/4.0},
        {1.0/9.0, -1.0/36.0, -2.0/36.0,  0.0/6.0,  0.0/12.0, -1.0/6.0,  2.0/12.0, -1.0/4.0,  0.0/4.0},
        {1.0/9.0,  2.0/36.0,  1.0/36.0,  1.0/6.0,  1.0/12.0,  1.0/6.0,  1.0/12.0,  0.0/4.0,  1.0/4.0},
        {1.0/9.0,  2.0/36.0,  1.0/36.0, -1.0/6.0, -1.0/12.0,  1.0/6.0,  1.0/12.0,  0.0/4.0, -1.0/4.0},
        {1.0/9.0,  2.0/36.0,  1.0/36.0, -1.0/6.0, -1.0/12.0, -1.0/6.0, -1.0/12.0,  0.0/4.0,  1.0/4.0},
        {1.0/9.0,  2.0/36.0,  1.0/36.0,  1.0/6.0,  1.0/12.0, -1.0/6.0, -1.0/12.0,  0.0/4.0, -1.0/4.0}
    };
};

struct latticeModelD3Q19
{

};

#endif // LATTICEMODEL_HXX_INCLUDED

#ifndef ZOUHENODE_HPP_INCLUDED
#define ZOUHENODE_HPP_INCLUDED

#include <vector>

#include "latticeBase.hpp"
#include "collisionBase.hxx"
#include "latticeNode.hxx"

#include "latticeModel.hxx"

#include "boundaryNode.hxx"

class ZouHeNode: public boundaryNode
{
    public:
        // Constructor: Creates Zou/He velocity boundary nodes
        // param lm lattice model which contains information on the number of rows,
        //       columns, dimensions, discrete directions and lattice velocity
        // param cm collision model which contains information on lattice density
        ZouHeNode
        (
            latticeBase &lb,
            collisionBase &cb,
            latticeModelD2Q9 &D2Q9,
            fluidField &field
        );
        // Destructor
        ~ZouHeNode() = default;
        // Adds a Zou/He velocity node to the nodes vector
        // param x x-coordinate of the node
        // param y y-coordinate of the node
        // param u_x x-velocity of the node
        // param u_y y-velocity of the node
        void addNode
        (
            std::size_t x,
            std::size_t y,
            double u_x,
            double u_y
        );
        // Updates the boundary nodes based on "On pressure and velocity boundary
        // conditions for the lattice Boltzmann"
        // param df lattice distribution functions stored row-wise in a 2D vector
        // param is_modify_stream boolean toggle for half-way bounceback nodes to
        //       perform functions during stream, set to FALSE for Zou/He velocity nodes
        void updateNode
        (
            std::vector<std::vector<double>> &df,
            bool is_modify_stream
        );
        // Updates the non-corner nodes
        // param df lattice distribution function stored row-wise in a 2D vector
        // param node Zou/He velocity node which contains information on the position
        //       of the boundary node and velocities of the node
        void updateEdge
        (
            std::vector<std::vector<double>> &df,
            latticeNode &node
        );
        // Updates the corner nodes, first-order expolation for node density
        // param df lattice distribution function stored row-wise in a 2D vector
        // param node Zou/He velocity node which contains information on the position
        //       of the boundary node and velocities of the node
        void updateCorner
        (
            std::vector<std::vector<double>> &df,
            latticeNode &node
        );
        // Toggles behaviour of Zou/He nodes when used as outlet, boundary node
        // velocity will be extrapolated (1st order) from the neighbouring nodes
        void toggleNormalFlow();
        // Boundary nodes stored in a 1D vector
        std::vector<latticeNode> nodes;
    protected:
        // Collision model which contains information on lattice density
        collisionBase &cb_;
        // Boolean toggle for open boundary condition (outlet)
        bool is_normal_flow_;
        // Additional constants beta1, beta2 and beta3 since the dx = dt = 1 condition
        // is not always maintained
        double beta1_;
        double beta2_;
        double beta3_;
//      std::vector<bool> is_corner_;
//      std::vector<bool> knowns_;
    private:
        // define fluid field;
        fluidField &field_;
        // define lattice model
        latticeModelD2Q9 &D2Q9_;
};

#endif // ZOUHENODE_HPP_INCLUDED

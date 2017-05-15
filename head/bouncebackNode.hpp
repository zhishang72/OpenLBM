#ifndef BOUNCEBACKNODE_HPP_INCLUDED
#define BOUNCEBACKNODE_HPP_INCLUDED

#include <vector>

#include "latticeBase.hpp"
#include "collisionBase.hxx"

#include "latticeModel.hxx"

#include "latticeNode.hxx"
#include "boundaryNode.hxx"

class bouncebackNode: public boundaryNode
{
    public:
        // Creates a full-way bounceback nodes according to "http://lbmworkshop.com/wp
        // -content/uploads/2011/08/Straight_boundaries.pdf"
        // param lb LatticeModel to provide information on number of rows, columns,
        //       dimensions, discrete directions and lattice velocity
        // param cb CollisionModel to indicate which nodes to be skipped during the
        //       collision step
        bouncebackNode
        (
            latticeBase &lb,
            collisionBase *cb,
            latticeModelD2Q9 &D2Q9,
            fluidField &field
        );
        // Creates a half-way bounceback node according to "http://lbmworkshop.com/wp-
        // content/uploads/2011/08/Straight_boundaries.pdf"
        // param lb LatticeModel to provide information on number of rows, columns,
        //       dimensions, discrete directions and lattice velocity
        // param sb StreamModel to copy the prestream node distribution functions
        //       so they can be bounced back in the same time step
        bouncebackNode
        (
            latticeBase &lb,
            streamBase *sb,
            latticeModelD2Q9 &D2Q9,
            fluidField &field
        );
        // Override copy constructor
        bouncebackNode(const bouncebackNode&) = default;
        // Override copy assignment
        bouncebackNode& operator= (const bouncebackNode&) = default;
        // Destructor
        ~bouncebackNode() = default;
        // Adds a bounceback node
        // param x x-coordinate of the node
        // param y y-coordinate of the node
        void addNode(std::size_t x, std::size_t y);
        // Adds a bounceback node
        // param x x-coordinate of the node
        // param y y-coordinate of the node
        // param z z-coordinate of the node
        void addNode(std::size_t x, std::size_t y, std::size_t z);
        // Performs the bounceback boundary condition on the boundary nodes based on
        // the type of bounceback nodes used.
        // Full-way bounceback: Reflects all the node distribution functions in the
        //       opposite direction (except center distribution function)
        // Half-way bounceback: Copies the prestream node distribution functions
        //       before streaming. Updates the post-stream unknown distribution functions
        //       with the prestream distribution functions in the opposite directions
        // param df lattice distribution functions store row-wise in a 2D vector
        // param is_modify_stream Boolean toggle for half-way bounceback as it has
        //       both pre-stream and post-stream functions. Used to fit in with how
        //       all boundary conditions are called
        void updateNode
        (
            std::vector<std::vector<double>> &df,
            bool is_modify_stream
        );
    protected:
        // Vector used to store information about the boundary nodes such as their
        // position at the lattice node
        std::vector<latticeNode> nodes;
        // Pointer to collision model as half-way bounceback nodes do not require
        // collision models and NULL references can't be declared
        collisionBase *cb_ = nullptr;
        // Pointer to stream model as full-way bounceback nodes do not require stream
        // models and NULL references can't be declared
        streamBase *sb_ = nullptr;
    private:
        // define fluid field;
        fluidField &field_;
        // define lattice model
        latticeModelD2Q9 &D2Q9_;
};

#endif // BOUNCEBACKNODE_HPP_INCLUDED

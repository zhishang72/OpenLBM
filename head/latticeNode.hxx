#ifndef LATTICENODE_HXX_INCLUDED
#define LATTICENODE_HXX_INCLUDED

#include <iostream>

class latticeNode
{
    public:
        // Constructor: Creates a node which stores information about its position in
        // the lattice
        // param x_position x-coordinate of the node
        // param y_position y-coordinate of the node
        // param n number of the lattice, for calculating index in the distribution function vector
        latticeNode
        (
            std::size_t x_position,
            std::size_t y_position,
            std::size_t n
        )
        : x_node {x_position},
          y_node {y_position},
          n_node {n},
          df_node {}
        {};
        // Constructor: Creates a node which stores information about its position in
        // the lattice
        // param x_position x-coordinate of the node
        // param y_position y-coordinate of the node
        // param z_position z-coordinate of the node
        // param n number of the lattice, for calculating index in the distribution function vector
        latticeNode
        (
            std::size_t x_position,
            std::size_t y_position,
            std::size_t z_position,
            std::size_t n
        )
        : x_node {x_position},
          y_node {y_position},
          z_node {z_position},
          n_node {n},
          df_node {}
        {};
        // Constructor: Creates a node which contains information on its position,
        // 2 double values, 1 boolean value and 1 integer value. To be used by
        // Zou/He velocity node
        // param x_position x-coordinate of the node
        // param y_position y-coordinate of the node
        // param nx number of columns of the lattice model
        // param v_x first double value to be stored in vector of double values, used
        //       as x-velocity of Zou/He velocity node
        // param v_y second double value to be stored in vector of double values,
        //       used as y-velocity of Zou/He velocity node
        // param index_i first integer value, used to indicate which side the node
        //       belongs to for non-corner nodes, and which corner the node belongs
        //       to for corner nodes in Zou/He velocity nodes
        latticeNode
        (
            std::size_t x_position,
            std::size_t y_position,
            std::size_t n,
            double v_x,
            double v_y,
            bool is_corner_node,
            int index
        )
        : x_node {x_position},
          y_node {y_position},
          n_node {n},
          u_node {v_x, v_y},
          corner {is_corner_node},
          index_i {index},
          df_node {}
        {};
        // Constructor: Creates a node which contains information on its position,
        // 1 double value, 1 boolean value and 1 integer value. To be used by
        // Zou/He pressure node
        // param x_position x-coordinate of the node
        // param y_position y-coordinate of the node
        // param nx number of columns of the lattice model
        // param pressure first double value, used as node pressure in Zou/He pressure
        //       nodes
        // param is_corner_node first boolean value, used to indicate if node is a corner node
        //       Zou/He pressure nodes
        // param index first integer value, used to indicate which side the node
        //       belongs to for non-corner nodes, and which corner the node belongs
        //       to for corner nodes in Zou/He pressure nodes
        latticeNode
        (
            std::size_t x_position,
            std::size_t y_position,
            std::size_t n,
            double pressure,
            bool is_corner_node,
            int index
        )
        : x_node {x_position},
          y_node {y_position},
          n_node {n},
          pressure_node {pressure},
          corner {is_corner_node},
          index_i {index},
          df_node {}
        {};
        // Virtual destructor since we are deriving from this class
        virtual ~latticeNode() = default;
        // x-coordinate
        std::size_t x_node;
        // y-coordinate
        std::size_t y_node;
        // z-coordinate
        std::size_t z_node;
        // Index in distribution function vector
        std::size_t n_node;
        // Double value, used as node pressure in Zou/He pressure nodes
        double pressure_node;
        // Vector containing double values, used as node velocity in Zou/He velocity nodes
        std::vector<double> u_node;
        // Boolean value, used to indicate if node is a corner node Zou/He velocity and pressure nodes
        bool corner;
        // Integer value, used to indicate which side the node belongs to for
        // non-corner and corner nodes, and which corner the node belongs to for corner nodes in
        // Zou/He velocity and pressure nodes
        int index_i;
        // Distribution functions of a single node, used by half-way bounceback nodes
        std::vector<double> df_node;
};

#endif // LATTICENODE_HXX_INCLUDED

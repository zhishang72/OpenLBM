#include <iostream>
#include <vector>

#include "latticeBase.hpp"
#include "collisionBase.hxx"
#include "streamBase.hxx"

#include "latticeModel.hxx"

#include "latticeNode.hxx"
#include "boundaryNode.hxx"

#include "bouncebackNode.hpp"

bouncebackNode::bouncebackNode
(
    latticeBase &lb,
    collisionBase *cb,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: boundaryNode(true, true, lb),
  nodes {},
  cb_ {cb},
  D2Q9_ (D2Q9),
  field_ (field)
{}

bouncebackNode::bouncebackNode
(
    latticeBase &lb,
    streamBase *sb,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: boundaryNode(true, true, lb),
  nodes {},
  sb_ {sb},
  D2Q9_ (D2Q9),
  field_ (field)
{}

void bouncebackNode::addNode
(
    std::size_t x,
    std::size_t y
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto n = y * nx + x;
    nodes.push_back(latticeNode(x, y, n));
    if (cb_) cb_->addNodeToSkip(n);
    // add node position to position vector
    position.push_back(n);
}

void bouncebackNode::addNode
(
    std::size_t x,
    std::size_t y,
    std::size_t z
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto n = y * nx + x;
    nodes.push_back(latticeNode(x, y, z, n));
    if (cb_) cb_->addNodeToSkip(n);
    // add node position to position vector
    position.push_back(n);
}

void bouncebackNode::updateNode
(
    std::vector<std::vector<double>> &df,
    bool is_modify_stream
)
{
    if (is_modify_stream)
    {
        const auto nx = lb_.getNumberOfNx();
        const auto ny = lb_.getNumberOfNy();
        for (auto &node : nodes)
        {
            const auto n = node.n_node;
            const auto left =   n % nx == 0;
            const auto right =  n % nx == nx - 1;
            const auto bottom = n / nx == 0;
            const auto top =    n / nx == ny - 1;
            if (bottom)          df[n][D2Q9_.N]  = node.df_node[D2Q9_.S];
            if (top)             df[n][D2Q9_.S]  = node.df_node[D2Q9_.N];
            if (left)            df[n][D2Q9_.E]  = node.df_node[D2Q9_.W];
            if (right)           df[n][D2Q9_.W]  = node.df_node[D2Q9_.E];
            if (bottom || left)  df[n][D2Q9_.NE] = node.df_node[D2Q9_.SW];
            if (bottom || right) df[n][D2Q9_.NW] = node.df_node[D2Q9_.SE];
            if (top || right)    df[n][D2Q9_.SW] = node.df_node[D2Q9_.NE];
            if (top || left)     df[n][D2Q9_.SE] = node.df_node[D2Q9_.NW];
        }  // node
    }
    else
    {
        if (cb_)
        {
            for (auto &node : nodes)
            {
                const auto n = node.n_node;
                auto temp_node = df[n];
                df[n][D2Q9_.E]  = temp_node[D2Q9_.W];
                df[n][D2Q9_.N]  = temp_node[D2Q9_.S];
                df[n][D2Q9_.W]  = temp_node[D2Q9_.E];
                df[n][D2Q9_.S]  = temp_node[D2Q9_.N];
                df[n][D2Q9_.NE] = temp_node[D2Q9_.SW];
                df[n][D2Q9_.NW] = temp_node[D2Q9_.SE];
                df[n][D2Q9_.SW] = temp_node[D2Q9_.NE];
                df[n][D2Q9_.SE] = temp_node[D2Q9_.NW];
                node.df_node = df[n];
            }  // node
        }
        if (sb_)
        {
            for (auto &node : nodes) node.df_node = df[node.n_node];
        }
    }
}

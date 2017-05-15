#include <iostream>
#include <stdexcept>
#include <vector>

#include "latticeBase.hpp"
#include "collisionBase.hxx"
#include "latticeNode.hxx"

#include "latticeModel.hxx"

#include "ZouHeNode.hpp"
#include "latticeNode.hxx"

ZouHeNode::ZouHeNode
(
    latticeBase &lb,
    collisionBase &cb,
    latticeModelD2Q9 &D2Q9,
    fluidField &field
)
: boundaryNode(false, false, lb),
  nodes {},
  cb_ (cb),
  is_normal_flow_ {false},
  beta1_ {},
  beta2_ {},
  beta3_ {},
  D2Q9_ (D2Q9),
  field_ (field)
{
    const auto c = lb_.getLatticeSpeed();
    const auto cs_sqr = c * c / 3.0;
    beta1_ = c / (9.0 *cs_sqr);
    beta2_ = 0.5 / c;
    beta3_ = beta2_ - beta1_;
}

void ZouHeNode::addNode
(
    std::size_t x,
    std::size_t y,
    double u_x,
    double u_y
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNy();
    const auto n  = y * nx + x;
    const auto left   = x == 0;
    const auto right  = x == nx - 1;
    const auto bottom = y == 0;
    const auto top    = y == ny - 1;

    auto edge_i = -1;
    if (right)  edge_i = 0;
    if (top)    edge_i = 1;
    if (left)   edge_i = 2;
    if (bottom) edge_i = 3;
    // adds a corner node
    if ((top || bottom) && (left || right))
    {
        auto corner_i = -1;
        if (bottom && left)  corner_i = 0;
        if (bottom && right) corner_i = 1;
        if (top && left)     corner_i = 2;
        if (top && right)    corner_i = 3;
        nodes.push_back(latticeNode(x, y, n, u_x, u_y, true, corner_i));
    }
    // adds a side node
    else
    {
        nodes.push_back(latticeNode(x, y, n, u_x, u_y, false, edge_i));
    }
}

void ZouHeNode::updateNode
(
    std::vector<std::vector<double>> &df,
    bool is_modify_stream
)
{
    if (!is_modify_stream)
    {
        for (auto node : nodes)
        {
            if (node.corner)
            {
                ZouHeNode::updateCorner(df, node);
            }
            else
            {
                ZouHeNode::updateEdge(df, node);
            }
        }  // n
    }
}

void ZouHeNode::updateEdge
(
    std::vector<std::vector<double>> &df,
    latticeNode &node
)
{
    const auto n = node.n_node;
    const auto nx = lb_.getNumberOfNx();
    const auto c = lb_.getLatticeSpeed();
    switch(node.index_i)
    {
        case 0:
        {  // right
            auto vel = is_normal_flow_ ? field_.u[n - 1] : node.u_node;
            const auto rho_node = (df[n][0] + df[n][D2Q9_.N] + df[n][D2Q9_.S] + 2.0 * (df[n][D2Q9_.E] +
                                   df[n][D2Q9_.NE] + df[n][D2Q9_.SE])) / (1.0 + vel[0] / c);
            const auto df_diff = 0.5 * (df[n][D2Q9_.S] - df[n][D2Q9_.N]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.W]  = df[n][D2Q9_.E] - 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.NW] = df[n][D2Q9_.SE] + df_diff - beta3_ * vel[0] + beta2_ * vel[1];
            df[n][D2Q9_.SW] = df[n][D2Q9_.NE] - df_diff - beta3_ * vel[0] - beta2_ * vel[1];
            break;
        }
        case 1:
        {  // top
            auto vel = is_normal_flow_ ? field_.u[n - nx] : node.u_node;
            const auto rho_node = (df[n][0] + df[n][D2Q9_.E] + df[n][D2Q9_.W] + 2.0 * (df[n][D2Q9_.N] +
                                   df[n][D2Q9_.NE] + df[n][D2Q9_.NW])) / (1.0 + vel[1] / c);
            const auto df_diff = 0.5 * (df[n][D2Q9_.E] - df[n][D2Q9_.W]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.S]  = df[n][D2Q9_.N] - 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.SW] = df[n][D2Q9_.NE] + df_diff - beta2_ * vel[0] - beta3_ * vel[1];
            df[n][D2Q9_.SE] = df[n][D2Q9_.NW] - df_diff + beta2_ * vel[0] - beta3_ * vel[1];
            break;
        }
        case 2:
        {  // left
            auto vel = is_normal_flow_ ? field_.u[n + 1] : node.u_node;
            const auto rho_node = (df[n][0] + df[n][D2Q9_.N] + df[n][D2Q9_.S] + 2.0 * (df[n][D2Q9_.W] +
                                   df[n][D2Q9_.NW] + df[n][D2Q9_.SW])) / (1.0 - vel[0] / c);
            const auto df_diff = 0.5 * (df[n][D2Q9_.S] - df[n][D2Q9_.N]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.E]  = df[n][D2Q9_.W] + 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.NE] = df[n][D2Q9_.SW] + df_diff + beta3_ * vel[0] + beta2_ * vel[1];
            df[n][D2Q9_.SE] = df[n][D2Q9_.NW] - df_diff + beta3_ * vel[0] - beta2_ * vel[1];
            break;
        }
        case 3:
        {  // bottom
            auto vel = is_normal_flow_ ? field_.u[n + nx] : node.u_node;
            const auto rho_node = (df[n][0] + df[n][D2Q9_.E] + df[n][D2Q9_.W] + 2.0 * (df[n][D2Q9_.S] +
                                   df[n][D2Q9_.SW] + df[n][D2Q9_.SE])) / (1.0 - vel[1] / c);
            const auto df_diff = 0.5 * (df[n][D2Q9_.W] - df[n][D2Q9_.E]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.N]  = df[n][D2Q9_.S] + 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.NE] = df[n][D2Q9_.SW] + df_diff + beta2_ * vel[0] + beta3_ * vel[1];
            df[n][D2Q9_.NW] = df[n][D2Q9_.SE] - df_diff - beta2_ * vel[0] + beta3_ * vel[1];
            break;
        }
        default:
        {
            throw std::runtime_error("Not a side");
        }
    }
}

void ZouHeNode::updateCorner
(
    std::vector<std::vector<double>> &df,
    latticeNode &node
)
{
    const auto n = node.n_node;
    auto vel = node.u_node;
    const auto nx = lb_.getNumberOfNx();
    const auto nc = lb_.getNumberOfDirections();
    switch (node.index_i)
    {
        case 0:
        {  // bottom-left
            auto rho_node = 0.5 * (cb_.rho_[n + nx] + cb_.rho_[n + 1]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.E]  = df[n][D2Q9_.W] + 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.N]  = df[n][D2Q9_.S] + 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.NE] = df[n][D2Q9_.SW] + 0.5 * beta1_ * vel[0] + 0.5 * beta1_ * vel[1];
            df[n][D2Q9_.NW] = -0.5 * beta3_ * vel[0] + 0.5 * beta3_ * vel[1];
            df[n][D2Q9_.SE] = 0.5 * beta3_ * vel[0] - 0.5 * beta3_ * vel[1];
            for (auto i = 1u; i < nc; ++i) rho_node -= df[n][i];
            df[n][0] = rho_node;
            break;
        }
        case 1:
        {  // bottom-right
            auto rho_node = 0.5 * (cb_.rho_[n + nx] + cb_.rho_[n - 1]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.W]  = df[n][D2Q9_.E] - 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.N]  = df[n][D2Q9_.S] + 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.NW] = df[n][D2Q9_.SE] - 0.5 * beta1_ * vel[0] + 0.5 * beta1_ * vel[1];
            df[n][D2Q9_.NE] = 0.5 * beta3_ * vel[0] + 0.5 * beta3_ * vel[1];
            df[n][D2Q9_.SW] = -0.5 * beta3_ * vel[0] - 0.5 * beta3_ * vel[1];
            for (auto i = 1u; i < nc; ++i) rho_node -= df[n][i];
            df[n][0] = rho_node;
            break;
        }
        case 2:
        {  // top-left
            auto rho_node = 0.5 * (cb_.rho_[n - nx] + cb_.rho_[n + 1]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.E]  = df[n][D2Q9_.W] + 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.S]  = df[n][D2Q9_.N] - 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.SE] = df[n][D2Q9_.NW] + 0.5 * beta1_ * vel[0] - 0.5 * beta1_ * vel[1];
            df[n][D2Q9_.NE] = 0.5 * beta3_ * vel[0] + 0.5 * beta3_ * vel[1];
            df[n][D2Q9_.SW] = -0.5 * beta3_ * vel[0] - 0.5 * beta3_ * vel[1];
            for (auto i = 1u; i < nc; ++i) rho_node -= df[n][i];
            df[n][0] = rho_node;
            break;
        }
        case 3:
        {  // top-right
            auto rho_node = 0.5 * (cb_.rho_[n - nx] + cb_.rho_[n - 1]);
            for (auto &u : vel) u *= rho_node;
            df[n][D2Q9_.W]  = df[n][D2Q9_.E] - 2.0 * beta1_ * vel[0];
            df[n][D2Q9_.S]  = df[n][D2Q9_.N] - 2.0 * beta1_ * vel[1];
            df[n][D2Q9_.SW] = df[n][D2Q9_.NE] - 0.5 * beta1_ * vel[0] - 0.5 * beta1_ * vel[1];
            df[n][D2Q9_.NW] = -0.5 * beta3_ * vel[0] + 0.5 * beta3_ * vel[1];
            df[n][D2Q9_.SE] = 0.5 * beta3_ * vel[0] - 0.5 * beta3_ * vel[1];
            for (auto i = 1u; i < nc; ++i) rho_node -= df[n][i];
            df[n][0] = rho_node;
            break;
        }
        default:
        {
            throw std::runtime_error("Not a corner");
        }
    }
}

void ZouHeNode::toggleNormalFlow()
{
    is_normal_flow_ = true;
}

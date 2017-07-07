#include <iostream>
#include <stdexcept>  // std::runtime_error
#include <vector>

#include "latticeBoltzmann.hpp"
#include "latticeBase.hpp"
#include "collisionBase.hxx"
#include "streamBase.hxx"
#include "boundaryNode.hxx"

latticeBoltzmann::latticeBoltzmann
(
    latticeBase &lb,
    collisionBase &cb,
    streamBase &sb
)
: lb_ (lb),
  cb_ (cb),
  sb_ (sb),
  df {},
  bn_ {}
{
    cb_.computefEq();
    df = cb_.eqdf;
}

void latticeBoltzmann::addBoundaryNode
(
    boundaryNode *bn
)
{
    bn_.push_back(bn);
}

void latticeBoltzmann::takeStep()
{
    cb_.computefEq();
    cb_.collide(df);
    for(auto bdr : bn_)
    {
        if(bdr->prestream) bdr->updateNode(df, false);
    }  // bdr
    df = sb_.stream(df);
    for(auto bdr : bn_)
    {
        if (bdr->streaming) bdr->updateNode(df, true);
        if (!bdr->prestream) bdr->updateNode(df, false);
    }  // bdr
    cb_.computeMacroscopicProperties(df);
}

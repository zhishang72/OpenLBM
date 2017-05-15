#include <iostream>
#include <vector>

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "latticeD2Q9.hpp"

latticeD2Q9::latticeD2Q9
(
    std::size_t num_nx,
    std::size_t num_ny,
    double dl,
    double dt,
    latticeModelD2Q9 &D2Q9
)
: latticeBase(num_nx, num_ny, 2, 9, dl, dt),
  D2Q9_ (D2Q9)
{
    auto c = latticeBase::getLatticeSpeed();
    for (auto &i : D2Q9_.e)
    {
        for (auto &d : i)
            d *= c;
    } // i
}

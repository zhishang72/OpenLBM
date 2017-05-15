#include <cmath>
#include <iostream>
#include <vector>

#include "latticeModel.hxx"

#include "latticeD2Q9.hpp"
#include "collisionD2Q9_BGK.hpp"
#include "streamD2Q9.hpp"
#include "bouncebackNode.hpp"
#include "ZouHeNode.hpp"
#include "latticeBoltzmann.hpp"
#include "boundaryNode.hxx"
#include "result.hpp"

int main()
{
    std::size_t ny = 256;
    std::size_t nx = 256;

    auto dt = 0.0001;
    auto dl = sqrt(dt);

    auto rho0_f = 1.0;
    auto visco_f = 1.0 / 18.0;

    std::vector<double> u0 = {0.0, 0.0};
    auto u_lid = 31.6;
    auto v_lid = 0.0;

    fluidField field
    (
        nx,
        ny,
        u0
    );
    latticeModelD2Q9 D2Q9;

    latticeD2Q9 lattice
    (
        nx,
        ny,
        dl,
        dt,
        D2Q9
    );

    collisionD2Q9_BGK collision
    (
        lattice,
        visco_f,
        rho0_f,
        D2Q9,
        field
    );

    streamD2Q9 stream
    (
        lattice,
        D2Q9
    );

    bouncebackNode bbnode
    (
        lattice,
        //&collision,
        &stream,
        D2Q9,
        field
    );
    ZouHeNode zhnode
    (
        lattice,
        collision,
        D2Q9,
        field
    );

    latticeBoltzmann run
    (
        lattice,
        collision,
        stream
    );

    for (auto y = 0u; y < ny; ++y)
    {
        bbnode.addNode(0, y);
        bbnode.addNode(nx - 1, y);
    }

    for (auto x = 0u; x < nx; ++x)
    {
        bbnode.addNode(x, 0);
        zhnode.addNode(x, ny - 1, u_lid, v_lid);
    }

    run.addBoundaryNode(&bbnode);
    run.addBoundaryNode(&zhnode);

    result results
    (
        lattice,
        field
    );

    for (auto t = 0u; t <= nx; ++t)
    {
        run.takeStep();
        if (t % (nx/8) == 0)
        {
            results.writeResultVTK(t);
        }
        std::cout << "t= " << t << std::endl;
    }
    return 0;
}

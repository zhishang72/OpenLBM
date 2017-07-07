#include <cmath>
#include <iostream>
#include <vector>

#include "latticeModel.hxx"

#include "latticeD2Q9.hpp"
#include "collisionD2Q9_BGK.hpp"
#include "collisionD2Q9_MRT.hpp"
#include "streamD2Q9.hpp"
#include "bouncebackNode.hpp"
#include "ZouHeNode.hpp"
#include "latticeBoltzmann.hpp"
#include "boundaryNode.hxx"
#include "momentComputing.h"
#include "result.hpp"

int main()
{
    std::size_t ny = 256;
    std::size_t nx = 256;
    auto tolerance = 1.0e-3;

    auto dt = 1.0;
    auto dl = sqrt(dt);

    auto rho0_f = 1.0;
    auto visco_f = 1.0 / 18.0;

    std::vector<double> u0 = {0.0, 0.0};
    auto u_lid = 0.3;
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

    //collisionD2Q9_BGK collision
    collisionD2Q9_MRT collision
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
        auto uprev = field.u;
        run.takeStep();
        auto error = checkError(uprev, field.u);
        std::cout << "t= " << t << "; error= " << error << std::endl;
        if (t % (nx/8) == 0)
        {
            results.writeResultVTK(t);
            if (checkSteadyState(uprev, field.u, tolerance))
            {
                results.writeResultVTK(t);
                break;
            }
        }
    }
    return 0;
}

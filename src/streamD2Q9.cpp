#include <vector>

#include "latticeModel.hxx"

#include "latticeBase.hpp"
#include "streamD2Q9.hpp"

streamD2Q9::streamD2Q9
(
    latticeBase &lb,
    latticeModelD2Q9 &D2Q9
)
: streamBase(lb),
  D2Q9_ (D2Q9)
{}

std::vector<std::vector<double>> streamD2Q9::stream
(
    const std::vector<std::vector<double>> &df
)
{
    const auto nx = lb_.getNumberOfNx();
    const auto ny = lb_.getNumberOfNx();
    auto temp_df = df;
    // Streaming
    for (auto n = 0u; n < nx * ny; ++n)
    {
        const auto left   = n % nx == 0;
        const auto right  = n % nx == nx - 1;
        const auto bottom = n / nx == 0;
        const auto top    = n / nx == ny - 1;
        if (!left)              temp_df[n][D2Q9_.E]  = df[n -  1][D2Q9_.E];
        if (!bottom)            temp_df[n][D2Q9_.N]  = df[n - nx][D2Q9_.N];
        if (!right)             temp_df[n][D2Q9_.W]  = df[n +  1][D2Q9_.W];
        if (!top)               temp_df[n][D2Q9_.S]  = df[n + nx][D2Q9_.S];
        if (!(bottom || left))  temp_df[n][D2Q9_.NE] = df[n - nx - 1][D2Q9_.NE];
        if (!(bottom || right)) temp_df[n][D2Q9_.NW] = df[n - nx + 1][D2Q9_.NW];
        if (!(top || right))    temp_df[n][D2Q9_.SW] = df[n + nx + 1][D2Q9_.SW];
        if (!(top || left))     temp_df[n][D2Q9_.SE] = df[n + nx - 1][D2Q9_.SE];
    }  // n
    return temp_df;
}

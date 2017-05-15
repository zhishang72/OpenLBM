#include "latticeBase.hpp"

latticeBase::latticeBase
(
    std::size_t nx,
    std::size_t ny,
    std::size_t num_dims,
    std::size_t num_dirs,
    double dl,
    double dt
)
: number_of_nx_ {nx},
  number_of_ny_ {ny},
  number_of_dimensions_ {num_dims},
  number_of_directions_ {num_dirs},
  space_step_ {dl},
  time_step_ {dt}
{}

latticeBase::latticeBase
(
    std::size_t nx,
    std::size_t ny,
    std::size_t nz,
    std::size_t num_dims,
    std::size_t num_dirs,
    double dl,
    double dt
)
: number_of_nx_ {nx},
  number_of_ny_ {ny},
  number_of_nz_ {nz},
  number_of_dimensions_ {num_dims},
  number_of_directions_ {num_dirs},
  space_step_ {dl},
  time_step_ {dt}
{}

std::size_t latticeBase::getNumberOfNx() const
{
    return number_of_nx_;
}

std::size_t latticeBase::getNumberOfNy() const
{
    return number_of_ny_;
}

std::size_t latticeBase::getNumberOfNz() const
{
    return number_of_nz_;
}

std::size_t latticeBase::getNumberOfDimensions() const
{
    return number_of_dimensions_;
}

std::size_t latticeBase::getNumberOfDirections() const
{
    return number_of_directions_;
}

double latticeBase::getSpaceStep() const
{
    return space_step_;
}

double latticeBase::getTimeStep() const
{
    return time_step_;
}

double latticeBase::getLatticeSpeed() const
{
    return c_;
}

bool latticeBase::checkInput()
{
    return number_of_dimensions_ == 0 || number_of_directions_ == 0 ||
           number_of_nx_ == 0 || number_of_ny_ == 0 || number_of_nz_ == 0;
}

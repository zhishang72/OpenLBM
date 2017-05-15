#include <vector>

#include "boundaryNodes.hpp"
#include "latticeBase.hpp"

// have to use parenthesis for reference initializing in initializer list due to
// bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=50025
boundaryNodes::boundaryNodes
(
    bool is_prestream,
    bool is_during_stream,
    latticeModel &lb
)
: prestream {is_prestream},
  during_stream {is_during_stream},
  position {},
  lb_ (lb)
{}

#include <stdbool.h>
#include "param.h"

const bool param_implicit_x = true;
const bool param_implicit_y = false;
#if NDIMS == 3
const bool param_implicit_z = false;
#endif


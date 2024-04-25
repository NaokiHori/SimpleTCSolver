#include "fluid.h"

double fluid_compute_momentum_diffusivity (
    const fluid_t * fluid
) {
  return 1. / fluid->Re;
}

double fluid_compute_scalar_diffusivity (
    const fluid_t * fluid
) {
  return 1. / fluid->Re / fluid->Pr;
}


#if !defined(FLUID_SOLVER_H)
#define FLUID_SOLVER_H

#include "array.h"
#include "domain.h"
#include "fluid.h"

// initialiser of fluid_t
extern int fluid_init(
    const char dirname_ic[],
    const domain_t * domain,
    fluid_t * fluid
);

// save flow field
extern int fluid_save(
    const char dirname[],
    const domain_t * domain,
    const fluid_t * fluid
);

extern double fluid_compute_momentum_diffusivity (
    const fluid_t * fluid
);

extern double fluid_compute_scalar_diffusivity (
    const fluid_t * fluid
);

// update fields using the previously-computed RK source terms
extern int fluid_predict_field(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// compute scalar potential by solving Poisson equation
extern int fluid_compute_potential(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// correct velocity field using scalar potential
extern int fluid_correct_velocity(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// update pressure
extern int fluid_update_pressure(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// exchange halos and impose boundary conditions

extern int fluid_update_boundaries_ux(
    const domain_t * domain,
    array_t * array
);

extern int fluid_update_boundaries_uy(
    const domain_t * domain,
    array_t * array
);

extern int fluid_update_boundaries_uz(
    const domain_t * domain,
    array_t * array
);

extern int fluid_update_boundaries_p(
    const domain_t * domain,
    array_t * array
);

extern int fluid_update_boundaries_psi(
    const domain_t * domain,
    array_t * array
);

extern int fluid_update_boundaries_t(
    const domain_t * domain,
    array_t * array
);

#endif // FLUID_SOLVER_H

#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"
#include "structure.h"

// definition of a structure fluid_t_
/**
 * @struct fluid_t_ (or fluid_t)
 * @brief struct storing fluid-related variables
 * @var ux, uy, uz             : velocity in each direction
 * @var p, psi                 : pressure, scalar potential
 * @var srcuxa, srcuxb, srcuxg : Runge-Kutta source terms, ux
 * @var srcuya, srcuyb, srcuyg : Runge-Kutta source terms, uy
 * @var srcuza, srcuzb, srcuzg : Runge-Kutta source terms, uz
 * @var diffusivity            : diffusivity (pre-factor in front of diffusive terms)
 */
struct fluid_t_ {
  array_t * restrict ux;
  array_t * restrict uy;
  array_t * restrict uz;
  array_t * restrict p;
  array_t * restrict psi;
  array_t * restrict srcuxa, * restrict srcuxb, * restrict srcuxg;
  array_t * restrict srcuya, * restrict srcuyb, * restrict srcuyg;
  array_t * restrict srcuza, * restrict srcuzb, * restrict srcuzg;
  double diffusivity;
};

// constructor and destructor
extern int fluid_init(const char dirname_ic[restrict], const domain_t * restrict domain, fluid_t * restrict * fluid);
extern int fluid_finalise(fluid_t * restrict fluid);

// compute right-hand-side of Runge-Kutta scheme
extern int fluid_compute_rhs(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid);

// update velocity field and its clean-up function
extern int fluid_update_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// compute scalar potential and its clean-up function
extern int fluid_compute_potential(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// correct velocity field using scalar potential
extern int fluid_correct_velocity(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// update pressure
extern int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

// exchange halos and impose boundary conditions
extern int fluid_update_boundaries_ux(const domain_t * restrict domain, double * restrict ux);
extern int fluid_update_boundaries_uy(const domain_t * restrict domain, double * restrict uy);
extern int fluid_update_boundaries_uz(const domain_t * restrict domain, double * restrict uz);
extern int fluid_update_boundaries_p (const domain_t * restrict domain, double * restrict  p);

// save flow field
extern int fluid_save(const char dirname[], const domain_t * restrict domain, const fluid_t * restrict fluid);

#endif // FLUID_H

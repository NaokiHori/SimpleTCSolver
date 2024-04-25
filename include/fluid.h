#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"

// definition of a structure fluid_t_
/**
 * @struct fluid_t
 * @brief struct storing fluid-related variables
 * @var u[x-z]      : velocity in each direction
 * @var p, psi      : pressure, scalar potential
 * @var t           : scalar field
 * @var srcu[x-z]   : Runge-Kutta source terms for velocities
 * @var l[x-z][x-z] : velocity-gradient tensor components
 * @var Re, Pr      : Reynolds and Prandtl number
 */
typedef struct {
  array_t ux;
  array_t uy;
#if NDIMS == 3
  array_t uz;
#endif
  array_t p;
  array_t psi;
  array_t t;
  array_t srcux[3];
  array_t srcuy[3];
#if NDIMS == 3
  array_t srcuz[3];
#endif
  array_t srct[3];
  array_t lxx;
  array_t lyx0, lyx1, lxy;
#if NDIMS == 3
  array_t lzx, lxz;
#endif
  array_t lyy0, lyy1;
#if NDIMS == 3
  array_t lzy, lyz;
  array_t lzz;
#endif
  double Re, Pr;
} fluid_t;

#endif // FLUID_H

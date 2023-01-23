#include <string.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


/**
 * @brief comute right-hand-side of Runge-Kutta scheme
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[inout] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return              : error code
 */
int fluid_compute_rhs(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
  compute_rhs_ux(domain, rkstep, fluid);
  compute_rhs_uy(domain, rkstep, fluid);
  compute_rhs_uz(domain, rkstep, fluid);
  return 0;
}


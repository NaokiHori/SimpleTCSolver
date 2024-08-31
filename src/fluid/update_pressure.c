#include "param.h"
#include "runge_kutta.h"
#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/psi.h"

#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++){
#define END \
      } \
    } \
  }

static inline int add_explicit(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // explicit contribution
    P(i, j, k) += PSI(i, j, k);
  END
  return 0;
}

static inline int add_implicit_x(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // x implicit contribution
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double l = jd_xm / hx_xm / hx_xm;
    const double u = jd_xp / hx_xp / hx_xp;
    const double c = - l - u;
    const double psi_xm = PSI(i-1, j  , k  );
    const double psi_x0 = PSI(i  , j  , k  );
    const double psi_xp = PSI(i+1, j  , k  );
    // NOTE: txx = "2" lxx
    P(i, j, k) -= 2. * prefactor / jd_x0 * (
        + l * psi_xm
        + c * psi_x0
        + u * psi_xp
    );
  END
  return 0;
}

static inline int add_implicit_y(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // y implicit contribution
    const double hy = HYXC(i  );
    const double jd = JDXC(i  );
    const double l = jd / hy / hy;
    const double u = jd / hy / hy;
    const double c = - l - u;
    const double psi_ym = PSI(i  , j-1, k  );
    const double psi_y0 = PSI(i  , j  , k  );
    const double psi_yp = PSI(i  , j+1, k  );
    // NOTE: tyy = "2" lyy
    P(i, j, k) -= 2. * prefactor / jd * (
        + l * psi_ym
        + c * psi_y0
        + u * psi_yp
    );
  END
  return 0;
}

static inline int add_implicit_z(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // z implicit contribution
    const double jd = JDXC(i  );
    const double l = jd / hz / hz;
    const double u = jd / hz / hz;
    const double c = - l - u;
    const double psi_zm = PSI(i  , j  , k-1);
    const double psi_z0 = PSI(i  , j  , k  );
    const double psi_zp = PSI(i  , j  , k+1);
    // NOTE: tzz = "2" lzz
    P(i, j, k) -= 2. * prefactor / jd * (
        + l * psi_zm
        + c * psi_z0
        + u * psi_zp
    );
  END
  return 0;
}

/**
 * @brief update pressure using scalar potential psi
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : scalar potential (in), pressure (out)
 * @return               : error code
 */
int fluid_update_pressure(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // explicit contribution, always present
  add_explicit(domain, fluid);
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid_compute_momentum_diffusivity(fluid);
  // additional corrections if diffusive terms
  //   in the direction is treated implicitly
  if(param_implicit_x){
    add_implicit_x(domain, prefactor, fluid);
  }
  if(param_implicit_y){
    add_implicit_y(domain, prefactor, fluid);
  }
  if(param_implicit_z){
    add_implicit_z(domain, prefactor, fluid);
  }
  // impose boundary conditions and communicate halo cells
  if(0 != fluid_update_boundaries_p(domain, &fluid->p)){
    return 1;
  }
  return 0;
}


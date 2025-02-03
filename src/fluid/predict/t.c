#include "param.h"
#include "memory.h"
#include "runge_kutta.h"
#include "linear_system.h"
#include "tdm.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/t.h"

static laplacians_t laplacians = {
  .is_initialised = false,
};

// [1 : isize]
#define LAPX(I) lapx[(I)-1]
#define LAPY(I) lapy[(I)-1]

static int init_laplacians(
    const domain_t * domain
){
  // laplacian x
  {
    const size_t isize = domain->glsizes[0];
    const double * hxxf = domain->hxxf;
    const double * jdxf = domain->jdxf;
    const double * jdxc = domain->jdxc;
    laplacians.lapx = memory_calloc(isize, sizeof(laplacian_t));
    // second-order derivative in x
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / JDXC(i  ) * JDXF(i  ) / HXXF(i  ) / HXXF(i  );
      const double u = 1. / JDXC(i  ) * JDXF(i+1) / HXXF(i+1) / HXXF(i+1);
      const double c = - l - u;
      laplacians.LAPX(i).l = l;
      laplacians.LAPX(i).c = c;
      laplacians.LAPX(i).u = u;
    }
  }
  // laplacian y
  {
    const size_t isize = domain->glsizes[0];
    const double * hyxc = domain->hyxc;
    laplacians.lapy = memory_calloc(isize, sizeof(laplacian_t));
    // second-order derivative in y
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / HYXC(i  ) / HYXC(i  );
      const double u = 1. / HYXC(i  ) / HYXC(i  );
      const double c = - l - u;
      laplacians.LAPY(i).l = l;
      laplacians.LAPY(i).c = c;
      laplacians.LAPY(i).u = u;
    }
  }
  laplacians.is_initialised = true;
  return 0;
}

#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++, cnt++){
#define END \
    } \
  }

static int advection_x(
    const domain_t * domain,
    const double * restrict t,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // t is advected in x
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double ux_xm = jd_xm / hx_xm * UX(i  , j  );
    const double ux_xp = jd_xp / hx_xp * UX(i+1, j  );
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
        + l * T(i-1, j  )
        + c * T(i  , j  )
        + u * T(i+1, j  )
    );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict t,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  // t is advected in y
  BEGIN
    const double jd = JDXC(i  );
    const double hy = HYXC(i  );
    const double uy_ym = jd / hy * UY(i  , j  );
    const double uy_yp = jd / hy * UY(i  , j+1);
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * T(i  , j-1)
        + c * T(i  , j  )
        + u * T(i  , j+1)
    );
  END
  return 0;
}

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * restrict lapx = laplacians.lapx;
  // diffusion in x, 0
  BEGIN
    src[cnt] += diffusivity * (
        + LAPX(i).l * T(i-1, j  )
        + LAPX(i).c * T(i  , j  )
        + LAPX(i).u * T(i+1, j  )
    );
  END
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * restrict lapy = laplacians.lapy;
  // diffusion in y, 0
  BEGIN
    src[cnt] += diffusivity * (
        + LAPY(i).l * T(i  , j-1)
        + LAPY(i).c * T(i  , j  )
        + LAPY(i).u * T(i  , j+1)
    );
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of scalar
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_t(
    const domain_t * domain,
    fluid_t * fluid
){
  if(!laplacians.is_initialised){
    if(0 != init_laplacians(domain)){
      return 1;
    }
  }
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  const double * restrict  t = fluid-> t.data;
  double * restrict srca = fluid->srct[rk_a].data;
  double * restrict srcg = fluid->srct[rk_g].data;
  const double diffusivity = fluid_compute_scalar_diffusivity(fluid);
  // advective contributions, always explicit
  advection_x(domain, t, ux, srca);
  advection_y(domain, t, uy, srca);
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, diffusivity, t, param_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, t, param_implicit_y ? srcg : srca);
  return 0;
}

/**
 * @brief update scalar
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int update_t(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  static linear_system_t linear_system = {
    .is_initialised = false,
  };
  if(!linear_system.is_initialised){
    // if not initialised yet, prepare linear solver
    //   for implicit diffusive term treatment
    const bool implicit[NDIMS] = {
      param_implicit_x,
      param_implicit_y,
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
    };
    if(0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)){
      return 1;
    }
  }
  // compute increments
  {
    const double coef_a = rkcoefs[rkstep][rk_a];
    const double coef_b = rkcoefs[rkstep][rk_b];
    const double coef_g = rkcoefs[rkstep][rk_g];
    const double * restrict srcta = fluid->srct[rk_a].data;
    const double * restrict srctb = fluid->srct[rk_b].data;
    const double * restrict srctg = fluid->srct[rk_g].data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * restrict delta = linear_system.x1pncl;
    const size_t nitems = isize * jsize;
    for(size_t n = 0; n < nitems; n++){
      delta[n] =
        + coef_a * dt * srcta[n]
        + coef_b * dt * srctb[n]
        + coef_g * dt * srctg[n];
    }
  }
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid_compute_scalar_diffusivity(fluid);
  // solve linear systems in x
  if(param_implicit_x){
    solve_in_x(
        prefactor,
        laplacians.lapx,
        &linear_system
    );
  }
  // solve linear systems in y
  if(param_implicit_y){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_y1,
        linear_system.x1pncl,
        linear_system.y1pncl
    );
    solve_in_y(
        // NOTE: tyy = "2" lyy
        2. * prefactor,
        laplacians.lapy,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_y1_to_x1,
        linear_system.y1pncl,
        linear_system.x1pncl
    );
  }
  // the field is actually updated here
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double * restrict delta = linear_system.x1pncl;
    double * restrict t = fluid->t.data;
    BEGIN
      T(i, j) += delta[cnt];
    END
    if(0 != fluid_update_boundaries_t(domain, &fluid->t)){
      return 1;
    }
  }
  return 0;
}


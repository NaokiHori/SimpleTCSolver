#if NDIMS == 3
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
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/lzx.h"
#include "array_macros/fluid/lzy.h"
#include "array_macros/fluid/p.h"

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
    // second-order derivative in x | 8
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
    // second-order derivative in y | 8
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / HYXC(i  ) / HYXC(i  );
      const double u = 1. / HYXC(i  ) / HYXC(i  );
      const double c = - l - u;
      laplacians.LAPY(i).l = l;
      laplacians.LAPY(i).c = c;
      laplacians.LAPY(i).u = u;
    }
  }
  // laplacian z
  {
    const double hz = domain->hz;
    // second-order derivative in z | 6
    const double l = 1. / hz / hz;
    const double u = 1. / hz / hz;
    const double c = - l - u;
    laplacians.lapz.l = l;
    laplacians.lapz.c = c;
    laplacians.lapz.u = u;
  }
  laplacians.is_initialised = true;
  return 0;
}

#define BEGIN \
  for(int cnt = 0, k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++, cnt++){
#define END \
      } \
    } \
  }

static int advection_x(
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // uz is advected in x | 19
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double ux_xm = + 0.5 * jd_xm / hx_xm * UX(i  , j  , k-1)
                         + 0.5 * jd_xm / hx_xm * UX(i  , j  , k  );
    const double ux_xp = + 0.5 * jd_xp / hx_xp * UX(i+1, j  , k-1)
                         + 0.5 * jd_xp / hx_xp * UX(i+1, j  , k  );
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
        + l * UZ(i-1, j  , k  )
        + c * UZ(i  , j  , k  )
        + u * UZ(i+1, j  , k  )
    );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  // uz is advected in y | 16
  BEGIN
    const double hy = HYXC(i  );
    const double jd = JDXC(i  );
    const double uy_ym = + 0.5 * jd / hy * UY(i  , j  , k-1)
                         + 0.5 * jd / hy * UY(i  , j  , k  );
    const double uy_yp = + 0.5 * jd / hy * UY(i  , j+1, k-1)
                         + 0.5 * jd / hy * UY(i  , j+1, k  );
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UZ(i  , j-1, k  )
        + c * UZ(i  , j  , k  )
        + u * UZ(i  , j+1, k  )
    );
  END
  return 0;
}

static int advection_z(
    const domain_t * domain,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  // uz is advected in z | 15
  BEGIN
    const double jd = JDXC(i  );
    const double uz_zm = + 0.5 * jd / hz * UZ(i  , j  , k-1)
                         + 0.5 * jd / hz * UZ(i  , j  , k  );
    const double uz_zp = + 0.5 * jd / hz * UZ(i  , j  , k  )
                         + 0.5 * jd / hz * UZ(i  , j  , k+1);
    const double l = - 0.5 * uz_zm;
    const double u = + 0.5 * uz_zp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UZ(i  , j  , k-1)
        + c * UZ(i  , j  , k  )
        + u * UZ(i  , j  , k+1)
    );
  END
  return 0;
}

static int diffusion_x0(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict lapx = laplacians.lapx;
  // diffusion in x, 0 | 7
  BEGIN
    src[cnt] += diffusivity * (
        + LAPX(i).l * UZ(i-1, j  , k  )
        + LAPX(i).c * UZ(i  , j  , k  )
        + LAPX(i).u * UZ(i+1, j  , k  )
    );
  END
  return 0;
}

static int diffusion_x1(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lzx,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  // diffusion in x, 1 | 13
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double lzx_xm = LZX(i  , j  , k  );
    const double lzx_xp = LZX(i+1, j  , k  );
    src[cnt] += diffusivity / jd_x0 * (
        - jd_xm / hx_xm * lzx_xm
        + jd_xp / hx_xp * lzx_xp
    );
  END
  return 0;
}

static int diffusion_y0(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict lapy = laplacians.lapy;
  // diffusion in y, 0 | 7
  BEGIN
    src[cnt] += diffusivity * (
        + LAPY(i).l * UZ(i  , j-1, k  )
        + LAPY(i).c * UZ(i  , j  , k  )
        + LAPY(i).u * UZ(i  , j+1, k  )
    );
  END
  return 0;
}

static int diffusion_y1(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict lzy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  // diffusion in y, 1 | 10
  BEGIN
    const double hy = HYXC(i  );
    const double jd = JDXC(i  );
    const double lzy_ym = LZY(i  , j  , k  );
    const double lzy_yp = LZY(i  , j+1, k  );
    src[cnt] += diffusivity / jd * (
        - jd / hy * lzy_ym
        + jd / hy * lzy_yp
    );
  END
  return 0;
}

static int diffusion_z(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict lapz = &laplacians.lapz;
  // diffusion in z | 8
  BEGIN
    // NOTE: tzz = "2" lzz
    src[cnt] += 2. * diffusivity * (
        + (*lapz).l * UZ(i  , j  , k-1)
        + (*lapz).c * UZ(i  , j  , k  )
        + (*lapz).u * UZ(i  , j  , k+1)
    );
  END
  return 0;
}

static int pressure(
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  // pressure gradient effect | 6
  BEGIN
    src[cnt] -= 1. / hz * (
        - P(i  , j  , k-1)
        + P(i  , j  , k  )
    );
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uz
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid
){
  if(!laplacians.is_initialised){
    if(0 != init_laplacians(domain)){
      return 1;
    }
  }
  const double * restrict  ux = fluid-> ux.data;
  const double * restrict  uy = fluid-> uy.data;
  const double * restrict  uz = fluid-> uz.data;
  const double * restrict   p = fluid->  p.data;
  const double * restrict lzx = fluid->lzx.data;
  const double * restrict lzy = fluid->lzy.data;
  double * restrict srca = fluid->srcuz[rk_a].data;
  double * restrict srcg = fluid->srcuz[rk_g].data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  // advective contributions, always explicit
  advection_x(domain, uz, ux, srca);
  advection_y(domain, uz, uy, srca);
  advection_z(domain, uz,     srca);
  // diagonal diffusive contributions, can be explicit or implicit
  diffusion_x0(domain, diffusivity, uz, param_implicit_x ? srcg : srca);
  diffusion_x1(domain, diffusivity, lzx,                          srca);
  diffusion_y0(domain, diffusivity, uz, param_implicit_y ? srcg : srca);
  diffusion_y1(domain, diffusivity, lzy,                          srca);
  diffusion_z (domain, diffusivity, uz, param_implicit_z ? srcg : srca);
  // pressure-gradient contribution, always implicit
  pressure(domain, p, srcg);
  return 0;
}

/**
 * @brief update uz
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int update_uz(
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
      param_implicit_z,
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
      domain->glsizes[2],
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
    const double * restrict srcuza = fluid->srcuz[rk_a].data;
    const double * restrict srcuzb = fluid->srcuz[rk_b].data;
    const double * restrict srcuzg = fluid->srcuz[rk_g].data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * restrict duz = linear_system.x1pncl;
    const size_t nitems = isize * jsize * ksize;
    for(size_t n = 0; n < nitems; n++){
      duz[n] =
        + coef_a * dt * srcuza[n]
        + coef_b * dt * srcuzb[n]
        + coef_g * dt * srcuzg[n];
    }
  }
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid_compute_momentum_diffusivity(fluid);
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
        prefactor,
        laplacians.lapy,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_y1_to_x1,
        linear_system.y1pncl,
        linear_system.x1pncl
    );
  }
  // solve linear systems in z
  if(param_implicit_z){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_z2,
        linear_system.x1pncl,
        linear_system.z2pncl
    );
    solve_in_z(
        // NOTE: tzz = "2" lzz
        2. * prefactor,
        &laplacians.lapz,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_z2_to_x1,
        linear_system.z2pncl,
        linear_system.x1pncl
    );
  }
  // the field is actually updated here
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double * restrict duz = linear_system.x1pncl;
    double * restrict uz = fluid->uz.data;
    BEGIN
      UZ(i, j, k) += duz[cnt];
    END
    if(0 != fluid_update_boundaries_uz(domain, &fluid->uz)){
      return 1;
    }
  }
  return 0;
}
#endif

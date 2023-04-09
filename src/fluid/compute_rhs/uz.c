#include <string.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "internal.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/uzlapx.h"
#include "arrays/domain/uzlapy.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"
#include "../arrays/srcuza.h"
#include "../arrays/srcuzb.h"
#include "../arrays/srcuzg.h"

static inline int advection_x(const domain_t * restrict domain, double * restrict srcuza, const double * restrict uz, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is transported by ux | 11
        // J = r dr dt dz
        const double jd = XC(i  ) * DXF(i  ) * dy * dz;
        const double cux_xm = 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k-1)
                            + 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  );
        const double cux_xp = 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k-1)
                            + 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k  );
        const double duz_xm = - UZ(i-1, j  , k  ) + UZ(i  , j  , k  );
        const double duz_xp = - UZ(i  , j  , k  ) + UZ(i+1, j  , k  );
        SRCUZA(i, j, k) +=
          - 1. / jd * 0.5 * cux_xm * duz_xm
          - 1. / jd * 0.5 * cux_xp * duz_xp;
      }
    }
  }
  return 0;
}

static inline int advection_y(const domain_t * restrict domain, double * restrict srcuza, const double * restrict uz, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double * restrict xc  = domain->xc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is transported by uy | 11
        // J = r dr dt dz
        const double jd = XC(i  ) * DXF(i  ) * dy * dz;
        const double cuy_ym = 0.5 * DXF(i  ) * dz * UY(i  , j  , k-1)
                            + 0.5 * DXF(i  ) * dz * UY(i  , j  , k  );
        const double cuy_yp = 0.5 * DXF(i  ) * dz * UY(i  , j+1, k-1)
                            + 0.5 * DXF(i  ) * dz * UY(i  , j+1, k  );
        const double duz_ym = - UZ(i  , j-1, k  ) + UZ(i  , j  , k  );
        const double duz_yp = - UZ(i  , j  , k  ) + UZ(i  , j+1, k  );
        SRCUZA(i, j, k) +=
          - 1. / jd * 0.5 * cuy_ym * duz_ym
          - 1. / jd * 0.5 * cuy_yp * duz_yp;
      }
    }
  }
  return 0;
}

static inline int advection_z(const domain_t * restrict domain, double * restrict srcuza, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double * restrict xc  = domain->xc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is transported by uz | 11
        // J = r dr dt dz
        const double jd = XC(i  ) * DXF(i  ) * dy * dz;
        const double cuz_zm = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k-1)
                            + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  );
        const double cuz_zp = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  )
                            + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k+1);
        const double duz_zm = - UZ(i  , j  , k-1) + UZ(i  , j  , k  );
        const double duz_zp = - UZ(i  , j  , k  ) + UZ(i  , j  , k+1);
        SRCUZA(i, j, k) +=
          - 1. / jd * 0.5 * cuz_zm * duz_zm
          - 1. / jd * 0.5 * cuz_zp * duz_zp;
      }
    }
  }
  return 0;
}

static inline int diffusion_x(const domain_t * restrict domain, double * restrict srcuza, const double diffusivity, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict uzlapx = domain->uzlapx;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is diffused in x | 5
        SRCUZA(i, j, k) += diffusivity * (
          + UZLAPX(i  ).l * UZ(i-1, j  , k  )
          + UZLAPX(i  ).c * UZ(i  , j  , k  )
          + UZLAPX(i  ).u * UZ(i+1, j  , k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_y(const domain_t * restrict domain, double * restrict srcuza, const double diffusivity, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict uzlapy = domain->uzlapy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is diffused in y | 5
        SRCUZA(i, j, k) += diffusivity * (
          + UZLAPY(i).l * UZ(i  , j-1, k  )
          + UZLAPY(i).c * UZ(i  , j  , k  )
          + UZLAPY(i).u * UZ(i  , j+1, k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_z(const domain_t * restrict domain, double * restrict srcuza, const double diffusivity, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t lap = domain->lapz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz is diffused in z | 5
        SRCUZA(i, j, k) += diffusivity * (
          + lap.l * UZ(i  , j  , k-1)
          + lap.c * UZ(i  , j  , k  )
          + lap.u * UZ(i  , j  , k+1)
        );
      }
    }
  }
  return 0;
}

static int pressure(const domain_t * restrict domain, double * restrict srcuzg, const double * restrict p){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        SRCUZG(i, j, k) -= 1. / dz * (
          - P(i  , j  , k-1)
          + P(i  , j  , k  )
        );
      }
    }
  }
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uz
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in,out] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return               : error code
 */
int fluid_compute_rhs_uz(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
  const int implicitx = config.get.implicitx();
  const int implicity = config.get.implicity();
  const int implicitz = config.get.implicitz();
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict ux = fluid->ux->data;
  const double * restrict uy = fluid->uy->data;
  const double * restrict uz = fluid->uz->data;
  const double * restrict  p = fluid-> p->data;
  double * restrict srcuza = fluid->srcuza->data;
  double * restrict srcuzb = fluid->srcuzb->data;
  double * restrict srcuzg = fluid->srcuzg->data;
  const double diffusivity = fluid->diffusivity;
  // copy previous k-step source term and reset
  if(rkstep != 0){
    memcpy(srcuzb, srcuza, (size_t)(SRCUZA_NITEMS_ALL) * sizeof(double));
  }
  memset(srcuza, 0, (size_t)(SRCUZA_NITEMS_ALL) * sizeof(double));
  memset(srcuzg, 0, (size_t)(SRCUZG_NITEMS_ALL) * sizeof(double));
  // advective contributions, always explicit
  advection_x(domain, srcuza, uz, ux);
  advection_y(domain, srcuza, uz, uy);
  advection_z(domain, srcuza, uz    );
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, implicitx ? srcuzg : srcuza, diffusivity, uz);
  diffusion_y(domain, implicity ? srcuzg : srcuza, diffusivity, uz);
  diffusion_z(domain, implicitz ? srcuzg : srcuza, diffusivity, uz);
  // pressure-gradient contribution, always implicit
  pressure(domain, srcuzg, p);
  return 0;
}


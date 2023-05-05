#include <string.h>
#include "param.h"
#include "domain.h"
#include "fluid.h"
#include "internal.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/domain/uylapx.h"
#include "arrays/domain/uylapy.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"
#include "../arrays/srcuya.h"
#include "../arrays/srcuyb.h"
#include "../arrays/srcuyg.h"

static inline int advection_x(const domain_t * restrict domain, double * restrict srcuya, const double * restrict uy, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxf = domain->dxf;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is transported by ux | 14
        const double pfactor = 0.5 / XC(i  ) / DXF(i  );
        const double l = + pfactor * (
            + 0.5 * XF(i  ) * UX(i  , j-1, k  )
            + 0.5 * XF(i  ) * UX(i  , j  , k  )
        );
        const double u = - pfactor * (
            + 0.5 * XF(i+1) * UX(i+1, j-1, k  )
            + 0.5 * XF(i+1) * UX(i+1, j  , k  )
        );
        const double c = - l - u;
        SRCUYA(i, j, k) +=
          + l * UY(i-1, j  , k  )
          + c * UY(i  , j  , k  )
          + u * UY(i+1, j  , k  );
      }
    }
  }
  return 0;
}

static inline int advection_y(const domain_t * restrict domain, double * restrict srcuya, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xc  = domain->xc;
  const double            dy  = domain->dy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is transported by uy | 14
        const double pfactor = 0.5 / XC(i  ) / dy;
        const double l = + pfactor * (
            + 0.5 * UY(i  , j-1, k  )
            + 0.5 * UY(i  , j  , k  )
        );
        const double u = - pfactor * (
            + 0.5 * UY(i  , j  , k  )
            + 0.5 * UY(i  , j+1, k  )
        );
        const double c = - l - u;
        SRCUYA(i, j, k) +=
          + l * UY(i  , j-1, k  )
          + c * UY(i  , j  , k  )
          + u * UY(i  , j+1, k  );
      }
    }
  }
  return 0;
}

static inline int advection_z(const domain_t * restrict domain, double * restrict srcuya, const double * restrict uy, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is transported by uz | 14
        const double pfactor = 0.5 / dz;
        const double l = + pfactor * (
            + 0.5 * UZ(i  , j-1, k  )
            + 0.5 * UZ(i  , j  , k  )
        );
        const double u = - pfactor * (
            + 0.5 * UZ(i  , j-1, k+1)
            + 0.5 * UZ(i  , j  , k+1)
        );
        const double c = - l - u;
        SRCUYA(i, j, k) +=
          + l * UY(i  , j  , k-1)
          + c * UY(i  , j  , k  )
          + u * UY(i  , j  , k+1);
      }
    }
  }
  return 0;
}

static inline int advection_c(const domain_t * restrict domain, double * restrict srcuya, const double * restrict uy, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xc  = domain->xc;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // Coriolis contribution | 13
        const double pfactor = 0.5 / XC(i  );
        const double uy_ym = 0.5 * UY(i  , j-1, k  )
                           + 0.5 * UY(i  , j  , k  );
        const double uy_yp = 0.5 * UY(i  , j  , k  )
                           + 0.5 * UY(i  , j+1, k  );
        const double ux_ym = 0.5 * UX(i  , j-1, k  )
                           + 0.5 * UX(i+1, j-1, k  );
        const double ux_yp = 0.5 * UX(i  , j  , k  )
                           + 0.5 * UX(i+1, j  , k  );
        SRCUYA(i, j, k) += - pfactor * (
            + uy_ym * ux_ym
            + uy_yp * ux_yp
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_x(const domain_t * restrict domain, double * restrict srcuya, const double diffusivity, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict uylapx = domain->uylapx;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is diffused in x | 5
        SRCUYA(i, j, k) += diffusivity * (
          + UYLAPX(i  ).l * UY(i-1, j  , k  )
          + UYLAPX(i  ).c * UY(i  , j  , k  )
          + UYLAPX(i  ).u * UY(i+1, j  , k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_y(const domain_t * restrict domain, double * restrict srcuya, const double diffusivity, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t *uylapy = domain->uylapy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is diffused in y | 5
        SRCUYA(i, j, k) += diffusivity * (
          + UYLAPY(i  ).l * UY(i  , j-1, k  )
          + UYLAPY(i  ).c * UY(i  , j  , k  )
          + UYLAPY(i  ).u * UY(i  , j+1, k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_z(const domain_t * restrict domain, double * restrict srcuya, const double diffusivity, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t lap = domain->lapz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy is diffused in z | 5
        SRCUYA(i, j, k) += diffusivity * (
          + lap.l * UY(i  , j  , k-1)
          + lap.c * UY(i  , j  , k  )
          + lap.u * UY(i  , j  , k+1)
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_c(const domain_t * restrict domain, double * restrict srcuya, const double diffusivity, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // diffusive effect from Christoffel symbol, 1 | 12
        {
          const double c_ym = 1. / XC(i  ) / DXF(i  );
          const double c_yp = 1. / XC(i  ) / DXF(i  );
          const double dxux_ym = - XF(i  ) * UX(i  , j-1, k  )
                                 + XF(i+1) * UX(i+1, j-1, k  );
          const double dxux_yp = - XF(i  ) * UX(i  , j  , k  )
                                 + XF(i+1) * UX(i+1, j  , k  );
          SRCUYA(i, j, k) += + diffusivity / XC(i  ) / dy * (
              - c_ym * dxux_ym
              + c_yp * dxux_yp
          );
        }
        // diffusive effect from Christoffel symbol, 2 | 12
        {
          const double c_ym = 1. / DXF(i  );
          const double c_yp = 1. / DXF(i  );
          const double duxx_ym = - 1. / XF(i  ) * UX(i  , j-1, k  )
                                 + 1. / XF(i+1) * UX(i+1, j-1, k  );
          const double duxx_yp = - 1. / XF(i  ) * UX(i  , j  , k  )
                                 + 1. / XF(i+1) * UX(i+1, j  , k  );
          SRCUYA(i, j, k) += - diffusivity / dy * (
              - c_ym * duxx_ym
              + c_yp * duxx_yp
          );
        }
      }
    }
  }
  return 0;
}

static int pressure(const domain_t * restrict domain, double * restrict srcuyg, const double * restrict p){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xc = domain->xc;
  const double            dy = domain->dy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        SRCUYG(i, j, k) +=
          - 1. / XC(i  ) / dy * (
            - P(i  , j-1, k  )
            + P(i  , j  , k  )
          );
      }
    }
  }
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in,out] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return               : error code
 */
int fluid_compute_rhs_uy(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict ux = fluid->ux->data;
  const double * restrict uy = fluid->uy->data;
  const double * restrict uz = fluid->uz->data;
  const double * restrict  p = fluid-> p->data;
  double * restrict srcuya = fluid->srcuya->data;
  double * restrict srcuyb = fluid->srcuyb->data;
  double * restrict srcuyg = fluid->srcuyg->data;
  const double diffusivity = fluid->diffusivity;
  // copy previous k-step source term and reset
  if(rkstep != 0){
    memcpy(srcuyb, srcuya, (size_t)(SRCUYA_NITEMS_ALL) * sizeof(double));
  }
  memset(srcuya, 0, (size_t)(SRCUYA_NITEMS_ALL) * sizeof(double));
  memset(srcuyg, 0, (size_t)(SRCUYG_NITEMS_ALL) * sizeof(double));
  // advective contributions, always explicit
  advection_x(domain, srcuya, uy, ux);
  advection_y(domain, srcuya, uy    );
  advection_z(domain, srcuya, uy, uz);
  advection_c(domain, srcuya, uy, ux);
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, param_implicit_x ? srcuyg : srcuya, diffusivity, uy);
  diffusion_y(domain, param_implicit_y ? srcuyg : srcuya, diffusivity, uy);
  diffusion_z(domain, param_implicit_z ? srcuyg : srcuya, diffusivity, uy);
  diffusion_c(domain, srcuya, diffusivity, ux);
  // pressure-gradient contribution, always implicit
  pressure(domain, srcuyg, p);
  return 0;
}


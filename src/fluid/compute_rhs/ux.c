#include <string.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "internal.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/domain/uxlapx.h"
#include "arrays/domain/uxlapy.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"
#include "../arrays/srcuxa.h"
#include "../arrays/srcuxb.h"
#include "../arrays/srcuxg.h"

static inline int advection_x(const domain_t * restrict domain, double * restrict srcuxa, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is transported by ux | 11
        // J = r dr dt dz
        const double jd = XF(i  ) * DXC(i  ) * dy * dz;
        const double cux_xm = 0.5 * XF(i-1) * dy * dz * UX(i-1, j  , k  )
                            + 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  );
        const double cux_xp = 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  )
                            + 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k  );
        const double dux_xm = - UX(i-1, j  , k  ) + UX(i  , j  , k  );
        const double dux_xp = - UX(i  , j  , k  ) + UX(i+1, j  , k  );
        SRCUXA(i, j, k) +=
          - 1. / jd * 0.5 * cux_xm * dux_xm
          - 1. / jd * 0.5 * cux_xp * dux_xp;
      }
    }
  }
  return 0;
}

static inline int advection_y(const domain_t * restrict domain, double * restrict srcuxa, const double * restrict ux, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is transported by uy | 11
        // J = r dr dt dz
        const double jd = XF(i  ) * DXC(i  ) * dy * dz;
        const double cuy_ym = 0.5 * DXF(i-1) * dz * UY(i-1, j  , k  )
                            + 0.5 * DXF(i  ) * dz * UY(i  , j  , k  );
        const double cuy_yp = 0.5 * DXF(i-1) * dz * UY(i-1, j+1, k  )
                            + 0.5 * DXF(i  ) * dz * UY(i  , j+1, k  );
        const double dux_ym = - UX(i  , j-1, k  ) + UX(i  , j  , k  );
        const double dux_yp = - UX(i  , j  , k  ) + UX(i  , j+1, k  );
        SRCUXA(i, j, k) +=
          - 1. / jd * 0.5 * cuy_ym * dux_ym
          - 1. / jd * 0.5 * cuy_yp * dux_yp;
      }
    }
  }
  return 0;
}

static inline int advection_z(const domain_t * restrict domain, double * restrict srcuxa, const double * restrict ux, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is transported by uz | 11
        // J = r dr dt dz
        const double jd = XF(i  ) * DXC(i  ) * dy * dz;
        const double cuz_zm = 0.5 * XC(i-1) * DXF(i-1) * dy * UZ(i-1, j  , k  )
                            + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  );
        const double cuz_zp = 0.5 * XC(i-1) * DXF(i-1) * dy * UZ(i-1, j  , k+1)
                            + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k+1);
        const double dux_zm = - UX(i  , j  , k-1) + UX(i  , j  , k  );
        const double dux_zp = - UX(i  , j  , k  ) + UX(i  , j  , k+1);
        SRCUXA(i, j, k) +=
          - 1. / jd * 0.5 * cuz_zm * dux_zm
          - 1. / jd * 0.5 * cuz_zp * dux_zp;
      }
    }
  }
  return 0;
}

static inline int advection_c(const domain_t * restrict domain, double * restrict srcuxa, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // centrifugal contribution | 11
        // J = r dr dt dz
        const double jd = XF(i  ) * DXC(i  ) * dy * dz;
        const double c_xm  = DXF(i-1) * dy * dz;
        const double c_xp  = DXF(i  ) * dy * dz;
        const double uy_xm = 0.5 * UY(i-1, j  , k  )
                           + 0.5 * UY(i-1, j+1, k  );
        const double uy_xp = 0.5 * UY(i  , j  , k  )
                           + 0.5 * UY(i  , j+1, k  );
        SRCUXA(i, j, k) +=
          + 1. / jd * 0.5 * c_xm * uy_xm * uy_xm
          + 1. / jd * 0.5 * c_xp * uy_xp * uy_xp;
      }
    }
  }
  return 0;
}

static inline int diffusion_x(const domain_t * restrict domain, double * restrict srcuxa, const double diffusivity, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict uxlapx = domain->uxlapx;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is diffused in x | 5
        SRCUXA(i, j, k) += diffusivity * (
          + UXLAPX(i  ).l * UX(i-1, j  , k  )
          + UXLAPX(i  ).c * UX(i  , j  , k  )
          + UXLAPX(i  ).u * UX(i+1, j  , k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_y(const domain_t * restrict domain, double * restrict srcuxa, const double diffusivity, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict uxlapy = domain->uxlapy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is diffused in y | 5
        SRCUXA(i, j, k) += diffusivity * (
          + UXLAPY(i  ).l * UX(i  , j-1, k  )
          + UXLAPY(i  ).c * UX(i  , j  , k  )
          + UXLAPY(i  ).u * UX(i  , j+1, k  )
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_z(const domain_t * restrict domain, double * restrict srcuxa, const double diffusivity, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t uxlapz = domain->lapz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // ux is diffused in z | 5
        SRCUXA(i, j, k) += diffusivity * (
          + uxlapz.l * UX(i  , j  , k-1)
          + uxlapz.c * UX(i  , j  , k  )
          + uxlapz.u * UX(i  , j  , k+1)
        );
      }
    }
  }
  return 0;
}

static inline int diffusion_c(const domain_t * restrict domain, double * restrict srcuxa, const double diffusivity, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxc = domain->dxc;
  const double            dy  = domain->dy;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // diffusive effect from Christoffel symbol, 1 | 12
        {
          const double c_ym = 1. / XF(i  ) / DXC(i  );
          const double c_yp = 1. / XF(i  ) / DXC(i  );
          const double dxuy_ym = - XC(i-1) * UY(i-1, j  , k  )
                                 + XC(i  ) * UY(i  , j  , k  );
          const double dxuy_yp = - XC(i-1) * UY(i-1, j+1, k  )
                                 + XC(i  ) * UY(i  , j+1, k  );
          SRCUXA(i, j, k) += - diffusivity / XF(i  ) / dy * (
              - c_ym * dxuy_ym
              + c_yp * dxuy_yp
          );
        }
        // diffusive effect from Christoffel symbol, 2 | 12
        {
          const double c_ym = 1. / DXC(i  );
          const double c_yp = 1. / DXC(i  );
          const double duyx_ym = - 1. / XC(i-1) * UY(i-1, j  , k  )
                                 + 1. / XC(i  ) * UY(i  , j  , k  );
          const double duyx_yp = - 1. / XC(i-1) * UY(i-1, j+1, k  )
                                 + 1. / XC(i  ) * UY(i  , j+1, k  );
          SRCUXA(i, j, k) += + diffusivity / dy * (
              - c_ym * duyx_ym
              + c_yp * duyx_yp
          );
        }
      }
    }
  }
  return 0;
}

static int pressure(const domain_t * restrict domain, double * restrict srcuxg, const double * restrict p){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxc = domain->dxc;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        SRCUXG(i, j, k) +=
          - 1. / DXC(i  ) * (
            - P(i-1, j  , k  )
            + P(i  , j  , k  )
          );
      }
    }
  }
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of ux
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in,out] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return               : error code
 */
int fluid_compute_rhs_ux(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
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
  double * restrict srcuxa = fluid->srcuxa->data;
  double * restrict srcuxb = fluid->srcuxb->data;
  double * restrict srcuxg = fluid->srcuxg->data;
  const double diffusivity = fluid->diffusivity;
  // copy previous k-step source term and reset
  if(rkstep != 0){
    memcpy(srcuxb, srcuxa, (size_t)(SRCUXA_NITEMS_ALL) * sizeof(double));
  }
  memset(srcuxa, 0, (size_t)(SRCUXA_NITEMS_ALL) * sizeof(double));
  memset(srcuxg, 0, (size_t)(SRCUXG_NITEMS_ALL) * sizeof(double));
  // advective contributions, always explicit
  advection_x(domain, srcuxa, ux    );
  advection_y(domain, srcuxa, ux, uy);
  advection_z(domain, srcuxa, ux, uz);
  advection_c(domain, srcuxa, uy    );
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, implicitx ? srcuxg : srcuxa, diffusivity, ux);
  diffusion_y(domain, implicity ? srcuxg : srcuxa, diffusivity, ux);
  diffusion_z(domain, implicitz ? srcuxg : srcuxa, diffusivity, ux);
  diffusion_c(domain, srcuxa, diffusivity, uy);
  // pressure-gradient contribution, always implicit
  pressure(domain, srcuxg, p);
  return 0;
}


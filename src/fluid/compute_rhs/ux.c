#include <string.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


/**
 * @brief comute right-hand-side of Runge-Kutta scheme of ux
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[inout] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return              : error code
 */
int compute_rhs_ux(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
  const int implicitx = config.get_bool("implicitx");
  const int implicity = config.get_bool("implicity");
  const int implicitz = config.get_bool("implicitz");
  const double Re = config.get_double("Re");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf = domain->xf;
  const double * restrict xc = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const laplace_t * restrict uxdifx = domain->uxdifx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict uz = fluid->uz;
  const double * restrict p = fluid->p;
  double * restrict srcuxa = fluid->srcuxa;
  double * restrict srcuxb = fluid->srcuxb;
  double * restrict srcuxg = fluid->srcuxg;
  /* previous k-step source term is copied */
  if(rkstep != 0){
    memcpy(srcuxb, srcuxa, SRCUXA_SIZE_0 * SRCUXA_SIZE_1 * SRCUXA_SIZE_2 * sizeof(double));
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        /* ! Jacobian determinant ! 2 ! */
        // J = r dr dt dz
        const double jd = XF(i  ) * DXC(i  ) * dy * dz;
        /* advection */
        /* ! ur transported by ur ! 13 ! */
        double advx;
        {
          const double cux_xm = 0.5 * XF(i-1) * dy * dz * UX(i-1, j  , k  )
                              + 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  );
          const double cux_xp = 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  )
                              + 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k  );
          const double dux_xm = - UX(i-1, j  , k  ) + UX(i  , j  , k  );
          const double dux_xp = - UX(i  , j  , k  ) + UX(i+1, j  , k  );
          advx = - 1. / jd * 0.5 * (
              + cux_xm * dux_xm
              + cux_xp * dux_xp
          );
        }
        /* ! ur transported by ut ! 13 ! */
        double advy;
        {
          const double cuy_ym = 0.5 * DXF(i-1) * dz * UY(i-1, j  , k  )
                              + 0.5 * DXF(i  ) * dz * UY(i  , j  , k  );
          const double cuy_yp = 0.5 * DXF(i-1) * dz * UY(i-1, j+1, k  )
                              + 0.5 * DXF(i  ) * dz * UY(i  , j+1, k  );
          const double dux_ym = - UX(i  , j-1, k  ) + UX(i  , j  , k  );
          const double dux_yp = - UX(i  , j  , k  ) + UX(i  , j+1, k  );
          advy = - 1. / jd * 0.5 * (
              + cuy_yp * dux_yp
              + cuy_ym * dux_ym
          );
        }
        /* ! ur transported by uz ! 13 ! */
        double advz;
        {
          const double cuz_zm = 0.5 * XC(i-1) * DXF(i-1) * dy * UZ(i-1, j  , k  )
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  );
          const double cuz_zp = 0.5 * XC(i-1) * DXF(i-1) * dy * UZ(i-1, j  , k+1)
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k+1);
          const double dux_zm = - UX(i  , j  , k-1) + UX(i  , j  , k  );
          const double dux_zp = - UX(i  , j  , k  ) + UX(i  , j  , k+1);
          advz = - 1. / jd * 0.5 * (
              + cuz_zp * dux_zp
              + cuz_zm * dux_zm
          );
        }
        /* ! centrifugal contribution ! 13 ! */
        double advc;
        {
          const double c_xm  = DXF(i-1) * dy * dz;
          const double c_xp  = DXF(i  ) * dy * dz;
          const double uy_xm = 0.5 * UY(i-1, j  , k  )
                             + 0.5 * UY(i-1, j+1, k  );
          const double uy_xp = 0.5 * UY(i  , j  , k  )
                             + 0.5 * UY(i  , j+1, k  );
          advc = + 1. / jd * 0.5 * (
              + c_xm * uy_xm * uy_xm
              + c_xp * uy_xp * uy_xp
          );
        }
        /* diffusion */
        /* ! diffused in r ! 8 ! */
        double difx;
        {
          difx = + 1. / Re * (
              + UXDIFX(i).l * UX(i-1, j  , k  )
              + UXDIFX(i).c * UX(i  , j  , k  )
              + UXDIFX(i).u * UX(i+1, j  , k  )
          );
        }
        /* ! diffused in t ! 8 ! */
        double dify;
        {
          dify = + 1. / Re / XF(i  ) / XF(i  ) / dy / dy * (
              + 1. * UX(i  , j-1, k  )
              - 2. * UX(i  , j  , k  )
              + 1. * UX(i  , j+1, k  )
          );
        }
        /* ! diffused in z ! 8 ! */
        double difz;
        {
          difz = + 1. / Re / dz / dz * (
              + 1. * UX(i  , j  , k-1)
              - 2. * UX(i  , j  , k  )
              + 1. * UX(i  , j  , k+1)
          );
        }
        /* ! additional diffusive term 0 ! 13 ! */
        double difa0;
        {
          const double c_ym = 1. / XF(i  ) / DXC(i  );
          const double c_yp = 1. / XF(i  ) / DXC(i  );
          const double dxuy_ym = - XC(i-1) * UY(i-1, j  , k  )
                                 + XC(i  ) * UY(i  , j  , k  );
          const double dxuy_yp = - XC(i-1) * UY(i-1, j+1, k  )
                                 + XC(i  ) * UY(i  , j+1, k  );
          difa0 = - 1. / Re / XF(i  ) / dy * (
              - c_ym * dxuy_ym
              + c_yp * dxuy_yp
          );
        }
        /* ! additional diffusive term 1 ! 13 ! */
        double difa1;
        {
          const double c_ym = 1. / DXC(i  );
          const double c_yp = 1. / DXC(i  );
          const double duyx_ym = - 1. / XC(i-1) * UY(i-1, j  , k  )
                                 + 1. / XC(i  ) * UY(i  , j  , k  );
          const double duyx_yp = - 1. / XC(i-1) * UY(i-1, j+1, k  )
                                 + 1. / XC(i  ) * UY(i  , j+1, k  );
          difa1 = + 1. / Re / dy * (
              - c_ym * duyx_ym
              + c_yp * duyx_yp
          );
        }
        double pre;
        {
          const double p_xm = P(i-1, j  , k  );
          const double p_xp = P(i  , j  , k  );
          pre = - 1. / DXC(i  ) * (
              + p_xp
              - p_xm
          );
        }
        /* sum of explicit terms */
        SRCUXA(i, j, k) =
          + advx
          + advy
          + advz
          + advc
          + (1. - 1. * implicitx) * difx
          + (1. - 1. * implicity) * dify
          + (1. - 1. * implicitz) * difz
          + difa0
          + difa1;
        /* sum of implicit terms */
        SRCUXG(i, j, k) =
          + (0. + 1. * implicitx) * difx
          + (0. + 1. * implicity) * dify
          + (0. + 1. * implicitz) * difz
          + pre;
      }
    }
  }
  return 0;
}


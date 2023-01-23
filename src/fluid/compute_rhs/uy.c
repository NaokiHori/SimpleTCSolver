#include <string.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[inout] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return              : error code
 */
int compute_rhs_uy(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
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
  const laplace_t * restrict uydifx = domain->uydifx;
  const laplace_t * restrict uydify = domain->uydify;
  const laplace_t            uydifz = domain->uydifz;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict uz = fluid->uz;
  const double * restrict p = fluid->p;
  double * restrict srcuya = fluid->srcuya;
  double * restrict srcuyb = fluid->srcuyb;
  double * restrict srcuyg = fluid->srcuyg;
  /* ! previous k-step source term is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srcuyb, srcuya, SRCUYA_SIZE_0 * SRCUYA_SIZE_1 * SRCUYA_SIZE_2 * sizeof(double));
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! Jacobian determinant ! 2 ! */
        // J = r dr dt dz
        const double jd = XC(i  ) * DXF(i  ) * dy * dz;
        /* advection */
        /* ! ut transported by ur ! 13 ! */
        double advx;
        {
          const double cux_xm = 0.5 * XF(i  ) * dy * dz * UX(i  , j-1, k  )
                              + 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  );
          const double cux_xp = 0.5 * XF(i+1) * dy * dz * UX(i+1, j-1, k  )
                              + 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k  );
          const double duy_xm = - UY(i-1, j  , k  ) + UY(i  , j  , k  );
          const double duy_xp = - UY(i  , j  , k  ) + UY(i+1, j  , k  );
          advx = - 1. / jd * 0.5 * (
              + cux_xm * duy_xm
              + cux_xp * duy_xp
          );
        }
        /* ! ut transported by ut ! 13 ! */
        double advy;
        {
          const double cuy_ym = 0.5 * DXF(i  ) * dz * UY(i  , j-1, k  )
                              + 0.5 * DXF(i  ) * dz * UY(i  , j  , k  );
          const double cuy_yp = 0.5 * DXF(i  ) * dz * UY(i  , j  , k  )
                              + 0.5 * DXF(i  ) * dz * UY(i  , j+1, k  );
          const double duy_ym = - UY(i  , j-1, k  ) + UY(i  , j  , k  );
          const double duy_yp = - UY(i  , j  , k  ) + UY(i  , j+1, k  );
          advy = - 1. / jd * 0.5 * (
              + cuy_ym * duy_ym
              + cuy_yp * duy_yp
          );
        }
        /* ! ut transported by uz ! 13 ! */
        double advz;
        {
          const double cuz_zm = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j-1, k  )
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  );
          const double cuz_zp = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j-1, k+1)
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k+1);
          const double duy_zm = - UY(i  , j  , k-1) + UY(i  , j  , k  );
          const double duy_zp = - UY(i  , j  , k  ) + UY(i  , j  , k+1);
          advz = - 1. / jd * 0.5 * (
              + cuz_zm * duy_zm
              + cuz_zp * duy_zp
          );
        }
        /* ! Coriolis contribution ! 17 ! */
        double advc;
        {
          const double c_ym = DXF(i  ) * dy * dz;
          const double c_yp = DXF(i  ) * dy * dz;
          const double uy_ym = 0.5 * UY(i  , j-1, k  )
                             + 0.5 * UY(i  , j  , k  );
          const double uy_yp = 0.5 * UY(i  , j  , k  )
                             + 0.5 * UY(i  , j+1, k  );
          const double ux_ym = 0.5 * UX(i  , j-1, k  )
                             + 0.5 * UX(i+1, j-1, k  );
          const double ux_yp = 0.5 * UX(i  , j  , k  )
                             + 0.5 * UX(i+1, j  , k  );
          advc = - 1. / jd * 0.5 * (
              + c_ym * uy_ym * ux_ym
              + c_yp * uy_yp * ux_yp
          );
        }
        /* diffusion */
        /* ! diffused in r ! 5 ! */
        const double difx = + 1. / Re * (
            + UYDIFX(i).l * UY(i-1, j  , k  )
            + UYDIFX(i).c * UY(i  , j  , k  )
            + UYDIFX(i).u * UY(i+1, j  , k  )
        );
        /* ! diffused in t ! 5 ! */
        const double dify = + 1. / Re * (
            + UYDIFY(i).l * UY(i  , j-1, k  )
            + UYDIFY(i).c * UY(i  , j  , k  )
            + UYDIFY(i).u * UY(i  , j+1, k  )
        );
        /* ! diffused in z ! 5 ! */
        const double difz = + 1. / Re * (
            + uydifz.l    * UY(i  , j  , k-1)
            + uydifz.c    * UY(i  , j  , k  )
            + uydifz.u    * UY(i  , j  , k+1)
        );
        /* ! additional diffusive term 0 ! 13 ! */
        double difa0;
        {
          const double c_ym = 1. / XC(i  ) / DXF(i  );
          const double c_yp = 1. / XC(i  ) / DXF(i  );
          const double dxux_ym = - XF(i  ) * UX(i  , j-1, k  )
                                 + XF(i+1) * UX(i+1, j-1, k  );
          const double dxux_yp = - XF(i  ) * UX(i  , j  , k  )
                                 + XF(i+1) * UX(i+1, j  , k  );
          difa0 = + 1. / Re / XC(i  ) / dy * (
              - c_ym * dxux_ym
              + c_yp * dxux_yp
          );
        }
        /* ! additional diffusive term 1 ! 13 ! */
        double difa1;
        {
          const double c_ym = 1. / DXF(i  );
          const double c_yp = 1. / DXF(i  );
          const double duxx_ym = - 1. / XF(i  ) * UX(i  , j-1, k  )
                                 + 1. / XF(i+1) * UX(i+1, j-1, k  );
          const double duxx_yp = - 1. / XF(i  ) * UX(i  , j  , k  )
                                 + 1. / XF(i+1) * UX(i+1, j  , k  );
          difa1 = - 1. / Re / dy * (
              - c_ym * duxx_ym
              + c_yp * duxx_yp
          );
        }
        double pre;
        {
          const double p_ym = P(i  , j-1, k  );
          const double p_yp = P(i  , j  , k  );
          pre = - 1. / XC(i  ) / dy * (
              - p_ym
              + p_yp
          );
        }
        /* sum of explicit terms */
        SRCUYA(i, j, k) =
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
        SRCUYG(i, j, k) =
          + (0. + 1. * implicitx) * difx
          + (0. + 1. * implicity) * dify
          + (0. + 1. * implicitz) * difz
          + pre;
      }
    }
  }
  return 0;
}


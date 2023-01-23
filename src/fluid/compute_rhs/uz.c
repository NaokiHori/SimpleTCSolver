#include <string.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "internal.h"


/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uz
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[inout] fluid  : velocity and pressure (in), RK source terms (inout)
 * @return              : error code
 */
int compute_rhs_uz(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid){
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
  const laplace_t * restrict uzdifx = domain->uzdifx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict ux = fluid->ux;
  const double * restrict uy = fluid->uy;
  const double * restrict uz = fluid->uz;
  const double * restrict p = fluid->p;
  double * restrict srcuza = fluid->srcuza;
  double * restrict srcuzb = fluid->srcuzb;
  double * restrict srcuzg = fluid->srcuzg;
  /* ! previous k-step source term is copied ! 3 ! */
  if(rkstep != 0){
    memcpy(srcuzb, srcuza, SRCUZA_SIZE_0 * SRCUZA_SIZE_1 * SRCUZA_SIZE_2 * sizeof(double));
  }
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        /* ! Jacobian determinant ! 2 ! */
        // J = r dr dt dz
        const double jd = XC(i  ) * DXF(i  ) * dy * dz;
        /* advection */
        /* ! uz transported by ur ! 13 ! */
        double advx;
        {
          const double cux_xm = 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k-1)
                              + 0.5 * XF(i  ) * dy * dz * UX(i  , j  , k  );
          const double cux_xp = 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k-1)
                              + 0.5 * XF(i+1) * dy * dz * UX(i+1, j  , k  );
          const double duz_xm = - UZ(i-1, j  , k  ) + UZ(i  , j  , k  );
          const double duz_xp = - UZ(i  , j  , k  ) + UZ(i+1, j  , k  );
          advx = - 0.5 / jd * (
              + cux_xm * duz_xm
              + cux_xp * duz_xp
          );
        }
        /* ! uz transported by ut ! 13 ! */
        double advy;
        {
          const double cuy_ym = 0.5 * DXF(i  ) * dz * UY(i  , j  , k-1)
                              + 0.5 * DXF(i  ) * dz * UY(i  , j  , k  );
          const double cuy_yp = 0.5 * DXF(i  ) * dz * UY(i  , j+1, k-1)
                              + 0.5 * DXF(i  ) * dz * UY(i  , j+1, k  );
          const double duz_ym = - UZ(i  , j-1, k  ) + UZ(i  , j  , k  );
          const double duz_yp = - UZ(i  , j  , k  ) + UZ(i  , j+1, k  );
          advy = - 0.5 / jd * (
              + cuy_ym * duz_ym
              + cuy_yp * duz_yp
          );
        }
        /* ! uz transported by uz ! 13 ! */
        double advz;
        {
          const double cuz_zm = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k-1)
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  );
          const double cuz_zp = 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k  )
                              + 0.5 * XC(i  ) * DXF(i  ) * dy * UZ(i  , j  , k+1);
          const double duz_zm = - UZ(i  , j  , k-1) + UZ(i  , j  , k  );
          const double duz_zp = - UZ(i  , j  , k  ) + UZ(i  , j  , k+1);
          advz = - 0.5 / jd * (
              + cuz_zm * duz_zm
              + cuz_zp * duz_zp
          );
        }
        /* diffusion */
        /* ! diffused in r ! 8 ! */
        double difx;
        {
          difx = + 1. / Re * (
              + UZDIFX(i).l * UZ(i-1, j  , k  )
              + UZDIFX(i).c * UZ(i  , j  , k  )
              + UZDIFX(i).u * UZ(i+1, j  , k  )
          );
        }
        /* ! diffused in t ! 8 ! */
        double dify;
        {
          dify = + 1. / Re / XC(i  ) / XC(i  ) / dy / dy * (
              + 1. * UZ(i  , j-1, k  )
              - 2. * UZ(i  , j  , k  )
              + 1. * UZ(i  , j+1, k  )
          );
        }
        /* ! diffused in z ! 8 ! */
        double difz;
        {
          difz = + 1. / Re / dz / dz * (
              + 1. * UZ(i  , j  , k-1)
              - 2. * UZ(i  , j  , k  )
              + 1. * UZ(i  , j  , k+1)
          );
        }
        double pre;
        {
          const double p_zm = P(i  , j  , k-1);
          const double p_zp = P(i  , j  , k  );
          pre = - 1. / dz * (
              + p_zp
              - p_zm
          );
        }
        /* sum of explicit terms */
        SRCUZA(i, j, k) =
          + advx
          + advy
          + advz
          + (1. - 1. * implicitx) * difx
          + (1. - 1. * implicity) * dify
          + (1. - 1. * implicitz) * difz;
        /* sum of implicit terms */
        SRCUZG(i, j, k) =
          + (0. + 1. * implicitx) * difx
          + (0. + 1. * implicity) * dify
          + (0. + 1. * implicitz) * difz
          + pre;
      }
    }
  }
  return 0;
}


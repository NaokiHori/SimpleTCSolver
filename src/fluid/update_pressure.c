#include <math.h>
#include "config.h"
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"



/**
 * @brief update pressure using scalar potential \psi
 * @param[in   ] domain : information related to MPI domain decomposition
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : scalar potential \psi (in), pressure (out)
 * @return              : error code
 */
int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict psi = fluid->psi;
  double * restrict p = fluid->p;
  /* ! add correction ! 7 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        P(i, j, k) += PSI(i, j, k);
      }
    }
  }
  /* ! correction for implicit treatment in x ! 23 ! */
  if(config.get_bool("implicitx")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const double * restrict xf = domain->xf;
    const double * restrict xc = domain->xc;
    const double * restrict dxf = domain->dxf;
    const double * restrict dxc = domain->dxc;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double c_xm = XF(i  ) / DXC(i  );
          const double c_xp = XF(i+1) / DXC(i+1);
          const double dpsi_xm = - PSI(i-1, j  , k  ) + PSI(i  , j  , k  );
          const double dpsi_xp = - PSI(i  , j  , k  ) + PSI(i+1, j  , k  );
          P(i, j, k) -= prefactor / XC(i  ) / DXF(i  ) * (
              - c_xm * dpsi_xm
              + c_xp * dpsi_xp
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in y ! 21 ! */
  if(config.get_bool("implicity")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const double * restrict xf = domain->xf;
    const double dy = domain->dy;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double c_ym = 1. / XF(i  ) / dy;
          const double c_yp = 1. / XF(i  ) / dy;
          const double dpsi_ym = - PSI(i  , j-1, k  ) + PSI(i  , j  , k  );
          const double dpsi_yp = - PSI(i  , j  , k  ) + PSI(i  , j+1, k  );
          P(i, j, k) -= prefactor / XF(i  ) / dy * (
              - c_ym * dpsi_ym
              + c_yp * dpsi_yp
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in z ! 20 ! */
  if(config.get_bool("implicitz")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const double dz = domain->dz;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double c_zm = 1. / dz;
          const double c_zp = 1. / dz;
          const double dpsi_zm = - PSI(i  , j  , k-1) + PSI(i  , j  , k  );
          const double dpsi_zp = - PSI(i  , j  , k  ) + PSI(i  , j  , k+1);
          P(i, j, k) -= prefactor / dz * (
              - c_zm * dpsi_zm
              + c_zp * dpsi_zp
          );
        }
      }
    }
  }
  /* ! boundary and halo values are updated ! 1 ! */
  fluid_update_boundaries_p(domain, p);
  return 0;
}


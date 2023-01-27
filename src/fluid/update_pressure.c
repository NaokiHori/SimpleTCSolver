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
  /* ! correction for implicit treatment in x ! 17 ! */
  if(config.get_bool("implicitx")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const laplace_t * restrict lappx = domain->lappx;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          P(i, j, k) -= prefactor * (
              + LAPPX(i).l * PSI(i-1, j  , k  )
              + LAPPX(i).c * PSI(i  , j  , k  )
              + LAPPX(i).u * PSI(i+1, j  , k  )
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in y ! 17 ! */
  if(config.get_bool("implicity")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const laplace_t * restrict lappy = domain->lappy;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          P(i, j, k) -= prefactor * (
              + LAPPY(i).l * PSI(i  , j-1, k  )
              + LAPPY(i).c * PSI(i  , j  , k  )
              + LAPPY(i).u * PSI(i  , j+1, k  )
          );
        }
      }
    }
  }
  /* ! correction for implicit treatment in z ! 17 ! */
  if(config.get_bool("implicitz")){
    const double Re = config.get_double("Re");
    const double gamma = RKCOEFS[rkstep].gamma;
    const double prefactor = (gamma * dt) / (2. * Re);
    const laplace_t lappz = domain->lappz;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          P(i, j, k) -= prefactor * (
              + lappz.l * PSI(i  , j  , k-1)
              + lappz.c * PSI(i  , j  , k  )
              + lappz.u * PSI(i  , j  , k+1)
          );
        }
      }
    }
  }
  /* ! boundary and halo values are updated ! 1 ! */
  fluid_update_boundaries_p(domain, p);
  return 0;
}


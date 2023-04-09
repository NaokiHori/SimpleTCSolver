#include "domain.h"
#include "fluid.h"
#include "arrays/domain/xc.h"
#include "arrays/fluid/uy.h"
#include "../arrays/psi.h"
#include "internal.h"

/**
 * @brief correct uy using scalar potential psi
 * @param[in]     domain    : information about domain decomposition and size
 * @param[in]     prefactor : pre-factor in front of grad psi
 * @param[in,out] fluid     : scalar potential psi (in), uy (out)
 * @return                  : error code
 */
int fluid_correct_velocity_uy(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xc = domain->xc;
  const double dy  = domain->dy;
  const double * restrict psi = fluid->psi->data;
  double * restrict uy = fluid->uy->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // correct y velocity 
        double psi_ym = PSI(i  , j-1, k  );
        double psi_yp = PSI(i  , j  , k  );
        UY(i, j, k) -= prefactor / XC(i  ) / dy * (
            + psi_yp
            - psi_ym
        );
      }
    }
  }
  // update boundary and halo cells 
  fluid_update_boundaries_uy(domain, uy);
  return 0;
}


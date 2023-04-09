#include "domain.h"
#include "fluid.h"
#include "arrays/domain/dxc.h"
#include "arrays/fluid/ux.h"
#include "../arrays/psi.h"
#include "internal.h"

/**
 * @brief correct ux using scalar potential psi
 * @param[in]     domain    : information about domain decomposition and size
 * @param[in]     prefactor : pre-factor in front of grad psi
 * @param[in,out] fluid     : scalar potential psi (in), ux (out)
 * @return                  : error code
 */
int fluid_correct_velocity_ux(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxc = domain->dxc;
  const double * restrict psi = fluid->psi->data;
  double * restrict ux = fluid->ux->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // correct x velocity 
        const double dx = DXC(i  );
        double psi_xm = PSI(i-1, j  , k  );
        double psi_xp = PSI(i  , j  , k  );
        UX(i, j, k) -= prefactor / dx * (
            + psi_xp
            - psi_xm
        );
      }
    }
  }
  // update boundary and halo cells 
  fluid_update_boundaries_ux(domain, ux);
  return 0;
}


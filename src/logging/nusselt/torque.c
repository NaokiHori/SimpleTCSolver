#include <stdio.h>
#include "param.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxc.h"
#include "arrays/fluid/uy.h"
#include "fileio.h"
#include "../internal.h"
#include "internal.h"

/**
 * @brief compute Nusselt number based on torque on the walls
 * @param[in] domain : information related to MPI domain decomposition
 * @return           : Nusselt number
 */
void logging_nusselt_compute_torque(const domain_t *domain, const fluid_t *fluid, double retvals[2]){
  const double diffusivity = fluid->diffusivity;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf = domain->xf;
  const double *xc = domain->xc;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *uy = fluid->uy->data;
  /* torque on the walls */
  // stress: diffusivity * r * d/dr ( ut / r )
  // torque = radius * stress
  retvals[0] = 0.;
  retvals[1] = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      // inner cylinder wall
      {
        // radius
        const double x = XF(1);
        // stress
        const int im = 0;
        const int ip = 1;
        const double stress = diffusivity * x * (
            - UY(im, j  , k  ) / XC(im)
            + UY(ip, j  , k  ) / XC(ip)
        ) / DXC(1);
        // surface element
        const double ds = x * dy * dz;
        retvals[0] -= x * stress * ds;
      }
      // outer cylinder wall
      {
        // radius
        const double x = XF(isize + 1);
        // stress
        const int im = isize + 0;
        const int ip = isize + 1;
        const double stress = diffusivity * x * (
            - UY(im, j  , k  ) / XC(im)
            + UY(ip, j  , k  ) / XC(ip)
        ) / DXC(isize + 1);
        // surface element
        const double ds = x * dy * dz;
        retvals[1] -= x * stress * ds;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, retvals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double ref = logging_nusselt_compute_torque_laminar(
      diffusivity,
      /* ri */ XF(      1),
      /* ro */ XF(isize+1),
      /* ui */ param_uy_xm,
      /* uo */ param_uy_xp,
      domain->lengths[1],
      domain->lengths[2]
  );
  retvals[0] /= ref;
  retvals[1] /= ref;
}


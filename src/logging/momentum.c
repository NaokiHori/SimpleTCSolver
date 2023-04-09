#include <stdio.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "internal.h"

/**
 * @brief compute total momenta
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_momentum(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf = domain->xf;
  const double * restrict xc = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict ux = fluid->ux->data;
  const double * restrict uy = fluid->uy->data;
  const double * restrict uz = fluid->uz->data;
  double moms[NDIMS] = {0.};
  // compute total x-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double cellsize = XF(i  ) * DXC(i  ) * dy * dz;
        moms[0] += UX(i, j, k) * cellsize;
      }
    }
  }
  // compute total y-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        moms[1] += UY(i, j, k) * cellsize;
      }
    }
  }
  // compute total z-momentum
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        moms[2] += UZ(i, j, k) * cellsize;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, moms, NDIMS, MPI_DOUBLE, MPI_SUM, comm_cart);
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    FILE *fp = fileio_fopen(fname, "a");
    if(NULL == fp) return 0;
    fprintf(fp, "%8.2f ", time);
    for(int n = 0; n < NDIMS; n++){
      fprintf(fp, "% 18.15e%c", moms[n], NDIMS - 1 == n ? '\n' : ' ');
    }
    fileio_fclose(fp);
  }
  return 0;
}


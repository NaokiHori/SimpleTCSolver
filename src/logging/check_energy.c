#include <stdio.h>
#include <math.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
#include "fileio.h"
#include "internal.h"



/**
 * @brief compute total kinetic and thermal energies
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int check_energy(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  /* compute kinetic and thermal energies */
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf = domain->xf;
  const double *xc = domain->xc;
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *uz = fluid->uz;
  // velocity in each dimension and thermal
  double quantities[NDIMS] = {0.};
  /* compute quadratic quantity in r direction */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        const double cellsize = XF(i  ) * DXC(i  ) * dy * dz;
        quantities[0] += 0.5 * pow(UX(i, j, k), 2.) * cellsize;
      }
    }
  }
  /* compute quadratic quantity in t direction */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        quantities[1] += 0.5 * pow(UY(i, j, k), 2.) * cellsize;
      }
    }
  }
  /* compute quadratic quantity in z direction */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        quantities[2] += 0.5 * pow(UZ(i, j, k), 2.) * cellsize;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, quantities, NDIMS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e % 18.15e\n",
          time,
          quantities[0],
          quantities[1],
          quantities[2],
          quantities[0] + quantities[1] + quantities[2]
      );
      fileio_fclose(fp);
    }
  }
  return 0;
}


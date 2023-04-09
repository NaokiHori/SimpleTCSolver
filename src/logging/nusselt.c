#include <stdio.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxc.h"
#include "arrays/fluid/uy.h"
#include "fileio.h"
#include "internal.h"

static double compute_torque_laminar(const double Re, const double ri, const double ro, const double ly, const double lz){
  // NOTE: assuming
  //   inner cylinder velocity is 1
  //   outer cylinder is at rest
  return 2. / Re * ly * lz * 1. / (1. / pow(ri, 2.) - 1. / pow(ro, 2.));
}

/**
 * @brief compute Nusselt number based on torque on the walls
 * @param[in] domain : information related to MPI domain decomposition
 * @return           : Nusselt number
 */
static void compute_nu_wall(const domain_t *domain, const fluid_t *fluid, double retvals[2]){
  const double Re = config.get.Re();
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
  // stress: 1 / Re * r * d/dr ( ut / r )
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
        const double stress = 1. / Re * x * (
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
        const double stress = 1. / Re * x * (
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
  const double ref = compute_torque_laminar(Re, XF(1), XF(isize+1), domain->lengths[1], domain->lengths[2]);
  retvals[0] /= ref;
  retvals[1] /= ref;
}

/**
 * @brief compute Nusselt number (normalised torque)
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_nusselt(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  double results[2] = {0.};
  compute_nu_wall(domain, fluid, results);
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e\n", time, results[0], results[1]);
      fileio_fclose(fp);
    }
  }
  return 0;
}


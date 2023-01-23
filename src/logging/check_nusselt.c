#include <stdio.h>
#include <math.h>
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"
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
static void compute_nu_wall(const domain_t *domain, const fluid_t *fluid, double *nu_walls){
  const double Re = config.get_double("Re");
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf = domain->xf;
  const double *xc = domain->xc;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double *uy = fluid->uy;
  /* torque on the walls */
  // stress: 1 / Re * r * d/dr ( ut / r )
  // torque = radius * stress
  nu_walls[0] = 0.;
  nu_walls[1] = 0.;
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
        nu_walls[0] -= x * stress * ds;
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
        nu_walls[1] -= x * stress * ds;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, nu_walls, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double ref = compute_torque_laminar(Re, XF(1), XF(isize+1), domain->lengths[1], domain->lengths[2]);
  nu_walls[0] /= ref;
  nu_walls[1] /= ref;
}

/**
 * @brief compute Nusselt number based on kinetic dissipation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : velocity
 * @return           : Nusselt number
 */
static double compute_nu_disk(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double lz = domain->lengths[2];
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double Re = config.get_double("Re");
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *uz = fluid->uz;
  double retval = 0.;
  /* ! x-momentum contribution ! 49 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        double cellsize = DXC(i)*dy*dz;
        // local dissipation at UX(i, j, k)
        double localval = 0.;
        {
          // x-negative contribution
          if(i != 1){
            double dux   = UX(i  , j  , k  )-UX(i-1, j  , k  );
            double duxdx = dux/DXF(i-1);
            localval += 0.5/DXC(i)*dux*duxdx;
          }
          // x-positive contribution
          if(i != isize+1){
            double dux   = UX(i+1, j  , k  )-UX(i  , j  , k  );
            double duxdx = dux/DXF(i  );
            localval += 0.5/DXC(i)*dux*duxdx;
          }
          // y-negative contribution
          {
            double dux   = UX(i  , j  , k  )-UX(i  , j-1, k  );
            double duxdy = dux/dy;
            localval += 0.5/dy*dux*duxdy;
          }
          // y-positive contribution
          {
            double dux   = UX(i  , j+1, k  )-UX(i  , j  , k  );
            double duxdy = dux/dy;
            localval += 0.5/dy*dux*duxdy;
          }
          // z-negative contribution
          {
            double dux   = UX(i  , j  , k  )-UX(i  , j  , k-1);
            double duxdz = dux/dz;
            localval += 0.5/dz*dux*duxdz;
          }
          // z-positive contribution
          {
            double dux   = UX(i  , j  , k+1)-UX(i  , j  , k  );
            double duxdz = dux/dz;
            localval += 0.5/dz*dux*duxdz;
          }
        }
        localval *= 1./Re;
        retval += Re*localval*cellsize;
      }
    }
  }
  /* ! y-momentum contribution ! 51 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        // local dissipation at UY(i, j, k)
        double localval = 0.;
        {
          // x-negative contribution
          {
            double c = i ==    1 ? 1. : 0.5;
            double duy   = UY(i  , j  , k  )-UY(i-1, j  , k  );
            double duydx = duy/DXC(i  );
            localval += c/DXF(i)*duy*duydx;
          }
          // x-positive contribution
          {
            double c = i == isize ? 1. : 0.5;
            double duy   = UY(i+1, j  , k  )-UY(i  , j  , k  );
            double duydx = duy/DXC(i+1);
            localval += c/DXF(i)*duy*duydx;
          }
          // y-negative contribution
          {
            double duy   = UY(i  , j  , k  )-UY(i  , j-1, k  );
            double duydy = duy/dy;
            localval += 0.5/dy*duy*duydy;
          }
          // y-positive contribution
          {
            double duy   = UY(i  , j+1, k  )-UY(i  , j  , k  );
            double duydy = duy/dy;
            localval += 0.5/dy*duy*duydy;
          }
          // z-negative contribution
          {
            double duy   = UY(i  , j  , k  )-UY(i  , j  , k-1);
            double duydz = duy/dz;
            localval += 0.5/dz*duy*duydz;
          }
          // z-positive contribution
          {
            double duy   = UY(i  , j  , k+1)-UY(i  , j  , k  );
            double duydz = duy/dz;
            localval += 0.5/dz*duy*duydz;
          }
        }
        localval *= 1./Re;
        retval += Re*localval*cellsize;
      }
    }
  }
  /* ! z-momentum contribution ! 51 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double cellsize = DXF(i)*dy*dz;
        // local dissipation at UZ(i, j, k)
        double localval = 0.;
        {
          // x-negative contribution
          {
            double c = i ==    1 ? 1. : 0.5;
            double duz   = UZ(i  , j  , k  )-UZ(i-1, j  , k  );
            double duzdx = duz/DXC(i  );
            localval += c/DXF(i)*duz*duzdx;
          }
          // x-positive contribution
          {
            double c = i == isize ? 1. : 0.5;
            double duz   = UZ(i+1, j  , k  )-UZ(i  , j  , k  );
            double duzdx = duz/DXC(i+1);
            localval += c/DXF(i)*duz*duzdx;
          }
          // y-negative contribution
          {
            double duz   = UZ(i  , j  , k  )-UZ(i  , j-1, k  );
            double duzdy = duz/dy;
            localval += 0.5/dy*duz*duzdy;
          }
          // y-positive contribution
          {
            double duz   = UZ(i  , j+1, k  )-UZ(i  , j  , k  );
            double duzdy = duz/dy;
            localval += 0.5/dy*duz*duzdy;
          }
          // z-negative contribution
          {
            double duz   = UZ(i  , j  , k  )-UZ(i  , j  , k-1);
            double duzdz = duz/dz;
            localval += 0.5/dz*duz*duzdz;
          }
          // z-positive contribution
          {
            double duz   = UZ(i  , j  , k+1)-UZ(i  , j  , k  );
            double duzdz = duz/dz;
            localval += 0.5/dz*duz*duzdz;
          }
        }
        localval *= 1./Re;
        retval += Re*localval*cellsize;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // average in the whole domain
  retval /= lx*ly*lz;
  retval += 1.;
  return retval;
}

/**
 * @brief compute Nusselt number based on various definitions
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int check_nusselt(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  /* compute Nusselt number based on several definintions */
  // compute Nu from torque on the walls
  double nu_walls[2] = {0.};
  compute_nu_wall(domain, fluid, nu_walls);
  // compute Nu from kinetic energy dissipation rate
  double nu_eps = compute_nu_disk(domain, fluid);
  int myrank;
  MPI_Comm_rank(domain->sdecomp->comm_cart, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, nu_walls[0], nu_walls[1], nu_eps);
      fileio_fclose(fp);
    }
  }
  return 0;
}


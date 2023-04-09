#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "internal.h"

/**
 * @brief check divergence and write the maximum value
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : domain information
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_divergence(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const double * restrict ux = fluid->ux->data;
  const double * restrict uy = fluid->uy->data;
  const double * restrict uz = fluid->uz->data;
  double divmax = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // compute local divergence
        const double dx    = DXF(i  );
        const double x     = XC(i  );
        const double x_xm  = XF(i  );
        const double x_xp  = XF(i+1);
        const double ux_xm = UX(i  , j  , k  );
        const double ux_xp = UX(i+1, j  , k  );
        const double uy_ym = UY(i  , j  , k  );
        const double uy_yp = UY(i  , j+1, k  );
        const double uz_zm = UZ(i  , j  , k  );
        const double uz_zp = UZ(i  , j  , k+1);
        const double div =
           + (x_xp * ux_xp - x_xm * ux_xm) / x / dx
           + (       uy_yp -        uy_ym) / x / dy
           + (       uz_zp -        uz_zm)     / dz;
        // check maximum
        divmax = fmax(divmax, fabs(div));
      }
    }
  }
  // collect information among all processes
  MPI_Allreduce(MPI_IN_PLACE, &divmax, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
  // result is written to a file from the main process
  int myrank = 0;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    FILE *fp = fileio_fopen(fname, "a");
    if(NULL == fp) return 0;
    fprintf(fp, "%8.2f % .1e\n", time, divmax);
    fileio_fclose(fp);
  }
  return 0;
}


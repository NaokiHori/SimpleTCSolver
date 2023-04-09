#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "fileio.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/domain/uxlapx.h"
#include "arrays/domain/uylapx.h"
#include "arrays/domain/uzlapx.h"
#include "arrays/domain/plapx.h"
#include "arrays/domain/uxlapy.h"
#include "arrays/domain/uylapy.h"
#include "arrays/domain/uzlapy.h"
#include "arrays/domain/plapy.h"
#include "internal.h"

/**
 * @brief define face-to-face distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return          : face-to-face distances in x direction
 */
static double *allocate_and_init_dxf(const int isize, const double *xf){
  // dxf: distance from cell face to cell face
  double *dxf = common_calloc((size_t)(DXF_NITEMS_0(isize)), sizeof(double));
  for(int i = 1; i <= isize; i++){
    DXF(i  ) = XF(i+1) - XF(i  );
  }
  return dxf;
}

/**
 * @brief define center-to-center distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xc    : cell-center positions in x direction
 * @return          : center-to-center distances in x direction
 */
static double *allocate_and_init_dxc(const int isize, const double *xc){
  // dxc: distance from cell center to cell center (generally)
  double *dxc = common_calloc((size_t)(DXC_NITEMS_0(isize)), sizeof(double));
  for(int i = 1; i <= isize + 1; i++){
    DXC(i  ) = XC(i  ) - XC(i-1);
  }
  return dxc;
}

/**
 * @brief compute discrete Laplacian in each direction for each scalar
 * @param[in,out] domain : grid positions (in), discrete Laplacian (out)
 * @return               : error code
 */
static int compute_laplacian(domain_t *domain){
  // Laplacian w.r.t. ux in x | 15
  {
    const size_t isize = domain->glsizes[0];
    const double *xf  = domain->xf;
    const double *xc  = domain->xc;
    const double *dxf = domain->dxf;
    const double *dxc = domain->dxc;
    laplacian_t *uxlapx = common_calloc(isize - 1, sizeof(laplacian_t));
    for(size_t i = 2; i <= isize; i++){
      UXLAPX(i).l = + XF(i-1) / DXF(i-1) / XC(i-1) / DXC(i  );
      UXLAPX(i).u = + XF(i+1) / DXF(i  ) / XC(i  ) / DXC(i  );
      UXLAPX(i).c = - XF(i  ) / DXF(i-1) / XC(i-1) / DXC(i  )
                    - XF(i  ) / DXF(i  ) / XC(i  ) / DXC(i  );
    }
    domain->uxlapx = uxlapx;
  }
  // Laplacian w.r.t. uy in x | 15
  {
    const size_t isize = domain->glsizes[0];
    const double *xf  = domain->xf;
    const double *xc  = domain->xc;
    const double *dxf = domain->dxf;
    const double *dxc = domain->dxc;
    laplacian_t *uylapx = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      UYLAPX(i).l = + 1. / XC(i-1) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
      UYLAPX(i).u = + 1. / XC(i+1) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
      UYLAPX(i).c = - 1. / XC(i  ) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.)
                    - 1. / XC(i  ) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
    }
    domain->uylapx = uylapx;
  }
  // Laplacian w.r.t. uz in x | 15
  {
    const size_t isize = domain->glsizes[0];
    const double *xf  = domain->xf;
    const double *xc  = domain->xc;
    const double *dxf = domain->dxf;
    const double *dxc = domain->dxc;
    laplacian_t *uzlapx = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      UZLAPX(i).l = + 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  );
      UZLAPX(i).u = + 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
      UZLAPX(i).c = - 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  )
                    - 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
    }
    domain->uzlapx = uzlapx;
  }
  // Laplacian w.r.t. p in x | 15
  {
    const size_t isize = domain->glsizes[0];
    const double *xf  = domain->xf;
    const double *xc  = domain->xc;
    const double *dxf = domain->dxf;
    const double *dxc = domain->dxc;
    laplacian_t *plapx = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      PLAPX(i).l = + 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  );
      PLAPX(i).u = + 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
      PLAPX(i).c = - 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  )
                   - 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
    }
    domain->plapx = plapx;
  }
  // Laplacian w.r.t. ux in y | 12
  {
    const size_t isize = domain->glsizes[0];
    const double *xf = domain->xf;
    const double dy  = domain->dy;
    laplacian_t *uxlapy = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 2; i <= isize; i++){
      UXLAPY(i).l = + 1. / XF(i  ) / XF(i  ) / dy / dy;
      UXLAPY(i).u = + 1. / XF(i  ) / XF(i  ) / dy / dy;
      UXLAPY(i).c = - 2. / XF(i  ) / XF(i  ) / dy / dy;
    }
    domain->uxlapy = uxlapy;
  }
  // Laplacian w.r.t. uy in y | 12
  {
    const size_t isize = domain->glsizes[0];
    const double *xc = domain->xc;
    const double dy  = domain->dy;
    laplacian_t *uylapy = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      UYLAPY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      UYLAPY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      UYLAPY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
    }
    domain->uylapy = uylapy;
  }
  // Laplacian w.r.t. uz in y | 12
  {
    const size_t isize = domain->glsizes[0];
    const double *xc = domain->xc;
    const double dy  = domain->dy;
    laplacian_t *uzlapy = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      UZLAPY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      UZLAPY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      UZLAPY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
    }
    domain->uzlapy = uzlapy;
  }
  // Laplacian w.r.t. p in y | 12
  {
    const size_t isize = domain->glsizes[0];
    const double *xc = domain->xc;
    const double dy  = domain->dy;
    laplacian_t *plapy = common_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      PLAPY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      PLAPY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
      PLAPY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
    }
    domain->plapy = plapy;
  }
  // Laplacian in z | 7
  {
    laplacian_t *lapz = &domain->lapz;
    const double dz = domain->dz;
    (*lapz).l = + 1. / dz / dz;
    (*lapz).u = + 1. / dz / dz;
    (*lapz).c = - 2. / dz / dz;
  }
  return 0;
}

/**
 * @brief constructor of the structure
 * @param[in] dirname_ic : name of directory in which initial conditions are stored
 * @return               : structure being allocated and initalised
 */
domain_t *domain_init(const char dirname_ic[]){
  domain_t *domain      = common_calloc(1, sizeof(domain_t));
  int    **glsizes      = &(domain->glsizes);
  int    **mysizes      = &(domain->mysizes);
  int    **offsets      = &(domain->offsets);
  double **lengths      = &(domain->lengths);
  sdecomp_info_t **info = &(domain->info);
  double **xf           = &(domain->xf);
  double **xc           = &(domain->xc);
  double **dxf          = &(domain->dxf);
  double **dxc          = &(domain->dxc);
  double *dy            = &(domain->dy);
  double *dz            = &(domain->dz);
  // load spatial information
  if(0 != domain_load(dirname_ic, domain)){
    return NULL;
  }
  // compute grid sizes
  // allocate and initialise x coordinates
  *dxf = allocate_and_init_dxf((*glsizes)[0], *xf);
  *dxc = allocate_and_init_dxc((*glsizes)[0], *xc);
  // grid sizes in homogeneous directions
  *dy = (*lengths)[1] / (*glsizes)[1];
  *dz = (*lengths)[2] / (*glsizes)[2];
  // compute discrete Laplace operators
  compute_laplacian(domain);
  // initialise sdecomp to distribute the domain
  *info = optimise_sdecomp_init(domain->uniformx, *glsizes);
  // local array sizes and offsets
  *mysizes = common_calloc(NDIMS, sizeof(int));
  *offsets = common_calloc(NDIMS, sizeof(int));
  for(int dim = 0; dim < NDIMS; dim++){
    size_t mysize = 0;
    size_t offset = 0;
    sdecomp.get_pencil_mysize(*info, SDECOMP_X1PENCIL, dim, (size_t)(*glsizes)[dim], &mysize);
    sdecomp.get_pencil_offset(*info, SDECOMP_X1PENCIL, dim, (size_t)(*glsizes)[dim], &offset);
    (*mysizes)[dim] = (int)mysize;
    (*offsets)[dim] = (int)offset;
  }
  return domain;
}


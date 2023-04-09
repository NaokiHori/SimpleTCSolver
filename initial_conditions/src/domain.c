#include "common.h"
#include "config.h"
#include "domain.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"

/**
 * @brief define cell-face positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] lx    : physical domain size in x direction
 *                      normally 1 since this is the reference length scale
 *                      to non-dimensionalise equations
 * @return : cell-face positions in x direction
 */
static double *allocate_and_init_xf(const int isize, const double lx, const bool uniformx){
  /* ! xf: cell face coordinates ! 28 ! */
  double *xf = common_calloc(XF_NITEMS_0(isize), sizeof(double));
  if(uniformx){
    // uniform grid, which enables us to use
    //   efficient DCT-based Poisson solver
    //   see src/fluid/compute_potential.c
    const double dx = lx / isize;
    for(int i = 1; i <= isize+1; i++){
      XF(i) = 1. * (i - 1) * dx;
    }
  }else{
    // stretched grid, clipped Chebyshev just as an example
    // number of grid points to be clipped at the edges
    const int nclip = 3;
    // changing from 1 to -1, dense in the vicinity of the edges
    for(int i = 1; i <= isize + 1; i++){
      XF(i) = cos(M_PI * (1. * (i - 1) + nclip) / (isize + 2. * nclip));
    }
    // normalse to force it changing from 0 to lx
    double min = XF(      1);
    double max = XF(isize+1);
    for(int i = 1; i <= isize+1; i++){
      XF(i) = lx * (XF(i) - min) / (max - min);
    }
  }
  // add 1 to make it from 1 to 1+lx
  for(int i = 1; i <= isize + 1; i++){
    XF(i) += 1.;
  }
  return xf;
}

/**
 * @brief define cell-center positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : cell-center positions in x direction
 */
static double *allocate_and_init_xc(const int isize, const double *xf){
  /* ! xc: cell center coordinates ! 8 ! */
  double *xc = common_calloc(XC_NITEMS_0(isize), sizeof(double));
  // center between two XFs
  for(int i = 1; i <= isize; i++){
    XC(i) = 0.5 * (XF(i  ) + XF(i+1));
  }
  // at boundaries, face positions are assigned
  XC(      0) = XF(      1);
  XC(isize+1) = XF(isize+1);
  return xc;
}

domain_t *domain_init(void){
  domain_t *domain = common_calloc(1, sizeof(domain_t));
  /* init */
  int *glisize = &(domain->glisize);
  int *gljsize = &(domain->gljsize);
  int *glksize = &(domain->glksize);
  double *lx = &(domain->lx);
  double *ly = &(domain->ly);
  double *lz = &(domain->lz);
  double **xf = &(domain->xf);
  double **xc = &(domain->xc);
  double *dy  = &(domain->dy);
  double *dz = &(domain->dz);
  const bool uniformx = config.get_bool("uniformx");
  *glisize = config.get_int("glisize");
  *gljsize = config.get_int("gljsize");
  *glksize = config.get_int("glksize");
  *lx = config.get_double("lx");
  *ly = config.get_double("ly");
  *lz = config.get_double("lz");
  *xf = allocate_and_init_xf(*glisize, *lx, uniformx);
  *xc = allocate_and_init_xc(*glisize, *xf);
  *dy = *ly / *gljsize;
  *dz = *lz / *glksize;
  /* save */
  common_save("uniformx.npy", 0, NULL, NPY_BOL, sizeof(bool), &uniformx);
  common_save("glisize.npy",  0, NULL, NPY_INT, sizeof(int),    glisize);
  common_save("gljsize.npy",  0, NULL, NPY_INT, sizeof(int),    gljsize);
  common_save("glksize.npy",  0, NULL, NPY_INT, sizeof(int),    glksize);
  common_save("lx.npy",       0, NULL, NPY_DBL, sizeof(double), lx     );
  common_save("ly.npy",       0, NULL, NPY_DBL, sizeof(double), ly     );
  common_save("lz.npy",       0, NULL, NPY_DBL, sizeof(double), lz     );
  size_t shape[1] = {0};
  shape[0] = *glisize + 1;
  common_save("xf.npy", 1, shape, NPY_DBL, sizeof(double) * shape[0], *xf);
  shape[0] = *glisize + 2;
  common_save("xc.npy", 1, shape, NPY_DBL, sizeof(double) * shape[0], *xc);
  return domain;
}

/**
 * @brief destruct a structure domain_t
 * @param[inout] domain : structure to be cleaned-up
 * @return              : error code
 */
int domain_finalise(domain_t *domain){
  common_free(domain->xf);
  common_free(domain->xc);
  common_free(domain);
  return 0;
}


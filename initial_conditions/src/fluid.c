#include <stdlib.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "snpyio.h"
#include "arrays/domain/xc.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"

static size_t get_nitems(const size_t *shape){
  size_t nitems = 1;
  for(size_t dim = 0; dim < NDIMS; dim++){
    nitems *= shape[dim];
  }
  return nitems;
}

static double myrandom(void){
  static const double fluc = 1.e-4;
  return fluc * (- 0.5 + 1. * rand() / RAND_MAX);
}

static void init_and_save_ux(const domain_t *domain){
  const int isize = domain->glisize;
  const int jsize = domain->gljsize;
  const int ksize = domain->glksize;
  double *ux = common_calloc(UX_NITEMS_ALL, sizeof(double));
  // bulk values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        UX(i, j, k) = myrandom();
      }
    }
  }
  // boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UX(      1, j, k) = 0.;
      UX(isize+1, j, k) = 0.;
    }
  }
  const size_t shape[NDIMS] = {UX_NITEMS_2(ksize), UX_NITEMS_1(jsize), UX_NITEMS_0(isize)};
  const size_t nitems = get_nitems(shape);
  common_save("ux.npy", NDIMS, shape, NPY_DBL, sizeof(double) * nitems, ux);
  common_free(ux);
}

static void init_and_save_uy(const domain_t *domain){
  const int isize = domain->glisize;
  const int jsize = domain->gljsize;
  const int ksize = domain->glksize;
  const double *xc = domain->xc;
  double *uy = common_calloc(UY_NITEMS_ALL, sizeof(double));
  // bulk values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UY(i, j, k) = 1. + (0. - 1.) / (XC(isize+1) - XC(0)) * (XC(i) - XC(0)) + myrandom();
      }
    }
  }
  // boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UY(      0, j, k) = 1.;
      UY(isize+1, j, k) = 0.;
    }
  }
  const size_t shape[NDIMS] = {UY_NITEMS_2(ksize), UY_NITEMS_1(jsize), UY_NITEMS_0(isize)};
  const size_t nitems = get_nitems(shape);
  common_save("uy.npy", NDIMS, shape, NPY_DBL, sizeof(double) * nitems, uy);
  common_free(uy);
}

static void init_and_save_uz(const domain_t *domain){
  const int isize = domain->glisize;
  const int jsize = domain->gljsize;
  const int ksize = domain->glksize;
  double *uz = common_calloc(UZ_NITEMS_ALL, sizeof(double));
  // bulk values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UZ(i, j, k) = myrandom();
      }
    }
  }
  // boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UZ(      0, j, k) = 0.;
      UZ(isize+1, j, k) = 0.;
    }
  }
  const size_t shape[NDIMS] = {UZ_NITEMS_2(ksize), UZ_NITEMS_1(jsize), UZ_NITEMS_0(isize)};
  const size_t nitems = get_nitems(shape);
  common_save("uz.npy", NDIMS, shape, NPY_DBL, sizeof(double) * nitems, uz);
  common_free(uz);
}

static void init_and_save_p(const domain_t *domain){
  const int isize = domain->glisize;
  const int jsize = domain->gljsize;
  const int ksize = domain->glksize;
  double *p = common_calloc(P_NITEMS_ALL, sizeof(double));
  // bulk values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        P(i, j, k) = 0.;
      }
    }
  }
  // boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      P(      0, j, k) = P(    1, j, k);
      P(isize+1, j, k) = P(isize, j, k);
    }
  }
  const size_t shape[NDIMS] = {P_NITEMS_2(ksize), P_NITEMS_1(jsize), P_NITEMS_0(isize)};
  const size_t nitems = get_nitems(shape);
  common_save("p.npy", NDIMS, shape, NPY_DBL, sizeof(double) * nitems, p);
  common_free(p);
}

void fluid_init(const domain_t * restrict domain){
  init_and_save_ux(domain);
  init_and_save_uy(domain);
  init_and_save_uz(domain);
  init_and_save_p(domain);
}


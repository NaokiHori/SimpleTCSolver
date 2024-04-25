#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/lyx0.h"
#include "array_macros/fluid/lyx1.h"
#include "internal.h"

#if NDIMS == 2
#define BEGIN \
  for(int j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize + 1; i++){
#define END \
    } \
  }
#else
#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize + 1; i++){
#define END \
      } \
    } \
  }
#endif

static int compute_lyx0(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hyxf = domain->hyxf;
  const double * restrict ux = fluid->ux.data;
  array_t * lyx0_array = &fluid->lyx0;
  double * restrict lyx0 = lyx0_array->data;
  // compute dominant term of lyx | 14
  BEGIN
    const double hy = HYXF(i  );
#if NDIMS == 2
    LYX0(i, j) = 1. / hy * (
        - UX(i  , j-1)
        + UX(i  , j  )
    );
#else
    LYX0(i, j, k) = 1. / hy * (
        - UX(i  , j-1, k  )
        + UX(i  , j  , k  )
    );
#endif
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lyx0_array)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lyx0_array)){
    return 1;
  }
#endif
  return 0;
}

static int compute_lyx1(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxc = domain->hxxc;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uy = fluid->uy.data;
  array_t * lyx1_array = &fluid->lyx1;
  double * restrict lyx1 = lyx1_array->data;
  // compute non-dominant term of lyx | 21
  BEGIN
    const double hx_xm = HXXC(i-1);
    const double hx_xp = HXXC(i  );
    const double jd_xm = JDXC(i-1);
    const double jd_x0 = JDXF(i  );
    const double jd_xp = JDXC(i  );
    const double jdhx_xm = jd_xm / hx_xm;
    const double jdhx_xp = jd_xp / hx_xp;
    const double djdhx = - jdhx_xm + jdhx_xp;
#if NDIMS == 2
    LYX1(i, j) = - 1. / jd_x0 * djdhx * (
        + (        1 == i ? 0. : 0.5) * UY(i-1, j  )
        + (isize + 1 == i ? 0. : 0.5) * UY(i  , j  )
    );
#else
    LYX1(i, j, k) = - 1. / jd_x0 * djdhx * (
        + (        1 == i ? 0. : 0.5) * UY(i-1, j  , k  )
        + (isize + 1 == i ? 0. : 0.5) * UY(i  , j  , k  )
    );
#endif
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lyx1_array)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lyx1_array)){
    return 1;
  }
#endif
  return 0;
}

int compute_lyx(
    const domain_t * domain,
    fluid_t * fluid
){
  compute_lyx0(domain, fluid);
  compute_lyx1(domain, fluid);
  return 0;
}


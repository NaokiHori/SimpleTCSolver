#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/lyy0.h"
#include "array_macros/fluid/lyy1.h"
#include "internal.h"

#define BEGIN \
  for(int j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++){
#define END \
    } \
  }

static int compute_lyy0(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict uy = fluid->uy.data;
  array_t * lyy0_array = &fluid->lyy0;
  double * restrict lyy0 = lyy0_array->data;
  // compute dominant term of lyy
  BEGIN
    const double hy = HYXC(i  );
    LYY0(i, j) = 1. / hy * (
        - UY(i  , j  )
        + UY(i  , j+1)
    );
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lyy0_array)){
    return 1;
  }
  return 0;
}

static int compute_lyy1(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  array_t * lyy1_array = &fluid->lyy1;
  double * restrict lyy1 = lyy1_array->data;
  // compute non-dominant term of lyy
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double jdhx_xm = jd_xm / hx_xm;
    const double jdhx_xp = jd_xp / hx_xp;
    const double djdhx = - jdhx_xm + jdhx_xp;
    LYY1(i, j) = 1. / jd_x0 * djdhx * (
        + 0.5 * UX(i  , j  )
        + 0.5 * UX(i+1, j  )
    );
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lyy1_array)){
    return 1;
  }
  return 0;
}

int compute_lyy(
    const domain_t * domain,
    fluid_t * fluid
){
  compute_lyy0(domain, fluid);
  compute_lyy1(domain, fluid);
  return 0;
}


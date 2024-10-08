#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/lxx.h"
#include "internal.h"

#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++){
#define END \
      } \
    } \
  }

int compute_lxx(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxc = domain->hxxc;
  const double * restrict ux = fluid->ux.data;
  array_t * lxx_array = &fluid->lxx;
  double * restrict lxx = lxx_array->data;
  // compute lxx
  BEGIN
    const double hx = HXXC(i  );
    LXX(i, j, k) = 1. / hx * (
        - UX(i  , j  , k  )
        + UX(i+1, j  , k  )
    );
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lxx_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lxx_array)){
    return 1;
  }
  return 0;
}


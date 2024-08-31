#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/lxy.h"
#include "internal.h"

#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize + 1; i++){
#define END \
      } \
    } \
  }

int compute_lxy(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict uy = fluid->uy.data;
  array_t * lxy_array = &fluid->lxy;
  double * restrict lxy = lxy_array->data;
  // compute lxy
  BEGIN
    const double hx = HXXF(i  );
    LXY(i, j, k) = 1. / hx * (
        - UY(i-1, j  , k  )
        + UY(i  , j  , k  )
    );
  END
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lxy_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lxy_array)){
    return 1;
  }
  return 0;
}


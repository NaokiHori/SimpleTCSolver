#if NDIMS == 3
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/lzy.h"
#include "internal.h"

int compute_lzy(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict uy = fluid->uy.data;
  array_t * lzy_array = &fluid->lzy;
  double * restrict lzy = lzy_array->data;
  // compute lzy | 10
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        LZY(i, j, k) = 1. / hz * (
            - UY(i  , j  , k-1)
            + UY(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lzy_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lzy_array)){
    return 1;
  }
  return 0;
}
#endif

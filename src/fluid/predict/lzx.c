#if NDIMS == 3
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/lzx.h"
#include "internal.h"

int compute_lzx(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict ux = fluid->ux.data;
  array_t * lzx_array = &fluid->lzx;
  double * restrict lzx = lzx_array->data;
  // compute lzx | 10
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        LZX(i, j, k) = 1. / hz * (
            - UX(i  , j  , k-1)
            + UX(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lzx_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lzx_array)){
    return 1;
  }
  return 0;
}
#endif

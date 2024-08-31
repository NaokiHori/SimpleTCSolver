#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/lzz.h"
#include "internal.h"

int compute_lzz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict uz = fluid->uz.data;
  array_t * lzz_array = &fluid->lzz;
  double * restrict lzz = lzz_array->data;
  // compute lzz
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        LZZ(i, j, k) = 1. / hz * (
            - UZ(i  , j  , k  )
            + UZ(i  , j  , k+1)
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lzz_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lzz_array)){
    return 1;
  }
  return 0;
}

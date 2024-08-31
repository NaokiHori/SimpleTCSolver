#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/lxz.h"
#include "internal.h"

int compute_lxz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict uz = fluid->uz.data;
  array_t * lxz_array = &fluid->lxz;
  double * restrict lxz = lxz_array->data;
  // compute lxz
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        const double hx = HXXF(i  );
        LXZ(i, j, k) = 1. / hx * (
            - UZ(i-1, j  , k  )
            + UZ(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lxz_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lxz_array)){
    return 1;
  }
  return 0;
}

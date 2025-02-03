#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/lyz.h"
#include "internal.h"

int compute_lyz(
    const domain_t * domain,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict uz = fluid->uz.data;
  array_t * lyz_array = &fluid->lyz;
  double * restrict lyz = lyz_array->data;
  // compute lyz
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double hy = HYXC(i  );
        LYZ(i, j, k) = 1. / hy * (
            - UZ(i  , j-1, k  )
            + UZ(i  , j  , k  )
        );
      }
    }
  }
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, lyz_array)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, lyz_array)){
    return 1;
  }
  return 0;
}

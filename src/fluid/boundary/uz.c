#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/uz.h"

/**
 * @brief update boundary values of z velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : z velocity
 * @return               : error code
 */
int fluid_update_boundaries_uz(
    const domain_t * domain,
    array_t * array
){
  // update boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * uz = array->data;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        // no-slip
        UZ(      0, j, k) = 0.;
        UZ(isize+1, j, k) = 0.;
      }
    }
  }
  // communicate
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
    if(0 != halo_communicate_in_z(domain, dtypes + 1, array)){
      return 1;
    }
  }
  return 0;
}

#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/t.h"

/**
 * @brief update boundary values of the scalar
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] t      : scalar
 * @return               : error code
 */
int fluid_update_boundaries_t(
    const domain_t * domain,
    array_t * array
){
  // update boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * t = array->data;
    for(int j = 1; j <= jsize; j++){
      T(      0, j) = param_t_xm;
      T(isize+1, j) = param_t_xp;
    }
  }
  // communicate
  {
    static MPI_Datatype dtypes[NDIMS - 1] = {
      MPI_DOUBLE,
    };
    if(0 != halo_communicate_in_y(domain, dtypes + 0, array)){
      return 1;
    }
  }
  return 0;
}


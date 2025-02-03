#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/p.h"

/**
 * @brief update boundary values of the pressure
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] array  : pressure
 * @return               : error code
 */
int fluid_update_boundaries_p(
    const domain_t * domain,
    array_t * array
){
  // update boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * p = array->data;
    for(int j = 1; j <= jsize; j++){
      // Neumann
      P(      0, j) = P(    1, j);
      P(isize+1, j) = P(isize, j);
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


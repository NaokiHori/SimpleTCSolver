#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/psi.h"

/**
 * @brief update boundary values of the scalar potential
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] psi    : scalar potential
 * @return               : error code
 */
int fluid_update_boundaries_psi(
    const domain_t * domain,
    array_t * array
){
  // update boundary values
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * psi = array->data;
    for(int j = 1; j <= jsize; j++){
      // Neumann
      PSI(      0, j) = PSI(    1, j);
      PSI(isize+1, j) = PSI(isize, j);
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


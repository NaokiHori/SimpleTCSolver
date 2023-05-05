#include <stdio.h>
#include "domain.h"
#include "fileio.h"
#include "../internal.h"
#include "internal.h"

/**
 * @brief compute Nusselt number (normalised torque)
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_nusselt(const char fname[], const domain_t *domain, const double time, const fluid_t *fluid){
  double results[3] = {0.};
  logging_nusselt_compute_torque(domain, fluid, results);
  logging_nusselt_compute_dissipation(domain, fluid, results + 2);
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  int myrank;
  MPI_Comm_rank(comm_cart, &myrank);
  if(0 == myrank){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      fprintf(fp, "%8.2f % 18.15e % 18.15e % 18.15e\n", time, results[0], results[1], results[2]);
      fileio_fclose(fp);
    }
  }
  return 0;
}


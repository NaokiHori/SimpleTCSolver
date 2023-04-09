#include <mpi.h>
#include "domain.h"
#include "sdecomp.h"
#include "fluid.h"
#include "arrays/fluid/uz.h"

static int communicate_in_y(const domain_t * restrict domain, double * restrict uz){
  // extract communicator
  const sdecomp_info_t * restrict info = domain->info;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(info, &comm_cart);
  // check negative / positive neighbour ranks
  int neighbours[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
  sdecomp.get_neighbours(info, SDECOMP_X1PENCIL, SDECOMP_YDIR, neighbours);
  // domain size
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  // fixed parameters
  // since data type is defined, number of items is 1
  const int nitems = 1;
  // same tag is fine since I use blocking communication
  const int tag = 0;
  // define datatype in y 
  MPI_Datatype dtype = MPI_DOUBLE;
  MPI_Type_vector(ksize, isize, (isize+2) * (jsize+2), dtype, &dtype);
  MPI_Type_commit(&dtype);
  // send to positive, receive from negative
  MPI_Sendrecv(
    &UZ(1,   jsize, 1), nitems, dtype, neighbours[1], tag,
    &UZ(1,       0, 1), nitems, dtype, neighbours[0], tag,
    comm_cart, MPI_STATUS_IGNORE
  );
  // send to negative, receive from positive
  MPI_Sendrecv(
    &UZ(1,       1, 1), nitems, dtype, neighbours[0], tag,
    &UZ(1, jsize+1, 1), nitems, dtype, neighbours[1], tag,
    comm_cart, MPI_STATUS_IGNORE
  );
  // clean-up used datatype
  MPI_Type_free(&dtype);
  return 0;
}

static int communicate_in_z(const domain_t * restrict domain, double * restrict uz){
  // extract communicator
  const sdecomp_info_t * restrict info = domain->info;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(info, &comm_cart);
  // check negative / positive neighbour ranks
  int neighbours[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
  sdecomp.get_neighbours(info, SDECOMP_X1PENCIL, SDECOMP_ZDIR, neighbours);
  // domain size
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  // define datatype in z 
  MPI_Datatype dtype = MPI_DOUBLE;
  MPI_Type_vector(jsize+2, isize, isize+2, MPI_DOUBLE, &dtype);
  MPI_Type_commit(&dtype);
  // fixed parameters
  // since data type is defined, number of items is 1
  const int nitems = 1;
  // same tag is fine since I use blocking communication
  const int tag = 0;
  // send to positive, receive from negative
  MPI_Sendrecv(
    &UZ(1, 0,   ksize), nitems, dtype, neighbours[1], tag,
    &UZ(1, 0,       0), nitems, dtype, neighbours[0], tag,
    comm_cart, MPI_STATUS_IGNORE
  );
  // send to negative, receive from positive
  MPI_Sendrecv(
    &UZ(1, 0,       1), nitems, dtype, neighbours[0], tag,
    &UZ(1, 0, ksize+1), nitems, dtype, neighbours[1], tag,
    comm_cart, MPI_STATUS_IGNORE
  );
  // clean-up used datatype
  MPI_Type_free(&dtype);
  return 0;
}

static int assign_boundary_conditions_in_x(const domain_t * restrict domain, double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  // set boundary values 
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UZ(      0, j, k) = 0.; // no-slip
      UZ(isize+1, j, k) = 0.; // no-slip
    }
  }
  return 0;
}

/**
 * @brief update boundary values of z velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] uz     : z velocity
 * @return               : error code
 */
int fluid_update_boundaries_uz(const domain_t * restrict domain, double * restrict uz){
  communicate_in_y(domain, uz);
  communicate_in_z(domain, uz);
  assign_boundary_conditions_in_x(domain, uz);
  return 0;
}


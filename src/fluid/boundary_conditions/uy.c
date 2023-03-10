#include "domain.h"
#include "fluid.h"
#include "arrays/fluid.h"



/**
 * @brief update boundary values of y velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[inout] uy     : y velocity
 * @return              : error code
 */
int fluid_update_boundaries_uy(const domain_t * restrict domain, double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update y halo values in 3D ! 23 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 1, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(ksize, isize, (isize+2)*(jsize+2), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &UY(1,   jsize, 1), 1, dtype, pos, 0,
      /* recv from neg. */ &UY(1,       0, 1), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &UY(1,       1, 1), 1, dtype, neg, 0,
      /* recv from pos. */ &UY(1, jsize+1, 1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! update z halo values in 3D ! 23 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int neg, pos;
    MPI_Cart_shift(comm, 2, 1, &neg, &pos);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(jsize+2, isize, isize+2, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    // communicate
    MPI_Sendrecv(
      /* send to   pos. */ &UY(1, 0,   ksize), 1, dtype, pos, 0,
      /* recv from neg. */ &UY(1, 0,       0), 1, dtype, neg, 0,
      comm, MPI_STATUS_IGNORE
    );
    MPI_Sendrecv(
      /* send to   neg. */ &UY(1, 0,       1), 1, dtype, neg, 0,
      /* recv from pos. */ &UY(1, 0, ksize+1), 1, dtype, pos, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! set boundary values ! 6 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UY(      0, j, k) = 1.; // no-slip
      UY(isize+1, j, k) = 0.; // no-slip
    }
  }
  return 0;
}


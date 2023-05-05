#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "common.h"

/**
 * @brief memory allocation with error handler
 * @param[in] count : number of elements
 * @param[in] size  : size of each element
 * @return          : pointer to the allocated buffer
 */
void *common_calloc(const size_t count, const size_t size){
  void *ptr = calloc(count, size);
  if(ptr == NULL){
    fprintf(stderr, "memory allocation error\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return ptr;
}

/**
 * @brief memory deallocation with error handler
 * @param[in] ptr : pointer to the allocated buffer
 */
void common_free(void *ptr){
  free(ptr);
}

/**
 * @brief get current wall time (NOT simulation time units)
 * @return : current time (synchronised among all processes)
 */
double common_get_wtime(void){
  double now = MPI_Wtime();
  // share value of main process
  MPI_Bcast(&now, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return now;
}


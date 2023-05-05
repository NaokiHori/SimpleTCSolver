#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include "config.h"

/**
 * @brief load environment variable and interpret it as an double-precision value
 * @param[in]  dsetname : name of the environment variable
 * @param[out] value    : resulting value
 * @return              : error code
 */
static int get_double(const char dsetname[], double *value){
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  char *string = getenv(dsetname);
  if(NULL == string){
    if(0 == myrank) printf("%s not found\n", dsetname);
    return 1;
  }
  errno = 0;
  *value = strtod(string, NULL);
  if(0 != errno){
    if(0 == myrank) printf("%s: invalid value as double\n", dsetname);
    return 1;
  }
  return 0;
}

const config_t config = {
  .get_double = get_double,
};


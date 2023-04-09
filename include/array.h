#if !defined(ARRAY_H)
#define ARRAY_H

#include "common.h"
#include "domain.h"

typedef struct {
  int glsizes[NDIMS];
  int mysizes[NDIMS];
  int offsets[NDIMS];
  int nadds[NDIMS][2];
  double * restrict data;
} array_t;

extern int array_prepare(const domain_t *domain, const int nadds[NDIMS][2], array_t * restrict * array);
extern int array_destroy(array_t * restrict array);

extern int array_load(const MPI_Comm comm_cart, const char dirname[], const char dsetname[],       array_t *array);
extern int array_dump(const MPI_Comm comm_cart, const char dirname[], const char dsetname[], const array_t *array);

#endif // ARRAY_H

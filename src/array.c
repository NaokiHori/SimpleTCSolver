#include "common.h"
#include "domain.h"
#include "array.h"
#include "fileio.h"

int array_prepare(const domain_t *domain, const int nadds[NDIMS][2], array_t * restrict * array){
  const int *glsizes = domain->glsizes;
  const int *mysizes = domain->mysizes;
  const int *offsets = domain->offsets;
  *array = common_calloc(1, sizeof(array_t));
  // get number of cells
  size_t nitems = 1;
  for(size_t dim = 0; dim < NDIMS; dim++){
    nitems *= mysizes[dim] + nadds[dim][0] + nadds[dim][1];
  }
  (*array)->data = common_calloc(nitems, sizeof(double));
  for(size_t dim = 0; dim < NDIMS; dim++){
    (*array)->glsizes[dim]  = glsizes[dim];
    (*array)->mysizes[dim]  = mysizes[dim];
    (*array)->offsets[dim]  = offsets[dim];
    (*array)->nadds[dim][0] = nadds[dim][0];
    (*array)->nadds[dim][1] = nadds[dim][1];
  }
  return 0;
}

int array_destroy(array_t * restrict array){
  common_free(array->data);
  common_free(array);
  return 0;
}

// array->data, including additional (boundary & halo) cells
//   [1 - nadds[0][0] : mysizes[0] + nadds[0][1]]
//   x
//   [1 - nadds[1][0] : mysizes[1] + nadds[1][1]]
//   x
//   [1 - nadds[2][0] : mysizes[2] + nadds[2][1]]
// buf, holding only data to be written / loaded
//   [1 - nadds[0][0] : mysizes[0] + nadds[0][1]]
//   x
//   [1               : mysizes[1]              ]
//   x
//   [1               : mysizes[2]              ]

static size_t get_index(const int mysizes[NDIMS], const int nadds[NDIMS][2], const size_t indices[NDIMS]){
  const size_t index =
    +  indices[0]
    + (indices[1] + nadds[1][0]) * (size_t)(mysizes[0] + nadds[0][0] + nadds[0][1])
    + (indices[2] + nadds[2][0]) * (size_t)(mysizes[0] + nadds[0][0] + nadds[0][1])
                                 * (size_t)(mysizes[1] + nadds[1][0] + nadds[1][1]);
  return index;
}

int array_load(const MPI_Comm comm_cart, const char dirname[], const char dsetname[], array_t *array){
  // prepare buffer
  double *buf = common_calloc(
      (array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1]) * array->mysizes[1] * array->mysizes[2],
      sizeof(double)
  );
  // read
  const int retval = fileio_r_nd_parallel(
      comm_cart,
      dirname,
      dsetname,
      NDIMS,
      (const int [NDIMS]){
        array->glsizes[2],
        array->glsizes[1],
        array->glsizes[0] + array->nadds[0][0] + array->nadds[0][1],
      },
      (const int [NDIMS]){
        array->mysizes[2],
        array->mysizes[1],
        array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1],
      },
      (const int [NDIMS]){
        array->offsets[2],
        array->offsets[1],
        array->offsets[0],
      },
      NPY_DBL,
      MPI_DOUBLE,
      buf
  );
  if(0 != retval){
    common_free(buf);
    return 1;
  }
  // copy
  double *data = array->data;
  for(size_t cnt = 0, k = 0; k < (size_t)(array->mysizes[2]); k++){
    for(size_t j = 0; j < (size_t)(array->mysizes[1]); j++){
      for(size_t i = 0; i < (size_t)(array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1]); i++){
        const size_t index = get_index(array->mysizes, array->nadds, (size_t [NDIMS]){i, j, k});
        data[index] = buf[cnt++];
      }
    }
  }
  common_free(buf);
  return 0;
}

int array_dump(const MPI_Comm comm_cart, const char dirname[], const char dsetname[], const array_t *array){
  // prepare buffer
  double *buf = common_calloc(
      (array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1]) * array->mysizes[1] * array->mysizes[2],
      sizeof(double)
  );
  // copy
  const double *data = array->data;
  for(size_t cnt = 0, k = 0; k < (size_t)(array->mysizes[2]); k++){
    for(size_t j = 0; j < (size_t)(array->mysizes[1]); j++){
      for(size_t i = 0; i < (size_t)(array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1]); i++){
        const size_t index = get_index(array->mysizes, array->nadds, (size_t [NDIMS]){i, j, k});
        buf[cnt++] = data[index];
      }
    }
  }
  // write
  fileio_w_nd_parallel(
      comm_cart,
      dirname,
      dsetname,
      NDIMS,
      (const int [NDIMS]){
        array->glsizes[2],
        array->glsizes[1],
        array->glsizes[0] + array->nadds[0][0] + array->nadds[0][1],
      },
      (const int [NDIMS]){
        array->mysizes[2],
        array->mysizes[1],
        array->mysizes[0] + array->nadds[0][0] + array->nadds[0][1],
      },
      (const int [NDIMS]){
        array->offsets[2],
        array->offsets[1],
        array->offsets[0],
      },
      NPY_DBL,
      MPI_DOUBLE,
      buf
  );
  common_free(buf);
  return 0;
}


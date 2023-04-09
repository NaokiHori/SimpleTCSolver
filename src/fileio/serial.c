#include <stdio.h>
#include <stdbool.h>
#include "common.h"
#include "fileio.h"
#include "internal.h"

/**
 * @brief read data from a npy file, by one process
 * @param[in]  dirname   : name of directory in which a target npy file is contained
 * @param[in]  dsetname  : name of dataset
 * @param[in]  ndims     : number of dimensions of dataset
 * @param[in]  shape     : shape of dataset
 * @param[in]  dtype     : datatype, e.g. '<f8'
 * @param[in]  totalsize : total datasize of dataset
 * @param[out] data      : pointer to the data to be loaded
 */
int fileio_r_serial(const char dirname[], const char dsetname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize, void *data){
  char *fname = fileio_internal_create_npyfname(dirname, dsetname);
  const size_t header_size = fileio_internal_r_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    fprintf(stderr, "%s: NPY header load failed\n", fname);
    common_free(fname);
    return 1;
  }
  FILE *fp = fileio_fopen(fname, "r");
  if(NULL == fp){
    common_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fileio_fclose(fp);
    common_free(fname);
    return 1;
  }
  const size_t nitems = fread(data, totalsize, 1, fp);
  if(1 != nitems){
    fprintf(stderr, "%s: fread failed\n", fname);
    fileio_fclose(fp);
    common_free(fname);
    return 1;
  }
  fileio_fclose(fp);
  common_free(fname);
  return 0;
}

/**
 * @brief write data to a npy file, by one process
 * @param[in] dirname   : name of directory in which a target npy file is contained
 * @param[in] dsetname  : name of dataset
 * @param[in] ndims     : number of dimensions of dataset
 * @param[in] shape     : shape of dataset
 * @param[in] dtype     : datatype, e.g. '<f8'
 * @param[in] totalsize : total datasize of dataset
 * @param[in] data      : pointer to the data to be written
 */
int fileio_w_serial(const char dirname[], const char dsetname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize, const void *data){
  char *fname = fileio_internal_create_npyfname(dirname, dsetname);
  const size_t header_size = fileio_internal_w_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    common_free(fname);
    return 1;
  }
  FILE *fp = fileio_fopen(fname, "a");
  if(NULL == fp){
    common_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fileio_fclose(fp);
    common_free(fname);
    return 1;
  }
  const size_t nitems = fwrite(data, totalsize, 1, fp);
  if(1 != nitems){
    fprintf(stderr, "%s: fwrite failed\n", fname);
    fileio_fclose(fp);
    common_free(fname);
    return 1;
  }
  fileio_fclose(fp);
  common_free(fname);
  return 0;
}


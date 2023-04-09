#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "common.h"
#include "fileio.h"
#include "snpyio.h"

// 1-byte boolean
const char NPY_BOL[] = {"'|b1'"};
// 4-byte little-endian integer
const char NPY_INT[] = {"'<i4'"};
// 8-byte little-endian floating point
const char NPY_DBL[] = {"'<f8'"};

/**
 * @brief allocate and initialise string having a NPY file name
 * @param[in] dirname  : name of directory, e.g. "output/save/step0000000000"
 * @param[in] dsetname : name of the dataset, e.g. "ux"
 * @return             : file name, e.g. "output/save/step0000000000/ux.npy"
 */
char *fileio_internal_create_npyfname(const char dirname[], const char dsetname[]){
  if(NULL == dirname){
    fprintf(stderr, "dirname is NULL\n");
    return NULL;
  }
  if(NULL == dsetname){
    fprintf(stderr, "dsetname is NULL\n");
    return NULL;
  }
  const char slash[] = {"/"};
  const char suffix[] = {".npy"};
  const size_t nchars =
    + strlen( dirname)
    + strlen(   slash)
    + strlen(dsetname)
    + strlen(  suffix);
  char *fname = common_calloc(nchars + 2, sizeof(char));
  snprintf(fname, nchars + 1, "%s%s%s%s", dirname, slash, dsetname, suffix);
  fname[nchars + 1] = '\0';
  return fname;
}

/**
 * @brief file opener
 * @param[in] path : pointer to the file name to be opened
 * @param[in] mode : file open mode
 * @return         : file pointer
 */
FILE *fileio_fopen(const char *path, const char *mode){
  errno = 0;
  FILE *stream = fopen(path, mode);
  if(NULL == stream){
    perror(path);
  }
  return stream;
}

/**
 * @brief file closer
 * @param[in] stream : file pointer to be closed
 * @return           : error code
 */
int fileio_fclose(FILE *stream){
  if(NULL == stream){
    return 1;
  }
  errno = 0;
  const int retval = fclose(stream);
  if(0 != retval){
    perror("");
    return 1;
  }
  return 0;
}

/**
 * @brief create directory name from prefix and current time step
 * @param[in] prefix : e.g. "output/save/step"
 * @param[in] step   : current time step, e.g. 0
 * @return           : resulting directory name, e.g. "output/save/step0000000000"
 */
char *fileio_create_dirname(const char prefix[], const int step){
  if(NULL == prefix){
    fprintf(stderr, "prefix is NULL\n");
    return NULL;
  }
  const int ndigits_step = 10;
  const size_t nchars = strlen(prefix) + ndigits_step + 1;
  char *dirname = common_calloc(nchars, sizeof(char));
  snprintf(dirname, nchars, "%s%0*d", prefix, ndigits_step, step);
  dirname[nchars-1] = '\0';
  return dirname;
}

/**
 * @brief create directory
 * @param[in] dirname : name of the directory to be created
 * @return            : error code
 */
int fileio_mkdir(const char dirname[]){
  // NOTE: call this function ONLY from the main process
  // NOTE: continue even if failed,
  //   since we want to override previous data (errorcode: EEXIST)
  // RWX masks for user, group, and others (0o777, ref. "man 2 chmod")
  const mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO;
  if(0 != mkdir(dirname, mode)){
    perror(dirname);
    return 1;
  }
  return 0;
}

/**
 * @brief wrapper function of snpyio_r_header with error handling
 * @param[in] fname            : name of the file from which data is loaded
 * @param[in] ndims            : number of dimensions of the data to be loaded
 * @param[in] shape            : sizes of the data in each dimension to be loaded
 * @param[in] dtype            : datatype of the data to be loaded
 * @param[in] is_fortran_order : memory contiguous direction, normally false
 * @return                     : (success) size of npy header in byte
 *                               (failure) 0
 */
size_t fileio_internal_r_npy_header(const char fname[], const size_t ndims, const size_t *shape, const char *dtype, const bool is_fortran_order){
  const char msg[] = {"NPY header read failed"};
  FILE *fp = fileio_fopen(fname, "r");
  if(NULL == fp){
    return 0;
  }
  // load header, return header size when succeeded, return 0 otherwise
  size_t ndims_ = 0;
  size_t *shape_ = NULL;
  char *dtype_ = NULL;
  bool is_fortran_order_ = false;
  size_t header_size = snpyio_r_header(&ndims_, &shape_, &dtype_, &is_fortran_order_, fp);
  fileio_fclose(fp);
  // check arguments, return loaded header size when all OK, return 0 otherwise
  // ndims
  if(ndims != ndims_){
    fprintf(stderr, "%s(%s), ndims: %zu expected, %zu obtained\n", msg, fname, ndims, ndims_);
    header_size = 0;
  }
  // shape (for each dimension)
  for(size_t n = 0; n < (ndims < ndims_ ? ndims : ndims_); n++){
    if(shape[n] != shape_[n]){
      fprintf(stderr, "%s(%s), shape[%zu]: %zu expected, %zu obtained\n", msg, fname, n, shape[n], shape_[n]);
      header_size = 0;
    }
  }
  // dtype
  if(0 != strcmp(dtype, dtype_)){
    fprintf(stderr, "%s(%s), dtype: %s expected, %s obtained\n", msg, fname, dtype, dtype_);
    header_size = 0;
  }
  // is_fortran_order
  if(is_fortran_order != is_fortran_order_){
    fprintf(stderr, "%s(%s), is_fortran_order: %s expected, %s obtained\n", msg, fname, is_fortran_order ? "true" : "false", is_fortran_order_ ? "true" : "false");
    header_size = 0;
  }
  common_free(shape_);
  common_free(dtype_);
  return header_size;
}

/**
 * @brief wrapper function of snpyio_w_header with error handling
 * @param[in] fname            : name of the file to which data is dumped
 * @param[in] ndims            : number of dimensions of the data to be written
 * @param[in] shape            : sizes of the data in each dimension to be written
 * @param[in] dtype            : datatype of the data to be written
 * @param[in] is_fortran_order : memory contiguous direction, normally false
 * @return                     : (success) size of npy header in byte
 *                               (failure) 0
 */
size_t fileio_internal_w_npy_header(const char fname[], const size_t ndims, const size_t *shape, const char dtype[], const bool is_fortran_order){
  const char msg[] = {"NPY header write failed"};
  FILE *fp = fileio_fopen(fname, "w");
  if(NULL == fp){
    return 0;
  }
  const size_t header_size = snpyio_w_header(ndims, shape, dtype, is_fortran_order, fp);
  if(0 == header_size){
    fprintf(stderr, "%s(%s)\n", msg, fname);
  }
  fileio_fclose(fp);
  return header_size;
}


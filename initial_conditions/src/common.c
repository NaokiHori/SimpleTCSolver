#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "common.h"
#include "config.h"
#include "snpyio.h"

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
    exit(EXIT_FAILURE);
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

// 1-byte (presumably  8-bit) boolean (byte-order does not make sense)
const char NPY_BOL[] = {"'|b1'"};
// 4-byte (presumably 32-bit) little-endian integer
const char NPY_INT[] = {"'<i4'"};
// 8-byte (presumably 64-bit) little-endian floating point
const char NPY_DBL[] = {"'<f8'"};

void common_save(const char fname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize, const void *data){
  // sanitise fname
  if(NULL == fname){
    fprintf(stderr, "ERROR: fname is NULL\n");
    return;
  }
  // get and sanitise dirname
  const char *dirname = config.get_string("dirname");
  if(NULL == dirname){
    fprintf(stderr, "ERROR: dirname is NULL\n");
    return;
  }
  if(0 == strcmp("", dirname)){
    fprintf(stderr, "ERROR: dirname is empty\n");
    fprintf(stderr, "ERROR:   if you intend to specify current directory,\n");
    fprintf(stderr, "ERROR:   give \".\" instead\n");
    return;
  }
  // concatenate names of directory and file
  const size_t nchars = strlen(dirname) + 1 + strlen(fname) + 1;
  char *filename = common_calloc(nchars, sizeof(char));
  snprintf(filename, nchars, "%s/%s", dirname, fname);
  filename[nchars-1] = '\0';
  // start to dump file
  errno = 0;
  FILE *fp = fopen(filename, "w");
  if(fp == NULL){
    // file open error, free memory and terminate
    perror(filename);
    goto terminate0;
  }
  size_t retval = 0;
  retval = snpyio_w_header(ndims, shape, dtype, false, fp);
  if(0 == retval){
    // failed to write header, close file, free memory and terminate
    fprintf(stderr, "%s:snpyio_w_header failed\n", fname);
    goto terminate1;
  }
  retval = fwrite(data, totalsize, 1, fp);
  if(1 != retval){
    // failed to write data, close file, free memory and terminate
    fprintf(stderr, "%s:fwrite failed\n", fname);
    goto terminate1;
  }
terminate1:
  fclose(fp);
terminate0:
  common_free(filename);
}


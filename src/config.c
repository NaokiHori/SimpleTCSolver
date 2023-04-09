#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include <mpi.h>
#include "common.h"
#include "fileio.h"
#include "config.h"


typedef enum {
  DTYPE_STR,
  DTYPE_BOL,
  DTYPE_DBL,
} dtype_t;

typedef struct {
  const char *label;
  const dtype_t dtype;
} key_t;

static void **values = NULL;

// available environment variables 
static const key_t keys[] = {
  // char *
  {.label = "dirname_ic",   .dtype = DTYPE_STR},
  // bool
  {.label = "implicitx",    .dtype = DTYPE_BOL},
  {.label = "implicity",    .dtype = DTYPE_BOL},
  {.label = "implicitz",    .dtype = DTYPE_BOL},
  // double
  {.label = "timemax",      .dtype = DTYPE_DBL},
  {.label = "wtimemax",     .dtype = DTYPE_DBL},
  {.label = "log_rate",     .dtype = DTYPE_DBL},
  {.label = "save_rate",    .dtype = DTYPE_DBL},
  {.label = "save_after",   .dtype = DTYPE_DBL},
  {.label = "stat_rate",    .dtype = DTYPE_DBL},
  {.label = "stat_after",   .dtype = DTYPE_DBL},
  {.label = "coef_dt_adv",  .dtype = DTYPE_DBL},
  {.label = "coef_dt_dif",  .dtype = DTYPE_DBL},
  {.label = "Re",           .dtype = DTYPE_DBL},
};

static size_t get_nitems(void){
  // get number of keys
  return sizeof(keys) / sizeof(keys[0]);
}

static int find_label_index(const char *label, size_t *index){
  // return index of the given key in the key list
  //   which is used to access the corresponding value
  const size_t nitems = get_nitems();
  for(size_t n = 0; n < nitems; n++){
    if(0 == strcmp(label, keys[n].label)){
      *index = n;
      return 0;
    }
  }
  return 1;
}

/*** getters ***/
// there should be a better way...

#define GETTER(type, name) \
  static type get_##name(void){ \
    size_t index = 0; \
    const int retval = find_label_index(#name, &index); \
    if(0 != retval){ \
      /* should not be here */ \
      return (type)(0); \
    }else{ \
      return *(type *)(values[index]); \
    } \
  }

static char * get_dirname_ic(void){
  size_t index = 0;
  const int retval = find_label_index("dirname_ic", &index);
  if(0 != retval){
    return NULL;
  }else{
    return (char *)(values[index]);
  }
}

GETTER(bool, implicitx)
GETTER(bool, implicity)
GETTER(bool, implicitz)
GETTER(double, timemax)
GETTER(double, wtimemax)
GETTER(double, log_rate)
GETTER(double, save_rate)
GETTER(double, save_after)
GETTER(double, stat_rate)
GETTER(double, stat_after)
GETTER(double, coef_dt_adv)
GETTER(double, coef_dt_dif)
GETTER(double, Re)

static int compute_nchars_max(void){
  // auxiliary function to compute
  //   maximum number of characters in "keys"
  //   just for pretty-print purpose
  const size_t nitems = get_nitems();
  int retval = 0;
  for(size_t n = 0; n < nitems; n++){
    const int nchars = (int)(strlen(keys[n].label));
    retval = retval < nchars ? nchars : retval;
  }
  return retval;
}

/**
 * @brief interpret string as a boolean
 * @param[in]  str   : string
 * @param[out] value : value after converted
 * @return           : (success) 0
 *                     (failure) 1
 */
static int convert_to_bol(const char str[], bool *value){
  if(
      0 == strcmp(str, "true")
   || 0 == strcmp(str, "True")
   || 0 == strcmp(str, "TRUE")
  ){
    *value = true;
    return 0;
  }
  if(
      0 == strcmp(str, "false")
   || 0 == strcmp(str, "False")
   || 0 == strcmp(str, "FALSE")
  ){
    *value = false;
    return 0;
  }
  return 1;
}

/**
 * @brief interpret string as a double-precision floating number
 * @param[in]  str   : string
 * @param[out] value : value after converted
 * @return           : (success) 0
 *                     (failure) 1
 */
static int convert_to_dbl(const char str[], double *value){
  errno = 0;
  *value = strtod(str, NULL);
  if(0 != errno){
    return 1;
  }
  return 0;
}

/**
 * @brief constructor to
 *          1. load environment variables,
 *          2. cnovert to proper types,
 *          3. and store them
 * @return : (success) 0
 *           (failure) 1
 */
static int construct(void){
  // allocate pointers to nitems key-value pairs
  const size_t nitems = get_nitems();
  values = common_calloc(nitems, sizeof(void *));
  // check parameters and print them
  int retval = 0;
  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const int nchars_max = compute_nchars_max();
  // for each key
  for(size_t n = 0; n < nitems; n++){
    // allocate space for one key-value
    const char *label = keys[n].label;
    char *value_in_str = getenv(label);
    if(NULL != value_in_str){
      if(0 == myrank) printf("#. ENV %*s is     found and ", nchars_max, label);
      // found, copy value
      const dtype_t dtype = keys[n].dtype;
      switch(dtype){
        case DTYPE_STR:
          {
            // type is string, copy and save
            const size_t nchars = strlen(value_in_str);
            // +1: NUL
            values[n] = common_calloc(nchars + 1, sizeof(char));
            memcpy(values[n], value_in_str, sizeof(char) * nchars);
            if(0 == myrank) printf("  valid");
            break;
          }
        case DTYPE_BOL:
          {
            // type is bool, convert and save
            values[n] = common_calloc(1, sizeof(bool));
            if(0 != convert_to_bol(value_in_str, values[n])){
              retval += 1;
              if(0 == myrank) printf("INVALID");
            }else{
              if(0 == myrank) printf("  valid");
            }
            break;
          }
        case DTYPE_DBL:
          {
            // type is double, convert and save
            values[n] = common_calloc(1, sizeof(double));
            if(0 != convert_to_dbl(value_in_str, values[n])){
              retval += 1;
              if(0 == myrank) printf("INVALID");
            }else{
              if(0 == myrank) printf("  valid");
            }
            break;
          }
      }
      if(0 == myrank) printf(": %s\n", value_in_str);
    }else{
      if(0 == myrank) printf("#. ENV %*s is NOT FOUND\n", nchars_max, label);
      retval += 1;
    }
  }
  // gather results just in case
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(0 != retval && 0 == myrank){
    printf("ERROR: Parameter(s) is missing / invalid\n");
    printf("ERROR: Check above table and give it correctly\n");
  }
  return retval;
}

/**
 * @brief destructor - clean-up memories used internally
 * @return : error code
 */
static int destruct(void){
  const size_t nitems = get_nitems();
  for(size_t n = 0; n < nitems; n++){
    common_free(values[n]);
  }
  common_free(values);
  return 0;
}

/**
 * @brief save parameters to NPY files
 * @param[in] dirname : name of directory to which data is saved
 * @return            : error code
 */
static int output(const char dirname[]){
  const double Re = get_Re();
  fileio_w_serial(dirname, "Re", 0, NULL, NPY_DBL, sizeof(double), &Re);
  return 0;
}

const config_t config = {
  .construct = construct,
  .destruct  = destruct,
  .output    = output,
  .get       = {
    // char *
    .dirname_ic   = get_dirname_ic,
    // bool
    .implicitx    = get_implicitx,
    .implicity    = get_implicity,
    .implicitz    = get_implicitz,
    // double
    .timemax      = get_timemax,
    .wtimemax     = get_wtimemax,
    .log_rate     = get_log_rate,
    .save_rate    = get_save_rate,
    .save_after   = get_save_after,
    .stat_rate    = get_stat_rate,
    .stat_after   = get_stat_after,
    .coef_dt_adv  = get_coef_dt_adv,
    .coef_dt_dif  = get_coef_dt_dif,
    .Re           = get_Re,
  },
};


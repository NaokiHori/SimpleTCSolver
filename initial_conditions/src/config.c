#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include "common.h"
#include "config.h"

typedef struct dict_t_ {
  const char *key;
  char *value;
} dict_t;

static dict_t **dict = NULL;

static const char *keys[] = {
  // name of environment variables
  "dirname",
  "lx",
  "ly",
  "lz",
  "glisize",
  "gljsize",
  "glksize",
  "uniformx",
};

static int get_nitems(void){
  // get number of keys,
  //   i.e. size of dictionary defined above
  return sizeof(keys)/sizeof(char *);
}

static int find_key_index(const char *key){
  // return index of the given key in the dictionary
  //   which is used to access the corresponding value
  const int nitems = get_nitems();
  for(int n = 0; n < nitems; n++){
    if(0 == strcmp(key, dict[n]->key)){
      return n;
    }
  }
  printf("unknown envname: %s\n", key);
  exit(EXIT_FAILURE);
}

/**
 * @brief load environment variables and create dictionary to store them
 * @return : error code
 */
static int load(void){
  // allocate pointers to nitems key-value pairs
  const int nitems = get_nitems();
  dict = common_calloc((size_t)nitems, sizeof(dict_t *));
  for(int n = 0; n < nitems; n++){
    // allocate space for one key-value
    dict[n] = common_calloc(1, sizeof(dict_t));
    // assign ptr to key
    const char *key = keys[n];
    dict[n]->key = key;
    // copy value instead of assign ptr,
    //   since ptr to env var might be changed by system
    char *value = getenv(key);
    if(value != NULL){
      size_t nchar = strlen(value);
      dict[n]->value = common_calloc(nchar+1, sizeof(char));
      memcpy(dict[n]->value, value, sizeof(char)*nchar);
    }else{
      dict[n]->value = NULL;
    }
  }
  // print the final dictionary
  {
    const int maxchar = 18;
    for(int n = 0; n < nitems; n++){
      const char *key = dict[n]->key;
      char *value = dict[n]->value;
      if(value == NULL){
        printf("#. ENV %*s is NOT found\n",     maxchar, key);
      }else{
        printf("#. ENV %*s is     found: %s\n", maxchar, key, value);
      }
    }
  }
  return 0;
}

/**
 * @brief clean-up memories used by storing dictionary
 * @return : error code
 */
static int unload(void){
  const int nitems = get_nitems();
  for(int n = 0; n < nitems; n++){
    common_free(dict[n]->value);
    common_free(dict[n]);
  }
  common_free(dict);
  return 0;
}

/**
 * @brief load ENV and return it
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static char *get_string(const char envname[]){
  /*
   * return "value" of key-value pair,
   *   where "key" is the given name of ENV
   */
  const int index = find_key_index(envname);
  char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    exit(EXIT_FAILURE);
  }
  return value;
}

/**
 * @brief load ENV and return it as a boolean value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static bool get_bool(const char envname[]){
  /*
   * return true  if true  is given
   * return false if false is given
   */
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    exit(EXIT_FAILURE);
  }
  if(0 == strcmp(value, "true")){
    return true;
  }
  if(0 == strcmp(value, "True")){
    return true;
  }
  if(0 == strcmp(value, "TRUE")){
    return true;
  }
  if(0 == strcmp(value, "false")){
    return false;
  }
  if(0 == strcmp(value, "False")){
    return false;
  }
  if(0 == strcmp(value, "FALSE")){
    return false;
  }
  printf("ERROR: %s cannot be interpreted as bool\n", value);
  exit(EXIT_FAILURE);
}

/**
 * @brief load ENV and return it as an integer value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static int get_int(const char envname[]){
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    exit(EXIT_FAILURE);
  }
  // try to convert value to an integer
  errno = 0;
  long retval_ = strtol(value, NULL, 10);
  if(errno != 0 || INT_MIN > retval_ || INT_MAX < retval_){
    printf("ERROR: over/underflow is detected: %s\n", value);
    exit(EXIT_FAILURE);
  }
  // conert long to int
  return (int)retval_;
}

/**
 * @brief load ENV and return it as a double value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static double get_double(const char envname[]){
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    exit(EXIT_FAILURE);
  }
  errno = 0;
  double retval = strtod(value, NULL);
  if(errno != 0){
    printf("ERROR: over/underflow is detected: %s\n", value);
    exit(EXIT_FAILURE);
  }
  return retval;
}

const config_t config = {
  .load   = load,
  .unload = unload,
  .get_string = get_string,
  .get_bool   = get_bool,
  .get_int    = get_int,
  .get_double = get_double,
};


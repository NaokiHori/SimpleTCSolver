/*
 * Copyright 2022-2023 Naoki Hori
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// https://github.com/NaokiHori/SimpleNpyIO

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "snpyio.h"


/* assumptions and definitions */
// check 1 byte is 8 bits
#if CHAR_BIT != 8
#error "CHAR_BIT is not 8"
#endif
/* ! all npy files should start from this magic string ! 1 ! */
static const char magic_string[] = {"\x93NUMPY"};

// null character
static const char null_char = '\0';

/* messages for logginig and error handlings */
#define SNPYIO_MESSAGE(stream, level, ...){    \
  fprintf(stream, "[NPYIO %s]\n", level);      \
  fprintf(stream, "    line: %d\n", __LINE__); \
  fprintf(stream, "    func: %s\n", __func__); \
  fprintf(stream, "    ");                     \
  fprintf(stream, __VA_ARGS__);                \
}

/* logger */
#if defined(SNPYIO_ENABLE_LOGGING)
#define SNPYIO_LOGGING(...){              \
  FILE *fp = fopen("snpyio.log", "a");    \
  SNPYIO_MESSAGE(fp, "LOG", __VA_ARGS__); \
  fclose(fp);                             \
}
#else
#define SNPYIO_LOGGING(...)
#endif

/* dump error message */
#define SNPYIO_ERROR(...){                      \
  SNPYIO_MESSAGE(stderr, "ERROR", __VA_ARGS__); \
}

/* fatal error, dump error message and abort */
#define SNPYIO_FATAL(...){                      \
  SNPYIO_MESSAGE(stderr, "FATAL", __VA_ARGS__); \
  exit(EXIT_FAILURE);                           \
}

/* take care of allocation / deallocation */
static void *kernel_calloc(const size_t count, const size_t size){
  if(count > SIZE_MAX / size){
    SNPYIO_FATAL("memory allocation failed: calloc(%zu, %zu)\n", count, size);
  }
  void *ptr = calloc(count, size);
  if(ptr == NULL){
    SNPYIO_FATAL("memory allocation failed: calloc(%zu, %zu)\n", count, size);
  }
  return ptr;
}

static void kernel_free(void *ptr){
  free(ptr);
}

/* ! singly-linked list ! 6 ! */
typedef struct node_t_ {
  // pointer to the data
  void *ptr;
  // pointer to the next node
  struct node_t_ *next;
} node_t;

/* ! memory pool storing all pointers to the allocated buffers in this library ! 1 ! */
static node_t *memory = NULL;

/* search pointer in the linked list whose root is "memory" */
static bool search_list(const void *ptr, node_t **node_prev){
  /*
   * return true if "ptr" is found,
   *   also node_prev is set to point the PREVIOUS node
   *   so that removing the node having "ptr" becomes easier
   *
   * return false if "ptr" is not found,
   *   also node_prev is set to NULL
   *
   * NOTE: being "node_prev == NULL" does not mean
   *   the pointer is not found (when "ptr" is found as the root node)
   */
  *node_prev = NULL;
  node_t *node_curr = memory;
  while(node_curr){
    if(node_curr->ptr == ptr){
      return true;
    }
    *node_prev = node_curr;
    node_curr = node_curr->next;
  }
  return false;
}

static int attach_list(void *ptr){
  if(ptr == NULL){
    return 0;
  }
  // check duplication,
  //   since we should NOT find the given pointer
  //   in the memory pool
  node_t *node_prev = NULL;
  const bool is_found = search_list(ptr, &node_prev);
  if(is_found){
    SNPYIO_FATAL("duplication is found in the pool\n");
  }
  // allocate a node to hold ptr information
  node_t *node_new = kernel_calloc(1, sizeof(node_t));
  node_new->ptr = ptr;
  // update root
  node_new->next = memory;
  memory = node_new;
  return 0;
}

static int detach_list(const void *ptr){
  if(ptr == NULL){
    return 0;
  }
  // check existence and location
  node_t *node_prev = NULL;
  const bool is_found = search_list(ptr, &node_prev);
  if(!is_found){
    SNPYIO_FATAL("memory to be removed is not found in the pool\n");
  }
  // deallocate node (NOT memory)
  node_t *node_curr = NULL;
  if(node_prev == NULL){
    // update root
    node_curr = memory;
    memory = memory->next;
  }else{
    // reconnect
    node_curr = node_prev->next;
    node_prev->next = node_curr->next;
  }
  kernel_free(node_curr);
  return 0;
}

static void *memory_calloc(const size_t count, const size_t size){
  void *ptr = kernel_calloc(count, size);
  attach_list(ptr);
  return ptr;
}

static void memory_free(void *ptr){
  detach_list(ptr);
  kernel_free(ptr);
}

static void error_handlings(void){
  // free all buffers attached to the memory pool
  while(memory){
    node_t *next = memory->next;
    kernel_free(memory->ptr);
    kernel_free(memory);
    memory = next;
  }
}

/* auxiliary functions which are used by writer and reader */

// https://gist.github.com/NaokiHori/91c560a59f4e4ef37eb33b8e1c055fbc
static bool is_big_endian(void){
  const uint16_t val = 1 << 8;
  return (bool)(((const uint8_t *)(&val))[0]);
}

// https://gist.github.com/NaokiHori/81ad6e1562e1ec23253246902c281cc2
static int convert_endian(void *val, const size_t size){
  // NULL check
  if(val == NULL){
    SNPYIO_ERROR("val is NULL\n");
    return -1;
  }
  // positive size check
  if(size <= 0){
    SNPYIO_ERROR("size of buffer should be positive: %zu\n", size);
    return -1;
  }
  // reject too large size
  // uint32_t (4 bytes) is the biggest datatype in this library
  if(size > 4){
    SNPYIO_ERROR("size of buffer is larger than 4: %zu\n", size);
    return -1;
  }
  const size_t n_bytes = size/sizeof(uint8_t);
  // swap and write back
  // use buffer for simplicity
  uint8_t *buf = memory_calloc(n_bytes, sizeof(uint8_t));
  for(size_t i = 0; i < n_bytes; i++){
    buf[i] = ((uint8_t *)val)[n_bytes-i-1];
  }
  memcpy(val, buf, size);
  memory_free(buf);
  return 0;
}

static int find_pattern(size_t *location, const char *p0, const size_t size_p0, const char *p1, const size_t size_p1){
  /*
   * try to find a pattern "p1" in "p0"
   *   and return its location IN BYTES
   * -1 is returned when error is detected
   *   if the pattern is not found
   * note that sizes of "p0" and "p1" are
   *   in BYTES, NOT number of elements
   * thus it is necessary to divide by the
   *   sizeof original datatype
   *   after the result is obtained
   */
  // NULL checks
  if(p0 == NULL){
    SNPYIO_ERROR("p0 is NULL\n");
    return -1;
  }
  if(p1 == NULL){
    SNPYIO_ERROR("p1 is NULL\n");
    return -1;
  }
  // p0 is shorter than p1, return not found
  if(size_p0 < size_p1){
    SNPYIO_ERROR("try to find a pattern p1 (size: %zu) in p0 (size: %zu)\n", size_p1, size_p0);
    return -1;
  }
  // e.g., size_p0 = 7, size_p1 = 3
  //     0 1 2 3 4 5 6
  // p0: a b c d e f g
  // p1: x y z
  //       x y z
  //         x y z
  //           x y z
  //             x y z
  //     ^       ^
  //    imin    imax
  size_t imin = 0;
  size_t imax = size_p0-size_p1;
  for(size_t i = imin; i <= imax; i++){
    if(memcmp(p0+i, p1, size_p1) == 0){
      *location = i;
      return 0;
    }
  }
  // not found
  // this is expected in some cases
  return -1;
}

/* reader */

static int load_magic_string(size_t *buf_size, FILE *fp){
  /*
   * all npy file should start with \x93NUMPY,
   *   which is checked here by comparing
   *   the fist 6 bytes of the file
   *   and the magic string given a priori
   */
  if(fp == NULL){
    SNPYIO_ERROR("fp is NULL\n");
    return -1;
  }
  const size_t nitems = strlen(magic_string);
  // allocate buffer and load from file
  // NOTE: file pointer is moved forward as well
  uint8_t *buf = memory_calloc(nitems, sizeof(uint8_t));
  {
    const size_t retval = fread(buf, sizeof(uint8_t), nitems, fp);
    if(retval != nitems){
      SNPYIO_ERROR("fread failed (%zu items loaded, while %zu items requested)\n", retval, nitems);
      return -1;
    }
  }
  *buf_size = sizeof(uint8_t)*nitems;
  // compare memories of
  //   1. "buf" (loaded from file)
  //   2. "magic_str" (answer)
  if(memcmp(buf, magic_string, *buf_size) != 0){
    SNPYIO_ERROR("file should start with magic string \"\\x93NUMPY\"\n");
    return -1;
  }
  memory_free(buf);
  return 0;
}

static int load_versions(uint8_t *major_version, uint8_t *minor_version, size_t *buf_size, FILE *fp){
  /*
   * check version of the file
   * for now 1.0, 2.0, and 3.0 are considered,
   *    and others are rejected
   */
  const size_t nitems_major_version = 1;
  const size_t nitems_minor_version = 1;
  const size_t buf_size_major_version = sizeof(uint8_t)*nitems_major_version;
  const size_t buf_size_minor_version = sizeof(uint8_t)*nitems_minor_version;
  // load from file
  // (file pointer is moved forward as well)
  {
    const size_t retval = fread(major_version, sizeof(uint8_t), nitems_major_version, fp);
    if(retval != nitems_major_version){
      SNPYIO_ERROR("fread failed (%zu items loaded, while %zu items requested)\n", retval, nitems_major_version);
      return -1;
    }
  }
  {
    const size_t retval = fread(minor_version, sizeof(uint8_t), nitems_minor_version, fp);
    if(retval != nitems_minor_version){
      SNPYIO_ERROR("fread failed (%zu items loaded, while %zu items requested)\n", retval, nitems_minor_version);
      return -1;
    }
  }
  // check version 1.x or 2.x or 3.x
  if(*major_version != 1 && *major_version != 2 && *major_version != 3){
    SNPYIO_ERROR("major version (%u) should be 1 or 2\n", *major_version);
    return -1;
  }
  // check version x.0
  if(*minor_version != 0){
    SNPYIO_ERROR("minor version (%u) should be 0\n", *minor_version);
    return -1;
  }
  SNPYIO_LOGGING("major version: %u\n", *major_version);
  SNPYIO_LOGGING("minor version: %u\n", *minor_version);
  *buf_size = buf_size_major_version + buf_size_minor_version;
  return 0;
}

static int load_header_len(size_t *header_len, size_t *buf_size, size_t major_version, FILE *fp){
  /*
   * check header size of the npy file
   * in particular HEADER_LEN = len(dict) + len(padding)
   *   is loaded
   * memory size of this variable depends on the major version
   *   of the npy file, 2 bytes for major_version = 1,
   *   while 4 bytes for major_version = 2
   */
  /* ! buffer size differs based on major_version ! 10 ! */
  size_t nitems;
  if(major_version == 1){
    *buf_size = sizeof(uint16_t);
    // usually 2
    nitems = *buf_size/sizeof(uint8_t);
  }else{
    *buf_size = sizeof(uint32_t);
    // usually 4
    nitems = *buf_size/sizeof(uint8_t);
  }
  /* ! allocate buffer and load corresponding memory size from file ! 8 ! */
  uint8_t *buf = memory_calloc(nitems, sizeof(uint8_t));
  {
    const size_t retval = fread(buf, sizeof(uint8_t), nitems, fp);
    if(retval != nitems){
      SNPYIO_ERROR("fread failed (%zu items loaded, while %zu items requested)\n", retval, nitems);
      return -1;
    }
  }
  /* ! convert endian of loaded buffer when needed ! 6 ! */
  if(is_big_endian()){
    if(convert_endian(header_len, sizeof(*header_len)) != 0){
      SNPYIO_ERROR("convert_endian failed\n");
      return -1;
    }
  }
  /* ! interpret buffer (sequence of uint8_t) as a value having corresponding datatype ! 11 ! */
  if(major_version == 1){
    // interpret as a 2-byte value
    uint16_t tmp = 0;
    memcpy(&tmp, buf, sizeof(uint8_t)*nitems);
    *header_len = (size_t)tmp;
  }else{
    // interpret as a 4-byte value
    uint32_t tmp = 0;
    memcpy(&tmp, buf, sizeof(uint8_t)*nitems);
    *header_len = (size_t)tmp;
  }
  memory_free(buf);
  SNPYIO_LOGGING("header_len: %zu\n", *header_len);
  return 0;
}

static int load_dict_and_padding(uint8_t **dict_and_padding, size_t *buf_size, size_t header_len, FILE *fp){
  /*
   * load dictionary and padding
   * loading padding is also necessary (or at least one of the easiest ways)
   *   to move file pointer "fp" forward
   */
  const size_t nitems = header_len/sizeof(uint8_t);
  *dict_and_padding = memory_calloc(nitems, sizeof(uint8_t));
  {
    const size_t retval = fread(*dict_and_padding, sizeof(uint8_t), nitems, fp);
    if(retval != nitems){
      SNPYIO_ERROR("fread failed (%zu items loaded, while %zu items requested)\n", retval, nitems);
      return -1;
    }
  }
  *buf_size = header_len;
  return 0;
}

static int extract_dict(char **dict, uint8_t *dict_and_padding, size_t header_len){
  /*
   * extract dictionary "dict" from "dict_and_padding",
   *   which contains dictionary and padding
   * also unnecessary spaces are removed from the original array
   *   to simplify the following procedures
   * note that spaces inside quotations (e.g. dictionary key might contain spaces)
   *   should NOT be eliminated, which are "necessary spaces" so to say
   */
  // look for "{" and "}" to find the range of dict in dict_and_padding,
  // e.g., dict_and_padding:
  //   <------------------------ dict ------------------------><- padding ->
  //   {'descr': VALUE, 'fortran_order': VALUE, 'shape': VALUE}           \n
  //   ^                                                      ^
  //   s                                                      e
  // s: start
  // this should be 0, because of NPY format definition,
  //   i.e., dict should start just after HEADER_LEN
  // confirm it just in case
  //   by checking whether the first byte is "{"
  const size_t s = 0;
  {
    // use char since it's "{"
    char p0 = (char)(dict_and_padding[0]);
    char p1 = '{';
    if(memcmp(&p0, &p1, sizeof(char)) != 0){
      SNPYIO_ERROR("dict_and_padding (%s) does not start with '{'\n", dict_and_padding);
      return -1;
    }
  }
  // e: end
  // checking from the last, since padding only contains
  //   space 0x20 and 0x0a, which is much safer than
  //   walking through all dicts which have much richer info
  size_t e = 0;
  {
    for(size_t i = header_len-1; i > 0; i--){
      // use uint8_t since padding is essentially binary
      //  rather than ascii
      const uint8_t p0 = dict_and_padding[i];
      const char p1 = '}';
      if(memcmp(&p0, &p1, sizeof(char)) == 0){
        e = i;
        break;
      }
      if(i == 1){
        SNPYIO_ERROR("dict_and_padding (%s) is empty\n", dict_and_padding);
        return -1;
      }
    }
  }
  // flag dict_and_padding to decide
  //   which part should be / should not be extracted
  // "meaningless spaces" (spaces outside quotations) are de-flagged
  // e.g., meaningless spaces are
  //   {'descr': VALUE, 'fortran_order': VALUE, 'shape': VALUE}
  //            ^      ^                ^      ^        ^
  size_t n_chars_dict = 0;
  bool *is_dict = memory_calloc(e-s+1, sizeof(bool));
  {
    bool is_inside_s_quotations = false;
    bool is_inside_d_quotations = false;
    for(size_t i = s; i <= e; i++){
      const uint8_t c = dict_and_padding[i];
      // check whether we are inside a pair of single quotations
      if(c == (uint8_t)('\'')){
        is_inside_s_quotations = !is_inside_s_quotations;
      }
      // check whether we are inside a pair of double quotations
      if(c == (uint8_t)('"')){
        is_inside_d_quotations = !is_inside_d_quotations;
      }
      if(c != (uint8_t)(' ')){
        // if "c" is not space, the information is meaningful
        //   as a member of the dictionary
        n_chars_dict++;
        is_dict[i-s] = true;
      }else{
        if(is_inside_s_quotations || is_inside_d_quotations){
          // key can contain spaces (not recommended though)
          // these spaces should NOT be removed
          n_chars_dict++;
          is_dict[i-s] = true;
        }else{
          // "c" is a space and outside pair of quotations,
          //   indicating this space is meaningless
          is_dict[i-s] = false;
        }
      }
    }
  }
  // copy flagged part to dict
  // +1: null character
  *dict = memory_calloc(n_chars_dict+1, sizeof(char));
  for(size_t i = s, j = 0; i <= e; i++){
    if(is_dict[i-s]){
      (*dict)[j] = (char)(dict_and_padding[i]);
      j++;
    }
  }
  memory_free(is_dict);
  SNPYIO_LOGGING("dict: %s\n", *dict);
  return 0;
}

static int find_dict_value(const char key[], char **val, const char *dict){
  /*
   * dictionary consists of pairs of "key" and "val"
   * this function extracts the "val" of the specified "key"
   */
  if(key == NULL){
    SNPYIO_ERROR("key is NULL\n");
    return -1;
  }
  if(dict == NULL){
    SNPYIO_ERROR("dict is NULL\n");
    return -1;
  }
  // number of characters
  const size_t n_chars_key  = strlen(key);
  const size_t n_chars_dict = strlen(dict);
  // do not accept zero-length strings
  if(n_chars_key == 0){
    SNPYIO_ERROR("key is empty\n");
    return -1;
  }
  if(n_chars_dict == 0){
    SNPYIO_ERROR("dict is empty\n");
    return -1;
  }
  // 1. find key locations (start and end)
  size_t key_s = 0;
  {
    size_t location;
    const int retval = find_pattern(
        &location,
        dict,
        sizeof(char)*n_chars_dict,
        key,
        sizeof(char)*n_chars_key
    );
    if(retval < 0){
      SNPYIO_ERROR("key (%s) not found in dict (%s)\n", key, dict);
      return -1;
    }
    key_s = location/sizeof(char);
  }
  // end is easy since we know the length of the key
  const size_t key_e = key_s + n_chars_key - 1;
  // 2. find val locations (start and end)
  // start: end of the key + ":",
  //   assuming no spaces
  //   ... key:value ...
  //         ^ ^
  const size_t val_s = key_e + 2;
  // parse dict to figure out the end location of value,
  //   which is based on the fact "python dict is delimited by ','"
  // we might not be able to find ","
  //   if the pair of key / value locates at the last of the given dict
  // in this case "}" is used to notice we are at the end of dict
  // note that "," might exist inside values
  // thus we need to regard
  //   "," ONLY outside pair of brackets as delimiters
  size_t val_e;
  int bracket_level_r = 0; // ( and )
  int bracket_level_s = 0; // [ and ]
  for(val_e = val_s; val_e < n_chars_dict; val_e++){
    const char c = dict[val_e];
    if(c == '('){
      bracket_level_r++;
    }
    if(c == ')'){
      bracket_level_r--;
    }
    if(c == '['){
      bracket_level_s++;
    }
    if(c == ']'){
      bracket_level_s--;
    }
    bool is_outside_brackets_r = false;
    bool is_outside_brackets_s = false;
    if(bracket_level_r == 0){
      is_outside_brackets_r = true;
    }else if(bracket_level_r < 0){
      // indicating ) is found but ( could not found in front of it
      SNPYIO_ERROR("strange dict (%s), ')' found but corresponding '(' not found in front of it\n", dict);
      return -1;
    }
    if(bracket_level_s == 0){
      is_outside_brackets_s = true;
    }else if(bracket_level_s < 0){
      // indicating ] is found but [ could not found in front of it
      SNPYIO_ERROR("strange dict (%s), ']' found but corresponding '[' not found in front of it\n", dict);
      return -1;
    }
    // we are at the end of val if "," is found outside all brackets
    if(c == ','){
      if(is_outside_brackets_r && is_outside_brackets_s){
        // end of val should be just before found ','
        val_e -= 1;
        break;
      }
    }
    // we are at the end of dict if "}" is found
    if(c == '}'){
      // end of val should be just before '}'
      val_e -= 1;
      break;
    }
  }
  // 3. now we know where val starts and terminates, so extract it
  const size_t n_chars_val = val_e-val_s+1;
  // +1: null character
  *val = memory_calloc(n_chars_val+1, sizeof(char));
  memcpy(*val, dict+val_s, sizeof(char)*n_chars_val);
  return 0;
}

static int extract_shape(size_t *ndim, size_t **shape, const char *val){
  /*
   * parse given python tuple "val" and obtain shape of data,
   *   which is necessary to be parsed to re-construct the data
   */
  if(val == NULL){
    SNPYIO_ERROR("val is NULL\n");
    return -1;
  }
  // load shape and number of dimensions
  *ndim = 0;
  *shape = NULL;
  {
    // copy "val" to a buffer "str" after removing parentheses
    const size_t n_chars_val = strlen(val);
    if(n_chars_val < 2){
      SNPYIO_ERROR("number of characters of val: %s, which is a tuple, is less than 2\n", val);
      return -1;
    }
    // with null character (+1), no parentheses (-2): -1
    char *str = memory_calloc(n_chars_val-1, sizeof(char));
    // +1: skip first "("
    memcpy(str, val+1, sizeof(char)*(n_chars_val-2));
    // parse "str" to know ndim and shape
    // e.g.,
    //   <empty> -> ndim = 0, shape = NULL
    //   314,    -> ndim = 1, shape[0] = 314
    //   31,4    -> ndim = 2, shape[0] = 31, shape[1] = 4
    //   3,1,4,  -> ndim = 3, shape[0] = 3, shape[1] = 1, shape[2] = 4
    // since ndim is unknown for now, store result as a linked list
    node_t *shape_ = NULL;
    for(size_t loc = 0;;){
      // tokenise
      // current location is already the end of string
      if(str[loc] == null_char){
        break;
      }
      // set pointer to the current location
      char *buf = str + loc;
      // find the next delimiter and replace it with null character
      for(char *c = buf; *c != null_char; c++, loc++){
        if(*c == ','){
          // replace delimiter with null character
          *c = null_char;
          // next starting point is one character ahead
          loc++;
          break;
        }
      }
      // try to interpret characters as a number
      const long long llnum = strtoll(buf, NULL, 10);
      if(llnum <= 0){
        SNPYIO_ERROR("non-positive shape: %lld\n", llnum);
        return -1;
      }
      const size_t num = (size_t)llnum;
      node_t *new_node = memory_calloc(1, sizeof(node_t));
      new_node->ptr = memory_calloc(1, sizeof(size_t));
      memcpy(new_node->ptr, &num, sizeof(size_t));
      new_node->next = shape_;
      shape_ = new_node;
      (*ndim)++;
    }
    // clean-up buffer
    memory_free(str);
    // convert "shape_" to normal pointer "shape"
    if(*ndim > 0){
      *shape = memory_calloc(*ndim, sizeof(size_t));
      for(size_t n = 0; n < *ndim; n++){
        (*shape)[*ndim-n-1] = *((size_t *)shape_->ptr);
        node_t *next = shape_->next;
        memory_free(shape_->ptr);
        memory_free(shape_);
        shape_ = next;
      }
    }
  }
  SNPYIO_LOGGING("ndim: %zu\n", *ndim);
  for(size_t i = 0; i < *ndim; i++){
    SNPYIO_LOGGING("shape[%zu]: %zu\n", i, (*shape)[i]);
  }
  return 0;
}

static int extract_dtype(char **dtype, const char *val){
  /*
   * find a key 'descr' and extract its value
   * return obtained value directly since it is enough
   */
  if(val == NULL){
    SNPYIO_ERROR("val is NULL\n");
    return -1;
  }
  const size_t n_chars_val = strlen(val);
  *dtype = memory_calloc(n_chars_val+1, sizeof(char));
  memcpy(*dtype, val, sizeof(char)*n_chars_val);
  (*dtype)[n_chars_val] = null_char;
  SNPYIO_LOGGING("dtype: %s\n", *dtype);
  return 0;
}

static int extract_is_fortran_order(bool *is_fortran_order, const char *val){
  /*
   * find a key 'fortran_order' and extract its value
   * check whether it is "True" or "False",
   *   convert it to boolean and return
   */
  if(val == NULL){
    SNPYIO_ERROR("val is NULL\n");
    return -1;
  }
  // try to find "True"
  bool true_is_found = false;
  {
    const char pattern[] = {"True"};
    size_t location;
    const size_t n_chars_val = strlen(val);
    const size_t n_chars_pattern = strlen(pattern);
    const int retval = find_pattern(
        &location,
        val,
        sizeof(char)*n_chars_val,
        pattern,
        sizeof(char)*n_chars_pattern
    );
    if(retval < 0){
      true_is_found = false;
    }else{
      true_is_found = true;
    }
  }
  // try to find "False"
  bool false_is_found = false;
  {
    const char pattern[] = {"False"};
    size_t location;
    const size_t n_chars_val = strlen(val);
    const size_t n_chars_pattern = strlen(pattern);
    const int retval = find_pattern(
        &location,
        val,
        sizeof(char)*n_chars_val,
        pattern,
        sizeof(char)*n_chars_pattern
    );
    if(retval < 0){
      false_is_found = false;
    }else{
      false_is_found = true;
    }
  }
  // check two results and decide final outcome
  if(true_is_found && false_is_found){
    SNPYIO_ERROR("both True and False are found: %s\n", val);
    return -1;
  }else if((!true_is_found) && (!false_is_found)){
    SNPYIO_ERROR("none of True and False are found: %s\n", val);
    return -1;
  }else if(true_is_found){
    *is_fortran_order = true;
  }else{
    *is_fortran_order = false;
  }
  SNPYIO_LOGGING("is_fortran_order: %u\n", *is_fortran_order);
  return 0;
}

size_t snpyio_r_header(size_t *ndim, size_t **shape, char **dtype, bool *is_fortran_order, FILE *fp){
  uint8_t major_version, minor_version;
  size_t header_len, header_size;
  uint8_t *dict_and_padding = NULL;
  // check file pointer is valid (at least non-NULL)
  if(fp == NULL){
    SNPYIO_ERROR("fp is NULL\n");
    SNPYIO_ERROR("check file is properly opened\n");
    error_handlings();
    return 0;
  }
  /* step 1: load header from file and move file pointer forward */
  // load header to get / sanitise input and move file pointer forward
  // also the total header size "header_size" is calculated
  //   by summing up the size of each data "buf_size"
  {
    header_size = 0;
    size_t buf_size;
    /* ! check magic string ! 5 ! */
    if(load_magic_string(&buf_size, fp) < 0){
      error_handlings();
      return 0;
    }
    header_size += buf_size;
    /* ! check versions ! 5 ! */
    if(load_versions(&major_version, &minor_version, &buf_size, fp) < 0){
      error_handlings();
      return 0;
    }
    header_size += buf_size;
    /* ! load HEADER_LEN ! 5 ! */
    if(load_header_len(&header_len, &buf_size, major_version, fp) < 0){
      error_handlings();
      return 0;
    }
    header_size += buf_size;
    /* ! load dictionary and padding ! 5 ! */
    if(load_dict_and_padding(&dict_and_padding, &buf_size, header_len, fp) < 0){
      error_handlings();
      return 0;
    }
    header_size += buf_size;
  }
  /* ! step 2: extract dictionary ! 10 ! */
  // extract dict from dict + padding
  // also non-crutial spaces (spaces outside quotations) are eliminated
  //   e.g., {'descr': '<i4','fortran_order': False,'shape': (3, 5, )}
  //      -> {'descr':'<i4','fortran_order':False,'shape':(3,5,)}
  char *dict = NULL;
  if(extract_dict(&dict, dict_and_padding, header_len) < 0){
    error_handlings();
    return 0;
  }
  memory_free(dict_and_padding);
  /* step 3: extract information which are needed to reconstruct array */
  /* in particular, shape, datatype, and memory order of the array */
  // shape
  {
    /* ! extract value of shape ! 10 ! */
    char *val = NULL;
    if(find_dict_value("'shape'", &val, dict) < 0){
      error_handlings();
      return 0;
    }
    if(extract_shape(ndim, shape, val) < 0){
      error_handlings();
      return 0;
    }
    memory_free(val);
  }
  // descr (data type)
  {
    /* ! extract value of descr ! 10 ! */
    char *val = NULL;
    if(find_dict_value("'descr'", &val, dict) < 0){
      error_handlings();
      return 0;
    }
    if(extract_dtype(dtype, val) < 0){
      error_handlings();
      return 0;
    }
    memory_free(val);
  }
  // fortran order (memory order)
  {
    /* ! extract value of fortran_order ! 10 ! */
    char *val = NULL;
    if(find_dict_value("'fortran_order'", &val, dict) < 0){
      error_handlings();
      return 0;
    }
    if(extract_is_fortran_order(is_fortran_order, val) < 0){
      error_handlings();
      return 0;
    }
    memory_free(val);
  }
  // clean-up buffer
  memory_free(dict);
  // detach memories from memory pool to use outside
  // user is responsible for deallocating them
  detach_list(*shape);
  detach_list(*dtype);
  return header_size;
}

/* writer */

static int create_descr_value(char **value, const char dtype[]){
  /*
   * create a value of a dictionary key: "descr",
   *   containing user-specified dtype
   * for now this function just copies the input
   *   after sanitising a bit
   * this function is kept here, however, for consistency
   *   and future extensions
   *
   * NOTE: user is responsible for giving a proper datatype
   * NOTE: this function allocates memory for value,
   *   which should be deallocated afterwards by the caller
   */
  if(dtype == NULL){
    SNPYIO_ERROR("dtype is NULL\n");
    return -1;
  }
  const size_t n_chars = strlen(dtype);
  // check empty string,
  //   since empty datatype is obviously strange
  if(n_chars == 0){
    SNPYIO_ERROR("given dtype is empty\n");
    return -1;
  }
  // +1: null character
  *value = memory_calloc(n_chars+1, sizeof(char));
  // copy dtype
  memcpy(*value, dtype, sizeof(char)*n_chars);
  (*value)[n_chars] = null_char;
  return 0;
}

static int create_fortran_order_value(char **value, const bool is_fortran_order){
  /*
   * Create a value of a dictionary key: "fortran_order",
   *   which is True or False
   *
   * NOTE: this function allocates memory for value,
   *   which should be deallocated afterwards by the caller
   */
  if(is_fortran_order){
    const char string[] = {"True"};
    const size_t n_chars = strlen(string);
    // "True" + null character
    *value = memory_calloc(n_chars+1, sizeof(char));
    memcpy(*value, string, sizeof(char)*n_chars);
    (*value)[n_chars] = null_char;
  }else{
    const char string[] = {"False"};
    const size_t n_chars = strlen(string);
    // "False" + null character
    *value = memory_calloc(n_chars+1, sizeof(char));
    memcpy(*value, string, sizeof(char)*n_chars);
    (*value)[n_chars] = null_char;
  }
  return 0;
}

static int create_shape_value(char **value, const size_t ndim, const size_t *shape){
  /*
   * Create a value of a dictionary key: "shape"
   * Examples:
   * 0D array: ndim = 0, *dims = NULL  -> ()
   * 1D array: ndim = 1, *dims = {5}   -> (5,)
   * 2D array: ndim = 2, *dims = {5,2} -> (5,2,)
   * from left to right,
   *   from outer (less contiguous) to inner (contiguous)
   *
   * NOTE: this function allocates memory for value,
   *   which should be deallocated afterwards by the caller
   *
   * WARNING: the user of this function is responsible for
   *   allocating memory for "shape" properly
   */
  // check whether shape contains only non-negative integers
  // numpy does accpect tuples containing 0 as shape,
  //   but they are not useful in most cases
  for(size_t i = 0; i < ndim; i++){
    if(shape[i] <= 0){
      SNPYIO_ERROR("shape[%zu] should be positive\n", i);
      return -1;
    }
  }
  // 1. count number of digits (e.g., 5: 1 digit, 15: 2 digits)
  //   of shape in each direction
  size_t *n_digits = memory_calloc(ndim, sizeof(size_t));
  for(size_t i = 0; i < ndim; i++){
    // simple way to compute digits,
    //   dividing by 10 as many times as possible
    size_t num = shape[i];
    n_digits[i] = 1;
    while(num /= 10){
      n_digits[i]++;
    }
  }
  // 2. compute total number of characters
  //   i.e., memory size to be allocated
  size_t n_chars = 3; // at least "(", ")", and null character exist
  for(size_t i = 0; i < ndim; i++){
    // number of digits in i-th direction
    // with comma (+1)
    n_chars += n_digits[i]+1;
  }
  // 3. allocate memory and assign values
  *value = memory_calloc(n_chars, sizeof(char));
  for(size_t i = 0, offset = 1; i < ndim; i++){
    // assign size of the array in each direction to "value"
    //   after converting the integer to characters, e.g., 128 -> "128"
    const size_t n_digit = n_digits[i];
    // + "," and null character
    char *buf = memory_calloc(n_digit+2, sizeof(char));
    // including ","
    {
      const int retval = snprintf(buf, n_digit+2, "%zu,", shape[i]);
      if(retval != (int)(n_digit+1)){
        SNPYIO_ERROR("snprintf failed (expected: %d, given: %d)\n", retval, (int)(n_digit+1));
        return -1;
      }
    }
    // copy result excluding null character
    memcpy((*value)+offset, buf, sizeof(char)*(n_digit+1));
    offset += n_digit+1;
    memory_free(buf);
  }
  // first character is a parenthesis
  (*value)[        0] = '(';
  // last-1 character is a parenthesis
  (*value)[n_chars-2] = ')';
  // last character is null character
  (*value)[n_chars-1] = null_char;
  // clean-up buffer
  memory_free(n_digits);
  return 0;
}

static int create_dict(char **dict, size_t *n_dict, const size_t ndim, const size_t *shape, const char dtype[], const bool is_fortran_order){
  /*
   * "dict" contains information which is necessary to recover the original array,
   *   1. datatype, 2. memory ordering, and 3. shape of the data
   * It is a python-like dictionary, whose structure is a pair of key and value:
   * --- -------------- -----------------------------------------------------------------
   *  1.  descr         Datatype, e.g., '<f8', 'float64'
   *  2.  fortran_order Memory order, True or False (usually False)
   *  3.  shape         A tuple having number of elements in each direction, e.g., (5,2,)
   * See corresponding function for how they are created
   * Also the number of elements of the dict is returned (to be consistent with "create_padding")
   */
  // keys, which are completely fixed
  char descr_key[]         = {"'descr'"};
  char fortran_order_key[] = {"'fortran_order'"};
  char shape_key[]         = {"'shape'"};
  // values, which depend on inputs
  char *descr_value         = NULL;
  char *fortran_order_value = NULL;
  char *shape_value         = NULL;
  /* 1. create dictionary values,
   *   in which inputs are evaluated and sanitised
   */
  /* ! create value of data type ! 4 ! */
  if(create_descr_value(&descr_value, dtype) != 0){
    SNPYIO_ERROR("create_descr_value failed\n");
    return -1;
  }
  /* ! create value of memory order ! 4 ! */
  if(create_fortran_order_value(&fortran_order_value, is_fortran_order) != 0){
    SNPYIO_ERROR("create_fortran_order_value failed\n");
    return -1;
  }
  /* ! create value of data sizes ! 4 ! */
  if(create_shape_value(&shape_value, ndim, shape) != 0){
    SNPYIO_ERROR("create_shape_value failed\n");
    return -1;
  }
  /* ! 2. assign all elements (strings) which compose dict ! 20 ! */
  const size_t n_elements_dict = 13;
  char **elements = memory_calloc(n_elements_dict, sizeof(char *));
  // initial wave bracket
  elements[ 0] = "{";
  // 'descr':descr_value
  elements[ 1] = (char *)descr_key;
  elements[ 2] = ":";
  elements[ 3] = (char *)descr_value;
  elements[ 4] = ",";
  // 'fortran_order':fortran_order_value
  elements[ 5] = (char *)fortran_order_key;
  elements[ 6] = ":";
  elements[ 7] = (char *)fortran_order_value;
  elements[ 8] = ",";
  // 'shape':shape_value
  elements[ 9] = (char *)shape_key;
  elements[10] = ":";
  elements[11] = (char *)shape_value;
  // final wave bracket
  elements[12] = "}";
  // 3. check total number of characters of
  //   {'descr':VALUE,'fortran_order':VALUE,'shape':VALUE}
  //   to allocate dict
  // NOTE: n_chars_dict is the number of characters of dict
  //   INCLUDING the last null character, while n_dict = strlen(dict),
  //   EXCLUDING the last null character.
  //   Thus n_dict = n_chars_dict - 1
  size_t n_chars_dict = 0;
  for(size_t i = 0; i < n_elements_dict; i++){
    // check each element and sum up its number of characters
    const size_t n_chars = strlen(elements[i]);
    n_chars_dict += n_chars;
  }
  // last null character
  n_chars_dict += 1;
  // 4. allocate dict and assign above "elements"
  *dict = memory_calloc(n_chars_dict, sizeof(char));
  for(size_t i = 0, offset = 0; i < n_elements_dict; i++){
    const size_t n_chars = strlen(elements[i]);
    memcpy((*dict)+offset, elements[i], sizeof(char)*n_chars);
    offset += n_chars;
  }
  (*dict)[n_chars_dict-1] = null_char;
  // clean-up all working memories
  memory_free(descr_value);
  memory_free(fortran_order_value);
  memory_free(shape_value);
  memory_free(elements);
  // as the length of "dict", use length WITHOUT null character,
  // i.e. strlen(*dict)
  *n_dict = strlen(*dict);
  SNPYIO_LOGGING("dict: %s\n", *dict);
  SNPYIO_LOGGING("size: %zu\n", *n_dict);
  return 0;
}

static int create_padding(uint8_t **padding, size_t *n_padding, uint8_t *major_version, const size_t n_dict){
  /*
   * The following relation holds for the header size
   *   size_header =
   *     + sizeof(magic string)      (= 6            bytes)
   *     + sizeof(major_version)     (= 1            byte )
   *     + sizeof(minor_version)     (= 1            byte )
   *     + sizeof(header_len)        (= 2 or 4       bytes)
   *     + sizeof(char)*strlen(dict) (= strlen(dict) bytes)
   *     + sizeof(uint8_t)*n_padding (= n_padding    bytes)
   *   is divisible by 64
   * Definitely this is not generally true, and we need some paddings
   *   consisting of some (0 or more) spaces ' ' and one newline '\n',
   *   whose length (number of elements) is returned
   */
  /* ! size of each element is computed ! 5 ! */
  const size_t n_magic_string = strlen(magic_string);
  const size_t size_magic_string  = sizeof(char)*n_magic_string;
  const size_t size_major_version = sizeof(uint8_t);
  const size_t size_minor_version = sizeof(uint8_t);
  const size_t size_dict          = sizeof(char)*n_dict;
  /* ! reject too large dict ! 4 ! */
  if(size_dict > UINT_MAX-64){
    SNPYIO_ERROR("size of dictionary is huge (%zu)\n", size_dict);
    return -1;
  }
  /* ! decide major version and datatype of HEADER_LEN ! 11 ! */
  // large portion of the header is occupied by dict
  // so check dict size, and if it is larger than USHRT_MAX-64,
  //   use major_version = 2
  size_t size_header_len;
  if(size_dict > USHRT_MAX-64){
    *major_version = 2;
    size_header_len = sizeof(uint32_t);
  }else{
    *major_version = 1;
    size_header_len = sizeof(uint16_t);
  }
  /* ! compute size of all data except padding ! 6 ! */
  const size_t size_except_padding =
    +size_magic_string
    +size_major_version
    +size_minor_version
    +size_header_len
    +size_dict;
  /* ! decide total size of the header, which should be 64 x N ! 8 ! */
  // increase total size by 64 until becoming larger than size_except_padding
  // NOTE: size_padding == 0 is NOT allowed since '\n' is necessary at the end
  //   thus the condition to continue loop is "<=", not "<"
  size_t size_header = 0;
  while(size_header <= size_except_padding){
    size_header += 64;
  }
  const size_t size_padding = size_header-size_except_padding;
  /* ! create padding ! 6 ! */
  *n_padding = size_padding/sizeof(uint8_t);
  *padding = memory_calloc(*n_padding, sizeof(uint8_t));
  // many ' 's: 0x20
  memset(*padding, 0x20, sizeof(uint8_t)*(*n_padding-1));
  // last '\n': 0x0a
  (*padding)[*n_padding-1] = 0x0a;
  SNPYIO_LOGGING("padding, size: %zu\n", *n_padding);
  return 0;
}

static int create_header_len(uint8_t **header_len, size_t *n_header_len, const uint8_t major_version, const size_t n_dict, const size_t n_padding){
  /*
   * In short, HEADER_LEN = n_dict + n_padding,
   * which should be written as a little-endian form
   *   (irrespective to the architecture)
   */
  /* ! reject too large dict / padding sizes ! 10 ! */
  // Here "too large" means header size (not data size)
  //   is larger than approx. 2GB, which would not happen normally
  if(n_dict >= UINT_MAX/2){
    SNPYIO_ERROR("dictionary size is huge (%zu)\n", n_dict);
    return -1;
  }
  if(n_padding >= UINT_MAX/2){
    SNPYIO_ERROR("padding size is huge (%zu)\n", n_padding);
    return -1;
  }
  /* ! compute header_len and store as an array of uint8_t ! 15 ! */
  if(major_version == 2){
    // major version 2, use uint32_t to store header_len
    const uint32_t header_len_uint32_t = (uint32_t)n_dict+(uint32_t)n_padding;
    *n_header_len = sizeof(uint32_t)/sizeof(uint8_t);
    *header_len = memory_calloc(*n_header_len, sizeof(uint8_t));
    memcpy(*header_len, &header_len_uint32_t, *n_header_len);
    SNPYIO_LOGGING("header_len (uint32_t): %u\n", header_len_uint32_t);
  }else{
    // major version 1, use uint16_t to store header_len
    const uint16_t header_len_uint16_t = (uint16_t)n_dict+(uint16_t)n_padding;
    *n_header_len = sizeof(uint16_t)/sizeof(uint8_t);
    *header_len = memory_calloc(*n_header_len, sizeof(uint8_t));
    memcpy(*header_len, &header_len_uint16_t, *n_header_len);
    SNPYIO_LOGGING("header_len (uint16_t): %hu\n", header_len_uint16_t);
  }
  /* ! convert endian of buffer which will be written if needed ! 6 ! */
  if(is_big_endian()){
    if(convert_endian(header_len, sizeof(*header_len)) != 0){
      SNPYIO_ERROR("convert_endian failed\n");
      return -1;
    }
  }
  return 0;
}

size_t snpyio_w_header(const size_t ndim, const size_t *shape, const char dtype[], const bool is_fortran_order, FILE *fp){
  /*
   * From input information (ndim, shape, dtype, is_fortran_order),
   *   create and output "header", which is enough and sufficient information consists of a *.npy file
   * NPY file format is defined here:
   *   https://numpy.org/devdocs/reference/generated/numpy.lib.format.html#format-version-1-0
   *
   * The size of the "header", defined as "header_size",
   *   is returned, which could be useful by the user
   *
   * Header of a npy file contains these 6 data (n_datasets):
   *       -----NAME----- --TYPE-- -# OF ELEMENTS- -------SIZE-------
   * No. 0 magic_string   char     6               6 bytes
   * No. 1 major_version  uint8_t  1               1 byte
   * No. 2 minor_version  uint8_t  1               1 byte
   * No. 3 header_len     uint8_t  n_header_len    n_header_len bytes
   * No. 4 dict           char     n_dict          n_dict bytes
   * No. 5 padding        uint8_t  n_padding       n_padding bytes
   *
   * See below and corresponding function for details of each member
   */
  // check file pointer is valid (at least non-NULL)
  if(fp == NULL){
    SNPYIO_ERROR("fp is NULL, check file is properly opened\n");
    return 0;
  }
  /* prepare all 6 datasets */
  /* ! magic_string ! */
  const size_t n_magic_string = strlen(magic_string);
  /* ! minor_version, always 0 ! 1 ! */
  const uint8_t minor_version = 0;
  /* ! dictionary (and its size) ! 6 ! */
  char *dict = NULL;
  size_t n_dict;
  if(create_dict(&dict, &n_dict, ndim, shape, dtype, is_fortran_order) != 0){
    error_handlings();
    return 0;
  }
  /* ! major_version and padding (and its size) ! 7 ! */
  uint8_t major_version;
  uint8_t *padding = NULL;
  size_t n_padding;
  if(create_padding(&padding, &n_padding, &major_version, n_dict) != 0){
    error_handlings();
    return 0;
  }
  /* ! comptue header_len ! 6 ! */
  uint8_t *header_len = NULL;
  size_t n_header_len;
  if(create_header_len(&header_len, &n_header_len, major_version, n_dict, n_padding) != 0){
    error_handlings();
    return 0;
  }
  /* dump all information to a buffer "header" and compute total size "header_size" */
  uint8_t *header = NULL;
  size_t header_size;
  size_t header_nitems;
  {
    const size_t n_datasets = 6;
    size_t *sizes   = NULL;
    size_t *offsets = NULL;
    sizes   = memory_calloc(n_datasets, sizeof(size_t));
    offsets = memory_calloc(n_datasets, sizeof(size_t));
    sizes[0] = sizeof(   char)*n_magic_string;
    sizes[1] = sizeof(uint8_t);
    sizes[2] = sizeof(uint8_t);
    sizes[3] = sizeof(uint8_t)*n_header_len;
    sizes[4] = sizeof(   char)*n_dict;
    sizes[5] = sizeof(uint8_t)*n_padding;
    // total size
    header_size = 0;
    for(uint8_t i = 0; i < n_datasets; i++){
      header_size += sizes[i];
    }
    // offsets
    offsets[0] = 0;
    for(uint8_t i = 1; i < n_datasets; i++){
      offsets[i] = offsets[i-1]+sizes[i-1];
    }
    // allocate buffer to store whole header
    header_nitems = header_size/sizeof(uint8_t);
    header = memory_calloc(header_nitems, sizeof(uint8_t));
    // write all information to a buffer "header"
    memcpy(header+offsets[0], magic_string,   sizes[0]);
    memcpy(header+offsets[1], &major_version, sizes[1]);
    memcpy(header+offsets[2], &minor_version, sizes[2]);
    memcpy(header+offsets[3], header_len,     sizes[3]);
    memcpy(header+offsets[4], dict,           sizes[4]);
    memcpy(header+offsets[5], padding,        sizes[5]);
    // clean-up buffers
    memory_free(sizes);
    memory_free(offsets);
  }
  SNPYIO_LOGGING("header_size: %zu\n", header_size);
  /* write to the given file stream */
  {
    const size_t retval = fwrite(header, sizeof(uint8_t), header_nitems, fp);
    if(retval != header_nitems){
      SNPYIO_ERROR("fwrite failed (%zu items written, while %zu items given)\n", retval, header_nitems);
      error_handlings();
      return 0;
    }
  }
  // clean-up all buffers
  memory_free(dict);
  memory_free(padding);
  memory_free(header_len);
  memory_free(header);
  return header_size;
}

#undef SNPYIO_MESSAGE
#undef SNPYIO_LOGGING
#undef SNPYIO_ERROR
#undef SNPYIO_FATAL


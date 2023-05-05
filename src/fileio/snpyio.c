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

// check 1 byte is 8 bits
#if 8 != CHAR_BIT
#error "CHAR_BIT is not 8"
#endif

// all npy files should start from this magic string | 1
static const char magic_string[] = {"\x93NUMPY"};

// header size is divisible by this number
static const size_t header_block_size = 64;

// null character
static const char null_char = '\0';

/* general messages for logginig and error handling */
#define SNPYIO_MESSAGE(stream, level, ...){    \
  fprintf(stream, "[NPYIO %s]\n", level);      \
  fprintf(stream, "    line: %d\n", __LINE__); \
  fprintf(stream, "    func: %s\n", __func__); \
  fprintf(stream, "    ");                     \
  fprintf(stream, __VA_ARGS__);                \
}

#if defined(SNPYIO_ENABLE_LOGGING)
#define SNPYIO_LOGGING(...){              \
  FILE *fp = fopen("snpyio.log", "a");    \
  SNPYIO_MESSAGE(fp, "LOG", __VA_ARGS__); \
  fclose(fp);                             \
}
#else
#define SNPYIO_LOGGING(...)
#endif

#define SNPYIO_ERROR(...){                      \
  SNPYIO_MESSAGE(stderr, "ERROR", __VA_ARGS__); \
}

#define SNPYIO_FATAL(...){                      \
  SNPYIO_MESSAGE(stderr, "FATAL", __VA_ARGS__); \
  exit(EXIT_FAILURE);                           \
}

// singly-linked list
typedef struct node_t_ {
  void *ptr;
  struct node_t_ *next;
} node_t;

/* memory manager */

// memory pool storing all pointers to the allocated buffers in this library
static node_t *memory = NULL;

// kernel memory allocator
static void *kernel_calloc(const size_t count, const size_t size){
  if(SIZE_MAX / size < count){
    SNPYIO_FATAL("memory allocation failed: calloc(%zu, %zu)\n", count, size);
  }
  void *ptr = calloc(count, size);
  if(NULL == ptr){
    SNPYIO_FATAL("memory allocation failed: calloc(%zu, %zu)\n", count, size);
  }
  return ptr;
}

// kernel memory deallocator
static void kernel_free(void *ptr){
  free(ptr);
}

// general memory allocator
static void *memory_calloc(const size_t count, const size_t size){
  void *ptr = kernel_calloc(count, size);
  node_t *node = kernel_calloc(1, sizeof(node_t));
  node->ptr = ptr;
  node->next = memory;
  memory = node;
  return ptr;
}

// deallocate a node holding the given pointer
static int detach_list(const void *ptr){
  // NOTE: ptr is not freed while the node is freed
  if(NULL == ptr) return 1;
  node_t **node = &memory;
  while(*node){
    if(ptr == (*node)->ptr){
      // keep next node
      node_t *node_next = (*node)->next;
      // clean-up current node
      kernel_free(*node);
      // update connection
      *node = node_next;
      return 0;
    }
    node = &((*node)->next);
  }
  return 1;
}

// general memory deallocator
static void memory_free(void *ptr){
  detach_list(ptr);
  kernel_free(ptr);
}

// free all buffers attached to the memory pool
static void error_handlings(void){
  while(memory){
    node_t *next = memory->next;
    kernel_free(memory->ptr);
    kernel_free(memory);
    memory = next;
  }
}

/* auxiliary functions */

static bool is_big_endian(void){
  const uint16_t val = 1 << 8;
  return (bool)(((const uint8_t *)(&val))[0]);
}

static int convert_endian(void *val, const size_t size){
  const size_t n_bytes = size / sizeof(uint8_t);
  // swap and write back
  // use buffer for simplicity
  uint8_t *buf = memory_calloc(n_bytes, sizeof(uint8_t));
  for(size_t i = 0; i < n_bytes; i++){
    buf[i] = ((uint8_t *)val)[n_bytes - i - 1];
  }
  memcpy(val, buf, size);
  memory_free(buf);
  return 0;
}

static int myfread(void *ptr, const size_t size, const size_t nitems, FILE *stream){
  if(nitems != fread(ptr, size, nitems, stream)){
    SNPYIO_ERROR("fread failed\n");
    return 1;
  }
  return 0;
}

static int myfwrite(const void *ptr, const size_t size, const size_t nitems, FILE *stream){
  if(nitems != fwrite(ptr, size, nitems, stream)){
    SNPYIO_ERROR("fwrite failed\n");
    return 1;
  }
  return 0;
}

static int sanitise_fp(const FILE *fp){
  if(NULL == fp){
    SNPYIO_ERROR("fp is NULL, check file is properly opened\n");
    return 1;
  }
  return 0;
}

/* reader */

static int load_magic_string(size_t *header_size, FILE *fp){
  /*
   * all npy file should start with a magic string,
   *   which is checked here by comparing
   *   the fist 6 bytes of the file
   */
  const size_t nitems = strlen(magic_string);
  // allocate buffer and load from file
  // NOTE: file pointer is moved forward as well
  uint8_t *buf = memory_calloc(nitems, sizeof(uint8_t));
  if(0 != myfread(buf, sizeof(uint8_t), nitems, fp)) return 1;
  const size_t buf_size = sizeof(uint8_t) * nitems;
  if(0 != memcmp(buf, magic_string, buf_size)){
    SNPYIO_ERROR("file does not begin with \"\\x93NUMPY\"\n");
    return 1;
  }
  memory_free(buf);
  *header_size += buf_size;
  return 0;
}

static int load_versions(uint8_t *major_version, uint8_t *minor_version, size_t *header_size, FILE *fp){
  /*
   * check version of the file
   * 1.0, 2.0, and 3.0 are accepted, while others are rejected
   */
  // load from file
  // NOTE: file pointer is moved forward as well
  if(0 != myfread(major_version, sizeof(uint8_t), 1, fp)) return 1;
  if(0 != myfread(minor_version, sizeof(uint8_t), 1, fp)) return 1;
  // check version 1.x or 2.x or 3.x
  if(1 != *major_version && 2 != *major_version && 3 != *major_version){
    SNPYIO_ERROR("major version (%hhu) should be 1, 2, or 3\n", *major_version);
    return 1;
  }
  // check version x.0
  if(0 != *minor_version){
    SNPYIO_ERROR("minor version (%hhu) should be 0\n", *minor_version);
    return 1;
  }
  SNPYIO_LOGGING("major version: %hhu\n", *major_version);
  SNPYIO_LOGGING("minor version: %hhu\n", *minor_version);
  *header_size += sizeof(uint8_t) + sizeof(uint8_t);
  return 0;
}

static int load_header_len(size_t *header_len, size_t *header_size, const size_t major_version, FILE *fp){
  /*
   * check header size of the npy file
   * in particular HEADER_LEN = len(dict) + len(padding)
   *   is loaded
   * memory size of this variable depends on the major version
   *   of the npy file, 2 bytes for major_version = 1,
   *   while 4 bytes for major_version = 2
   */
  // buffer size differs based on major_version | 1
  const size_t buf_size = 1 == major_version ? sizeof(uint16_t) : sizeof(uint32_t);
  const size_t nitems = buf_size / sizeof(uint8_t);
  // allocate buffer and load corresponding memory size from file | 2
  uint8_t *buf = memory_calloc(nitems, sizeof(uint8_t));
  if(0 != myfread(buf, sizeof(uint8_t), nitems, fp)) return 1;
  // convert endian of loaded buffer when needed | 4
  if(is_big_endian() && 0 != convert_endian(header_len, sizeof(*header_len))){
    SNPYIO_ERROR("convert_endian failed\n");
    return 1;
  }
  // interpret buffer (sequence of uint8_t) as a value having corresponding datatype | 11
  if(1 == major_version){
    // interpret as a 2-byte value
    uint16_t tmp = 0;
    memcpy(&tmp, buf, sizeof(uint8_t) * nitems);
    *header_len = (size_t)tmp;
  }else{
    // interpret as a 4-byte value
    uint32_t tmp = 0;
    memcpy(&tmp, buf, sizeof(uint8_t) * nitems);
    *header_len = (size_t)tmp;
  }
  memory_free(buf);
  SNPYIO_LOGGING("header_len: %zu\n", *header_len);
  *header_size += buf_size;
  return 0;
}

static int load_dict_and_padding(uint8_t **dict_and_padding, size_t *header_size, const size_t header_len, FILE *fp){
  /*
   * load dictionary and padding
   * loading padding is also necessary (or at least one of the easiest ways)
   *   to move file pointer "fp" forward
   */
  const size_t nitems = header_len / sizeof(uint8_t);
  *dict_and_padding = memory_calloc(nitems, sizeof(uint8_t));
  if(0 != myfread(*dict_and_padding, sizeof(uint8_t), nitems, fp)) return 1;
  *header_size += header_len;
  return 0;
}

static int extract_dict(char **dict, const uint8_t *dict_and_padding, const size_t header_len){
  /*
   * extract "dict" from "dict_and_padding",
   *   which contains dictionary and padding
   * also unnecessary spaces are removed from the original array
   *   to simplify the following procedures
   * NOTE: spaces inside quotations (e.g. dictionary key might contain spaces)
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
  const size_t s = 0;
  // confirm it just in case
  //   by checking whether the first byte is "{"
  if(0 != memcmp(dict_and_padding, (char []){"{"}, sizeof(char))){
    SNPYIO_ERROR("dict_and_padding (%s) does not start with '{'\n", dict_and_padding);
    return 1;
  }
  // e: end
  // checking from the last, since padding only contains
  //   space 0x20 and 0x0a, which is safer than
  //   walking through all dicts from the head
  //   which have much richer info
  size_t e = 0;
  for(e = header_len - 1; ; e--){
    if(0 == memcmp(dict_and_padding + e, (char []){"}"}, sizeof(char))){
      // end of wavy bracket is found
      break;
    }
    if(1 == e){
      SNPYIO_ERROR("dict_and_padding (%s) is empty\n", dict_and_padding);
      return 1;
    }
  }
  // flag dict_and_padding to decide
  //   which part should be / should not be extracted
  // "meaningless spaces" (spaces outside quotations) are de-flagged
  // e.g., the following spaces are removed here
  //   {'descr': VALUE, 'fortran_order': VALUE, 'shape': VALUE}
  //            ^      ^                ^      ^        ^
  size_t n_chars_dict = e - s + 1;
  bool *is_dict = memory_calloc(n_chars_dict, sizeof(bool));
  // default: true
  for(size_t i = 0; i < n_chars_dict; i++){
    is_dict[i] = true;
  }
  // flags to check whether inside single/double quotations
  const uint8_t quots[] = {'\'', '"'};
  bool is_inside[] = {false, false};
  for(size_t i = s; i <= e; i++){
    const uint8_t c = dict_and_padding[i];
    // check whether we are inside a pair of quotations
    for(size_t j = 0; j < 2; j++){
      if(c == quots[j]) is_inside[j] = !is_inside[j];
    }
    // if "c" is not a space, the information is meaningful
    //   as a member of the dictionary
    if(' ' != c) continue;
    // this character is a space but inside key or value
    // these spaces should NOT be removed
    if(is_inside[0] || is_inside[1]) continue;
    // this character is a space and outside pair of quotations,
    //   indicating this space is meaningless and deflagged
    n_chars_dict--;
    is_dict[i - s] = false;
  }
  // copy flagged part to dict
  // +1: null character
  *dict = memory_calloc(n_chars_dict + 1, sizeof(char));
  for(size_t i = s, j = 0; i <= e; i++){
    if(is_dict[i - s]){
      (*dict)[j++] = (char)(dict_and_padding[i]);
    }
  }
  (*dict)[n_chars_dict] = null_char;
  memory_free(is_dict);
  SNPYIO_LOGGING("dict: %s\n", *dict);
  return 0;
}

static int find_dict_value(const char key[], char **val, const char *dict){
  /*
   * dictionary consists of pairs of "key" and "val"
   * this function extracts the "val" of the specified "key"
   */
  // find the starting point of value
  // first find key
  char *val_s = strstr(dict, key);
  if(NULL == val_s){
    SNPYIO_ERROR("key (%s) not found in dict (%s)\n", key, dict);
    return 1;
  }
  // move pointer to the head of value
  // NOTE: +1 is to account for ":"
  val_s += strlen(key) + 1;
  // find ",", which is the deliminator of python dict,
  //   or "}", which is the terminator of python dict
  // NOTE: "," is used as a deliminator of python tuple
  //   so it is important not to take it as a deliminator
  //   if the parser is inside a tuple
  // initialise end of value as the starting point of it
  char *val_e = val_s;
  const uint8_t quots[] = {'\'', '"'};
  bool is_inside[] = {false, false};
  for(int bracket_level = 0; ; val_e++){
    // check whether we are inside a pair of quotations
    for(size_t j = 0; j < 2; j++){
      if(*val_e == quots[j]) is_inside[j] = !is_inside[j];
    }
    // do not check when the tokeniser is inside quotations
    if(is_inside[0] || is_inside[1]) continue;
    // otherwise check tuple
    if('(' == *val_e) bracket_level++;
    if(')' == *val_e) bracket_level--;
    if(0 > bracket_level){
      SNPYIO_ERROR("strange dict (%s), ')' found but corresponding '(' not found in front of it\n", dict);
      return 1;
    }
    if(1 < bracket_level){
      SNPYIO_ERROR("nested tuple is found, which is not supported: %s\n", dict);
      return 1;
    }
    // we are at the end of val if "," is found and outside all brackets
    if(',' == *val_e && 0 == bracket_level) break;
    // we are at the end of dict if "}" is found
    if('}' == *val_e) break;
    // should not reach, null_char is found
    if(null_char == *val_e){
      SNPYIO_ERROR("strange dict (%s), no terminator and null_char is found\n", dict);
      return 1;
    }
  }
  // end of val should be one-character-ahead of the deliminator / terminator
  val_e -= 1;
  // now we know where val starts and terminates, so extract it
  const size_t n_chars_val = (val_e + 1 - val_s) / sizeof(char);
  // +1: null character
  *val = memory_calloc(n_chars_val + 1, sizeof(char));
  memcpy(*val, val_s, sizeof(char) * n_chars_val);
  (*val)[n_chars_val] = null_char;
  return 0;
}

static int extract_shape(size_t *ndim, size_t **shape, const char *dict){
  /*
   * parse given python tuple "val" and obtain shape of data,
   *   which is necessary to be parsed to re-construct the data
   */
  char *val = NULL;
  if(0 != find_dict_value("'shape'", &val, dict)) return 1;
  // load shape and number of dimensions
  *ndim = 0;
  *shape = NULL;
  // copy "val" to a buffer "str" after removing parentheses
  const size_t n_chars_val = strlen(val);
  if(2 > n_chars_val){
    SNPYIO_ERROR("number of characters of val: %s, which is a tuple, is less than 2\n", val);
    return 1;
  }
  // with null character (+1), no parentheses (-2): -1
  char *str = memory_calloc(n_chars_val - 1, sizeof(char));
  // +1: skip first "("
  memcpy(str, val + 1, sizeof(char) * (n_chars_val - 2));
  memory_free(val);
  // parse "str" to know ndim and shape
  // e.g.,
  //   <empty> -> ndim = 0, shape = NULL
  //   314,    -> ndim = 1, shape[0] = 314
  //   31,4    -> ndim = 2, shape[0] = 31, shape[1] = 4
  //   3,1,4,  -> ndim = 3, shape[0] = 3, shape[1] = 1, shape[2] = 4
  // since ndim is unknown for now, store result as a linked list
  node_t *shape_ = NULL;
  for(size_t loc = 0; ; ){
    // tokenise
    // set pointer to the current location
    char *buf = str + loc;
    // current location is already the end of string
    if(null_char == *buf) break;
    // find the next delimiter and replace it with null character
    for(char *c = buf; null_char != *c; c++, loc++){
      if(',' == *c){
        // replace delimiter with null character
        *c = null_char;
        // next investigation starts one character ahead
        loc++;
        break;
      }
    }
    // try to interpret the NULL-terminated string as a number
    const long long llnum = strtoll(buf, NULL, 10);
    if(0 >= llnum){
      SNPYIO_ERROR("non-positive shape: %lld\n", llnum);
      return 1;
    }
    // valid shape, store result to linked list
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
  // convert "shape_" (linked list) to normal pointer "shape"
  if(0 < *ndim){
    *shape = memory_calloc(*ndim, sizeof(size_t));
    for(size_t n = 0; n < *ndim; n++){
      (*shape)[*ndim - n - 1] = *((size_t *)shape_->ptr);
      node_t *next = shape_->next;
      memory_free(shape_->ptr);
      memory_free(shape_);
      shape_ = next;
    }
  }
  SNPYIO_LOGGING("ndim: %zu\n", *ndim);
  for(size_t i = 0; i < *ndim; i++){
    SNPYIO_LOGGING("shape[%zu]: %zu\n", i, (*shape)[i]);
  }
  return 0;
}

static int extract_dtype(char **dtype, const char *dict){
  /*
   * find a key 'descr' and extract its value
   * return obtained value directly since it is enough
   */
  char *val = NULL;
  if(0 != find_dict_value("'descr'", &val, dict)) return 1;
  const size_t n_chars_val = strlen(val);
  *dtype = memory_calloc(n_chars_val + 1, sizeof(char));
  memcpy(*dtype, val, sizeof(char) * n_chars_val);
  (*dtype)[n_chars_val] = null_char;
  memory_free(val);
  SNPYIO_LOGGING("dtype: %s\n", *dtype);
  return 0;
}

static int extract_is_fortran_order(bool *is_fortran_order, const char *dict){
  /*
   * find a key 'fortran_order' and extract its value
   * check whether it is "True" or "False",
   *   convert it to a boolean value and return
   */
  char *val = NULL;
  if(0 != find_dict_value("'fortran_order'", &val, dict)) return 1;
  const bool true_is_found  = (NULL != strstr(val,  "True"));
  const bool false_is_found = (NULL != strstr(val, "False"));
  if(true_is_found && false_is_found){
    SNPYIO_ERROR("both True and False are found: %s\n", val);
    return 1;
  }
  if((!true_is_found) && (!false_is_found)){
    SNPYIO_ERROR("neither True nor False was found: %s\n", val);
    return 1;
  }
  *is_fortran_order = true_is_found ? true : false;
  memory_free(val);
  SNPYIO_LOGGING("is_fortran_order: %s\n", *is_fortran_order ? "True" : "False");
  return 0;
}

/**
 * @brief read NPY header
 * @param[out] ndim             : number of dimensions of the data set, e.g. 2
 * @param[out] shape            : number of points of the data set in each dimension, e.g. [3, 4]
 * @param[out] dtype            : data type, e.g. "'<f8'"
 * @param[out] is_fortran_order : row-major order (false) or column-major order (true)
 * @param[in]  fp               : file stream to which the header is loaded
 * @return                      : (success) loaded header size (in bytes)
 *                                (failure) 0
 */
size_t snpyio_r_header(size_t *ndim, size_t **shape, char **dtype, bool *is_fortran_order, FILE *fp){
  if(0 != sanitise_fp(fp)) goto r_err_hndl;
  // load header from file and move file pointer forward
  // also the total header size "header_size" is calculated
  //   by summing up the size of each data
  // check magic string | 4
  size_t header_size = 0;
  if(0 != load_magic_string(&header_size, fp)){
    goto r_err_hndl;
  }
  // check versions | 5
  uint8_t major_version = 1;
  uint8_t minor_version = 0;
  if(0 != load_versions(&major_version, &minor_version, &header_size, fp)){
    goto r_err_hndl;
  }
  // load HEADER_LEN | 4
  size_t header_len  = 0;
  if(0 != load_header_len(&header_len, &header_size, major_version, fp)){
    goto r_err_hndl;
  }
  // load dictionary and padding | 4
  uint8_t *dict_and_padding = NULL;
  if(0 != load_dict_and_padding(&dict_and_padding, &header_size, header_len, fp)){
    goto r_err_hndl;
  }
  // extract dictionary | 9
  // extract dict from dict + padding
  // also non-crutial spaces (spaces outside quotations) are eliminated
  //   e.g., {'descr': '<i4','fortran_order': False,'shape': (3, 5, )}
  //      -> {'descr':'<i4','fortran_order':False,'shape':(3,5,)}
  char *dict = NULL;
  if(0 != extract_dict(&dict, dict_and_padding, header_len)){
    goto r_err_hndl;
  }
  memory_free(dict_and_padding);
  // extract information which are needed to reconstruct array
  // in particular, shape, datatype, and memory order of the array
  // extract value of shape | 3
  if(0 != extract_shape(ndim, shape, dict)){
    goto r_err_hndl;
  }
  // extract value of descr | 3
  if(0 != extract_dtype(dtype, dict)){
    goto r_err_hndl;
  }
  // extract value of fortran_order | 3
  if(0 != extract_is_fortran_order(is_fortran_order, dict)){
    goto r_err_hndl;
  }
  // clean-up buffer storing dictionary
  memory_free(dict);
  // detach memories from memory pool
  // user is responsible for deallocating them
  detach_list(*shape);
  detach_list(*dtype);
  return header_size;
r_err_hndl:
  error_handlings();
  return 0;
}

/* writer */

static int sanitise_shape(const size_t ndim, const size_t *shape){
  for(size_t i = 0; i < ndim; i++){
    if(0 >= shape[i]){
      SNPYIO_ERROR("shape[%zu] should be positive\n", i);
      return 1;
    }
  }
  return 0;
}

static int sanitise_dtype(const char dtype[]){
  if(NULL == dtype){
    SNPYIO_ERROR("dtype is NULL, give a proper address which can be dereferenced\n");
    return 1;
  }
  const size_t n_chars = strlen(dtype);
  if(2 > n_chars){
    SNPYIO_ERROR("invalid dtype: %s\n", dtype);
    return 1;
  }
  const uint8_t quots[] = {'\'', '"'};
  const char *squots[] = {"single", "double"};
  for(size_t i = 0; i < 2; i++){
    if(quots[i] == dtype[0]){
      if(quots[i] == dtype[n_chars - 1]) return 0;
      SNPYIO_ERROR("dtype starts with but does not end with a %s quotation\n", squots[i]);
      return 1;
    }
  }
  SNPYIO_ERROR("dtype does not start with a single quotation nor a double quotation\n");
  return 1;
}

static int create_descr_value(char **value, const char dtype[]){
  /*
   * create a value of a dictionary key: "descr",
   *   containing user-specified dtype
   */
  const size_t n_chars = strlen(dtype);
  // +1: null character
  *value = memory_calloc(n_chars + 1, sizeof(char));
  // copy dtype
  memcpy(*value, dtype, sizeof(char) * n_chars);
  (*value)[n_chars] = null_char;
  return 0;
}

static int create_fortran_order_value(char **value, const bool is_fortran_order){
  /*
   * Create a value of a dictionary key: "fortran_order",
   *   which is True or False
   */
  const char *string = is_fortran_order ? "True" : "False";
  const size_t n_chars = strlen(string);
  *value = memory_calloc(n_chars + 1, sizeof(char));
  memcpy(*value, string, sizeof(char) * n_chars);
  (*value)[n_chars] = null_char;
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
   */
  // 1. count number of digits (e.g., 5: 1 digit, 15: 2 digits)
  //   of shape in each direction
  size_t *n_digits = memory_calloc(ndim, sizeof(size_t));
  for(size_t i = 0; i < ndim; i++){
    // simple way to compute digits,
    //   dividing by 10 as many times as possible
    size_t num = shape[i];
    n_digits[i] = 1;
    while(num /= 10) n_digits[i]++;
  }
  // 2. compute total number of characters
  //   i.e., memory size to be allocated
  size_t n_chars = 2; // at least "(", ")"
  for(size_t i = 0; i < ndim; i++){
    // number of digits in i-th direction
    // with comma (+1)
    n_chars += n_digits[i] + 1;
  }
  // 3. allocate memory and assign values
  *value = memory_calloc(n_chars + 1, sizeof(char));
  for(size_t i = 0, offset = 1; i < ndim; i++){
    // assign size of the array in each direction to "value"
    //   after converting the integer to characters, e.g., 128 -> "128"
    const int n_digit = n_digits[i];
    // + "," and null character
    char *buf = memory_calloc(n_digit + 2, sizeof(char));
    // including ","
    if(n_digit + 1 != snprintf(buf, n_digit + 2, "%zu,", shape[i])){
      SNPYIO_ERROR("snprintf failed\n");
      return 1;
    }
    // copy result excluding null character
    memcpy((*value) + offset, buf, sizeof(char) * (n_digit + 1));
    offset += n_digit + 1;
    memory_free(buf);
  }
  // first character is a parenthesis
  (*value)[          0] = '(';
  // last-1 character is a parenthesis
  (*value)[n_chars - 1] = ')';
  // last character is null character
  (*value)[    n_chars] = null_char;
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
  // buffers to store values
  char *descr_value         = NULL;
  char *fortran_order_value = NULL;
  char *shape_value         = NULL;
  // create dictionary values,
  //   in which inputs are evaluated and sanitised
  // create value of data type | 4
  if(0 != create_descr_value(&descr_value, dtype)){
    SNPYIO_ERROR("create_descr_value failed\n");
    return 1;
  }
  // create value of memory order | 4
  if(0 != create_fortran_order_value(&fortran_order_value, is_fortran_order)){
    SNPYIO_ERROR("create_fortran_order_value failed\n");
    return 1;
  }
  // create value of data sizes | 4
  if(0 != create_shape_value(&shape_value, ndim, shape)){
    SNPYIO_ERROR("create_shape_value failed\n");
    return 1;
  }
  // assign all elements (strings) which compose dict | 7
  const char *elements[] = {
    "{",
    "'descr'",         ":", (char *)descr_value,         ",",
    "'fortran_order'", ":", (char *)fortran_order_value, ",",
    "'shape'",         ":", (char *)shape_value,         ",",
    "}",
  };
  // check total number of characters of
  //   {'descr':VALUE,'fortran_order':VALUE,'shape':VALUE}
  //   to allocate dict
  // NOTE: n_chars_dict is the number of characters of dict
  //   INCLUDING the last null character, while n_dict = strlen(dict),
  //   EXCLUDING the last null character.
  //   Thus n_dict = n_chars_dict - 1
  size_t n_chars_dict = 0;
  for(size_t i = 0; i < sizeof(elements) / sizeof(char *); i++){
    n_chars_dict += strlen(elements[i]);
  }
  // allocate dict and assign above "elements"
  // +1: null character
  *dict = memory_calloc(n_chars_dict + 1, sizeof(char));
  for(size_t i = 0, offset = 0; i < sizeof(elements) / sizeof(char *); i++){
    const size_t n_chars = strlen(elements[i]);
    memcpy((*dict) + offset, elements[i], sizeof(char) * n_chars);
    offset += n_chars;
  }
  (*dict)[n_chars_dict] = null_char;
  // clean-up all working memories
  memory_free(descr_value);
  memory_free(fortran_order_value);
  memory_free(shape_value);
  // as the length of "dict", use length WITHOUT null character,
  // i.e. strlen(*dict)
  *n_dict = strlen(*dict);
  SNPYIO_LOGGING("dict: %s\n", *dict);
  SNPYIO_LOGGING("size: %zu\n", *n_dict);
  return 0;
}

static int create_padding(uint8_t **padding, size_t *n_padding, uint8_t *major_version, size_t *header_size, const size_t n_dict){
  /*
   * The following relation holds for the header size
   *   header_size =
   *     + sizeof(magic string)      (= 6            bytes)
   *     + sizeof(major_version)     (= 1            byte )
   *     + sizeof(minor_version)     (= 1            byte )
   *     + sizeof(header_len)        (= 2 or 4       bytes)
   *     + sizeof(char)*strlen(dict) (= strlen(dict) bytes)
   *     + sizeof(uint8_t)*n_padding (= n_padding    bytes)
   *   is divisible by header_block_size
   * Definitely this is not generally true, and we need some paddings
   *   consisting of some (0 or more) spaces ' ' and one newline '\n',
   *   whose length (number of elements) is returned
   */
  // size of each element is computed | 5
  const size_t n_magic_string = strlen(magic_string);
  const size_t size_magic_string  = n_magic_string * sizeof(   char);
  const size_t size_major_version = 1              * sizeof(uint8_t);
  const size_t size_minor_version = 1              * sizeof(uint8_t);
  const size_t size_dict          = n_dict         * sizeof(   char);
  // reject too large dict | 4
  if(UINT_MAX - header_block_size < size_dict){
    SNPYIO_ERROR("size of dictionary is huge (%zu)\n", size_dict);
    return 1;
  }
  // decide major version and datatype of HEADER_LEN | 11
  // large portion of the header is occupied by dict
  // so check dict size, and if it is larger than USHRT_MAX - header_block_size,
  //   use major_version = 2
  size_t size_header_len = 0;
  if(USHRT_MAX - header_block_size < size_dict){
    *major_version = 2;
    size_header_len = sizeof(uint32_t);
  }else{
    *major_version = 1;
    size_header_len = sizeof(uint16_t);
  }
  // compute size of all data except padding | 6
  const size_t size_except_padding =
    + size_magic_string
    + size_major_version
    + size_minor_version
    + size_header_len
    + size_dict;
  // decide total size of the header, which should be header_block_size x N | 8
  // increase total size by header_block_size until becoming larger than size_except_padding
  // NOTE: size_padding == 0 is NOT allowed since '\n' is necessary at the end
  //   thus the condition to continue loop is "<=", not "<"
  *header_size = 0;
  while(*header_size <= size_except_padding){
    *header_size += header_block_size;
  }
  const size_t size_padding = *header_size - size_except_padding;
  // create padding | 6
  *n_padding = size_padding / sizeof(uint8_t);
  *padding = memory_calloc(*n_padding, sizeof(uint8_t));
  // many ' 's: 0x20
  memset(*padding, 0x20, sizeof(uint8_t) * (*n_padding - 1));
  // last '\n': 0x0a
  (*padding)[*n_padding - 1] = 0x0a;
  SNPYIO_LOGGING("padding, size: %zu\n", *n_padding);
  return 0;
}

static int create_header_len(uint8_t **header_len, size_t *n_header_len, const uint8_t major_version, const size_t n_dict, const size_t n_padding){
  /*
   * HEADER_LEN = n_dict + n_padding,
   *   which is written in little-endian
   *   (irrespective to the architecture)
   */
  // reject too large dict / padding sizes | 12
  // Here "too large" means header size (not data size)
  //   is larger than approx. 2GB, which would not happen normally
  if(UINT_MAX / 2 <= n_dict){
    SNPYIO_ERROR("dictionary size is huge (%zu)\n", n_dict);
    return 1;
  }
  // padding is to make header size header_block_size x N
  // so it should not exceed header_block_size
  if(header_block_size < n_padding){
    SNPYIO_ERROR("padding size is huge (%zu)\n", n_padding);
    return 1;
  }
  // compute header_len and store as an array of uint8_t | 15
  if(1 == major_version){
    // major version 1, use uint16_t to store header_len
    const uint16_t header_len_ = (uint16_t)n_dict + (uint16_t)n_padding;
    *n_header_len = sizeof(uint16_t) / sizeof(uint8_t);
    *header_len = memory_calloc(*n_header_len, sizeof(uint8_t));
    memcpy(*header_len, &header_len_, *n_header_len);
    SNPYIO_LOGGING("header_len (uint16_t): %hu\n", header_len_);
  }else{
    // major version 2, use uint32_t to store header_len
    const uint32_t header_len_ = (uint32_t)n_dict + (uint32_t)n_padding;
    *n_header_len = sizeof(uint32_t) / sizeof(uint8_t);
    *header_len = memory_calloc(*n_header_len, sizeof(uint8_t));
    memcpy(*header_len, &header_len_, *n_header_len);
    SNPYIO_LOGGING("header_len (uint32_t): %u\n", header_len_);
  }
  // convert endian of buffer which will be written if needed | 4
  if(is_big_endian() && 0 != convert_endian(header_len, sizeof(*header_len))){
    SNPYIO_ERROR("convert_endian failed\n");
    return 1;
  }
  return 0;
}

/**
 * @brief write NPY header
 * @param[in] ndim             : number of dimensions of the data set, e.g. 2
 * @param[in] shape            : number of points of the data set in each dimension, e.g. [3, 4]
 * @param[in] dtype            : data type, e.g. "'<f8'"
 * @param[in] is_fortran_order : row-major order (false) or column-major order (true)
 * @param[in] fp               : file stream to which the header is written
 * @return                     : (success) written header size (in bytes)
 *                               (failure) 0
 */
size_t snpyio_w_header(const size_t ndim, const size_t *shape, const char dtype[], const bool is_fortran_order, FILE *fp){
  if(0 != sanitise_shape(ndim, shape)) goto w_err_hndl;
  if(0 != sanitise_dtype(dtype)) goto w_err_hndl;
  if(0 != sanitise_fp(fp)) goto w_err_hndl;
  // magic_string
  const size_t n_magic_string = strlen(magic_string);
  // minor_version, always 0 | 1
  const uint8_t minor_version = 0;
  // dictionary (and its size) | 5
  char *dict = NULL;
  size_t n_dict = 0;
  if(0 != create_dict(&dict, &n_dict, ndim, shape, dtype, is_fortran_order)){
    goto w_err_hndl;
  }
  // major_version and padding (and its size) | 7
  uint8_t *padding = NULL;
  size_t n_padding = 0;
  uint8_t major_version = 1;
  size_t header_size = 0;
  if(0 != create_padding(&padding, &n_padding, &major_version, &header_size, n_dict)){
    goto w_err_hndl;
  }
  // comptue header_len | 5
  uint8_t *header_len = NULL;
  size_t n_header_len = 0;
  if(0 != create_header_len(&header_len, &n_header_len, major_version, n_dict, n_padding)){
    goto w_err_hndl;
  }
  // dump all information to a buffer "header" and compute total size "header_size"
  const size_t sizes[6] = {
    n_magic_string * sizeof(   char),
    1              * sizeof(uint8_t),
    1              * sizeof(uint8_t),
    n_header_len   * sizeof(uint8_t),
    n_dict         * sizeof(   char),
    n_padding      * sizeof(uint8_t),
  };
  if(0 != myfwrite(magic_string,   sizes[0], 1, fp)) goto w_err_hndl;
  if(0 != myfwrite(&major_version, sizes[1], 1, fp)) goto w_err_hndl;
  if(0 != myfwrite(&minor_version, sizes[2], 1, fp)) goto w_err_hndl;
  if(0 != myfwrite(header_len,     sizes[3], 1, fp)) goto w_err_hndl;
  if(0 != myfwrite(dict,           sizes[4], 1, fp)) goto w_err_hndl;
  if(0 != myfwrite(padding,        sizes[5], 1, fp)) goto w_err_hndl;
  // clean-up all buffers
  memory_free(dict);
  memory_free(padding);
  memory_free(header_len);
  SNPYIO_LOGGING("header_size: %zu\n", header_size);
  return header_size;
w_err_hndl:
  error_handlings();
  return 0;
}

#undef SNPYIO_MESSAGE
#undef SNPYIO_LOGGING
#undef SNPYIO_ERROR
#undef SNPYIO_FATAL


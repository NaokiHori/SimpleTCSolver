/*
 * Copyright 2022 Naoki Hori
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

// https://github.com/NaokiHori/SimpleDecomp

#include "sdecomp.h"
#include "../internal.h"


/* process distributions */

typedef enum {
  TYPE_NPROCS = 0,
  TYPE_MYRANK = 1
} type_t;

static int get_process_config_2d(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const type_t type, int *result){
#define SDECOMP_INTERNAL_NDIMS 2
  // compute
  //   type == TYPE_NPROCS: number of processes
  // or
  //   type == TYPE_MYRANK: my location
  // of the given decomposition "info"
  //   in the given dimension "dim"
  // First, get info about x1pencil using MPI_Cart_get
  int x1pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int         periods[SDECOMP_INTERNAL_NDIMS] = {0};
  int x1pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  MPI_Cart_get(
      info->comm_cart,
      SDECOMP_INTERNAL_NDIMS,
      x1pencil_nprocs,
      periods,
      x1pencil_myrank
  );
  if(SDECOMP_X1PENCIL == pencil){
    *result = TYPE_NPROCS == type
      ? x1pencil_nprocs[dir]
      : x1pencil_myrank[dir];
    return 0;
  }
  // y1 pencil, determined by x1 pencil
  int y1pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int y1pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  y1pencil_nprocs[0] = x1pencil_nprocs[1];
  y1pencil_nprocs[1] = x1pencil_nprocs[0];
  y1pencil_myrank[0] = x1pencil_myrank[1];
  y1pencil_myrank[1] = x1pencil_myrank[0];
  if(SDECOMP_Y1PENCIL == pencil){
    *result = TYPE_NPROCS == type
      ? y1pencil_nprocs[dir]
      : y1pencil_myrank[dir];
    return 0;
  }
  SDECOMP_ERROR(
      "should not reach here (%s:%d)\n",
      error_label, __FILE__, __LINE__
  );
  return 1;
#undef SDECOMP_INTERNAL_NDIMS
}

static int get_process_config_3d(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const type_t type, int *result){
#define SDECOMP_INTERNAL_NDIMS 3
  // compute
  //   type == TYPE_NPROCS: number of processes
  // or
  //   type == TYPE_MYRANK: my location
  // of the given decomposition "info"
  //   in the given dimension "dim"
  // First, get info about x1pencil using MPI_Cart_get
  int x1pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int         periods[SDECOMP_INTERNAL_NDIMS] = {0};
  int x1pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  MPI_Cart_get(
      info->comm_cart,
      SDECOMP_INTERNAL_NDIMS,
      x1pencil_nprocs,
      periods,
      x1pencil_myrank
  );
  if(SDECOMP_X1PENCIL == pencil){
    *result = TYPE_NPROCS == type
      ? x1pencil_nprocs[dir]
      : x1pencil_myrank[dir];
    return 0;
  }
  // y1 pencil, determined by x1 pencil
  int y1pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int y1pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  y1pencil_nprocs[0] = x1pencil_nprocs[1];
  y1pencil_nprocs[1] = x1pencil_nprocs[0];
  y1pencil_nprocs[2] = x1pencil_nprocs[2];
  y1pencil_myrank[0] = x1pencil_myrank[1];
  y1pencil_myrank[1] = x1pencil_myrank[0];
  y1pencil_myrank[2] = x1pencil_myrank[2];
  if(SDECOMP_Y1PENCIL == pencil){
    *result = TYPE_NPROCS == type
      ? y1pencil_nprocs[dir]
      : y1pencil_myrank[dir];
    return 0;
  }
  // z1 pencil, determined by y1 pencil
  int z1pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int z1pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  z1pencil_nprocs[0] = y1pencil_nprocs[0];
  z1pencil_nprocs[1] = y1pencil_nprocs[2];
  z1pencil_nprocs[2] = y1pencil_nprocs[1];
  z1pencil_myrank[0] = y1pencil_myrank[0];
  z1pencil_myrank[1] = y1pencil_myrank[2];
  z1pencil_myrank[2] = y1pencil_myrank[1];
  if(SDECOMP_Z1PENCIL == pencil){
    *result = TYPE_NPROCS == type
      ? z1pencil_nprocs[dir]
      : z1pencil_myrank[dir];
    return 0;
  }
  // x2 pencil, determined by z1 pencil
  int x2pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int x2pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  x2pencil_nprocs[0] = z1pencil_nprocs[2];
  x2pencil_nprocs[1] = z1pencil_nprocs[1];
  x2pencil_nprocs[2] = z1pencil_nprocs[0];
  x2pencil_myrank[0] = z1pencil_myrank[2];
  x2pencil_myrank[1] = z1pencil_myrank[1];
  x2pencil_myrank[2] = z1pencil_myrank[0];
  if(SDECOMP_X2PENCIL == pencil){
    *result = type == TYPE_NPROCS
      ? x2pencil_nprocs[dir]
      : x2pencil_myrank[dir];
    return 0;
  }
  // y2 pencil, determined by x2 pencil
  int y2pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int y2pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  y2pencil_nprocs[0] = x2pencil_nprocs[1];
  y2pencil_nprocs[1] = x2pencil_nprocs[0];
  y2pencil_nprocs[2] = x2pencil_nprocs[2];
  y2pencil_myrank[0] = x2pencil_myrank[1];
  y2pencil_myrank[1] = x2pencil_myrank[0];
  y2pencil_myrank[2] = x2pencil_myrank[2];
  if(SDECOMP_Y2PENCIL == pencil){
    *result = type == TYPE_NPROCS
      ? y2pencil_nprocs[dir]
      : y2pencil_myrank[dir];
    return 0;
  }
  // z2 pencil, determined by y2 pencil
  int z2pencil_nprocs[SDECOMP_INTERNAL_NDIMS] = {0};
  int z2pencil_myrank[SDECOMP_INTERNAL_NDIMS] = {0};
  z2pencil_nprocs[0] = y2pencil_nprocs[0];
  z2pencil_nprocs[1] = y2pencil_nprocs[2];
  z2pencil_nprocs[2] = y2pencil_nprocs[1];
  z2pencil_myrank[0] = y2pencil_myrank[0];
  z2pencil_myrank[1] = y2pencil_myrank[2];
  z2pencil_myrank[2] = y2pencil_myrank[1];
  if(SDECOMP_Z2PENCIL == pencil){
    *result = type == TYPE_NPROCS
      ? z2pencil_nprocs[dir]
      : z2pencil_myrank[dir];
    return 0;
  }
  SDECOMP_ERROR(
      "should not reach here (%s:%d)\n",
      error_label, __FILE__, __LINE__
  );
  return 1;
#undef SDECOMP_INTERNAL_NDIMS
}

/**
 * @brief get number of processes in the given dimension
 * @param[in]  info   : struct containing information of process distribution
 * @param[in]  pencil : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in]  dir    : direction which I am interested in
 * @param[out] nprocs : (success) number of processes in the given dimension
 *                      (failure) undefined
 * @return            : (success) 0
 *                      (failure) non-zero value
 */
int sdecomp_internal_get_nprocs(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int *nprocs){
  const char error_label[] = {"sdecomp.get_nprocs"};
  // NULL check
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "nprocs", nprocs)) return 1;
  // sanitise
  if(0 != sdecomp_internal_sanitise_pencil(error_label, info->ndims, pencil)) return 1;
  if(0 != sdecomp_internal_sanitise_dir   (error_label, info->ndims,    dir)) return 1;
  // main
  if(2 == info->ndims){
    if(0 != get_process_config_2d(error_label, info, pencil, dir, TYPE_NPROCS, nprocs)) return 1;
  }else{
    if(0 != get_process_config_3d(error_label, info, pencil, dir, TYPE_NPROCS, nprocs)) return 1;
  }
  return 0;
}

/**
 * @brief get my process position in the whole domain
 * @param[in]  info   : struct containing information of process distribution
 * @param[in]  pencil : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in]  dir    : direction which I am interested in
 * @param[out] myrank : (success) my position in the given dimension
 *                      (failure) undefined
 * @return            : (success) 0
 *                      (failure) non-zero value
 */
int sdecomp_internal_get_myrank(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int *myrank){
  const char error_label[] = {"sdecomp.get_myrank"};
  // NULL check
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "myrank", myrank)) return 1;
  // sanitise
  if(0 != sdecomp_internal_sanitise_pencil(error_label, info->ndims, pencil)) return 1;
  if(0 != sdecomp_internal_sanitise_dir   (error_label, info->ndims,    dir)) return 1;
  // main
  if(2 == info->ndims){
    if(0 != get_process_config_2d(error_label, info, pencil, dir, TYPE_MYRANK, myrank)) return 1;
  }else{
    if(0 != get_process_config_3d(error_label, info, pencil, dir, TYPE_MYRANK, myrank)) return 1;
  }
  return 0;
}

/**
 * @brief get neighbour process ranks
 * @param[in]  info       : struct containing information of process distribution
 * @param[in]  pencil     : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in]  dir        : direction which I am interested in
 * @param[out] neighbours : (success) neighbour ranks in sdecomp->comm_cart
 *                          (failure) undefined
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int sdecomp_internal_get_neighbours(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int neighbours[2]){
  const char error_label[] = {"sdecomp.get_neighbours"};
  if(0 != sdecomp_internal_sanitise_null(error_label,       "info",       info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "neighbours", neighbours)) return 1;
  // reject all invalid pencils and directions to simplify the following conditions
  if(0 != sdecomp_internal_sanitise_pencil(error_label, info->ndims, pencil)) return 1;
  if(0 != sdecomp_internal_sanitise_dir   (error_label, info->ndims,    dir)) return 1;
  // initialise with "no neighbour"
  neighbours[0] = MPI_PROC_NULL;
  neighbours[1] = MPI_PROC_NULL;
  // check direction, I need direction w.r.t. comm_cart (x1 pencil)
  // NOTE: see figure to see the correspondence
  int direction = -1;
  if(2 == info->ndims){
    // 2D
    const int alldirs[2][2] = {
      // x1 pencil
      {(int)SDECOMP_XDIR, (int)SDECOMP_YDIR,},
      // y1 pencil
      {(int)SDECOMP_YDIR, (int)SDECOMP_XDIR,},
    };
    direction = alldirs[pencil][dir];
  }else{
    // 3D
    const int alldirs[6][3] = {
      // x1 pencil
      {(int)SDECOMP_XDIR, (int)SDECOMP_YDIR, (int)SDECOMP_ZDIR,},
      // y1 pencil
      {(int)SDECOMP_YDIR, (int)SDECOMP_XDIR, (int)SDECOMP_ZDIR,},
      // z1 pencil
      {(int)SDECOMP_YDIR, (int)SDECOMP_ZDIR, (int)SDECOMP_XDIR,},
      // x2 pencil
      {(int)SDECOMP_XDIR, (int)SDECOMP_ZDIR, (int)SDECOMP_YDIR,},
      // y2 pencil
      {(int)SDECOMP_ZDIR, (int)SDECOMP_XDIR, (int)SDECOMP_YDIR,},
      // z2 pencil
      {(int)SDECOMP_ZDIR, (int)SDECOMP_YDIR, (int)SDECOMP_XDIR,},
    };
    direction = alldirs[pencil][dir];
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(info, &comm_cart);
  // always interested in one-process away
  const int disp = 1;
  MPI_Cart_shift(
      comm_cart,
      direction,
      disp,
      neighbours + 0,
      neighbours + 1
  );
  return 0;
}

/* pencil information */

static int get_pencil_config(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const size_t glsize, int (*kernel)(const char error_label[], const size_t glsize, const int nprocs, const int myrank, size_t *result), size_t *result){
  // wrapper function to get mysize / offset of pencil
  // NOTE: since this is only called by specific entrypoints,
  //   input sanitisers are omitted
  // compute number of process in the given direction
  int nprocs = 0;
  if(0 != sdecomp_internal_get_nprocs(info, pencil, dir, &nprocs)) return 1;
  // compute my location in the given direction
  int myrank = 0;
  if(0 != sdecomp_internal_get_myrank(info, pencil, dir, &myrank)) return 1;
  // call one of "number of grid calculator" and "offset calculator"
  if(0 != kernel(error_label, glsize, nprocs, myrank, result)) return 1;
  return 0;
}

/**
 * @brief get number of grid points of the given pencil in the given direction
 * @param[in]  info   : struct containing information of process distribution
 * @param[in]  pencil : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in]  dir    : direction which I am interested in
 * @param[in]  glsize : number of total grid points in the direction
 * @param[out] mysize : (success) number of local grid points in the direction
 *                      (failure) undefined
 * @return            : (success) 0
 *                      (failure) non-zero value
 */
int sdecomp_internal_get_pencil_mysize(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const size_t glsize, size_t *mysize){
  const char error_label[] = {"sdecomp.get_pencil_mysize"};
  // NULL check
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "mysize", mysize)) return 1;
  // sanitise
  if(0 != sdecomp_internal_sanitise_pencil(error_label, info->ndims, pencil)) return 1;
  if(0 != sdecomp_internal_sanitise_dir   (error_label, info->ndims,    dir)) return 1;
  // main
  if(0 != get_pencil_config(error_label, info, pencil, dir, glsize, sdecomp_internal_kernel_get_mysize, mysize)) return 1;
  return 0;
}

/**
 * @brief get offset of grid points of the given pencil in the given direction
 * @param[in]  info   : struct containing information of process distribution
 * @param[in]  pencil : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in]  dir    : direction which I am interested in
 * @param[in]  glsize : number of global grid points in the direction
 * @param[out] offset : (success) offset in the give direction
 *                      (failure) undefined
 * @return            : (success) 0
 *                      (failure) non-zero value
 */
int sdecomp_internal_get_pencil_offset(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const size_t glsize, size_t *offset){
  const char error_label[] = {"sdecomp.get_pencil_mysize"};
  // NULL check
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "offset", offset)) return 1;
  // sanitise
  if(0 != sdecomp_internal_sanitise_pencil(error_label, info->ndims, pencil)) return 1;
  if(0 != sdecomp_internal_sanitise_dir   (error_label, info->ndims,    dir)) return 1;
  // main
  if(0 != get_pencil_config(error_label, info, pencil, dir, glsize, sdecomp_internal_kernel_get_offset, offset)) return 1;
  return 0;
}


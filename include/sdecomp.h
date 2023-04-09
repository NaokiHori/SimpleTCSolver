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

#if !defined(SDECOMP_H)
#define SDECOMP_H

/*** public structures and APIs ***/

#include <stddef.h>  // size_t
#include <stdbool.h> // bool
#include <mpi.h>     // MPI_Comm

// struct storing domain decomposition information
typedef struct sdecomp_info_t_ sdecomp_info_t;
// struct storing pencil transpose plan
typedef struct sdecomp_transpose_plan_t_ sdecomp_transpose_plan_t;

// spatial directions
typedef unsigned short sdecomp_dir_t;
extern const sdecomp_dir_t SDECOMP_XDIR; // 0
extern const sdecomp_dir_t SDECOMP_YDIR; // 1
extern const sdecomp_dir_t SDECOMP_ZDIR; // 2

// pencil types
typedef unsigned short sdecomp_pencil_t;
extern const sdecomp_pencil_t SDECOMP_X1PENCIL; // 0
extern const sdecomp_pencil_t SDECOMP_Y1PENCIL; // 1
extern const sdecomp_pencil_t SDECOMP_Z1PENCIL; // 2
extern const sdecomp_pencil_t SDECOMP_X2PENCIL; // 3
extern const sdecomp_pencil_t SDECOMP_Y2PENCIL; // 4
extern const sdecomp_pencil_t SDECOMP_Z2PENCIL; // 5

/* APIs of sdecomp_transpose_t */
// accessed by sdecomp.transpose.xxx
typedef struct {
  // constructor of sdecomp_transpose_plan_t
  int (* const construct)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil_bef,
      const sdecomp_pencil_t pencil_aft,
      const size_t *glsizes,
      const size_t size_of_element,
      sdecomp_transpose_plan_t **plan // out
  );
  // transpose runner
  int (* const execute)(
      sdecomp_transpose_plan_t * restrict plan,
      const void * restrict sendbuf,
      void * restrict recvbuf
  );
  // destructor of sdecomp_transpose_plan_t
  int (* const destruct)(
      sdecomp_transpose_plan_t *plan
  );
} sdecomp_transpose_t;

/* APIs of sdecomp_t */
// accessed by sdecomp.xxx
typedef struct {
  // constructor of sdecomp_info_t
  int (* const construct)(
      const MPI_Comm comm_default,
      const size_t ndims,
      const size_t *dims,
      const bool *periods,
      sdecomp_info_t **info // out
  );
  // destructor of sdecomp_info_t
  int (* const destruct)(
      sdecomp_info_t *sdecomp
  );
  // getter, number of dimensions
  int (* const get_ndims)(
      const sdecomp_info_t *info,
      size_t *ndims // out
  );
  // getter, number of total processes in comm_cart
  int (* const get_comm_size)(
      const sdecomp_info_t *info,
      int *nprocs // out
  );
  // getter, my ID in comm_cart
  int (* const get_comm_rank)(
      const sdecomp_info_t *info,
      int *myrank // out
  );
  // getter, number of processes in one dimension
  int (* const get_nprocs)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil,
      const sdecomp_dir_t dir,
      int *nprocs // out
  );
  // getter, position of my process
  int (* const get_myrank)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil,
      const sdecomp_dir_t dir,
      int *myrank // out
  );
  // getter, ranks of my neighbours
  int (* const get_neighbours)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil,
      const sdecomp_dir_t dir,
      int neighbours[2] // out
  );
  // getter, local size of my pencil
  int (* const get_pencil_mysize)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil,
      const sdecomp_dir_t dir,
      const size_t glsize,
      size_t *mysize // out
  );
  // getter, offset of my pencil
  int (* const get_pencil_offset)(
      const sdecomp_info_t *info,
      const sdecomp_pencil_t pencil,
      const sdecomp_dir_t dir,
      const size_t glsize,
      size_t *offset // out
  );
  // getter, default communicator
  int (* const get_comm_cart)(
      const sdecomp_info_t *info,
      MPI_Comm *comm // out
  );
  // transpose functions sdecomp.transpose
  const sdecomp_transpose_t transpose;
} sdecomp_t;

extern const sdecomp_t sdecomp;

#endif // SDECOMP_H

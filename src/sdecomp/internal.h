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

#if !defined(SDECOMP_INTERNAL_H)
#define SDECOMP_INTERNAL_H

/*
 * FOR INTERNAL USE
 * DO NOT INCLUDE THIS HEADER
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "sdecomp.h"

// error messages
#define SDECOMP_ERROR(...){ \
  FILE *stream = stderr; \
  fprintf(stream, "[SDECOMP ERROR (%s)] " __VA_ARGS__); \
  fflush(stream); \
}

struct sdecomp_info_t_ {
  MPI_Comm comm_cart;
  size_t ndims;
};

// memory management
extern void *sdecomp_internal_calloc(const char error_label[], const size_t count, const size_t size);
extern void sdecomp_internal_free(void *ptr);

// kernel functions to decide local size / offset of pencils
extern int sdecomp_internal_kernel_get_mysize(const char error_label[], const size_t glsize, const int nprocs, const int myrank, size_t *mysize);
extern int sdecomp_internal_kernel_get_offset(const char error_label[], const size_t glsize, const int nprocs, const int myrank, size_t *offset);

// constructor and destructor
extern int sdecomp_internal_construct(const MPI_Comm comm_default, const size_t ndims, const size_t *dims, const bool *periods, sdecomp_info_t **info);
extern int sdecomp_internal_destruct(sdecomp_info_t *info);

// getters, just return members
extern int sdecomp_internal_get_ndims(const sdecomp_info_t *info, size_t *ndims);
extern int sdecomp_internal_get_comm_size(const sdecomp_info_t *info, int *nprocs);
extern int sdecomp_internal_get_comm_rank(const sdecomp_info_t *info, int *myrank);
extern int sdecomp_internal_get_comm_cart(const sdecomp_info_t *info, MPI_Comm *comm_cart);
// getters which require certain amount of operation
extern int sdecomp_internal_get_nprocs(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int *nprocs);
extern int sdecomp_internal_get_myrank(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int *myrank);
extern int sdecomp_internal_get_neighbours(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, int neighbours[2]);
extern int sdecomp_internal_get_pencil_mysize(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const size_t glsize, size_t *mysize);
extern int sdecomp_internal_get_pencil_offset(const sdecomp_info_t *info, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const size_t glsize, size_t *offset);

// functions which take care of transpose operations
extern int sdecomp_internal_transpose_construct(const sdecomp_info_t *info, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft, const size_t *glsizes, const size_t size_of_element, sdecomp_transpose_plan_t **plan);
extern int sdecomp_internal_transpose_execute(sdecomp_transpose_plan_t *plan, const void * restrict sendbuf, void * restrict recvbuf);
extern int sdecomp_internal_transpose_destruct(sdecomp_transpose_plan_t *plan);

// sanitisers
extern int sdecomp_internal_sanitise_null(const char error_label[], const char ptr_name[], const void *ptr);
extern int sdecomp_internal_sanitise_ndims(const char error_label[], const size_t ndims);
extern int sdecomp_internal_sanitise_comm(const char error_label[], const MPI_Comm comm);
extern int sdecomp_internal_sanitise_pencil(const char error_label[], const size_t ndims, const sdecomp_pencil_t pencil);
extern int sdecomp_internal_sanitise_dir(const char error_label[], const size_t ndims, const sdecomp_dir_t dir);
extern int sdecomp_internal_sanitise_nprocs(const char error_label[], const int nprocs);
extern int sdecomp_internal_sanitise_myrank(const char error_label[], const int myrank);
extern int sdecomp_internal_sanitise_glsize(const char error_label[], const size_t glsize);
extern int sdecomp_internal_sanitise_size_of_element(const char error_label[], const size_t size_of_element);
extern int sdecomp_internal_sanitise_pencil_pair_2d(const char error_label[], const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft);
extern int sdecomp_internal_sanitise_pencil_pair_3d(const char error_label[], const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft, bool *is_forward);

#endif // SDECOMP_INTERNAL_H

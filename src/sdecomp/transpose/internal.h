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

#if !defined(SDECOMP_INTERNAL_TRANSPOSE_H)
#define SDECOMP_INTERNAL_TRANSPOSE_H

/*
 * FOR INTERNAL USE
 * DO NOT INCLUDE THIS HEADER
 */

struct sdecomp_transpose_plan_t_ {
  int * restrict scounts;
  int * restrict rcounts;
  int * restrict sdispls;
  int * restrict rdispls;
  MPI_Datatype * restrict stypes;
  MPI_Datatype * restrict rtypes;
  MPI_Comm comm_2d;
};

extern int sdecomp_internal_transpose_allocate(const char error_label[], const int nprocs_2d, sdecomp_transpose_plan_t **plan);
extern int sdecomp_internal_transpose_deallocate(sdecomp_transpose_plan_t *plan);

extern int sdecomp_internal_transpose_init_2d(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil_bef, sdecomp_pencil_t pencil_aft, const size_t *glsizes, const size_t size_of_element, sdecomp_transpose_plan_t **plan);
extern int sdecomp_internal_transpose_init_3d(const char error_label[], const sdecomp_info_t *info, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft, const size_t *glsizes, const size_t size_of_element, sdecomp_transpose_plan_t **plan);

extern int sdecomp_internal_execute(sdecomp_transpose_plan_t * restrict plan, const void * restrict sendbuf, void * restrict recvbuf);

#endif // SDECOMP_INTERNAL_TRANSPOSE_H

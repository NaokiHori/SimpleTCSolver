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
#include "internal.h"

// spatial directions
const sdecomp_dir_t SDECOMP_XDIR = 0;
const sdecomp_dir_t SDECOMP_YDIR = 1;
const sdecomp_dir_t SDECOMP_ZDIR = 2;

// pencil types
const sdecomp_pencil_t SDECOMP_X1PENCIL = 0;
const sdecomp_pencil_t SDECOMP_Y1PENCIL = 1;
const sdecomp_pencil_t SDECOMP_Z1PENCIL = 2;
const sdecomp_pencil_t SDECOMP_X2PENCIL = 3;
const sdecomp_pencil_t SDECOMP_Y2PENCIL = 4;
const sdecomp_pencil_t SDECOMP_Z2PENCIL = 5;

/* assign all "methods" (pointers to all internal functions) */
const sdecomp_t sdecomp = {
  .construct         = sdecomp_internal_construct,
  .destruct          = sdecomp_internal_destruct,
  .get_ndims         = sdecomp_internal_get_ndims,
  .get_comm_size     = sdecomp_internal_get_comm_size,
  .get_comm_rank     = sdecomp_internal_get_comm_rank,
  .get_nprocs        = sdecomp_internal_get_nprocs,
  .get_myrank        = sdecomp_internal_get_myrank,
  .get_neighbours    = sdecomp_internal_get_neighbours,
  .get_pencil_mysize = sdecomp_internal_get_pencil_mysize,
  .get_pencil_offset = sdecomp_internal_get_pencil_offset,
  .get_comm_cart     = sdecomp_internal_get_comm_cart,
  .transpose         = {
    .construct = sdecomp_internal_transpose_construct,
    .execute   = sdecomp_internal_transpose_execute,
    .destruct  = sdecomp_internal_transpose_destruct,
  },
};


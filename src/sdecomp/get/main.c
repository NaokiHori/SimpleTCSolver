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


/* just return members */

int sdecomp_internal_get_ndims(const sdecomp_info_t *info, size_t *ndims){
  const char error_label[] = {"sdecomp.get_ndims"};
  if(0 != sdecomp_internal_sanitise_null(error_label,  "info",  info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "ndims", ndims)) return 1;
  *ndims = info->ndims;
  return 0;
}

int sdecomp_internal_get_comm_cart(const sdecomp_info_t *info, MPI_Comm *comm_cart){
  const char error_label[] = {"sdecomp.get_comm_cart"};
  if(0 != sdecomp_internal_sanitise_null(error_label,      "info",      info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "comm_cart", comm_cart)) return 1;
  *comm_cart = info->comm_cart;
  return 0;
}

/* simple MPI wrappers */

int sdecomp_internal_get_comm_size(const sdecomp_info_t *info, int *nprocs){
  const char error_label[] = {"sdecomp.get_comm_size"};
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "nprocs", nprocs)) return 1;
  MPI_Comm_size(info->comm_cart, nprocs);
  return 0;
}

int sdecomp_internal_get_comm_rank(const sdecomp_info_t *info, int *myrank){
  const char error_label[] = {"sdecomp.get_comm_rank"};
  if(0 != sdecomp_internal_sanitise_null(error_label,   "info",   info)) return 1;
  if(0 != sdecomp_internal_sanitise_null(error_label, "myrank", myrank)) return 1;
  MPI_Comm_rank(info->comm_cart, myrank);
  return 0;
}


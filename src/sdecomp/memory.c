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

// functions which take care of memory management

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "sdecomp.h"
#include "internal.h"


void *sdecomp_internal_calloc(const char error_label[], const size_t count, const size_t size){
  errno = 0;
  void *ptr = calloc(count, size);
  if(NULL == ptr){
    const char *message = strerror(errno);
    SDECOMP_ERROR(
        "%s",
        error_label, message
    );
  }
  return ptr;
}

void sdecomp_internal_free(void *ptr){
  // no NULL check needes since free is NULL-safe
  free(ptr);
}


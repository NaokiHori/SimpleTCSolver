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

#if !defined(SNPYIO_H)
#define SNPYIO_H

#include <stdio.h>   // size_t, FILE
#include <stdbool.h> // bool

extern size_t snpyio_w_header(const size_t  ndim, const size_t  *shape, const char dtype[], const bool  is_fortran_order, FILE *fp);
extern size_t snpyio_r_header(      size_t *ndim,       size_t **shape,       char **dtype,       bool *is_fortran_order, FILE *fp);

#endif // SNPYIO_H

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
extern size_t snpyio_w_header(const size_t  ndim, const size_t  *shape, const char dtype[], const bool  is_fortran_order, FILE *fp);

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
extern size_t snpyio_r_header(      size_t *ndim,       size_t **shape,       char **dtype,       bool *is_fortran_order, FILE *fp);

#endif // SNPYIO_H

#if !defined(COMMON_H)
#define COMMON_H

#define NDIMS 3

#include <stdlib.h> // size_t
#include <math.h>
// C99 does not specify M_PI in math.h
#if !defined(M_PI)
#define M_PI 3.1415926535897932
#endif

extern void *common_calloc(const size_t count, const size_t size);
extern void common_free(void *ptr);

/* save files */
// ref: https://numpy.org/doc/stable/reference/generated/numpy.dtype.html
// dtype, boolean
extern const char NPY_BOL[];
// dtype, integer
extern const char NPY_INT[];
// dtype, double
extern const char NPY_DBL[];
// general wrapper to save NPY files
extern void common_save(const char fname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize, const void *data);

#endif // COMMON_H

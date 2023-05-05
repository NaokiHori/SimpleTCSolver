#if !defined(COMMON_H)
#define COMMON_H

#define NDIMS 3

#include <stdlib.h> // size_t

extern void *common_calloc(const size_t count, const size_t size);
extern void common_free(void *ptr);
extern double common_get_wtime(void);

#endif // COMMON_H

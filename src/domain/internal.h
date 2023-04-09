#if !defined(DOMAIN_INTERNAL_H)
#define DOMAIN_INTERNAL_H

#include "common.h"

extern int domain_load(const char dirname[], domain_t *domain);

extern sdecomp_info_t *optimise_sdecomp_init(const bool uniformx, const int *glsizes);

#endif // DOMAIN_INTERNAL_H

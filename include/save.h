#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"

// next time to trigger output
extern double save_next;

extern int save(const domain_t *domain, const int step, const double time, const fluid_t *fluid);

#endif // SAVE_H

#if !defined(LOGGING_H)
#define LOGGING_H

#include "domain.h"
#include "fluid.h"

// next time to trigger output
extern double logging_next;

extern int logging(const domain_t *domain, const int step, const double time, const double dt, const double wtime, const fluid_t *fluid);

#endif // LOGGING_H

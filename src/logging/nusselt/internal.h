#if !defined(LOGGING_NUSSELT_INTERNAL_H)
#define LOGGING_NUSSELT_INTERNAL_H

#include "domain.h"
#include "fluid.h"

extern void logging_nusselt_compute_torque(const domain_t *domain, const fluid_t *fluid, double retvals[2]);
extern void logging_nusselt_compute_dissipation(const domain_t *domain, const fluid_t *fluid, double *retvals);

extern double logging_nusselt_compute_torque_laminar(const double diffusivity, const double ri, const double ro, const double ui, const double uo, const double ly, const double lz);

#endif // LOGGING_NUSSELT_INTERNAL_H

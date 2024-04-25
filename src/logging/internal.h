#if !defined(LOGGING_INTERNAL_H)
#define LOGGING_INTERNAL_H

extern int logging_check_divergence(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_injection(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_dissipation(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

#endif // LOGGING_INTERNAL_H

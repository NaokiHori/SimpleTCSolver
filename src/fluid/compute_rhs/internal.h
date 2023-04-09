#if !defined(FLUID_COMPUTE_RHS_INTERNAL)
#define FLUID_COMPUTE_RHS_INTERNAL

extern int fluid_compute_rhs_ux(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid);
extern int fluid_compute_rhs_uy(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid);
extern int fluid_compute_rhs_uz(const domain_t * restrict domain, const int rkstep, fluid_t * restrict fluid);

#endif // FLUID_COMPUTE_RHS_INTERNAL

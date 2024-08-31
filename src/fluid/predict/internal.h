#if !defined(FLUID_INTEGRATE_INTERNAL)
#define FLUID_INTEGRATE_INTERNAL

#include "linear_system.h"
#include "fluid.h"

// store approximation of laplacian
typedef struct {
  double l;
  double c;
  double u;
} laplacian_t;

// store Laplacian for each directoin
typedef struct {
  bool is_initialised;
  laplacian_t * restrict lapx;
  laplacian_t * restrict lapy;
  laplacian_t lapz;
} laplacians_t;

extern int compute_lxx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lxy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lxz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lyx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lyy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lyz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lzx(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lzy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_lzz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_ux(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_t(
    const domain_t * domain,
    fluid_t * fluid
);

extern int update_ux(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_uz(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_t(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int solve_in_x(
    const double prefactor,
    const laplacian_t * lapx,
    linear_system_t * linear_system
);

extern int solve_in_y(
    const double prefactor,
    const laplacian_t * lapy,
    linear_system_t * linear_system
);

extern int solve_in_z(
    const double prefactor,
    const laplacian_t * lapz,
    linear_system_t * linear_system
);

#endif // FLUID_INTEGRATE_INTERNAL

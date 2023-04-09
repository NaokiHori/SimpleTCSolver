#if !defined(COMPUTE_POTENTIAL_INTERNAL_H)
#define COMPUTE_POTENTIAL_INTERNAL_H

/* FOR INTERNAL USE */
/* DO NOT INCLUDE THIS HEADER OUTSIDE comnpute_potential */

#include <stdbool.h>
#include "domain.h"
#include "tdm.h"

// structure only used to solve Poisson equation
// NOTE: this is shared among normal and efficient solvers
//   thus some variables may not be used
typedef struct {
  void *buf0;
  void *buf1;
  fftw_plan fftw_plan_x[2];
  fftw_plan fftw_plan_y[2];
  fftw_plan fftw_plan_z[2];
  size_t tdm_sizes[NDIMS];
  tdm_info_t *tdm_info;
  double *evalys, *evalzs;
  sdecomp_transpose_plan_t *r_transposer_x1_to_y1;
  sdecomp_transpose_plan_t *r_transposer_y1_to_x1;
  sdecomp_transpose_plan_t *c_transposer_x1_to_y1;
  sdecomp_transpose_plan_t *c_transposer_y1_to_x1;
  sdecomp_transpose_plan_t *c_transposer_y1_to_z1;
  sdecomp_transpose_plan_t *c_transposer_z1_to_y1;
  sdecomp_transpose_plan_t *c_transposer_z1_to_x2;
  sdecomp_transpose_plan_t *c_transposer_x2_to_z1;
} poisson_solver_t;

extern poisson_solver_t *poisson_solver;

extern int fluid_compute_potential_init(const domain_t *domain);
extern int fluid_compute_potential_finalise(void);

#endif // COMPUTE_POTENTIAL_INTERNAL_H

#include <fftw3.h>
#include "common.h"
#include "sdecomp.h"
#include "tdm.h"
#include "internal.h"


static void wrap_sdecomp_transpose_destruct(sdecomp_transpose_plan_t *plan){
  if(NULL == plan) return;
  sdecomp.transpose.destruct(plan);
}

/**
 * @brief destruct a structure poisson_solver_t
 * @return : error code
 */
int fluid_compute_potential_finalise(void){
  // allocated via fftw_malloc, should be freed correspondingly
  fftw_free(poisson_solver->buf0);
  fftw_free(poisson_solver->buf1);
  tdm.destruct(poisson_solver->tdm_info);
  common_free(poisson_solver->evalys);
  common_free(poisson_solver->evalzs);
  wrap_sdecomp_transpose_destruct(poisson_solver->r_transposer_x1_to_y1);
  wrap_sdecomp_transpose_destruct(poisson_solver->r_transposer_y1_to_x1);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_x1_to_y1);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_y1_to_x1);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_y1_to_z1);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_z1_to_y1);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_z1_to_x2);
  wrap_sdecomp_transpose_destruct(poisson_solver->c_transposer_x2_to_z1);
  fftw_destroy_plan(poisson_solver->fftw_plan_x[0]);
  fftw_destroy_plan(poisson_solver->fftw_plan_x[1]);
  fftw_destroy_plan(poisson_solver->fftw_plan_y[0]);
  fftw_destroy_plan(poisson_solver->fftw_plan_y[1]);
  fftw_destroy_plan(poisson_solver->fftw_plan_z[0]);
  fftw_destroy_plan(poisson_solver->fftw_plan_z[1]);
  fftw_cleanup();
  common_free(poisson_solver);
  return 0;
}


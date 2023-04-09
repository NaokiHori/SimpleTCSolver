#include "sdecomp.h"
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "linear_system.h"
#include "internal.h"
#include "arrays/domain/uxlapx.h"
#include "arrays/domain/uxlapy.h"
#include "arrays/fluid/ux.h"
#include "../arrays/srcuxa.h"
#include "../arrays/srcuxb.h"
#include "../arrays/srcuxg.h"

static linear_system_t * restrict linear_system = NULL;

static int compute_delta(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const double alpha = RKCOEFS[rkstep].alpha;
  const double beta  = RKCOEFS[rkstep].beta;
  const double gamma = RKCOEFS[rkstep].gamma;
  const double adt = alpha * dt;
  const double bdt = beta  * dt;
  const double gdt = gamma * dt;
  const double * restrict srcuxa = fluid->srcuxa->data;
  const double * restrict srcuxb = fluid->srcuxb->data;
  const double * restrict srcuxg = fluid->srcuxg->data;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict dux = linear_system->x1pncl;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        *(dux + (cnt++)) =
          + adt * SRCUXA(i, j, k)
          + bdt * SRCUXB(i, j, k)
          + gdt * SRCUXG(i, j, k);
      }
    }
  }
  return 0;
}

static int solve_in_x(const domain_t * restrict domain, const double prefactor){
  tdm_info_t * restrict tdm_info = linear_system->tdm_x;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t * restrict uxlapx = domain->uxlapx;
  const int isize = linear_system->x1pncl_mysizes[0];
  for(int i = 2; i <= isize; i++){
    tdm_l[i-2] =    - prefactor * UXLAPX(i).l;
    tdm_c[i-2] = 1. - prefactor * UXLAPX(i).c;
    tdm_u[i-2] =    - prefactor * UXLAPX(i).u;
  }
  tdm.solve(tdm_info, linear_system->x1pncl);
  return 0;
}

static int solve_in_y(const domain_t * restrict domain, const double prefactor){
  tdm_info_t * restrict tdm_info = linear_system->tdm_y;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t * restrict uxlapy = domain->uxlapy;
  const int isize = linear_system->y1pncl_mysizes[0];
  const int jsize = linear_system->y1pncl_mysizes[1];
  const int ioffs = linear_system->y1pncl_offsets[0];
  for(int i = ioffs + 2; i <= isize + ioffs; i++){
    for(int j = 0; j < jsize; j++){
      tdm_l[j] =    - prefactor * UXLAPY(i).l;
      tdm_c[j] = 1. - prefactor * UXLAPY(i).c;
      tdm_u[j] =    - prefactor * UXLAPY(i).u;
    }
    tdm.solve(tdm_info, linear_system->y1pncl);
  }
  return 0;
}

static int solve_in_z(const domain_t * restrict domain, const double prefactor){
  tdm_info_t * restrict tdm_info = linear_system->tdm_z;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t uxlapz = domain->lapz;
  const int ksize = linear_system->z2pncl_mysizes[2];
  for(int k = 0; k < ksize; k++){
    tdm_l[k] =    - prefactor * uxlapz.l;
    tdm_c[k] = 1. - prefactor * uxlapz.c;
    tdm_u[k] =    - prefactor * uxlapz.u;
  }
  tdm.solve(tdm_info, linear_system->z2pncl);
  return 0;
}

static int update_field(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dux = linear_system->x1pncl;
  double * restrict ux = fluid->ux->data;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        UX(i, j, k) += *(dux + (cnt++));
      }
    }
  }
  fluid_update_boundaries_ux(domain, ux);
  return 0;
}

/**
 * @brief update ux
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int fluid_update_velocity_ux(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  if(linear_system == NULL){
    const size_t glsizes[NDIMS] = {
      (size_t)(domain->glsizes[0]-1),
      (size_t)(domain->glsizes[1]  ),
      (size_t)(domain->glsizes[2]  ),
    };
    linear_system = init_linear_system(domain->info, glsizes);
  }
  // compute increments
  compute_delta(domain, rkstep, dt, fluid);
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * RKCOEFS[rkstep].gamma * dt * fluid->diffusivity;
  if(config.get.implicitx()){
    solve_in_x(domain, prefactor);
  }
  if(config.get.implicity()){
    sdecomp.transpose.execute(
        linear_system->transposer_x1_to_y1,
        linear_system->x1pncl,
        linear_system->y1pncl
    );
    solve_in_y(domain, prefactor);
    sdecomp.transpose.execute(
        linear_system->transposer_y1_to_x1,
        linear_system->y1pncl,
        linear_system->x1pncl
    );
  }
  if(config.get.implicitz()){
    sdecomp.transpose.execute(
        linear_system->transposer_x1_to_z2,
        linear_system->x1pncl,
        linear_system->z2pncl
    );
    solve_in_z(domain, prefactor);
    sdecomp.transpose.execute(
        linear_system->transposer_z2_to_x1,
        linear_system->z2pncl,
        linear_system->x1pncl
    );
  }
  // the field is actually updated here
  update_field(domain, fluid);
  return 0;
}

int fluid_update_velocity_finalise_ux(void){
  linear_system_finalise(linear_system);
  return 0;
}


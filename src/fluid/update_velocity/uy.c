#include "sdecomp.h"
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "linear_system.h"
#include "internal.h"
#include "arrays/domain/uylapx.h"
#include "arrays/domain/uylapy.h"
#include "arrays/fluid/uy.h"
#include "../arrays/srcuya.h"
#include "../arrays/srcuyb.h"
#include "../arrays/srcuyg.h"

static linear_system_t *linear_system = NULL;

static int compute_delta(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const double alpha = RKCOEFS[rkstep].alpha;
  const double beta  = RKCOEFS[rkstep].beta;
  const double gamma = RKCOEFS[rkstep].gamma;
  const double adt = alpha * dt;
  const double bdt = beta  * dt;
  const double gdt = gamma * dt;
  const double * restrict srcuya = fluid->srcuya->data;
  const double * restrict srcuyb = fluid->srcuyb->data;
  const double * restrict srcuyg = fluid->srcuyg->data;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict duy = linear_system->x1pncl;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        *(duy + (cnt++)) =
          + adt * SRCUYA(i, j, k)
          + bdt * SRCUYB(i, j, k)
          + gdt * SRCUYG(i, j, k);
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
  const laplacian_t * restrict uylapx = domain->uylapx;
  const int isize = linear_system->x1pncl_mysizes[0];
  for(int i = 1; i <= isize; i++){
    tdm_l[i-1] =    - prefactor * UYLAPX(i).l;
    tdm_c[i-1] = 1. - prefactor * UYLAPX(i).c;
    tdm_u[i-1] =    - prefactor * UYLAPX(i).u;
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
  const laplacian_t * restrict uylapy = domain->uylapy;
  const int isize = linear_system->y1pncl_mysizes[0];
  const int jsize = linear_system->y1pncl_mysizes[1];
  const int ioffs = linear_system->y1pncl_offsets[0];
  for(int i = ioffs + 1; i <= isize + ioffs; i++){
    for(int j = 0; j < jsize; j++){
      tdm_l[j] =    - prefactor * UYLAPY(i).l;
      tdm_c[j] = 1. - prefactor * UYLAPY(i).c;
      tdm_u[j] =    - prefactor * UYLAPY(i).u;
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
  const laplacian_t uylapz = domain->lapz;
  const int ksize = linear_system->z2pncl_mysizes[2];
  for(int k = 0; k < ksize; k++){
    tdm_l[k] =    - prefactor * uylapz.l;
    tdm_c[k] = 1. - prefactor * uylapz.c;
    tdm_u[k] =    - prefactor * uylapz.u;
  }
  tdm.solve(tdm_info, linear_system->z2pncl);
  return 0;
}

static int update_field(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict duy = linear_system->x1pncl;
  double * restrict uy = fluid->uy->data;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UY(i, j, k) += *(duy + (cnt++));
      }
    }
  }
  fluid_update_boundaries_uy(domain, uy);
  return 0;
}

/**
 * @brief update uy
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int fluid_update_velocity_uy(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  if(linear_system == NULL){
    const size_t glsizes[NDIMS] = {
      (size_t)(domain->glsizes[0]),
      (size_t)(domain->glsizes[1]),
      (size_t)(domain->glsizes[2]),
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

int fluid_update_velocity_finalise_uy(void){
  linear_system_finalise(linear_system);
  return 0;
}


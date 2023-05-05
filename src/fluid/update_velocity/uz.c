#include "param.h"
#include "sdecomp.h"
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "linear_system.h"
#include "internal.h"
#include "arrays/domain/uzlapx.h"
#include "arrays/domain/uzlapy.h"
#include "arrays/fluid/uz.h"
#include "../arrays/srcuza.h"
#include "../arrays/srcuzb.h"
#include "../arrays/srcuzg.h"

static linear_system_t *linear_system = NULL;

static int compute_delta(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const double alpha = param_rkcoefs[rkstep].alpha;
  const double beta  = param_rkcoefs[rkstep].beta;
  const double gamma = param_rkcoefs[rkstep].gamma;
  const double adt = alpha * dt;
  const double bdt = beta  * dt;
  const double gdt = gamma * dt;
  const double * restrict srcuza = fluid->srcuza->data;
  const double * restrict srcuzb = fluid->srcuzb->data;
  const double * restrict srcuzg = fluid->srcuzg->data;
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict duz = linear_system->x1pncl;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        *(duz + (cnt++)) =
          + adt * SRCUZA(i, j, k)
          + bdt * SRCUZB(i, j, k)
          + gdt * SRCUZG(i, j, k);
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
  const laplacian_t * restrict uzlapx = domain->uzlapx;
  const int isize = linear_system->x1pncl_mysizes[0];
  for(int i = 1; i <= isize; i++){
    tdm_l[i-1] =    - prefactor * UZLAPX(i).l;
    tdm_c[i-1] = 1. - prefactor * UZLAPX(i).c;
    tdm_u[i-1] =    - prefactor * UZLAPX(i).u;
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
  const laplacian_t * restrict uzlapy = domain->uzlapy;
  const int isize = linear_system->y1pncl_mysizes[0];
  const int jsize = linear_system->y1pncl_mysizes[1];
  const int ioffs = linear_system->y1pncl_offsets[0];
  for(int i = ioffs + 1; i <= isize + ioffs; i++){
    for(int j = 0; j < jsize; j++){
      tdm_l[j] =    - prefactor * UZLAPY(i).l;
      tdm_c[j] = 1. - prefactor * UZLAPY(i).c;
      tdm_u[j] =    - prefactor * UZLAPY(i).u;
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
  const laplacian_t uzlapz = domain->lapz;
  const int ksize = linear_system->z2pncl_mysizes[2];
  for(int k = 0; k < ksize; k++){
    tdm_l[k] =    - prefactor * uzlapz.l;
    tdm_c[k] = 1. - prefactor * uzlapz.c;
    tdm_u[k] =    - prefactor * uzlapz.u;
  }
  tdm.solve(tdm_info, linear_system->z2pncl);
  return 0;
}

static int update_field(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict duz = linear_system->x1pncl;
  double * restrict uz = fluid->uz->data;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UZ(i, j, k) += *(duz + (cnt++));
      }
    }
  }
  fluid_update_boundaries_uz(domain, uz);
  return 0;
}

/**
 * @brief update uz
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int fluid_update_velocity_uz(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
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
    0.5 * param_rkcoefs[rkstep].gamma * dt * fluid->diffusivity;
  if(param_implicit_x){
    solve_in_x(domain, prefactor);
  }
  if(param_implicit_y){
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
  if(param_implicit_z){
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

int fluid_update_velocity_finalise_uz(void){
  linear_system_finalise(linear_system);
  return 0;
}


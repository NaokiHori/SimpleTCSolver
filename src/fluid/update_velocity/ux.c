#include <stdbool.h>
#include <math.h>
#include "sdecomp.h"
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "linear_system.h"
#include "tdm.h"
#include "arrays/fluid.h"
#include "internal.h"


// local variable to store buffers and transpose plans
//   to solve linear systems in each direction
static linear_system_t *linear_system = NULL;

static void solve_linear_systems_in_x(const domain_t * restrict domain, const double prefactor){
  // x1 pencil
  const int isize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, linear_system->glsizes[0]);
  const int jsize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, linear_system->glsizes[1]);
  const int ksize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, linear_system->glsizes[2]);
  // solve linear systems
  double * restrict tdm_l = linear_system->tdm_l;
  double * restrict tdm_c = linear_system->tdm_c;
  double * restrict tdm_u = linear_system->tdm_u;
  const laplace_t * restrict uxdifx = domain->uxdifx;
  // set diagonal components of linear system
  for(int i = 0; i < isize; i++){
    // N.B. UXDIFX(2) <-> i = 0
    tdm_l[i] =    - prefactor * UXDIFX(i+2).l;
    tdm_c[i] = 1. - prefactor * UXDIFX(i+2).c;
    tdm_u[i] =    - prefactor * UXDIFX(i+2).u;
  }
  tdm_solve_double(
      /* size of system  */ isize,
      /* number of rhs   */ jsize * ksize,
      /* periodic        */ false,
      /* lower- diagonal */ tdm_l,
      /* center-diagonal */ tdm_c,
      /* upper- diagonal */ tdm_u,
      /* input / output  */ linear_system->x1pncl
  );
}

static void solve_linear_systems_in_y(const domain_t * restrict domain, const double prefactor){
  // y1 pencil
  const int isize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, linear_system->glsizes[0]);
  const int jsize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, linear_system->glsizes[1]);
  const int ksize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_ZDIR, linear_system->glsizes[2]);
  const int ioffs = sdecomp_get_pencil_offset(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, linear_system->glsizes[0]);
  // solve linear systems
  double * restrict tdm_l = linear_system->tdm_l;
  double * restrict tdm_c = linear_system->tdm_c;
  double * restrict tdm_u = linear_system->tdm_u;
  const laplace_t * restrict uxdify = domain->uxdify;
  // for each x (this is needed since Laplacian varies in x)
  for(int i = ioffs; i < isize + ioffs; i++){
    // set diagonal components of linear system
    for(int j = 0; j < jsize; j++){
      // N.B. UXDIFY(2) <-> i = 0
      tdm_l[j] =    - prefactor * UXDIFY(i+2).l;
      tdm_c[j] = 1. - prefactor * UXDIFY(i+2).c;
      tdm_u[j] =    - prefactor * UXDIFY(i+2).u;
    }
    tdm_solve_double(
        /* size of system  */ jsize,
        /* number of rhs   */ ksize,
        /* periodic        */ true,
        /* lower- diagonal */ tdm_l,
        /* center-diagonal */ tdm_c,
        /* upper- diagonal */ tdm_u,
        /* input / output  */ linear_system->y1pncl
    );
  }
}

static void solve_linear_systems_in_z(const domain_t * restrict domain, const double prefactor){
  // z1 pencil
  const int isize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Z1PENCIL, SDECOMP_XDIR, linear_system->glsizes[0]);
  const int jsize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Z1PENCIL, SDECOMP_YDIR, linear_system->glsizes[1]);
  const int ksize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Z1PENCIL, SDECOMP_ZDIR, linear_system->glsizes[2]);
  // solve linear systems
  double * restrict tdm_l = linear_system->tdm_l;
  double * restrict tdm_c = linear_system->tdm_c;
  double * restrict tdm_u = linear_system->tdm_u;
  const laplace_t uxdifz = domain->uxdifz;
  // set diagonal components of linear system
  for(int k = 0; k < ksize; k++){
    tdm_l[k] =    - prefactor * uxdifz.l;
    tdm_c[k] = 1. - prefactor * uxdifz.c;
    tdm_u[k] =    - prefactor * uxdifz.u;
  }
  tdm_solve_double(
      /* size of system  */ ksize,
      /* number of rhs   */ isize * jsize,
      /* periodic        */ true,
      /* lower- diagonal */ tdm_l,
      /* center-diagonal */ tdm_c,
      /* upper- diagonal */ tdm_u,
      /* input / output  */ linear_system->z1pncl
  );
}

/**
 * @brief update ux
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return              : error code
 */
int fluid_update_velocity_ux(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // buffer to store \delta u_x
#define DUX(I, J, K) (dux[ ((K)-1) * (jsize) * (isize-1) + ((J)-1) * (isize-1) + ((I)-2) ])
  /* ! initalise linear solver ! 8 ! */
  if(linear_system == NULL){
    // size of linear system to be solved for "ux": (glisize-1) x gljsize
    const int glisize = domain->glsizes[0];
    const int gljsize = domain->glsizes[1];
    const int glksize = domain->glsizes[2];
    const int glsizes_ux[NDIMS] = {glisize-1, gljsize, glksize};
    linear_system = init_linear_system(domain->sdecomp, glsizes_ux);
  }
  /* ! compute increments ! 22 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double alpha = RKCOEFS[rkstep].alpha;
    const double beta  = RKCOEFS[rkstep].beta;
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict srcuxa = fluid->srcuxa;
    const double * restrict srcuxb = fluid->srcuxb;
    const double * restrict srcuxg = fluid->srcuxg;
    double * restrict dux = linear_system->x1pncl;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          DUX(i, j, k) =
            + alpha * dt * SRCUXA(i, j, k)
            + beta  * dt * SRCUXB(i, j, k)
            + gamma * dt * SRCUXG(i, j, k);
        }
      }
    }
  }
  /* ! pre-factor in front of coefficients ! 6 ! */
  double prefactor;
  {
    const double gamma = RKCOEFS[rkstep].gamma;
    const double Re = config.get_double("Re");
    prefactor = (gamma * dt) / (2. * Re);
  }
  /* ! solve linear system in x direction ! 3 ! */
  if(config.get_bool("implicitx")){
    solve_linear_systems_in_x(domain, prefactor);
  }
  /* ! transpose x1pencil to y1pencil ! 8 ! */
  const bool needs_x_y_transpose = config.get_bool("implicity") || config.get_bool("implicitz");
  if(needs_x_y_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_x1_to_y1,
        linear_system->x1pncl,
        linear_system->y1pncl
    );
  }
  /* ! solve linear system in y direction ! 3 ! */
  if(config.get_bool("implicity")){
    solve_linear_systems_in_y(domain, prefactor);
  }
  /* ! transpose y1pencil to z1pencil ! 8 ! */
  const bool needs_y_z_transpose = config.get_bool("implicitz");
  if(needs_y_z_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_y1_to_z1,
        linear_system->y1pncl,
        linear_system->z1pncl
    );
  }
  /* ! solve linear system in z direction ! 3 ! */
  if(config.get_bool("implicitz")){
    solve_linear_systems_in_z(domain, prefactor);
  }
  /* ! transpose z1pencil to y1pencil ! 7 ! */
  if(needs_y_z_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_z1_to_y1,
        linear_system->z1pncl,
        linear_system->y1pncl
    );
  }
  /* ! transpose y1pencil to x1pencil ! 7 ! */
  if(needs_x_y_transpose){
    sdecomp_transpose_execute(
        linear_system->transposer_y1_to_x1,
        linear_system->y1pncl,
        linear_system->x1pncl
    );
  }
  /* ! update ! 15 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double * restrict dux = linear_system->x1pncl;
    double * restrict ux = fluid->ux;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          UX(i, j, k) += DUX(i, j, k);
        }
      }
    }
    fluid_update_boundaries_ux(domain, ux);
  }
#undef DUX
  return 0;
}

int fluid_update_velocity_finalise_ux(void){
  linear_system_finalise(linear_system);
  return 0;
}


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
  const laplace_t * restrict uzdifx = domain->uzdifx;
  // set diagonal components of linear system
  for(int i = 0; i < isize; i++){
    // N.B. UZDIFX(1) <-> i = 0
    tdm_l[i] =    - prefactor * UZDIFX(i+1).l;
    tdm_c[i] = 1. - prefactor * UZDIFX(i+1).c;
    tdm_u[i] =    - prefactor * UZDIFX(i+1).u;
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
  const laplace_t * restrict uzdify = domain->uzdify;
  // for each x (this is needed since Laplacian varies in x)
  for(int i = ioffs; i < isize + ioffs; i++){
    // set diagonal components of linear system
    for(int j = 0; j < jsize; j++){
      // N.B. UZDIFY(1) <-> i = 0
      tdm_l[j] =    - prefactor * UZDIFY(i+1).l;
      tdm_c[j] = 1. - prefactor * UZDIFY(i+1).c;
      tdm_u[j] =    - prefactor * UZDIFY(i+1).u;
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
  const laplace_t uzdifz = domain->uzdifz;
  // set diagonal components of linear system
  for(int k = 0; k < ksize; k++){
    tdm_l[k] =    - prefactor * uzdifz.l;
    tdm_c[k] = 1. - prefactor * uzdifz.c;
    tdm_u[k] =    - prefactor * uzdifz.u;
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
 * @brief update uz
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return              : error code
 */
int fluid_update_velocity_uz(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // buffer to store \delta u_x
#define DUZ(I, J, K) (duz[ ((K)-1) * (jsize) * (isize) + ((J)-1) * (isize) + ((I)-1) ])
  /* ! initalise linear solver ! 8 ! */
  if(linear_system == NULL){
    // size of linear system to be solved for "uz": glisize x gljsize x glksize
    const int glisize = domain->glsizes[0];
    const int gljsize = domain->glsizes[1];
    const int glksize = domain->glsizes[2];
    const int glsizes_uz[NDIMS] = {glisize, gljsize, glksize};
    linear_system = init_linear_system(domain->sdecomp, glsizes_uz);
  }
  /* ! compute increments ! 22 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double alpha = RKCOEFS[rkstep].alpha;
    const double beta  = RKCOEFS[rkstep].beta;
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict srcuza = fluid->srcuza;
    const double * restrict srcuzb = fluid->srcuzb;
    const double * restrict srcuzg = fluid->srcuzg;
    double * restrict duz = linear_system->x1pncl;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          DUZ(i, j, k) =
            + alpha * dt * SRCUZA(i, j, k)
            + beta  * dt * SRCUZB(i, j, k)
            + gamma * dt * SRCUZG(i, j, k);
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
    const double * restrict duz = linear_system->x1pncl;
    double * restrict uz = fluid->uz;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          UZ(i, j, k) += DUZ(i, j, k);
        }
      }
    }
    fluid_update_boundaries_uz(domain, uz);
  }
#undef DUZ
  return 0;
}

int fluid_update_velocity_finalise_uz(void){
  linear_system_finalise(linear_system);
  return 0;
}


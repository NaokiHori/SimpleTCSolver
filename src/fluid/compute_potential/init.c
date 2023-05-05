#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <fftw3.h>
#include "common.h"
#include "domain.h"
#include "sdecomp.h"
#include "tdm.h"
#include "arrays/domain/plapx.h"
#include "internal.h"

// C99 does not specify M_PI in math.h
#if !defined(M_PI)
#define M_PI 3.1415926535897932
#endif

/* initialse Poisson solver */
// several pencils for different data types are treated
//   and thus this source is very complicated
// used prefixes are as follows:
//   r_: real    (double)       type
//   c_: complex (fftw_complex) type
//
//   gl_    : global array size (not pencils)
//   x1pncl_: each pencl (x1, y1, ...)

/* size of domain and pencils */
// NOTE: define globally to reduce the number of arguments of static functions
// global domain size in real space
static size_t r_gl_sizes[NDIMS] = {0};
// global domain size in complex space (only used by normal solver)
static size_t c_gl_sizes[NDIMS] = {0};
// local domain size (x1 pencil) in real space
static size_t r_x1pncl_sizes[NDIMS] = {0};
static size_t r_x1pncl_bytes = sizeof(double);
// local domain size (y1 pencil) in real space
static size_t r_y1pncl_sizes[NDIMS] = {0};
static size_t r_y1pncl_bytes = sizeof(double);
// local domain size (y1 pencil) in complex space (only used by normal solver)
static size_t c_y1pncl_sizes[NDIMS] = {0};
static size_t c_y1pncl_bytes = sizeof(fftw_complex);
// local domain size (z1 pencil) in complex space
static size_t c_z1pncl_sizes[NDIMS] = {0};
static size_t c_z1pncl_bytes = sizeof(fftw_complex);
// local domain size (x2 pencil) in complex space (only used by normal solver)
static size_t c_x2pncl_sizes[NDIMS] = {0};
static size_t c_x2pncl_bytes = sizeof(fftw_complex);

static int report_failure(const char type[]){
  // function to just dump error message and abort
  fprintf(stderr, "Poisson solver, initialisation failed: %s\n", type);
  fprintf(stderr, "  FFTW:    A possible reason is you link Intel-MKL lib\n");
  fprintf(stderr, "           Make sure you use FFTW3 directly,\n");
  fprintf(stderr, "           NOT its crazy wrapper offered by MKL\n");
  fprintf(stderr, "  SDECOMP: Check sdecomp.log and check arguments\n");
  fprintf(stderr, "           If they are all correct, PLEASE CONTACT ME\n");
  fflush(stderr);
  return 0;
}

static size_t larger_of(const size_t val0, const size_t val1){
  if(val0 > val1){
    return val0;
  }else{
    return val1;
  }
}

static int compute_pencil_sizes(const domain_t *domain){
  // NOTE: those variables are defined globally at the top of this file
  //   to reduce the nhumber of arguments which functions take
  // global domain size in real space
  const sdecomp_info_t *info = domain->info;
  r_gl_sizes[0] = (size_t)(domain->glsizes[0]        );
  r_gl_sizes[1] = (size_t)(domain->glsizes[1]        );
  r_gl_sizes[2] = (size_t)(domain->glsizes[2]        );
  // global domain size in complex space (only used by normal solver)
  c_gl_sizes[0] = (size_t)(domain->glsizes[0]        );
  c_gl_sizes[1] = (size_t)(domain->glsizes[1] / 2 + 1); // Hermite symmetry
  c_gl_sizes[2] = (size_t)(domain->glsizes[2]        );
  // local domain size (x1 pencil) in real space
  for(int dim = 0; dim < NDIMS; dim++){
    size_t size = 0;
    sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_gl_sizes[dim], &size);
    r_x1pncl_sizes[dim] = size;
    r_x1pncl_bytes *= size;
  }
  // local domain size (y1 pencil) in real space
  for(int dim = 0; dim < NDIMS; dim++){
    size_t size = 0;
    sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_gl_sizes[dim], &size);
    r_y1pncl_sizes[dim] = size;
    r_y1pncl_bytes *= size;
  }
  // local domain size (y1 pencil) in complex space (only used by normal solver)
  for(int dim = 0; dim < NDIMS; dim++){
    size_t size = 0;
    sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_gl_sizes[dim], &size);
    c_y1pncl_sizes[dim] = size;
    c_y1pncl_bytes *= size;
  }
  // local domain size (z1 pencil) in complex space
  for(int dim = 0; dim < NDIMS; dim++){
    size_t size = 0;
    sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, c_gl_sizes[dim], &size);
    c_z1pncl_sizes[dim] = size;
    c_z1pncl_bytes *= size;
  }
  // local domain size (x2 pencil) in complex space (only used by normal solver)
  for(int dim = 0; dim < NDIMS; dim++){
    size_t size = 0;
    sdecomp.get_pencil_mysize(info, SDECOMP_X2PENCIL, dim, c_gl_sizes[dim], &size);
    c_x2pncl_sizes[dim] = size;
    c_x2pncl_bytes *= size;
  }
  return 0;
}

static int allocate_buffers(void){
  // two buffers are enough to do the job,
  //   which are allocated here
  void **buf0 = &poisson_solver->buf0;
  void **buf1 = &poisson_solver->buf1;
  size_t buf0_bytes = 0;
  size_t buf1_bytes = 0;
  // r_x1pncl -> rotate -> r_y1pncl -> FFT -> c_y1pncl -> rotate -> c_z1pncl -> FFT -> c_z1pncl -> rotate -> c_x2pncl
  // buffer0               buffer1            buffer0               buffer1            buffer0               buffer1
  buf0_bytes = larger_of(larger_of(r_x1pncl_bytes, c_y1pncl_bytes), c_z1pncl_bytes);
  buf1_bytes = larger_of(larger_of(r_y1pncl_bytes, c_z1pncl_bytes), c_x2pncl_bytes);
  // allocate them using fftw_malloc to enforce them 16bit-aligned for SIMD
  *buf0 = fftw_malloc(buf0_bytes);
  *buf1 = fftw_malloc(buf1_bytes);
  if(NULL == *buf0 || NULL == *buf1){
    fprintf(stderr, "FATAL: fftw_malloc failed\n");
    fflush(stderr);
    return 1;
  }
  return 0;
}

static int init_tri_diagonal_solver(const domain_t *domain){
  // N x N tri-diagonal matrix,
  //   which are solved for My x Mz times
  // tdm_sizes[0] = N, tdm_sizes[1] = My, tdm_sizes[2] = Mz
  // since lower- and upper-diagonal components are
  //   independent to y, z directions and time,
  //   we compute here and re-use them
  // center-diagonal components are, on the other hand,
  //   dependent on time and thus needs to compute everytime
  //   in the solver
  size_t *tdm_sizes = poisson_solver->tdm_sizes;
  tdm_info_t **tdm_info = &poisson_solver->tdm_info;
  // in x: d^2p / dx^2 = q
  tdm_sizes[0] = c_x2pncl_sizes[0];
  tdm_sizes[1] = c_x2pncl_sizes[1];
  tdm_sizes[2] = c_x2pncl_sizes[2];
  tdm.construct(
    /* size of system */ (int)tdm_sizes[0],
    /* number of rhs  */ 1,
    /* is periodic    */ false,
    /* is complex     */ true,
    /* output         */ tdm_info
  );
  // initialise tri-diagonal matrix in x direction
  double *tdm_l = NULL;
  double *tdm_u = NULL;
  tdm.get_l(*tdm_info, &tdm_l);
  tdm.get_u(*tdm_info, &tdm_u);
  const laplacian_t *plapx = domain->plapx;
  for(size_t i = 1; i <= tdm_sizes[0]; i++){
    tdm_l[i-1] = PLAPX(i  ).l;
    tdm_u[i-1] = PLAPX(i  ).u;
  }
  return 0;
}

static int init_pencil_rotations(const domain_t *domain){
  const sdecomp_info_t *info = domain->info;
  const size_t r_dsize = sizeof(      double);
  const size_t c_dsize = sizeof(fftw_complex);
  // between x1 and y1 (real)
  {
    sdecomp_transpose_plan_t **plan0 = &poisson_solver->r_transposer_x1_to_y1;
    sdecomp_transpose_plan_t **plan1 = &poisson_solver->r_transposer_y1_to_x1;
    sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_gl_sizes, r_dsize, plan0);
    sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_gl_sizes, r_dsize, plan1);
    if(NULL == *plan0){
      report_failure("SDECOMP x1 to y1 for real");
      return 1;
    }
    if(NULL == *plan1){
      report_failure("SDECOMP y1 to x1 for real");
      return 1;
    }
  }
  // between y1 and z1 (complex)
  {
    sdecomp_transpose_plan_t **plan0 = &poisson_solver->c_transposer_y1_to_z1;
    sdecomp_transpose_plan_t **plan1 = &poisson_solver->c_transposer_z1_to_y1;
    sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, c_dsize, plan0);
    sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, c_dsize, plan1);
    if(NULL == *plan0){
      report_failure("SDECOMP y1 to z1 for complex");
      return 1;
    }
    if(NULL == *plan1){
      report_failure("SDECOMP z1 to y1 for complex");
      return 1;
    }
  }
  // between z1 and x2 (complex)
  {
    sdecomp_transpose_plan_t **plan0 = &poisson_solver->c_transposer_z1_to_x2;
    sdecomp_transpose_plan_t **plan1 = &poisson_solver->c_transposer_x2_to_z1;
    sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_X2PENCIL, c_gl_sizes, c_dsize, plan0);
    sdecomp.transpose.construct(info, SDECOMP_X2PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, c_dsize, plan1);
    if(NULL == *plan0){
      report_failure("SDECOMP z1 to x2 for complex");
      return 1;
    }
    if(NULL == *plan1){
      report_failure("SDECOMP x2 to z1 for complex");
      return 1;
    }
  }
  return 0;
}

static int init_ffts(void){
  const unsigned flags = FFTW_PATIENT | FFTW_DESTROY_INPUT;
  // NOTE: two buffers should be properly given
  //   see "allocate_buffers" above
  // y, real / complex
  {
    fftw_plan *fplan = &poisson_solver->fftw_plan_y[0];
    fftw_plan *bplan = &poisson_solver->fftw_plan_y[1];
    const int r_signal_length = (int)(r_y1pncl_sizes[SDECOMP_YDIR]);
    const int c_signal_length = (int)(c_y1pncl_sizes[SDECOMP_YDIR]);
    const int repeat_for = (int)(r_y1pncl_sizes[SDECOMP_ZDIR] * r_y1pncl_sizes[SDECOMP_XDIR]);
    *fplan = fftw_plan_many_dft_r2c(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        flags
    );
    *bplan = fftw_plan_many_dft_c2r(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        flags
    );
    if(NULL == *fplan){
      report_failure("FFTW y-forward");
      return 1;
    }
    if(NULL == *bplan){
      report_failure("FFTW y-backward");
      return 1;
    }
  }
  // z, complex to complex
  {
    const int signal_length = (int)(c_z1pncl_sizes[SDECOMP_ZDIR]);
    const int repeat_for = (int)(c_z1pncl_sizes[SDECOMP_XDIR] * c_z1pncl_sizes[SDECOMP_YDIR]);
    fftw_plan *fplan = &poisson_solver->fftw_plan_z[0];
    fftw_plan *bplan = &poisson_solver->fftw_plan_z[1];
    *fplan = fftw_plan_many_dft(
        1, &signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, signal_length,
        poisson_solver->buf0, NULL, 1, signal_length,
        FFTW_FORWARD, flags
    );
    *bplan = fftw_plan_many_dft(
        1, &signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, signal_length,
        poisson_solver->buf1, NULL, 1, signal_length,
        FFTW_BACKWARD, flags
    );
    if(NULL == *fplan){
      report_failure("FFTW z-forward");
      return 1;
    }
    if(NULL == *bplan){
      report_failure("FFTW z-backward");
      return 1;
    }
  }
  return 0;
}

static int init_eigenvalues(const domain_t *domain){
  const sdecomp_info_t *info = domain->info;
  double **evalys = &poisson_solver->evalys;
  double **evalzs = &poisson_solver->evalzs;
  // x2 pencil, DFTs in y and z directions
  const sdecomp_pencil_t pencil = SDECOMP_X2PENCIL;
  const double y_signal_length = 1. * domain->glsizes[1];
  const double z_signal_length = 1. * domain->glsizes[2];
  size_t myjsize = 0;
  size_t myksize = 0;
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], &myjsize);
  sdecomp.get_pencil_mysize(info, pencil, SDECOMP_ZDIR, c_gl_sizes[SDECOMP_ZDIR], &myksize);
  size_t joffset = 0;
  size_t koffset = 0;
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], &joffset);
  sdecomp.get_pencil_offset(info, pencil, SDECOMP_ZDIR, c_gl_sizes[SDECOMP_ZDIR], &koffset);
  const double dy = domain->dy;
  const double dz = domain->dz;
  // initialise eigenvalues in homogeneous directions
  *evalys = common_calloc(myjsize, sizeof(double));
  *evalzs = common_calloc(myksize, sizeof(double));
  for(size_t cnt = 0, j = joffset; j < myjsize + joffset; cnt++, j++){
    (*evalys)[cnt] =
      - 4. / pow(dy, 2.) * pow(
        sin( M_PI * j / y_signal_length ),
        2.
      );
  }
  for(size_t cnt = 0, k = koffset; k < myksize + koffset; cnt++, k++){
    (*evalzs)[cnt] =
      - 4. / pow(dz, 2.) * pow(
        sin( M_PI * k / z_signal_length ),
        2.
      );
  }
  return 0;
}

/**
 * @brief initialise variables which are used to solve Poisson equation
 * @param[in] domain : information about domain decomposition and domain size
 * @return           : (success) 0
 *                   : (failure) 1
 */
int fluid_compute_potential_init(const domain_t *domain){
  poisson_solver = common_calloc(1, sizeof(poisson_solver_t));
  // check domain size (global, local, pencils)
  compute_pencil_sizes(domain);
  // initialise each part of poisson_solver_t
  if(0 != allocate_buffers())               return 1;
  if(0 != init_tri_diagonal_solver(domain)) return 1;
  if(0 != init_pencil_rotations(domain))    return 1;
  if(0 != init_ffts())                      return 1;
  if(0 != init_eigenvalues(domain))         return 1;
  return 0;
}


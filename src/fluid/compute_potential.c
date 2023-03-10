#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "sdecomp.h"
#include "fluid.h"
#include "tdm.h"
#include "arrays/fluid.h"


static int tell_fftw_plan_creation_failure(const int line){
  // function to just dump error message and abort
  printf("Line %d: FFTW plan creation failed\n", line);
  printf("A possible reason is you link Intel-MKL lib\n");
  printf("Make sure you use FFTW3, NOT its wrapper by MKL\n");
  MPI_Abort(MPI_COMM_WORLD, 0);
  return 0;
}

// structure only used in this source,
//   which contains variables to solve Poisson equation
typedef struct {
  double *x1pncl_r, *y1pncl_r;
  fftw_complex *y1pncl_c, *z1pncl_c, *x2pncl_c;
  fftw_plan fftw_plan_fwrd_x, fftw_plan_bwrd_x;
  fftw_plan fftw_plan_fwrd_y, fftw_plan_bwrd_y;
  fftw_plan fftw_plan_fwrd_z, fftw_plan_bwrd_z;
  double *tdm_l, *tdm_c, *tdm_u;
  double *eigenvalues_x, *eigenvalues_y, *eigenvalues_z;
  sdecomp_transpose_t *transposer_x1_to_y1_r, *transposer_y1_to_x1_r;
  sdecomp_transpose_t *transposer_y1_to_z1_c, *transposer_z1_to_y1_c;
  sdecomp_transpose_t *transposer_z1_to_x2_c, *transposer_x2_to_z1_c;
} poisson_solver_t;

static poisson_solver_t *poisson_solver = NULL;

/**
 * @brief initialise utility variables which are used to solve Poisson equation
 * @param[in] domain : information about domain decomposition and size
 * @return           : error code
 */
static int init_poisson_solver(const domain_t * restrict domain){
  poisson_solver = common_calloc(1, sizeof(poisson_solver_t));
  const sdecomp_t *sdecomp = domain->sdecomp;
  /* ! compute size of each pencil in each direction ! 23 ! */
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  // x1 pencil, real
  const int x1pncl_isize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize    );
  const int x1pncl_jsize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize    );
  const int x1pncl_ksize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, glksize    );
  // y1 pencil, real
  const int y1pncl_isize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize    );
  const int y1pncl_jsize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, gljsize    );
  const int y1pncl_ksize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_ZDIR, glksize    );
  // y1 pencil, complex
  const int y1pncl_isize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize    );
  const int y1pncl_jsize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, gljsize/2+1);
  const int y1pncl_ksize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_ZDIR, glksize    );
  // z1 pencil, complex
  const int z1pncl_isize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_XDIR, glisize    );
  const int z1pncl_jsize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_YDIR, gljsize/2+1);
  const int z1pncl_ksize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_ZDIR, glksize    );
  // x2 pencil, complex
  const int x2pncl_isize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X2PENCIL, SDECOMP_XDIR, glisize    );
  const int x2pncl_jsize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X2PENCIL, SDECOMP_YDIR, gljsize/2+1);
  const int x2pncl_ksize_c = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X2PENCIL, SDECOMP_ZDIR, glksize    );
  // tri-diagonal matrix solver in x
  {
    // although right-hand-side is fftw_complex, matrix itself is in double
    poisson_solver->tdm_l = common_calloc(x2pncl_isize_c, sizeof(double));
    poisson_solver->tdm_c = common_calloc(x2pncl_isize_c, sizeof(double));
    poisson_solver->tdm_u = common_calloc(x2pncl_isize_c, sizeof(double));
    const laplace_t * restrict lappx = domain->lappx;
    for(int i = 1; i <= x2pncl_isize_c; i++){
      poisson_solver->tdm_l[i-1] = LAPPX(i).l;
      poisson_solver->tdm_u[i-1] = LAPPX(i).u;
    }
  }
  /* parallel matrix transpose */
  {
    int sizes[NDIMS];
    // x1 to y1 (real)
    sizes[0] = glisize;
    sizes[1] = gljsize;
    sizes[2] = glksize;
    poisson_solver->transposer_x1_to_y1_r = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_X1PENCIL, sizes, sizeof(      double),           MPI_DOUBLE);
    poisson_solver->transposer_y1_to_x1_r = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Y1PENCIL, sizes, sizeof(      double),           MPI_DOUBLE);
    // y1 to z1 (complex)
    sizes[0] = glisize;
    sizes[1] = gljsize/2+1;
    sizes[2] = glksize;
    poisson_solver->transposer_y1_to_z1_c = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_Y1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
    poisson_solver->transposer_z1_to_y1_c = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Z1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
    // z1 to x2 (complex)
    sizes[0] = glisize;
    sizes[1] = gljsize/2+1;
    sizes[2] = glksize;
    poisson_solver->transposer_z1_to_x2_c = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_Z1PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
    poisson_solver->transposer_x2_to_z1_c = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_X2PENCIL, sizes, sizeof(fftw_complex), MPI_C_DOUBLE_COMPLEX);
  }
  /* ! allocate pencils ! 28 ! */
  {
    int sizes[NDIMS];
    // real x1 pencil
    sizes[0] = x1pncl_isize_r;
    sizes[1] = x1pncl_jsize_r;
    sizes[2] = x1pncl_ksize_r;
    poisson_solver->x1pncl_r = fftw_malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(double));
    // real y1 pencil
    sizes[0] = y1pncl_isize_r;
    sizes[1] = y1pncl_jsize_r;
    sizes[2] = y1pncl_ksize_r;
    poisson_solver->y1pncl_r = fftw_malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(double));
    // complex y1 pencil
    sizes[0] = y1pncl_isize_c;
    sizes[1] = y1pncl_jsize_c;
    sizes[2] = y1pncl_ksize_c;
    poisson_solver->y1pncl_c = fftw_malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(fftw_complex));
    // complex z1 pencil
    sizes[0] = z1pncl_isize_c;
    sizes[1] = z1pncl_jsize_c;
    sizes[2] = z1pncl_ksize_c;
    poisson_solver->z1pncl_c = fftw_malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(fftw_complex));
    // complex x2 pencil
    sizes[0] = x2pncl_isize_c;
    sizes[1] = x2pncl_jsize_c;
    sizes[2] = x2pncl_ksize_c;
    poisson_solver->x2pncl_c = fftw_malloc(sizes[0]*sizes[1]*sizes[2]*sizeof(fftw_complex));
  }
  /* ! create y forward fftw plan ! 21 ! */
  {
    // forward transform in y direction, real to complex
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = gljsize, // size of FFT: number of elements in y
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1  // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = y1pncl_ksize_r * y1pncl_isize_r, // how many times do we want to repeat the same FFT?
      .is = y1pncl_jsize_r, // distance between each FFT input
      .os = y1pncl_jsize_c  // distance between each FFT output
    }};
    // create plan
    poisson_solver->fftw_plan_fwrd_y = fftw_plan_guru_dft_r2c(1, dims, 1, hdims, poisson_solver->y1pncl_r, poisson_solver->y1pncl_c, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_fwrd_y == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* ! create y backward fftw plan ! 21 ! */
  {
    // backward transform in y direction, complex to real
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = gljsize, // size of FFT: number of elements in y (NOTE: consider in real space, not gljsize/2+1)
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1  // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = y1pncl_ksize_c * y1pncl_isize_c, // how many times do we want to repeat the same FFT?
      .is = y1pncl_jsize_c, // distance between each FFT input
      .os = y1pncl_jsize_r // distance between each FFT output
    }};
    // create plan
    poisson_solver->fftw_plan_bwrd_y = fftw_plan_guru_dft_c2r(1, dims, 1, hdims, poisson_solver->y1pncl_c, poisson_solver->y1pncl_r, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_bwrd_y == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* ! create z forward fftw plan ! 21 ! */
  {
    // forward transform in z direction, complex to complex
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = glksize, // size of FFT: number of elements in z
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1 // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = z1pncl_isize_c * z1pncl_jsize_c, // how many times do we want to repeat the same FFT?
      .is = z1pncl_ksize_c, // distance between each FFT input
      .os = z1pncl_ksize_c // distance between each FFT output
    }};
    // create plan
    poisson_solver->fftw_plan_fwrd_z = fftw_plan_guru_dft(1, dims, 1, hdims, poisson_solver->z1pncl_c, poisson_solver->z1pncl_c, FFTW_FORWARD, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_fwrd_z == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* ! create z backward fftw plan ! 21 ! */
  {
    // forward transform in z direction, complex to complex
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = glksize, // size of FFT: number of elements in z
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1 // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = z1pncl_isize_c * z1pncl_jsize_c, // how many times do we want to repeat the same FFT?
      .is = z1pncl_ksize_c, // distance between each FFT input
      .os = z1pncl_ksize_c // distance between each FFT output
    }};
    // create plan
    poisson_solver->fftw_plan_bwrd_z = fftw_plan_guru_dft(1, dims, 1, hdims, poisson_solver->z1pncl_c, poisson_solver->z1pncl_c, FFTW_BACKWARD, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_bwrd_z == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* ! compute eigenvalues in y ! 14 ! */
  {
    const double dy = domain->dy;
    const int myjsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X2PENCIL, SDECOMP_YDIR, gljsize/2+1);
    const int joffset = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X2PENCIL, SDECOMP_YDIR, gljsize/2+1);
    poisson_solver->eigenvalues_y = common_calloc(myjsize, sizeof(double));
    for(int local_j = 0; local_j < myjsize; local_j++){
      const int global_j = local_j + joffset;
      double val = -4./pow(dy, 2.)*pow(
          sin( M_PI * global_j / gljsize ),
          2.
      );
      poisson_solver->eigenvalues_y[local_j] = val;
    }
  }
  /* ! compute eigenvalues in z ! 14 ! */
  {
    const double dz = domain->dz;
    const int myksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X2PENCIL, SDECOMP_ZDIR, glksize    );
    const int koffset = sdecomp_get_pencil_offset(sdecomp, SDECOMP_X2PENCIL, SDECOMP_ZDIR, glksize    );
    poisson_solver->eigenvalues_z = common_calloc(myksize, sizeof(double));
    for(int local_k = 0; local_k < myksize; local_k++){
      const int global_k = local_k + koffset;
      double val = -4./pow(dz, 2.)*pow(
          sin( M_PI * global_k / glksize ),
          2.
      );
      poisson_solver->eigenvalues_z[local_k] = val;
    }
  }
  return 0;
}

#define X1PNCL_R(I, J, K) (x1pncl_r[((K)-1)*(jsize)*(isize)+((J)-1)*(isize)+((I)-1)])

/**
 * @brief compute scalar potential \psi to correct velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : velocity (in), scalar potential \psi (out)
 * @return              : error code
 */
int fluid_compute_potential(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // initalise solver, whose buffer "x1pncl_r" is necessary below
  if(poisson_solver == NULL){
    init_poisson_solver(domain);
  }
  /* ! compute right-hand-side ! 30 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict xf = domain->xf;
    const double * restrict xc = domain->xc;
    const double * restrict dxf = domain->dxf;
    const double            dy  = domain->dy;
    const double            dz  = domain->dz;
    const double * restrict ux = fluid->ux;
    const double * restrict uy = fluid->uy;
    const double * restrict uz = fluid->uz;
    double * restrict x1pncl_r = poisson_solver->x1pncl_r;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          const double x_xm  = XF(i  );
          const double x_xp  = XF(i+1);
          const double ux_xm = UX(i  , j  , k  );
          const double ux_xp = UX(i+1, j  , k  );
          const double uy_ym = UY(i  , j  , k  );
          const double uy_yp = UY(i  , j+1, k  );
          const double uz_zm = UZ(i  , j  , k  );
          const double uz_zp = UZ(i  , j  , k+1);
          X1PNCL_R(i, j, k) = 1. / (gamma * dt) * (
              + (x_xp * ux_xp - x_xm * ux_xm) / XC(i  ) / DXF(i  )
              + (       uy_yp -        uy_ym) / XC(i  ) / dy
              + (       uz_zp -        uz_zm) / dz
          );
        }
      }
    }
  }
  /* ! transpose real x1pencil to y1pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_x1_to_y1_r,
      poisson_solver->x1pncl_r,
      poisson_solver->y1pncl_r
  );
  /* ! project to wave space (y) ! 1 ! */
  fftw_execute(poisson_solver->fftw_plan_fwrd_y);
  /* ! transpose complex y1pencil to z1pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_y1_to_z1_c,
      poisson_solver->y1pncl_c,
      poisson_solver->z1pncl_c
  );
  /* ! project to wave space (z) ! 1 ! */
  fftw_execute(poisson_solver->fftw_plan_fwrd_z);
  /* ! transpose complex z1pencil to x2pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_z1_to_x2_c,
      poisson_solver->z1pncl_c,
      poisson_solver->x2pncl_c
  );
  /* solve linear systems */
  {
    const int glisize = domain->glsizes[0];
    const int gljsize = domain->glsizes[1];
    const int glksize = domain->glsizes[2];
    const int isize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X2PENCIL, SDECOMP_XDIR, glisize    );
    const int jsize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X2PENCIL, SDECOMP_YDIR, gljsize/2+1);
    const int ksize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_X2PENCIL, SDECOMP_ZDIR, glksize    );
    const double * restrict xc = domain->xc;
    const laplace_t * restrict lappx = domain->lappx;
    const double * restrict tdm_l = poisson_solver->tdm_l;
    const double * restrict tdm_u = poisson_solver->tdm_u;
    double * restrict tdm_c = poisson_solver->tdm_c;
    for(int k = 0; k < ksize; k++){
      const double eval_z = poisson_solver->eigenvalues_z[k];
      for(int j = 0; j < jsize; j++){
        const double eval_y = poisson_solver->eigenvalues_y[j];
        /* set center diagonal components */
        for(int i = 1; i <= isize; i++){
          tdm_c[i-1] = LAPPX(i).c
            + eval_y / XC(i  ) / XC(i  )
            + eval_z;
        }
        fftw_complex * restrict rhs = poisson_solver->x2pncl_c + k * jsize * isize + j * isize;
        /* ! boundary treatment (Neumann boundary condition) ! 2 ! */
        tdm_c[      0] += tdm_l[      0];
        tdm_c[isize-1] += tdm_u[isize-1];
        /* ! solve linear system ! 9 ! */
        tdm_solve_fftw_complex(
            /* size of system  */ isize,
            /* number of rhs   */ 1,
            /* periodic        */ false,
            /* lower- diagonal */ tdm_l,
            /* center-diagonal */ tdm_c,
            /* upper- diagonal */ tdm_u,
            /* input / output  */ rhs
        );
        /* ! normalise FFT ! 3 ! */
        for(int i = 0; i < isize; i++){
          rhs[i] /= 1. * gljsize * glksize;
        }
      }
    }
  }
  /* ! transpose complex x2pencil to z1pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_x2_to_z1_c,
      poisson_solver->x2pncl_c,
      poisson_solver->z1pncl_c
  );
  /* ! project to physical space (z) ! 1 ! */
  fftw_execute(poisson_solver->fftw_plan_bwrd_z);
  /* ! transpose complex z1pencil to y1pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_z1_to_y1_c,
      poisson_solver->z1pncl_c,
      poisson_solver->y1pncl_c
  );
  /* ! project to physical space (y) ! 1 ! */
  fftw_execute(poisson_solver->fftw_plan_bwrd_y);
  /* ! transpose real y1pencil to x1pencil ! 5 ! */
  sdecomp_transpose_execute(
      poisson_solver->transposer_y1_to_x1_r,
      poisson_solver->y1pncl_r,
      poisson_solver->x1pncl_r
  );
  /* ! store result ! 16 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double * restrict x1pncl_r = poisson_solver->x1pncl_r;
    double * restrict psi = fluid->psi;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          PSI(i, j, k) = X1PNCL_R(i, j, k);
        }
      }
    }
    // NOTE: psi and p are assumed to have the same array shape
    fluid_update_boundaries_p(domain, psi);
  }
  return 0;
}

#undef X1PNCL_R

/**
 * @brief destruct a structure poisson_solver_t
 * @return : error code
 */
int fluid_compute_potential_finalise(void){
  // use fftw_free since they are allocated using fftw_malloc
  fftw_free(poisson_solver->x1pncl_r);
  fftw_free(poisson_solver->y1pncl_r);
  fftw_free(poisson_solver->y1pncl_c);
  fftw_free(poisson_solver->z1pncl_c);
  fftw_free(poisson_solver->x2pncl_c);
  fftw_destroy_plan(poisson_solver->fftw_plan_fwrd_y);
  fftw_destroy_plan(poisson_solver->fftw_plan_bwrd_y);
  fftw_destroy_plan(poisson_solver->fftw_plan_fwrd_z);
  fftw_destroy_plan(poisson_solver->fftw_plan_bwrd_z);
  fftw_cleanup();
  common_free(poisson_solver->eigenvalues_y);
  common_free(poisson_solver->eigenvalues_z);
  sdecomp_transpose_finalise(poisson_solver->transposer_x1_to_y1_r);
  sdecomp_transpose_finalise(poisson_solver->transposer_y1_to_x1_r);
  sdecomp_transpose_finalise(poisson_solver->transposer_y1_to_z1_c);
  sdecomp_transpose_finalise(poisson_solver->transposer_z1_to_y1_c);
  sdecomp_transpose_finalise(poisson_solver->transposer_z1_to_x2_c);
  sdecomp_transpose_finalise(poisson_solver->transposer_x2_to_z1_c);
  common_free(poisson_solver);
  return 0;
}


#include <fftw3.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "sdecomp.h"
#include "fluid.h"
#include "tdm.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "../arrays/psi.h"
#include "internal.h"

// definition of solver
poisson_solver_t *poisson_solver = NULL;

static int assign_input(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double gamma = RKCOEFS[rkstep].gamma;
  const double * restrict xf  = domain->xf;
  const double * restrict xc  = domain->xc;
  const double * restrict dxf = domain->dxf;
  const double            dy  = domain->dy;
  const double            dz  = domain->dz;
  const double * restrict ux = fluid->ux->data;
  const double * restrict uy = fluid->uy->data;
  const double * restrict uz = fluid->uz->data;
  double * restrict rhs = poisson_solver->buf0;
  // normalise FFT beforehand
  const double norm = 1. * domain->glsizes[1] * domain->glsizes[2];
  const double prefactor = 1. / (gamma * dt) / norm;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dx    = DXF(i  );
        const double x     = XC(i  );
        const double x_xm  = XF(i  );
        const double x_xp  = XF(i+1);
        const double ux_xm = UX(i  , j  , k  );
        const double ux_xp = UX(i+1, j  , k  );
        const double uy_ym = UY(i  , j  , k  );
        const double uy_yp = UY(i  , j+1, k  );
        const double uz_zm = UZ(i  , j  , k  );
        const double uz_zp = UZ(i  , j  , k+1);
        rhs[cnt++] = prefactor * (
           + (x_xp * ux_xp - x_xm * ux_xm) / x / dx
           + (       uy_yp -        uy_ym) / x / dy
           + (       uz_zp -        uz_zm)     / dz
        );
      }
    }
  }
  return 0;
}

static int extract_output(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict rhs = poisson_solver->buf0;
  double * restrict psi = fluid->psi->data;
  for(int cnt = 0, k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        PSI(i, j, k) = *(rhs + (cnt++));
      }
    }
  }
  // NOTE: psi and p are assumed to have the same array shape (halos)
  fluid_update_boundaries_p(domain, psi);
  return 0;
}

static int solve_linear_systems(const domain_t * restrict domain){
  const double * restrict xc = domain->xc;
  // tri-diagonal matrix
  tdm_info_t * restrict tdm_info = poisson_solver->tdm_info;
  const size_t * restrict tdm_sizes = poisson_solver->tdm_sizes;
  double * restrict tdm_l = NULL;
  double * restrict tdm_u = NULL;
  double * restrict tdm_c = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_u(tdm_info, &tdm_u);
  tdm.get_c(tdm_info, &tdm_c);
  // eigenvalues coming from Fourier projection
  const double * restrict evalys = poisson_solver->evalys;
  const double * restrict evalzs = poisson_solver->evalzs;
  fftw_complex * restrict rhs = poisson_solver->buf1;
  for(size_t cnt = 0, k = 0; k < tdm_sizes[2]; k++){
    const double evalz = evalzs[k];
    for(size_t j = 0; j < tdm_sizes[1]; j++){
      const double evaly = evalys[j];
      // set center diagonal components
      for(size_t i = 1; i <= tdm_sizes[0]; i++){
        tdm_c[i-1] =
          - tdm_l[i-1]
          - tdm_u[i-1]
          + evaly / XC(i  ) / XC(i  )
          + evalz;
      }
      // boundary treatment (Neumann boundary condition)
      tdm_c[             0] += tdm_l[             0];
      tdm_c[tdm_sizes[0]-1] += tdm_u[tdm_sizes[0]-1];
      tdm.solve(tdm_info, rhs + (cnt++) * tdm_sizes[0]);
    }
  }
  return 0;
}

/**
 * @brief compute scalar potential psi to correct velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : velocity (in), scalar potential psi (out)
 * @return               : (success) 0
 *                       : (failure) 1
 */
int fluid_compute_potential(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // initialise Poisson solver 
  if(NULL == poisson_solver){
    const int retval = fluid_compute_potential_init(domain);
    if(0 != retval) return 1;
    MPI_Comm comm_cart = MPI_COMM_NULL;
    sdecomp.get_comm_cart(domain->info, &comm_cart);
    int myrank = 0;
    MPI_Comm_rank(comm_cart, &myrank);
    if(0 == myrank){
      printf("Poisson solver initialised\n");
    }
  }
  // compute right-hand side of Poisson equation 
  // assigned to buf0
  assign_input(domain, rkstep, dt, fluid);
  // solve the equation
  /* normal solver with DFT */
  // dft solver, transpose real x1pencil to y1pencil 
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver->r_transposer_x1_to_y1,
      poisson_solver->buf0,
      poisson_solver->buf1
  );
  // dft solver, project y to wave space 
  // f(x, y)    -> f(x, k_y)
  // f(x, y, z) -> f(x, k_y, z)
  // from buf1 to buf0
  fftw_execute(poisson_solver->fftw_plan_y[0]);
  // dft solver, transpose complex y1pencil to z1pencil 
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver->c_transposer_y1_to_z1,
      poisson_solver->buf0,
      poisson_solver->buf1
  );
  // dft solver, project z to wave space 
  // f(x, k_y, z) -> f(x, k_y, k_z)
  // from buf1 to buf0
  fftw_execute(poisson_solver->fftw_plan_z[0]);
  // dft solver, transpose complex z1pencil to x2pencil 
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver->c_transposer_z1_to_x2,
      poisson_solver->buf0,
      poisson_solver->buf1
  );
  // dft solver, solve linear systems 
  solve_linear_systems(domain);
  // dft solver, transpose complex x2pencil to z1pencil 
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver->c_transposer_x2_to_z1,
      poisson_solver->buf1,
      poisson_solver->buf0
  );
  // dft solver, project z to physical space 
  // f(x, k_y, k_z) -> f(x, k_y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver->fftw_plan_z[1]);
  // dft solver, transpose complex z1pencil to y1pencil 
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver->c_transposer_z1_to_y1,
      poisson_solver->buf1,
      poisson_solver->buf0
  );
  // dft solver, project y to physical space 
  // f(x, k_y)    -> f(x, y)
  // f(x, k_y, z) -> f(x, y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver->fftw_plan_y[1]);
  // dft solver, transpose real y1pencil to x1pencil 
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver->r_transposer_y1_to_x1,
      poisson_solver->buf1,
      poisson_solver->buf0
  );
  extract_output(domain, fluid);
  return 0;
}


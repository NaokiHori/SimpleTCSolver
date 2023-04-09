#include "config.h"
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/domain/plapx.h"
#include "arrays/domain/plapy.h"
#include "arrays/fluid/p.h"
#include "arrays/psi.h"

static inline int add_explicit(const domain_t * restrict domain, const fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict psi = fluid->psi->data;
  double * restrict p = fluid->p->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // explicit contribution 
        P(i, j, k) += PSI(i, j, k);
      }
    }
  }
  return 0;
}

static inline int add_implicit_x(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict plapx = domain->plapx;
  const double * restrict psi = fluid->psi->data;
  double * restrict p = fluid->p->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // x implicit contribution 
        P(i, j, k) -= prefactor * (
            + PLAPX(i  ).l * PSI(i-1, j  , k  )
            + PLAPX(i  ).c * PSI(i  , j  , k  )
            + PLAPX(i  ).u * PSI(i+1, j  , k  )
        );
      }
    }
  }
  return 0;
}

static inline int add_implicit_y(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict plapy = domain->plapy;
  const double * restrict psi = fluid->psi->data;
  double * restrict p = fluid->p->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // y implicit contribution 
        P(i, j, k) -= prefactor * (
            + PLAPY(i  ).l * PSI(i  , j-1, k  )
            + PLAPY(i  ).c * PSI(i  , j  , k  )
            + PLAPY(i  ).u * PSI(i  , j+1, k  )
        );
      }
    }
  }
  return 0;
}

static inline int add_implicit_z(const domain_t * restrict domain, const double prefactor, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t plapz = domain->lapz;
  const double * restrict psi = fluid->psi->data;
  double * restrict p = fluid->p->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // z implicit contribution 
        P(i, j, k) -= prefactor * (
            + plapz.l * PSI(i  , j  , k-1)
            + plapz.c * PSI(i  , j  , k  )
            + plapz.u * PSI(i  , j  , k+1)
        );
      }
    }
  }
  return 0;
}

/**
 * @brief update pressure using scalar potential psi
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : scalar potential (in), pressure (out)
 * @return               : error code
 */
int fluid_update_pressure(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // explicit contribution, always present
  add_explicit(domain, fluid);
  // gamma dt diffusivity / 2 
  const double prefactor =
    0.5 * RKCOEFS[rkstep].gamma * dt * fluid->diffusivity;
  // additional corrections if diffusive terms
  //   in the direction is treated implicitly
  if(config.get.implicitx()){
    add_implicit_x(domain, prefactor, fluid);
  }
  if(config.get.implicity()){
    add_implicit_y(domain, prefactor, fluid);
  }
  if(config.get.implicitz()){
    add_implicit_z(domain, prefactor, fluid);
  }
  // impose boundary conditions and communicate halo cells 
  fluid_update_boundaries_p(domain, fluid->p->data);
  return 0;
}


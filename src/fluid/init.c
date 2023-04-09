#include <math.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"
#include "arrays/psi.h"
#include "arrays/srcuxa.h"
#include "arrays/srcuxb.h"
#include "arrays/srcuxg.h"
#include "arrays/srcuya.h"
#include "arrays/srcuyb.h"
#include "arrays/srcuyg.h"
#include "arrays/srcuza.h"
#include "arrays/srcuzb.h"
#include "arrays/srcuzg.h"

#define CONCAT(a, b) a##b

#define ARRAY_PREPARE(dsetname, DSETNAME) \
  array_prepare( \
      domain, \
      (int [NDIMS][2]){ \
        {CONCAT(DSETNAME, _LNADD_0), CONCAT(DSETNAME, _UNADD_0)}, \
        {CONCAT(DSETNAME, _LNADD_1), CONCAT(DSETNAME, _UNADD_1)}, \
        {CONCAT(DSETNAME, _LNADD_2), CONCAT(DSETNAME, _UNADD_2)}, \
      }, \
      &fluid->dsetname \
  )

/**
 * @brief allocate fluid_t
 * @param[in] domain : information about domain decomposition and size
 * @return           : structure being allocated
 */
static fluid_t *allocate(const domain_t * restrict domain){
  // main structure
  fluid_t * restrict fluid = common_calloc(1, sizeof(fluid_t));
  // velocity
  ARRAY_PREPARE(    ux,     UX);
  ARRAY_PREPARE(    uy,     UY);
  ARRAY_PREPARE(    uz,     UZ);
  // pressure
  ARRAY_PREPARE(     p,      P);
  // scalar potential
  ARRAY_PREPARE(   psi,    PSI);
  // Runge-Kutta source terms
  ARRAY_PREPARE(srcuxa, SRCUXA);
  ARRAY_PREPARE(srcuxb, SRCUXB);
  ARRAY_PREPARE(srcuxg, SRCUXG);
  ARRAY_PREPARE(srcuya, SRCUYA);
  ARRAY_PREPARE(srcuyb, SRCUYB);
  ARRAY_PREPARE(srcuyg, SRCUYG);
  ARRAY_PREPARE(srcuza, SRCUZA);
  ARRAY_PREPARE(srcuzb, SRCUZB);
  ARRAY_PREPARE(srcuzg, SRCUZG);
  return fluid;
}

/**
 * @brief constructor of the structure
 * @param[in] dirname_ic : name of directory in which initial flow fields are stored
 * @param[in] domain     : information about domain decomposition and size
 * @return               : structure being allocated and initalised
 */
fluid_t *fluid_init(const char dirname_ic[], const domain_t * restrict domain){
  // allocate structure and its members 
  fluid_t * restrict fluid = allocate(domain);
  // load flow fields 
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  if(0 != array_load(comm_cart, dirname_ic, "ux", fluid->ux)) return NULL;
  if(0 != array_load(comm_cart, dirname_ic, "uy", fluid->uy)) return NULL;
  if(0 != array_load(comm_cart, dirname_ic, "uz", fluid->uz)) return NULL;
  if(0 != array_load(comm_cart, dirname_ic,  "p", fluid-> p)) return NULL;
  // impose boundary conditions and communicate halo cells 
  fluid_update_boundaries_ux(domain, fluid->ux->data);
  fluid_update_boundaries_uy(domain, fluid->uy->data);
  fluid_update_boundaries_uz(domain, fluid->uz->data);
  fluid_update_boundaries_p (domain, fluid-> p->data);
  // compute diffusivity 
  const double Re = config.get.Re();
  fluid->diffusivity = 1. / Re;
  return fluid;
}


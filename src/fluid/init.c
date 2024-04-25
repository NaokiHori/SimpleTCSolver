#include <math.h>
#include "param.h"
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "fileio.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/psi.h"
#include "array_macros/fluid/t.h"
#include "array_macros/fluid/lxx.h"
#include "array_macros/fluid/lyx0.h"
#include "array_macros/fluid/lyx1.h"
#if NDIMS == 3
#include "array_macros/fluid/lzx.h"
#endif
#include "array_macros/fluid/lxy.h"
#include "array_macros/fluid/lyy0.h"
#include "array_macros/fluid/lyy1.h"
#if NDIMS == 3
#include "array_macros/fluid/lzy.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lxz.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lyz.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lzz.h"
#endif
#include "array_macros/fluid/srcux.h"
#include "array_macros/fluid/srcuy.h"
#if NDIMS == 3
#include "array_macros/fluid/srcuz.h"
#endif
#include "array_macros/fluid/srct.h"

/**
 * @brief allocate members
 * @param[in]  domain : information about domain decomposition and size
 * @param[out] fluid  : structure storing flow fields and auxiliary buffers
 * @return            : error code
 */
static int allocate(
    const domain_t * domain,
    fluid_t * fluid
){
  // velocity
  if(0 != array.prepare(domain, UX_NADDS, sizeof(double), &fluid->ux )) return 1;
  if(0 != array.prepare(domain, UY_NADDS, sizeof(double), &fluid->uy )) return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, UZ_NADDS, sizeof(double), &fluid->uz )) return 1;
#endif
  // pressure and scalar potential
  if(0 != array.prepare(domain, P_NADDS,   sizeof(double), &fluid->p  )) return 1;
  if(0 != array.prepare(domain, PSI_NADDS, sizeof(double), &fluid->psi)) return 1;
  // scalar field
  if(0 != array.prepare(domain, T_NADDS, sizeof(double), &fluid->t)) return 1;
  // velocity-gradient tensor components
  if(0 != array.prepare(domain, LXX_NADDS,  sizeof(double), &fluid->lxx )) return 1;
  if(0 != array.prepare(domain, LYX0_NADDS, sizeof(double), &fluid->lyx0)) return 1;
  if(0 != array.prepare(domain, LYX1_NADDS, sizeof(double), &fluid->lyx1)) return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, LZX_NADDS,  sizeof(double), &fluid->lzx )) return 1;
#endif
  if(0 != array.prepare(domain, LXY_NADDS,  sizeof(double), &fluid->lxy )) return 1;
  if(0 != array.prepare(domain, LYY0_NADDS, sizeof(double), &fluid->lyy0)) return 1;
  if(0 != array.prepare(domain, LYY1_NADDS, sizeof(double), &fluid->lyy1)) return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, LZY_NADDS,  sizeof(double), &fluid->lzy )) return 1;
#endif
#if NDIMS == 3
  if(0 != array.prepare(domain, LXZ_NADDS,  sizeof(double), &fluid->lxz )) return 1;
#endif
#if NDIMS == 3
  if(0 != array.prepare(domain, LYZ_NADDS,  sizeof(double), &fluid->lyz )) return 1;
#endif
#if NDIMS == 3
  if(0 != array.prepare(domain, LZZ_NADDS,  sizeof(double), &fluid->lzz )) return 1;
#endif
  // Runge-Kutta source terms
  for(size_t n = 0; n < 3; n++){
    if(0 != array.prepare(domain, SRCUX_NADDS, sizeof(double), &fluid->srcux[n])) return 1;
    if(0 != array.prepare(domain, SRCUY_NADDS, sizeof(double), &fluid->srcuy[n])) return 1;
#if NDIMS == 3
    if(0 != array.prepare(domain, SRCUZ_NADDS, sizeof(double), &fluid->srcuz[n])) return 1;
#endif
    if(0 != array.prepare(domain, SRCT_NADDS,  sizeof(double), &fluid->srct[n]))  return 1;
  }
  return 0;
}

static void report(
    const sdecomp_info_t * info,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(info, &myrank);
  if(root == myrank){
    printf("FLUID\n");
    printf("\tRe: % .7e\n", fluid->Re);
    printf("\tPr: % .7e\n", fluid->Pr);
    printf("\tdiffusive treatment in x: %s\n", param_implicit_x ? "implicit" : "explicit");
    printf("\tdiffusive treatment in y: %s\n", param_implicit_y ? "implicit" : "explicit");
#if NDIMS == 3
    printf("\tdiffusive treatment in z: %s\n", param_implicit_z ? "implicit" : "explicit");
#endif
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial flow fields are stored
 * @param[in]  domain     : information about domain decomposition and size
 * @param[out]            : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int fluid_init(
    const char dirname_ic[],
    const domain_t * domain,
    fluid_t * fluid
){
  // allocate arrays
  if(0 != allocate(domain, fluid)) return 1;
  // load flow fields
  if(0 != array.load(domain, dirname_ic, "ux", fileio.npy_double, &fluid->ux)) return 1;
  if(0 != array.load(domain, dirname_ic, "uy", fileio.npy_double, &fluid->uy)) return 1;
#if NDIMS == 3
  if(0 != array.load(domain, dirname_ic, "uz", fileio.npy_double, &fluid->uz)) return 1;
#endif
  if(0 != array.load(domain, dirname_ic,  "p", fileio.npy_double, &fluid-> p)) return 1;
  if(0 != array.load(domain, dirname_ic,  "t", fileio.npy_double, &fluid-> t)) return 1;
  // impose boundary conditions and communicate halo cells
  if(0 != fluid_update_boundaries_ux(domain, &fluid->ux)) return 1;
  if(0 != fluid_update_boundaries_uy(domain, &fluid->uy)) return 1;
#if NDIMS == 3
  if(0 != fluid_update_boundaries_uz(domain, &fluid->uz)) return 1;
#endif
  if(0 != fluid_update_boundaries_p(domain, &fluid->p)) return 1;
  if(0 != fluid_update_boundaries_t(domain, &fluid->t)) return 1;
  // get non-dimensional numbers
  if(0 != config.get_double("Re", &fluid->Re)) return 1;
  if(0 != config.get_double("Pr", &fluid->Pr)) return 1;
  report(domain->info, fluid);
  return 0;
}


#include "common.h"
#include "array.h"
#include "fluid.h"
#include "internal.h"

/**
 * @brief destruct a structure fluid_t
 * @param[in,out] fluid : structure to be cleaned-up
 * @return              : error code
 */
int fluid_finalise(fluid_t * restrict fluid){
  fluid_update_velocity_finalise();
  fluid_compute_potential_finalise();
  array_destroy(fluid->ux);
  array_destroy(fluid->uy);
  array_destroy(fluid->uz);
  array_destroy(fluid->p);
  array_destroy(fluid->psi);
  array_destroy(fluid->srcuxa);
  array_destroy(fluid->srcuxb);
  array_destroy(fluid->srcuxg);
  array_destroy(fluid->srcuya);
  array_destroy(fluid->srcuyb);
  array_destroy(fluid->srcuyg);
  array_destroy(fluid->srcuza);
  array_destroy(fluid->srcuzb);
  array_destroy(fluid->srcuzg);
  common_free(fluid);
  return 0;
}


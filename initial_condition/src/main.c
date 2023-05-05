#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"

int main(void){
  // load env
  config.load();
  // save step.npy, time.npy
  const    int z_int = 0;
  const double z_dbl = 0.;
  common_save("step.npy", 0, NULL, NPY_INT, sizeof(   int), &z_int);
  common_save("time.npy", 0, NULL, NPY_DBL, sizeof(double), &z_dbl);
  // init domain and save
  domain_t *domain = domain_init();
  // init fluid and save
  fluid_init(domain);
  // free domain
  domain_finalise(domain);
  // release env
  config.unload();
  return 0;
}


#include "sdecomp.h"
#include "array.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/fluid/p.h"

int fluid_save(const char dirname[], const domain_t * restrict domain, const fluid_t * restrict fluid){
  // serial
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    fileio_w_serial(dirname, "f_diffusivity", 0, NULL, NPY_DBL, sizeof(double), &fluid->diffusivity);
  }
  // collective
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  // ux
  array_dump(comm_cart, dirname, "ux", fluid->ux);
  // uy
  array_dump(comm_cart, dirname, "uy", fluid->uy);
  // uz
  array_dump(comm_cart, dirname, "uz", fluid->uz);
  // p
  array_dump(comm_cart, dirname,  "p",  fluid->p);
  return 0;
}


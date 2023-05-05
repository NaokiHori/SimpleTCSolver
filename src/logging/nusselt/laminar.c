#include "internal.h"

double logging_nusselt_compute_torque_laminar(const double diffusivity, const double ri, const double ro, const double ui, const double uo, const double ly, const double lz){
  const double num = ri * ro * (ui * ro - uo * ri);
  const double den = ro * ro - ri * ri;
  return 2. * diffusivity * num / den * ly * lz;
}


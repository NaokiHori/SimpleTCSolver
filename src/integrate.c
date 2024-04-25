#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "integrate.h"
#include "decide_dt.h"

// integrate the equations for one time step
int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    double * dt
){
  // decide time step size
  if(0 != decide_dt(domain, fluid, dt)){
    return 1;
  }
  // Runge-Kutta iterations
  // max iteration, should be three
  for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    if(0 != fluid_predict_field(domain, rkstep, *dt, fluid)){
      return 1;
    }
    if(0 != fluid_compute_potential(domain, rkstep, *dt, fluid)){
      return 1;
    }
    if(0 != fluid_correct_velocity(domain, rkstep, *dt, fluid)){
      return 1;
    }
    if(0 != fluid_update_pressure(domain, rkstep, *dt, fluid)){
      return 1;
    }
  }
  return 0;
}

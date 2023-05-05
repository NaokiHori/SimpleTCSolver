#include "param.h"

// coefficients of three-stage Runge-Kutta scheme
const rkcoef_t param_rkcoefs[3] = {
  {.alpha = 32. / 60., .beta  =   0. / 60., .gamma = 32. / 60.},
  {.alpha = 25. / 60., .beta  = -17. / 60., .gamma =  8. / 60.},
  {.alpha = 45. / 60., .beta  = -25. / 60., .gamma = 20. / 60.}
};


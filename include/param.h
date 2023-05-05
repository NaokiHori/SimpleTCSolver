#pragma once

#include <stdbool.h>

/* runge-kutta.c */
// Runge-Kutta configurations
typedef struct {
  double alpha;
  double beta;
  double gamma;
} rkcoef_t;
extern const rkcoef_t param_rkcoefs[3];

// implicit.c
extern const bool param_implicit_x;
extern const bool param_implicit_y;
extern const bool param_implicit_z;

// boundary-condition.c
extern const double param_uy_xm;
extern const double param_uy_xp;
extern const double param_uz_xm;
extern const double param_uz_xp;


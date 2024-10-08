#if !defined(PARAM_H)
#define PARAM_H

// fixed parameters, which are usually fixed
//   but still user can easily control, are declared
// they are defined under src/param/xxx.c

#include <stdbool.h>

/* implicit.c */
// flags to specify the diffusive treatment of the momentum equations
extern const bool param_implicit_x;
extern const bool param_implicit_y;
extern const bool param_implicit_z;

/* boundary-condition.c */
// NOTE: impermeable walls and Neumann BC for the pressure are unchangeable
// negative-x-wall velocity in y direction
extern const double param_uy_xm;
// positive-x-wall velocity in y direction
extern const double param_uy_xp;
// scalar value on the negative wall
extern const double param_t_xm;
// scalar value on the positive wall
extern const double param_t_xp;

#endif // PARAM_H

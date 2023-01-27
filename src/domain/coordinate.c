#include <math.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "internal.h"


/**
 * @brief define cell-face positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] lx    : physical domain size in x direction
 *                      normally 1 since this is the reference length scale
 *                      to non-dimensionalise equations
 * @return : cell-face positions in x direction
 */
double *allocate_and_init_xf(const int isize, const double lx){
  /* ! xf: cell face coordinates ! 28 ! */
  double *xf = common_calloc(XF_SIZE, sizeof(double));
  const bool use_stretched_grid = config.get_bool("use_stretched_grid");
  if(use_stretched_grid){
    // stretched grid, clipped Chebyshev just as an example
    // number of grid points to be clipped at the edges
    const int nclip = 3;
    for(int i = 1; i <= isize+1; i++){
      XF(i  ) = cos( M_PI * (1. * (i - 1) + nclip) / (isize + 2. * nclip) );
    }
    // map XF(1) and XF(isize+1) to r_i and r_o
    const double r_min = 1.;
    const double r_max = 1. + lx;
    const double xf_min = XF(        1);
    const double xf_max = XF(isize + 1);
    for(int i = 1; i <= isize+1; i++){
      XF(i  ) = r_min
              + (r_max - r_min) / (xf_max - xf_min)
              * (XF(i) - xf_min);
    }
  }else{
    // uniform grid, which enables us to use
    //   efficient DCT-based Poisson solver
    //   see src/fluid/compute_potential.c
    const double dx = lx / isize;
    for(int i = 1; i <= isize+1; i++){
      XF(i  ) = 1. + 1. * (i - 1) * dx;
    }
  }
  return xf;
}

/**
 * @brief define cell-center positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : cell-center positions in x direction
 */
double *allocate_and_init_xc(const int isize, const double *xf){
  /* ! xc: cell center coordinates ! 8 ! */
  double *xc = common_calloc(XC_SIZE, sizeof(double));
  // center between two XFs
  for(int i = 1; i <= isize; i++){
    XC(i) = 0.5 * (XF(i  ) + XF(i+1));
  }
  // at boundaries, face positions are assigned
  XC(      0) = XF(      1);
  XC(isize+1) = XF(isize+1);
  return xc;
}

/**
 * @brief define face-to-face distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : face-to-face distances in x direction
 */
double *allocate_and_init_dxf(const int isize, const double *xf){
  /* ! dxf: distance from cell face to cell face ! 4 ! */
  double *dxf = common_calloc(DXF_SIZE, sizeof(double));
  for(int i = 1; i <= isize; i++){
    DXF(i) = XF(i+1) - XF(i  );
  }
  return dxf;
}

/**
 * @brief define center-to-center distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : center-to-center distances in x direction
 */
double *allocate_and_init_dxc(const int isize, const double *xc){
  /* ! dxc: distance from cell center to cell center (generally) ! 4 ! */
  double *dxc = common_calloc(DXC_SIZE, sizeof(double));
  for(int i = 1; i <= isize+1; i++){
    DXC(i) = XC(i  ) - XC(i-1);
  }
  return dxc;
}

// discrete Laplace operator in x direction with respect to ux
laplace_t *allocate_and_init_lapuxx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *lapuxx = common_calloc(LAPUXX_SIZE, sizeof(laplace_t));
  /* ! lapuxx ! 6 ! */
  for(int i = 2; i <= isize; i++){
    LAPUXX(i).l = + XF(i-1) / DXF(i-1) / XC(i-1) / DXC(i  );
    LAPUXX(i).u = + XF(i+1) / DXF(i  ) / XC(i  ) / DXC(i  );
    LAPUXX(i).c = - XF(i  ) / DXF(i-1) / XC(i-1) / DXC(i  )
                  - XF(i  ) / DXF(i  ) / XC(i  ) / DXC(i  );
  }
  return lapuxx;
}

// discrete Laplace operator in x direction with respect to uy
laplace_t *allocate_and_init_lapuyx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *lapuyx = common_calloc(LAPUYX_SIZE, sizeof(laplace_t));
  /* ! lapuyx ! 6 ! */
  for(int i = 1; i <= isize; i++){
    LAPUYX(i).l = + 1. / XC(i-1) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
    LAPUYX(i).u = + 1. / XC(i+1) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
    LAPUYX(i).c = - 1. / XC(i  ) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.)
                  - 1. / XC(i  ) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
  }
  return lapuyx;
}

// discrete Laplace operator in x direction with respect to uz
laplace_t *allocate_and_init_lapuzx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *lapuzx = common_calloc(LAPUZX_SIZE, sizeof(laplace_t));
  /* ! lapuzx ! 6 ! */
  for(int i = 1; i <= isize; i++){
    LAPUZX(i).l = + 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  );
    LAPUZX(i).u = + 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
    LAPUZX(i).c = - 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  )
                  - 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
  }
  return lapuzx;
}

// discrete Laplace operator in x direction with respect to p and psi
laplace_t *allocate_and_init_lappx (const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *lappx = common_calloc(LAPPX_SIZE, sizeof(laplace_t));
  for(int i = 1; i <= isize; i++){
    LAPPX(i).l = + 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  );
    LAPPX(i).u = + 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
    LAPPX(i).c = - 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  )
                 - 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
  }
  return lappx;
}

// discrete Laplace operator in y direction with respect to ux
laplace_t *allocate_and_init_lapuxy(const int isize, const double *xf, const double dy){
  laplace_t *lapuxy = common_calloc(LAPUXY_SIZE, sizeof(laplace_t));
  /* ! lapuxy ! 6 ! */
  for(int i = 2; i <= isize; i++){
    // for each x
    LAPUXY(i).l = + 1. / XF(i  ) / XF(i  ) / dy / dy;
    LAPUXY(i).u = + 1. / XF(i  ) / XF(i  ) / dy / dy;
    LAPUXY(i).c = - 2. / XF(i  ) / XF(i  ) / dy / dy;
  }
  return lapuxy;
}

// discrete Laplace operator in y direction with respect to uy
laplace_t *allocate_and_init_lapuyy(const int isize, const double *xc, const double dy){
  laplace_t *lapuyy = common_calloc(LAPUYY_SIZE, sizeof(laplace_t));
  /* ! lapuyy ! 6 ! */
  for(int i = 1; i <= isize; i++){
    // for each x
    LAPUYY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPUYY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPUYY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
  }
  return lapuyy;
}

// discrete Laplace operator in y direction with respect to uz
laplace_t *allocate_and_init_lapuzy(const int isize, const double *xc, const double dy){
  laplace_t *lapuzy = common_calloc(LAPUZY_SIZE, sizeof(laplace_t));
  /* ! lapuzy ! 6 ! */
  for(int i = 1; i <= isize; i++){
    // for each x
    LAPUZY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPUZY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPUZY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
  }
  return lapuzy;
}

// discrete Laplace operator in y direction with respect to p
laplace_t *allocate_and_init_lappy (const int isize, const double *xc, const double dy){
  laplace_t *lappy = common_calloc(LAPPY_SIZE, sizeof(laplace_t));
  for(int i = 1; i <= isize; i++){
    // for each x
    LAPPY(i).l = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPPY(i).u = + 1. / XC(i  ) / XC(i  ) / dy / dy;
    LAPPY(i).c = - 2. / XC(i  ) / XC(i  ) / dy / dy;
  }
  return lappy;
}

// discrete Laplace operator in z direction with respect to ux
laplace_t init_lapuxz(const double dz){
  /* ! lapuxz ! 6 ! */
  const laplace_t lapuxz = {
    // constant for all x and y directions
    .l = + 1. / dz / dz,
    .u = + 1. / dz / dz,
    .c = - 2. / dz / dz
  };
  return lapuxz;
}

// discrete Laplace operator in z direction with respect to uy
laplace_t init_lapuyz(const double dz){
  /* ! lapuyz ! 6 ! */
  const laplace_t lapuyz = {
    // constant for all x and y directions
    .l = + 1. / dz / dz,
    .u = + 1. / dz / dz,
    .c = - 2. / dz / dz
  };
  return lapuyz;
}

// discrete Laplace operator in z direction with respect to uz
laplace_t init_lapuzz(const double dz){
  /* ! lapuzz ! 6 ! */
  const laplace_t lapuzz = {
    // constant for all x and y directions
    .l = + 1. / dz / dz,
    .u = + 1. / dz / dz,
    .c = - 2. / dz / dz
  };
  return lapuzz;
}

// discrete Laplace operator in z direction with respect to p
laplace_t init_lappz(const double dz){
  const laplace_t lappz = {
    // constant for all x and y directions
    .l = + 1. / dz / dz,
    .u = + 1. / dz / dz,
    .c = - 2. / dz / dz
  };
  return lappz;
}


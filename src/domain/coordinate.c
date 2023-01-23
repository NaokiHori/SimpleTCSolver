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

laplace_t *allocate_and_init_uxdifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *uxdifx = common_calloc(UXDIFX_SIZE, sizeof(laplace_t));
  /* ! uxdifx ! 6 ! */
  for(int i = 2; i <= isize; i++){
    UXDIFX(i).l = + XF(i-1) / DXF(i-1) / XC(i-1) / DXC(i  );
    UXDIFX(i).u = + XF(i+1) / DXF(i  ) / XC(i  ) / DXC(i  );
    UXDIFX(i).c = - XF(i  ) / DXF(i-1) / XC(i-1) / DXC(i  )
                  - XF(i  ) / DXF(i  ) / XC(i  ) / DXC(i  );
  }
  return uxdifx;
}

laplace_t *allocate_and_init_uydifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *uydifx = common_calloc(UYDIFX_SIZE, sizeof(laplace_t));
  /* ! uydifx ! 6 ! */
  for(int i = 1; i <= isize; i++){
    UYDIFX(i).l = + 1. / XC(i-1) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
    UYDIFX(i).u = + 1. / XC(i+1) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
    UYDIFX(i).c = - 1. / XC(i  ) / DXC(i  ) * pow(XF(i  ), 3.) / DXF(i  ) / pow(XC(i  ), 2.)
                  - 1. / XC(i  ) / DXC(i+1) * pow(XF(i+1), 3.) / DXF(i  ) / pow(XC(i  ), 2.);
  }
  return uydifx;
}

laplace_t *allocate_and_init_uzdifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc){
  laplace_t *uzdifx = common_calloc(UZDIFX_SIZE, sizeof(laplace_t));
  /* ! uzdifx ! 6 ! */
  for(int i = 1; i <= isize; i++){
    UZDIFX(i).l = + 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  );
    UZDIFX(i).u = + 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
    UZDIFX(i).c = - 1. / DXC(i  ) * XF(i  ) / DXF(i  ) / XC(i  )
                  - 1. / DXC(i+1) * XF(i+1) / DXF(i  ) / XC(i  );
  }
  return uzdifx;
}


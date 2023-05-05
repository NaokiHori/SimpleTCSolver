#include <stdio.h>
#include <math.h>
#include "param.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "arrays/domain/xf.h"
#include "arrays/domain/xc.h"
#include "arrays/domain/dxf.h"
#include "arrays/domain/dxc.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "fileio.h"
#include "../internal.h"
#include "internal.h"

static double get_ux_x_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf  = domain->xf;
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *ux = fluid->ux->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        // ux-x contribution
        const double cellsize = XF(i  ) * DXC(i  ) * dy * dz;
        // negative direction
        if(1 != i){
          const double duxdx = 1. / DXF(i-1) * (UX(i  , j  , k  ) - UX(i-1, j  , k  ));
          const double w = XC(i-1) * DXF(i-1) / XF(i  ) / DXC(i  );
          dissipation += diffusivity * 0.5 * w * pow(duxdx, 2.) * cellsize;
        }
        // positive direction
        if(isize + 1 != i){
          const double duxdx = 1. / DXF(i  ) * (UX(i+1, j  , k  ) - UX(i  , j  , k  ));
          const double w = XC(i  ) * DXF(i  ) / XF(i  ) / DXC(i  );
          dissipation += diffusivity * 0.5 * w * pow(duxdx, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_ux_y_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf  = domain->xf;
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *ux = fluid->ux->data;
  const double *uy = fluid->uy->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        // ux-y contribution
        const double cellsize = XF(i  ) * DXC(i  ) * dy * dz;
        // negative direction
        {
          const double duxdy = 1. / XF(i  ) / dy * (UX(i  , j  , k  ) - UX(i  , j-1, k  ));
          const double uyx =
            1 == i ? 1. / XF(i  ) * (      UY(i-1, j  , k  )                          )
                   : 1. / XF(i  ) * (0.5 * UY(i-1, j  , k  ) + 0.5 * UY(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duxdy - uyx, 2.) * cellsize;
        }
        // positive direction
        {
          const double duxdy = 1. / XF(i  ) / dy * (UX(i  , j+1, k  ) - UX(i  , j  , k  ));
          const double uyx =
            isize + 1 == i ? 1. / XF(i  ) * (                                UY(i  , j+1, k  ))
                           : 1. / XF(i  ) * (0.5 * UY(i-1, j+1, k  ) + 0.5 * UY(i  , j+1, k  ));
          dissipation += diffusivity * 0.5 * pow(duxdy - uyx, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_ux_z_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf  = domain->xf;
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *ux = fluid->ux->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        // ux-z contribution
        const double cellsize = XF(i  ) * DXC(i  ) * dy * dz;
        // negative direction
        {
          const double duxdz = 1. / dz * (UX(i  , j  , k  ) - UX(i  , j  , k-1));
          dissipation += diffusivity * 0.5 * pow(duxdz, 2.) * cellsize;
        }
        // positive direction
        {
          const double duxdz = 1. / dz * (UX(i  , j  , k+1) - UX(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duxdz, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uy_x_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf  = domain->xf;
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *uy = fluid->uy->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy-x contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        // negative direction
        {
          const double duydx = 1. / DXC(i  ) * (UY(i  , j  , k  ) - UY(i-1, j  , k  ));
          const double w0 = XF(i  ) * DXC(i  ) / XC(i  ) / DXF(i  );
          const double w1 = 1 == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(duydx, 2.) * cellsize;
        }
        // positive direction
        {
          const double duydx = 1. / DXC(i+1) * (UY(i+1, j  , k  ) - UY(i  , j  , k  ));
          const double w0 = XF(i+1) * DXC(i+1) / XC(i  ) / DXF(i  );
          const double w1 = isize == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(duydx, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uy_y_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *ux = fluid->ux->data;
  const double *uy = fluid->uy->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy-y contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        // negative direction
        {
          const double duydy = 1. / XC(i  ) / dy * (UY(i  , j  , k  ) - UY(i  , j-1, k  ));
          const double uxx = 1. / XC(i  ) * (0.5 * UX(i  , j-1, k  ) + 0.5 * UX(i+1, j-1, k  ));
          dissipation += diffusivity * 0.5 * pow(duydy + uxx, 2.) * cellsize;
        }
        // positive direction
        {
          const double duydy = 1. / XC(i  ) / dy * (UY(i  , j+1, k  ) - UY(i  , j  , k  ));
          const double uxx = 1. / XC(i  ) * (0.5 * UX(i  , j  , k  ) + 0.5 * UX(i+1, j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duydy + uxx, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uy_z_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *uy = fluid->uy->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uy-z contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        // negative direction
        {
          const double duydz = 1. / dz * (UY(i  , j  , k  ) - UY(i  , j  , k-1));
          dissipation += diffusivity * 0.5 * pow(duydz, 2.) * cellsize;
        }
        // positive direction
        {
          const double duydz = 1. / dz * (UY(i  , j  , k+1) - UY(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duydz, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uz_x_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xf  = domain->xf;
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double *dxc = domain->dxc;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *uz = fluid->uz->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz-x contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        // negative direction
        {
          const double duzdx = 1. / DXC(i  ) * (UZ(i  , j  , k  ) - UZ(i-1, j  , k  ));
          const double w0 = XF(i  ) * DXC(i  ) / XC(i  ) / DXF(i  );
          const double w1 = 1 == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(duzdx, 2.) * cellsize;
        }
        // positive direction
        {
          const double duzdx = 1. / DXC(i+1) * (UZ(i+1, j  , k  ) - UZ(i  , j  , k  ));
          const double w0 = XF(i+1) * DXC(i+1) / XC(i  ) / DXF(i  );
          const double w1 = isize == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(duzdx, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uz_y_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *uz = fluid->uz->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz-y contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        // negative direction
        {
          const double duzdy = 1. / XC(i  ) / dy * (UZ(i  , j  , k  ) - UZ(i  , j-1, k  ));
          dissipation += diffusivity * 0.5 * pow(duzdy, 2.) * cellsize;
        }
        // positive direction
        {
          const double duzdy = 1. / XC(i  ) / dy * (UZ(i  , j+1, k  ) - UZ(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duzdy, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

static double get_uz_z_contribution(const domain_t *domain, const fluid_t *fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *xc  = domain->xc;
  const double *dxf = domain->dxf;
  const double dy   = domain->dy;
  const double dz   = domain->dz;
  const double *uz = fluid->uz->data;
  const double diffusivity = fluid->diffusivity;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // uz-z contribution
        const double cellsize = XC(i  ) * DXF(i  ) * dy * dz;
        {
          const double duzdz = 1. / dz * (UZ(i  , j  , k  ) - UZ(i  , j  , k-1));
          dissipation += diffusivity * 0.5 * pow(duzdz, 2.) * cellsize;
        }
        // positive direction
        {
          const double duzdz = 1. / dz * (UZ(i  , j  , k+1) - UZ(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(duzdz, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}

void logging_nusselt_compute_dissipation(const domain_t *domain, const fluid_t *fluid, double *retval){
  const double diffusivity = fluid->diffusivity;
  const int isize = domain->mysizes[0];
  const double *xf = domain->xf;
  *retval = 0.;
  // ux contribution
  *retval += get_ux_x_contribution(domain, fluid);
  *retval += get_ux_y_contribution(domain, fluid);
  *retval += get_ux_z_contribution(domain, fluid);
  // uy contribution
  *retval += get_uy_x_contribution(domain, fluid);
  *retval += get_uy_y_contribution(domain, fluid);
  *retval += get_uy_z_contribution(domain, fluid);
  // uz contribution
  *retval += get_uz_x_contribution(domain, fluid);
  *retval += get_uz_y_contribution(domain, fluid);
  *retval += get_uz_z_contribution(domain, fluid);
  MPI_Allreduce(MPI_IN_PLACE, retval, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  const double ref = logging_nusselt_compute_torque_laminar(
      diffusivity,
      /* ri */ XF(      1),
      /* ro */ XF(isize+1),
      /* ui */ param_uy_xm,
      /* uo */ param_uy_xp,
      domain->lengths[1],
      domain->lengths[2]
  );
  *retval = (*retval + diffusivity * (pow(param_uy_xm, 2.) - pow(param_uy_xp, 2.)) * domain->lengths[1] * domain->lengths[2]) / ref / (param_uy_xm / XF(1) - param_uy_xp / XF(isize+1));
}


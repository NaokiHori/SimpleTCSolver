#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "fileio.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/hyxf.h"
#include "array_macros/domain/hyxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/t.h"
#include "array_macros/fluid/lxx.h"
#include "array_macros/fluid/lyx0.h"
#include "array_macros/fluid/lyx1.h"
#if NDIMS == 3
#include "array_macros/fluid/lzx.h"
#endif
#include "array_macros/fluid/lxy.h"
#include "array_macros/fluid/lyy0.h"
#include "array_macros/fluid/lyy1.h"
#if NDIMS == 3
#include "array_macros/fluid/lzy.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lxz.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lyz.h"
#endif
#if NDIMS == 3
#include "array_macros/fluid/lzz.h"
#endif
#include "internal.h"

static int compute_lxx_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lxx = fluid->lxx.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i  );
      const double lij = LXX(i  , j  );
      const double lji = LXX(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double lij = LXX(i  , j  , k  );
        const double lji = LXX(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
#endif
  return 0;
}

static int compute_lyx_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i  );
      const double lij = LYX0(i  , j  ) + LYX1(i  , j  );
      const double lji = LXY(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i  );
        const double lij = LYX0(i  , j  , k  ) + LYX1(i  , j  , k  );
        const double lji = LXY(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
static int compute_lzx_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lzx = fluid->lzx.data;
  const double * restrict lxz = fluid->lxz.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i  );
        const double lij = LZX(i  , j  , k  );
        const double lji = LXZ(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
  return 0;
}
#endif

static int compute_lxy_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i  );
      const double lij = LXY(i  , j  );
      const double lji = LYX0(i  , j  ) + LYX1(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i  );
        const double lij = LXY(i  , j  , k  );
        const double lji = LYX0(i  , j  , k  ) + LYX1(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
#endif
  return 0;
}

static int compute_lyy_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lyy0 = fluid->lyy0.data;
  const double * restrict lyy1 = fluid->lyy1.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i  );
      const double lij = LYY0(i  , j  ) + LYY1(i  , j  );
      const double lji = LYY0(i  , j  ) + LYY1(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double lij = LYY0(i  , j  , k  ) + LYY1(i  , j  , k  );
        const double lji = LYY0(i  , j  , k  ) + LYY1(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
static int compute_lzy_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lzy = fluid->lzy.data;
  const double * restrict lyz = fluid->lyz.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double lij = LZY(i  , j  , k  );
        const double lji = LYZ(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
  return 0;
}
#endif

#if NDIMS == 3
static int compute_lxz_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lzx = fluid->lzx.data;
  const double * restrict lxz = fluid->lxz.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i  );
        const double lij = LXZ(i  , j  , k  );
        const double lji = LZX(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
  return 0;
}
#endif

#if NDIMS == 3
static int compute_lyz_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lzy = fluid->lzy.data;
  const double * restrict lyz = fluid->lyz.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double lij = LYZ(i  , j  , k  );
        const double lji = LZY(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
  return 0;
}
#endif

#if NDIMS == 3
static int compute_lzz_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lzz = fluid->lzz.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double lij = LZZ(i  , j  , k  );
        const double lji = LZZ(i  , j  , k  );
        *value += jd * lij * (lij + lji);
      }
    }
  }
  return 0;
}
#endif

static int compute_x_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict t = fluid->t.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i  );
      const double hx = HXXF(i  );
      const double t_xm = T(i-1, j  );
      const double t_xp = T(i  , j  );
      const double dt = - t_xm + t_xp;
      *value += jd * pow(dt / hx, 2.);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i  );
        const double hx = HXXF(i  );
        const double t_xm = T(i-1, j  , k  );
        const double t_xp = T(i  , j  , k  );
        const double dt = - t_xm + t_xp;
        *value += jd * pow(dt / hx, 2.);
      }
    }
  }
#endif
  return 0;
}

static int compute_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict t = fluid->t.data;
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i  );
      const double hy = HYXC(i  );
      const double t_ym = T(i  , j-1);
      const double t_yp = T(i  , j  );
      const double dt = - t_ym + t_yp;
      *value += jd * pow(dt / hy, 2.);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double hy = HYXC(i  );
        const double t_ym = T(i  , j-1, k  );
        const double t_yp = T(i  , j  , k  );
        const double dt = - t_ym + t_yp;
        *value += jd * pow(dt / hy, 2.);
      }
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
static int compute_z_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict t = fluid->t.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i  );
        const double t_zm = T(i  , j  , k-1);
        const double t_zp = T(i  , j  , k  );
        const double dt = - t_zm + t_zp;
        *value += jd * pow(dt / hz, 2.);
      }
    }
  }
  return 0;
}
#endif

/**
 * @brief compute dissipated energy
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_dissipation(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  // kinetic energy (0) and squared scalar (1)
  double values[2] = {0., 0.};
  // kinetic energy
  compute_lxx_contribution(domain, fluid, values + 0);
  compute_lyx_contribution(domain, fluid, values + 0);
#if NDIMS == 3
  compute_lzx_contribution(domain, fluid, values + 0);
#endif
  compute_lxy_contribution(domain, fluid, values + 0);
  compute_lyy_contribution(domain, fluid, values + 0);
#if NDIMS == 3
  compute_lzy_contribution(domain, fluid, values + 0);
#endif
#if NDIMS == 3
  compute_lxz_contribution(domain, fluid, values + 0);
#endif
#if NDIMS == 3
  compute_lyz_contribution(domain, fluid, values + 0);
#endif
#if NDIMS == 3
  compute_lzz_contribution(domain, fluid, values + 0);
#endif
  // squared scalar
  compute_x_contribution(domain, fluid, values + 1);
  compute_y_contribution(domain, fluid, values + 1);
#if NDIMS == 3
  compute_z_contribution(domain, fluid, values + 1);
#endif
  // communicate among all processes and save
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : values;
  void * recvbuf = values;
  MPI_Reduce(sendbuf, recvbuf, sizeof(values) / sizeof(values[0]), MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 1;
    }
    fprintf(fp, "%8.2f ", time);
    fprintf(fp, "% 18.15e ", fluid_compute_momentum_diffusivity(fluid) * values[0]);
    fprintf(fp, "% 18.15e ", fluid_compute_scalar_diffusivity  (fluid) * values[1]);
    fprintf(fp, "\n");
    fileio.fclose(fp);
  }
  return 0;
}


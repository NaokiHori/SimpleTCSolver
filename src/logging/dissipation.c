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
#include "array_macros/fluid/lxy.h"
#include "array_macros/fluid/lyy0.h"
#include "array_macros/fluid/lyy1.h"
#include "internal.h"

static int compute_lxx_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lxx = fluid->lxx.data;
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i  );
      const double lij = LXX(i  , j  );
      const double lji = LXX(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
  return 0;
}

static int compute_lyx_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy.data;
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i  );
      const double lij = LYX0(i  , j  ) + LYX1(i  , j  );
      const double lji = LXY(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
  return 0;
}

static int compute_lxy_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxf = domain->jdxf;
  const double * restrict lyx0 = fluid->lyx0.data;
  const double * restrict lyx1 = fluid->lyx1.data;
  const double * restrict lxy  = fluid->lxy.data;
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i  );
      const double lij = LXY(i  , j  );
      const double lji = LYX0(i  , j  ) + LYX1(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
  return 0;
}

static int compute_lyy_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict jdxc = domain->jdxc;
  const double * restrict lyy0 = fluid->lyy0.data;
  const double * restrict lyy1 = fluid->lyy1.data;
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i  );
      const double lij = LYY0(i  , j  ) + LYY1(i  , j  );
      const double lji = LYY0(i  , j  ) + LYY1(i  , j  );
      *value += jd * lij * (lij + lji);
    }
  }
  return 0;
}

static int compute_x_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict t = fluid->t.data;
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
  return 0;
}

static int compute_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid,
    double * value
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hyxc = domain->hyxc;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict t = fluid->t.data;
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
  return 0;
}

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
  compute_lxy_contribution(domain, fluid, values + 0);
  compute_lyy_contribution(domain, fluid, values + 0);
  // squared scalar
  compute_x_contribution(domain, fluid, values + 1);
  compute_y_contribution(domain, fluid, values + 1);
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


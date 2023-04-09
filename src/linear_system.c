#include "common.h"
#include "config.h"
#include "sdecomp.h"
#include "domain.h"
#include "linear_system.h"
#include "tdm.h"

static size_t get_nitems(const size_t sizes[NDIMS]){
  // helper function just to multiply sizes and return
  size_t nitems = 1;
  for(int dim = 0; dim < NDIMS; dim++){
    nitems *= sizes[dim];
  }
  return nitems;
}

/**
 * @brief initialise linear solver to update field implicitly
 * @param[in] info    : information about domain decomposition
 * @param[in] glsizes : GLOBAL size of array
 * @return            : structure storing buffers and plans to solve linear systems in each direction
 */
linear_system_t *init_linear_system(sdecomp_info_t * restrict info, const size_t glsizes[NDIMS]){
  const bool implicitx = config.get.implicitx();
  const bool implicity = config.get.implicity();
  const bool implicitz = config.get.implicitz();
  // allocate main structure
  linear_system_t *linear_system = common_calloc(1, sizeof(linear_system_t));
  // pencils (and their sizes) to store input and output of linear systems
  double **x1pncl = &linear_system->x1pncl;
  double **y1pncl = &linear_system->y1pncl;
  double **z2pncl = &linear_system->z2pncl;
  // check size first
  size_t *x1pncl_mysizes = linear_system->x1pncl_mysizes;
  size_t *y1pncl_mysizes = linear_system->y1pncl_mysizes;
  size_t *z2pncl_mysizes = linear_system->z2pncl_mysizes;
  size_t *x1pncl_offsets = linear_system->x1pncl_offsets;
  size_t *y1pncl_offsets = linear_system->y1pncl_offsets;
  size_t *z2pncl_offsets = linear_system->z2pncl_offsets;
  for(int dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, glsizes[dim], &x1pncl_mysizes[dim]);
    sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, glsizes[dim], &y1pncl_mysizes[dim]);
    sdecomp.get_pencil_mysize(info, SDECOMP_Z2PENCIL, dim, glsizes[dim], &z2pncl_mysizes[dim]);
    sdecomp.get_pencil_offset(info, SDECOMP_X1PENCIL, dim, glsizes[dim], &x1pncl_offsets[dim]);
    sdecomp.get_pencil_offset(info, SDECOMP_Y1PENCIL, dim, glsizes[dim], &y1pncl_offsets[dim]);
    sdecomp.get_pencil_offset(info, SDECOMP_Z2PENCIL, dim, glsizes[dim], &z2pncl_offsets[dim]);
  }
  // allocate pencils if needed
  // NOTE: x1pncl is not needed for fully-explicit case,
  //   but I always allocate it here for simplicity (to store delta values)
  *x1pncl =             common_calloc(get_nitems(x1pncl_mysizes), sizeof(double));
  *y1pncl = implicity ? common_calloc(get_nitems(y1pncl_mysizes), sizeof(double)) : NULL;
  *z2pncl = implicitz ? common_calloc(get_nitems(z2pncl_mysizes), sizeof(double)) : NULL;
  // initialise parallel matrix transpose if needed
  // between x1 and y1
  if(implicity){
    sdecomp_transpose_plan_t **plan_f = &linear_system->transposer_x1_to_y1;
    sdecomp_transpose_plan_t **plan_b = &linear_system->transposer_y1_to_x1;
    sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, glsizes, sizeof(double), plan_f);
    sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, glsizes, sizeof(double), plan_b);
  }
  // between x1 and z2
  if(implicitz){
    sdecomp_transpose_plan_t **plan_f = &linear_system->transposer_x1_to_z2;
    sdecomp_transpose_plan_t **plan_b = &linear_system->transposer_z2_to_x1;
    sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Z2PENCIL, glsizes, sizeof(double), plan_f);
    sdecomp.transpose.construct(info, SDECOMP_Z2PENCIL, SDECOMP_X1PENCIL, glsizes, sizeof(double), plan_b);
  }
  // initialise tri-diagonal matrix solvers
  // Thomas algorithm in x direction
  if(implicitx){
    tdm.construct(
        /* size of system */ (int)(x1pncl_mysizes[0]),
        /* number of rhs  */ (int)(x1pncl_mysizes[1] * x1pncl_mysizes[2]),
        /* is periodic    */ false,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_x
    );
  }
  // Thomas algorithm in y direction
  if(implicity){
    tdm.construct(
        /* size of system */ (int)(y1pncl_mysizes[1]),
        /* number of rhs  */ (int)(y1pncl_mysizes[2]), // coefficients vary for each x
        /* is periodic    */ true,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_y
    );
  }
  // Thomas algorithm in z direction
  if(implicitz){
    tdm.construct(
        /* size of system */ (int)(z2pncl_mysizes[2]),
        /* number of rhs  */ (int)(z2pncl_mysizes[0] * z2pncl_mysizes[1]),
        /* is periodic    */ true,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_z
    );
  }
  return linear_system;
}

int linear_system_finalise(linear_system_t * restrict linear_system){
  const bool implicitx = config.get.implicitx();
  const bool implicity = config.get.implicity();
  const bool implicitz = config.get.implicitz();
  if(NULL != linear_system){
    common_free(linear_system->x1pncl);
    if(implicitx){
      tdm.destruct(linear_system->tdm_x);
    }
    if(implicity){
      common_free(linear_system->y1pncl);
      sdecomp.transpose.destruct(linear_system->transposer_x1_to_y1);
      sdecomp.transpose.destruct(linear_system->transposer_y1_to_x1);
      tdm.destruct(linear_system->tdm_y);
    }
    if(implicitz){
      common_free(linear_system->z2pncl);
      sdecomp.transpose.destruct(linear_system->transposer_x1_to_z2);
      sdecomp.transpose.destruct(linear_system->transposer_z2_to_x1);
      tdm.destruct(linear_system->tdm_z);
    }
  }
  common_free(linear_system);
  return 0;
}


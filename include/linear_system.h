#if !defined(LINEAR_SYSTEM_H)
#define LINEAR_SYSTEM_H

#include "sdecomp.h"
#include "domain.h"
#include "tdm.h"

/**
 * @struct linear_system_t
 * @brief structure storing buffers and plans to solve tri-diagonal linear systems in each dimension A x = b
 * @var x1pncl              : buffers to store x1-pencil
 * @var y1pncl              : buffers to store y1-pencil
 * @var z2pncl              : buffers to store z2-pencil
 * @var x1pncl_mysizes      : mysize of x1pencil
 * @var y1pncl_mysizes      : mysize of y1pencil
 * @var z2pncl_mysizes      : mysize of z2pencil
 * @var x1pncl_offsets      : offset of x1pencil
 * @var y1pncl_offsets      : offset of y1pencil
 * @var z2pncl_offsets      : offset of z2pencil
 * @var tdm_[x-z]           : thomas algorithm solvers in all directions
 * @var transposer_xx_to_xx : plans to transpose between two pencils
 */
typedef struct {
  double *x1pncl;
  double *y1pncl;
  double *z2pncl;
  size_t x1pncl_mysizes[NDIMS];
  size_t y1pncl_mysizes[NDIMS];
  size_t z2pncl_mysizes[NDIMS];
  size_t x1pncl_offsets[NDIMS];
  size_t y1pncl_offsets[NDIMS];
  size_t z2pncl_offsets[NDIMS];
  tdm_info_t *tdm_x;
  tdm_info_t *tdm_y;
  tdm_info_t *tdm_z;
  sdecomp_transpose_plan_t *transposer_x1_to_y1;
  sdecomp_transpose_plan_t *transposer_y1_to_x1;
  sdecomp_transpose_plan_t *transposer_x1_to_z2;
  sdecomp_transpose_plan_t *transposer_z2_to_x1;
} linear_system_t;

extern linear_system_t *init_linear_system(sdecomp_info_t * restrict info, const size_t glsizes[NDIMS]);
extern int linear_system_finalise(linear_system_t * restrict linear_system);

#endif // LINEAR_SYSTEM_H

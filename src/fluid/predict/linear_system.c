#include <assert.h>
#include "linear_system.h"
#include "tdm.h"
#include "internal.h"

int solve_in_x(
    const double prefactor,
    const laplacian_t * lapx,
    linear_system_t * linear_system
){
  const size_t * mysizes = linear_system->x1pncl_mysizes;
  const size_t isize = mysizes[0];
  const size_t jsize = mysizes[1];
#if NDIMS == 3
  const size_t ksize = mysizes[2];
#endif
  tdm_info_t * tdm_info = linear_system->tdm_x;
  int tdm_size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &tdm_size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  assert((size_t)tdm_size == isize);
  double * restrict x1pncl = linear_system->x1pncl;
#if NDIMS == 2
  for(size_t j = 0; j < jsize; j++){
    for(size_t i = 0; i < isize; i++){
      tdm_l[i] =    - prefactor * lapx[i].l;
      tdm_c[i] = 1. - prefactor * lapx[i].c;
      tdm_u[i] =    - prefactor * lapx[i].u;
    }
    tdm.solve(tdm_info, x1pncl + j * isize);
  }
#else
  for(size_t k = 0; k < ksize; k++){
    for(size_t j = 0; j < jsize; j++){
      for(size_t i = 0; i < isize; i++){
        tdm_l[i] =    - prefactor * lapx[i].l;
        tdm_c[i] = 1. - prefactor * lapx[i].c;
        tdm_u[i] =    - prefactor * lapx[i].u;
      }
      tdm.solve(tdm_info, x1pncl + (k * jsize + j) * isize);
    }
  }
#endif
  return 0;
}

int solve_in_y(
    const double prefactor,
    const laplacian_t * lapy,
    linear_system_t * linear_system
){
  const size_t * mysizes = linear_system->y1pncl_mysizes;
  const size_t isize = mysizes[0];
  const size_t jsize = mysizes[1];
#if NDIMS == 3
  const size_t ksize = mysizes[2];
#endif
  tdm_info_t * tdm_info = linear_system->tdm_y;
  int tdm_size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &tdm_size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  assert((size_t)tdm_size == jsize);
  double * restrict y1pncl = linear_system->y1pncl;
#if NDIMS == 2
  for(size_t i = 0; i < isize; i++){
    for(size_t j = 0; j < jsize; j++){
      tdm_l[j] =    - prefactor * lapy[i].l;
      tdm_c[j] = 1. - prefactor * lapy[i].c;
      tdm_u[j] =    - prefactor * lapy[i].u;
    }
    tdm.solve(tdm_info, y1pncl + i * jsize);
  }
#else
  for(size_t i = 0; i < isize; i++){
    for(size_t k = 0; k < ksize; k++){
      for(size_t j = 0; j < jsize; j++){
        tdm_l[j] =    - prefactor * lapy[i].l;
        tdm_c[j] = 1. - prefactor * lapy[i].c;
        tdm_u[j] =    - prefactor * lapy[i].u;
      }
      tdm.solve(tdm_info, y1pncl + (i * ksize + k) * jsize);
    }
  }
#endif
  return 0;
}

#if NDIMS == 3
int solve_in_z(
    const double prefactor,
    const laplacian_t * lapz,
    linear_system_t * linear_system
){
  const size_t * mysizes = linear_system->z2pncl_mysizes;
  const size_t isize = mysizes[0];
  const size_t jsize = mysizes[1];
  const size_t ksize = mysizes[2];
  tdm_info_t * tdm_info = linear_system->tdm_z;
  int tdm_size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &tdm_size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  assert((size_t)tdm_size == ksize);
  double * restrict z2pncl = linear_system->z2pncl;
  for(size_t j = 0; j < jsize; j++){
    for(size_t i = 0; i < isize; i++){
      for(size_t k = 0; k < ksize; k++){
        tdm_l[k] =    - prefactor * (*lapz).l;
        tdm_c[k] = 1. - prefactor * (*lapz).c;
        tdm_u[k] =    - prefactor * (*lapz).u;
      }
      tdm.solve(tdm_info, z2pncl + (j * isize + i) * ksize);
    }
  }
  return 0;
}
#endif


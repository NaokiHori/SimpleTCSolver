#include <math.h>
#include "fluid.h"
#include "statistics.h"
#include "array.h"
#include "arrays/fluid/ux.h"
#include "arrays/fluid/uy.h"
#include "arrays/fluid/uz.h"
#include "arrays/ux1.h"
#include "arrays/ux2.h"
#include "arrays/uy1.h"
#include "arrays/uy2.h"
#include "arrays/uz1.h"
#include "arrays/uz2.h"
#include "internal.h"

/**
 * @brief compute ux^1 and ux^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] ux     : x velocity
 */
static void collect_mean_ux(const domain_t * restrict domain, const double * restrict ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict ux1 = g_st_int->ux1->data;
  double * restrict ux2 = g_st_int->ux2->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        UX1(i, j, k) += pow(UX(i, j, k), 1.);
        UX2(i, j, k) += pow(UX(i, j, k), 2.);
      }
    }
  }
}

/**
 * @brief compute uy^1 and uy^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uy     : y velocity
 */
static void collect_mean_uy(const domain_t * restrict domain, const double * restrict uy){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uy1 = g_st_int->uy1->data;
  double * restrict uy2 = g_st_int->uy2->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        UY1(i, j, k) += pow(UY(i, j, k), 1.);
        UY2(i, j, k) += pow(UY(i, j, k), 2.);
      }
    }
  }
}

/**
 * @brief compute uz^1 and uz^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uz     : z velocity
 */
static void collect_mean_uz(const domain_t * restrict domain, const double * restrict uz){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uz1 = g_st_int->uz1->data;
  double * restrict uz2 = g_st_int->uz2->data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        UZ1(i, j, k) += pow(UZ(i, j, k), 1.);
        UZ2(i, j, k) += pow(UZ(i, j, k), 2.);
      }
    }
  }
}

/**
 * @brief accumulate statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : velocity
 */
void collect(const domain_t * restrict domain, const fluid_t * restrict fluid){
  // collect temporally-averaged quantities
  collect_mean_ux(domain, fluid->ux->data);
  collect_mean_uy(domain, fluid->uy->data);
  collect_mean_uz(domain, fluid->uz->data);
  // increment number of samples
  g_st_int->num += 1;
  // schedule next event
  g_st_int->next += g_st_int->rate;
}


#include <stdio.h>
#include <math.h>
#include <float.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "statistics.h"
#include "arrays/ux1.h"
#include "arrays/ux2.h"
#include "arrays/uy1.h"
#include "arrays/uy2.h"
#include "arrays/uz1.h"
#include "arrays/uz2.h"
#include "internal.h"

statistics_internal_t * restrict g_st_int = NULL;

#define CONCAT(a, b) a##b

#define ARRAY_PREPARE(dsetname, DSETNAME) \
  array_prepare( \
      domain, \
      (int [NDIMS][2]){ \
        {CONCAT(DSETNAME, _LNADD_0), CONCAT(DSETNAME, _UNADD_0)}, \
        {CONCAT(DSETNAME, _LNADD_1), CONCAT(DSETNAME, _UNADD_1)}, \
        {CONCAT(DSETNAME, _LNADD_2), CONCAT(DSETNAME, _UNADD_2)}, \
      }, \
      &g_st_int->dsetname \
  )

/**
 * @brief constructor - initialise and allocate internal buffers, schedule collection
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current time (hereafter in free-fall time units)
 */
static int init(const domain_t *domain, const double time){
  g_st_int = common_calloc(1, sizeof(statistics_internal_t));
  ARRAY_PREPARE(ux1, UX1);
  ARRAY_PREPARE(ux2, UX2);
  ARRAY_PREPARE(uy1, UY1);
  ARRAY_PREPARE(uy2, UY2);
  ARRAY_PREPARE(uz1, UZ1);
  ARRAY_PREPARE(uz2, UZ2);
  double after = 0.;
  if(0 != config.get_double("stat_rate", &g_st_int->rate)) return 1;
  if(0 != config.get_double("stat_after",         &after)) return 1;
  g_st_int->next = g_st_int->rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / g_st_int->rate
  );
  // report
  int myrank = 0;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(0 == myrank){
    printf("STATISTICS\n");
    printf("\tnext: % .3e\n", g_st_int->next);
    printf("\trate: % .3e\n", g_st_int->rate);
    fflush(stdout);
  }
  return 0;
}

/**
 * @brief destructor
 */
static void finalise(void){
  array_destroy(g_st_int->ux1);
  array_destroy(g_st_int->ux2);
  array_destroy(g_st_int->uy1);
  array_destroy(g_st_int->uy2);
  array_destroy(g_st_int->uz1);
  array_destroy(g_st_int->uz2);
  common_free(g_st_int);
}

/**
 * @brief getter of a member of g_st_int: next
 * @return : g_st_int->next
 */
static double get_next_time(void){
  return g_st_int->next;
}

// definition
const statistics_t statistics = {
  .init          = init,
  .finalise      = finalise,
  .collect       = collect,
  .output        = output,
  .get_next_time = get_next_time,
};


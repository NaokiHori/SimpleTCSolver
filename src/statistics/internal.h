#if !defined(STATISTICS_INTERNAL_H)
#define STATISTICS_INTERNAL_H

// definition of statistics_internal_t 
/**
 * @struct statistics_internal_t
 * @brief struct storing statistics-related variables
 * @var num        : number of samples which have been collected
 * @var ux1, ux2   : mean and squared ux
 * @var uy1, uy2   : mean and squared uy
 * @var uz1, uz2   : mean and squared uz
 * @var next, rate : next time for collection and collection rate
 */
typedef struct {
  int num;
  array_t * restrict ux1;
  array_t * restrict ux2;
  array_t * restrict uy1;
  array_t * restrict uy2;
  array_t * restrict uz1;
  array_t * restrict uz2;
  double next, rate;
} statistics_internal_t;

extern statistics_internal_t * restrict g_st_int;

extern void collect(const domain_t * restrict domain, const fluid_t * restrict fluid);
extern void output(const domain_t * restrict domain, const int step, const double time);

#endif // STATISTICS_INTERNAL_H

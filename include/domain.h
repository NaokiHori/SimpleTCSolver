#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stdbool.h>
#include "sdecomp.h"

// definition of a structure domain_t
/**
 * @struct domain_t
 * @brief struct storing parameters relevant to spatial domain
 * @var info       : MPI domain decomposition
 * @var is_curved  : cylindrical (true) or Cartesian (false)
 * @var glsizes    : global     number of grid points in each direction
 * @var mysizes    : local (my) number of grid points in each direction
 * @var offsets    : offsets to my starting index in each direction
 * @var lengths    : domain size in each direction
 * @var xf, xc     : cell-face and cell-center locations in x direction
 * @var hxxf, hxxc : wall-normal scale factors at x faces and centers
 * @var hyxf, hyxc : stream-wise scale factors at x faces and centers
 * @var hz         : span-wise scale factor
 * @var jdxf, jdxc : Jacobian determinants at x faces and centers
 */
typedef struct {
  sdecomp_info_t * info;
  bool is_curved;
  size_t glsizes[NDIMS];
  size_t mysizes[NDIMS];
  size_t offsets[NDIMS];
  double lengths[NDIMS];
  double * restrict xf, * restrict xc;
  double * restrict jdxf, * restrict jdxc;
  double * restrict hxxf, * restrict hxxc;
  double * restrict hyxf, * restrict hyxc;
#if NDIMS == 3
  double hz;
#endif
} domain_t;

// constructor
extern int domain_init(
    const char dirname_ic[],
    domain_t * domain
);

// save members which are necessary to restart
extern int domain_save(
    const char dirname[],
    const domain_t * domain
);

#endif // DOMAIN_H

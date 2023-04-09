#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stdbool.h>
#include "sdecomp.h"

typedef struct{
  double l;
  double c;
  double u;
} laplacian_t;

// definition of a structure domain_t 
/**
 * @struct domain_t
 * @brief struct storing parameters relevant to spatial domain
 * @var info     : MPI domain decomposition
 * @var glsizes  : global     number of grid points in each direction
 * @var mysizes  : local (my) number of grid points in each direction
 * @var offsets  : offsets to my starting index in each direction
 * @var lengths  : domain size in each direction
 * @var xf, xc   : cell-face and cell-center locations in x direction
 * @var dxf, dxc : face-to-face and center-to-center distances in x direction
 * @var dy, dz   : grid sizes in homogeneous directions
 * @var uniformx : grid is uniformly distributed (true) or not (false)
 * @var lap      : Laplace operator in each direction for each scalar
 */
typedef struct {
  sdecomp_info_t *info;
  int *glsizes;
  int *mysizes;
  int *offsets;
  double *lengths;
  double *xf, *xc;
  double *dxf, *dxc;
  double dy;
  double dz;
  bool uniformx;
  laplacian_t *uxlapx, *uylapx, *uzlapx, *plapx;
  laplacian_t *uxlapy, *uylapy, *uzlapy, *plapy;
  laplacian_t lapz;
} domain_t;

// constructor and destructor
extern domain_t *domain_init(const char dirname_ic[]);
extern int domain_finalise(domain_t *domain);

// save members which are necessary to restart
extern int domain_save(const char dirname[], const domain_t *domain);

#endif // DOMAIN_H

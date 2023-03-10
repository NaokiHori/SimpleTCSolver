#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "arrays/domain.h"
#include "sdecomp.h"

// struct to store discrete Laplace operators
// d^2q / dx^2
//   \approx "u" q_{n+3/2} + "c" q_{n+1/2} + "l" q_{n-1/2}
//   or
//   \approx "u" q_{n+1  } + "c" q_{n    } + "l" q_{n-1  }
//   (depending on where q are positioned, face or center)
typedef struct {
  double l;
  double c;
  double u;
} laplace_t;

/* ! definition of a structure domain_t ! 21 !*/
/** @struct domain_t
 *  @brief struct storing parameters relevant to spatial domain
 *  @var sdecomp  : MPI domain decomposition
 *  @var glsizes  : global     number of grid points in each direction
 *  @var mysizes  : local (my) number of grid points in each direction
 *  @var offsets  : offsets to my starting index in each direction
 *  @var lengths  : domain size in each direction
 *  @var xf, xc   : cell-face and cell-center locations in x direction
 *  @var dxf, dxc : face-to-face and center-to-center distances in x direction
 *  @var d[xyz]   : representative grid sizes (dx is dummy unless uniform in x)
 */
typedef struct {
  sdecomp_t *sdecomp;
  int *glsizes;
  int *mysizes;
  int *offsets;
  double *lengths;
  double *xf, *xc;
  double *dxf, *dxc;
  double dx, dy, dz;
  laplace_t *lapuxx, *lapuyx, *lapuzx, *lappx;
  laplace_t *lapuxy, *lapuyy, *lapuzy, *lappy;
  laplace_t  lapuxz,  lapuyz,  lapuzz,  lappz;
} domain_t;


// constructor and destructor
extern domain_t *domain_init(void);
extern int domain_finalise(domain_t *domain);

// save members which will be used for post-processing
extern int domain_save(const char dirname[], const domain_t *domain);

#endif // DOMAIN_H

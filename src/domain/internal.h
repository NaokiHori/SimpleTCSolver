#if !defined(DOMAIN_INTERNAL_H)
#define DOMAIN_INTERNAL_H

#include "domain.h"

// coordinate.c
extern double *allocate_and_init_xf(const int isize, const double lx);
extern double *allocate_and_init_xc(const int isize, const double *xf);
extern double *allocate_and_init_dxf(const int isize, const double *xf);
extern double *allocate_and_init_dxc(const int isize, const double *xc);

extern laplace_t *allocate_and_init_lapuxx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);
extern laplace_t *allocate_and_init_lapuyx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);
extern laplace_t *allocate_and_init_lapuzx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);
extern laplace_t *allocate_and_init_lappx (const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);

extern laplace_t *allocate_and_init_lapuxy(const int isize, const double *xf, const double dy);
extern laplace_t *allocate_and_init_lapuyy(const int isize, const double *xc, const double dy);
extern laplace_t *allocate_and_init_lapuzy(const int isize, const double *xc, const double dy);
extern laplace_t *allocate_and_init_lappy (const int isize, const double *xc, const double dy);

extern laplace_t init_lapuxz(const double dz);
extern laplace_t init_lapuyz(const double dz);
extern laplace_t init_lapuzz(const double dz);
extern laplace_t init_lappz (const double dz);

// init.c
extern domain_t *domain_init(void);

// finalise.c
extern int domain_finalise(domain_t *domain);

// save.c
extern int domain_save(const char dirname[], const domain_t *domain);

#endif // DOMAIN_INTERNAL_H

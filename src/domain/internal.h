#if !defined(DOMAIN_INTERNAL_H)
#define DOMAIN_INTERNAL_H

#include "domain.h"

// coordinate.c
extern double *allocate_and_init_xf(const int isize, const double lx);
extern double *allocate_and_init_xc(const int isize, const double *xf);
extern double *allocate_and_init_dxf(const int isize, const double *xf);
extern double *allocate_and_init_dxc(const int isize, const double *xc);

extern laplace_t *allocate_and_init_uxdifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);
extern laplace_t *allocate_and_init_uydifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);
extern laplace_t *allocate_and_init_uzdifx(const int isize, const double *xf, const double *xc, const double *dxf, const double *dxc);

extern laplace_t *allocate_and_init_uxdify(const int isize, const double *xf, const double dy);
extern laplace_t *allocate_and_init_uydify(const int isize, const double *xc, const double dy);
extern laplace_t *allocate_and_init_uzdify(const int isize, const double *xc, const double dy);

extern laplace_t init_uxdifz(const double dz);
extern laplace_t init_uydifz(const double dz);
extern laplace_t init_uzdifz(const double dz);

// init.c
extern domain_t *domain_init(void);

// finalise.c
extern int domain_finalise(domain_t *domain);

// save.c
extern int domain_save(const char dirname[], const domain_t *domain);

#endif // DOMAIN_INTERNAL_H

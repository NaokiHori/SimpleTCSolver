#if !defined(DOMAIN_H)
#define DOMAIN_H

typedef struct {
  int glisize, gljsize, glksize;
  double lx, ly, lz;
  double *xf, *xc;
  double dy, dz;
} domain_t;

// constructor and destructor
extern domain_t *domain_init(void);
extern int domain_finalise(domain_t *domain);

#endif // DOMAIN_H

#include "common.h"
#include "domain.h"
#include "internal.h"


/**
 * @brief destruct a structure domain_t
 * @param[inout] domain : structure to be cleaned-up
 * @return              : error code
 */
int domain_finalise(domain_t *domain){
  common_free(domain->xf);
  common_free(domain->xc);
  common_free(domain->dxf);
  common_free(domain->dxc);
  common_free(domain->glsizes);
  common_free(domain->mysizes);
  common_free(domain->offsets);
  common_free(domain->lengths);
  common_free(domain->lapuxx);
  common_free(domain->lapuyx);
  common_free(domain->lapuzx);
  common_free(domain->lappx );
  common_free(domain->lapuxy);
  common_free(domain->lapuyy);
  common_free(domain->lapuzy);
  common_free(domain->lappy );
  sdecomp_finalise(domain->sdecomp);
  common_free(domain);
  return 0;
}


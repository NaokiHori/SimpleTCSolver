#include "common.h"
#include "domain.h"


/**
 * @brief destruct a structure domain_t
 * @param[in,out] domain : structure to be cleaned-up
 * @return               : error code
 */
int domain_finalise(domain_t *domain){
  common_free(domain->xf);
  common_free(domain->xc);
  common_free(domain->dxf);
  common_free(domain->dxc);
  common_free(domain->uxlapx);
  common_free(domain->uylapx);
  common_free(domain->uzlapx);
  common_free(domain->plapx);
  common_free(domain->uxlapy);
  common_free(domain->uylapy);
  common_free(domain->uzlapy);
  common_free(domain->plapy);
  common_free(domain->glsizes);
  common_free(domain->mysizes);
  common_free(domain->offsets);
  common_free(domain->lengths);
  sdecomp.destruct(domain->info);
  common_free(domain);
  return 0;
}


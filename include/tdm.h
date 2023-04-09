#if !defined(TDM_H)
#define TDM_H

#include <stdbool.h>

typedef struct tdm_info_t_ tdm_info_t;

typedef struct {
  int (* const construct)(
      const int size,
      const int nrhs,
      const bool is_periodic,
      const bool is_complex,
      tdm_info_t **info
  );
  int (* const get_l)(
      const tdm_info_t * restrict info,
      double * restrict *l
  );
  int (* const get_c)(
      const tdm_info_t * restrict info,
      double * restrict *c
  );
  int (* const get_u)(
      const tdm_info_t * restrict info,
      double * restrict *u
  );
  int (* const get_size)(
      const tdm_info_t * restrict info,
      int * restrict size
  );
  int (* const get_nrhs)(
      const tdm_info_t * restrict info,
      int * restrict nrhs
  );
  int (* const solve)(
      const tdm_info_t * restrict info,
      void * restrict data
  );
  int (* const destruct)(
      tdm_info_t * restrict info
  );
} tdm_t;

extern const tdm_t tdm;

#endif // TDM_H

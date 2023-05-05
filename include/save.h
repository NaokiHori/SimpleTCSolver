#pragma once

#include "domain.h"
#include "fluid.h"

typedef struct save_t_ {
  // constructor
  int (* const init)(
    const domain_t *domain,
    const double time
  );
  // destructor
  void (* const finalise)(
      void
  );
  // make space to save flow fields, save some scalars
  int (* const prepare)(
      const domain_t * restrict domain,
      const int step,
      const double time,
      char * restrict * dirname
  );
  // getter, next timing to call "output"
  double (* const get_next_time)(
      void
  );
} save_t;

extern const save_t save;


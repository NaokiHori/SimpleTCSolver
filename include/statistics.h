#pragma once

#include "domain.h"
#include "fluid.h"

typedef struct {
  // constructor
  int (* const init)(
      const domain_t * restrict domain,
      const double time
  );
  // destructor
  void (* const finalise)(
      void
  );
  // collecting statistics
  void (* const collect)(
      const domain_t * restrict domain,
      const fluid_t * restrict fluid
  );
  // save statistics to files
  void (* const output)(
      const domain_t * restrict domain,
      const int step,
      const double time
  );
  // getter, next timing to call "collect"
  double (* const get_next_time)(
      void
  );
} statistics_t;

extern const statistics_t statistics;


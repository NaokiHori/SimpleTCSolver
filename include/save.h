#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"

typedef struct save_t_ {
  // constructor
  void (*init)(
    const domain_t *domain,
    const double time,
    const double rate,
    const double after
  );
  // destructor
  void (*finalise)(
      void
  );
  // save flow fields etc. to files
  void (*output)(
      const domain_t *domain,
      const int step,
      const double time,
      const fluid_t *fluid
  );
  // getter, next timing to call "output"
  double (*get_next_time)(
      void
  );
} save_t;

extern const save_t save;

#endif // SAVE_H

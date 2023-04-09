#if !defined(STATISTICS_H)
#define STATISTICS_H

#include "domain.h"
#include "fluid.h"

typedef struct {
  // constructor
  void (* const init)(
      const domain_t *domain,
      const double time,
      const double rate,
      const double after
  );
  // destructor
  void (*finalise)(
      void
  );
  // collecting statistics
  void (*collect)(
      const domain_t *domain,
      const fluid_t *fluid
  );
  // save statistics to files
  void (*output)(
      const domain_t *domain,
      const int step,
      const double time
  );
  // getter, next timing to call "collect"
  double (*get_next_time)(
      void
  );
} statistics_t;

extern const statistics_t statistics;

#endif // STATISTICS_H

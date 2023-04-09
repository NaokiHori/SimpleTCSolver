#if !defined(CONFIG_H)
#define CONFIG_H

#include <stdbool.h>

typedef struct {
  char * (*dirname_ic)(void);
  bool   (*implicitx)(void);
  bool   (*implicity)(void);
  bool   (*implicitz)(void);
  double (*timemax)(void);
  double (*wtimemax)(void);
  double (*log_rate)(void);
  double (*save_rate)(void);
  double (*save_after)(void);
  double (*stat_rate)(void);
  double (*stat_after)(void);
  double (*visu_rate)(void);
  double (*visu_after)(void);
  double (*coef_dt_adv)(void);
  double (*coef_dt_dif)(void);
  double (*Re)(void);
} config_getter_t;

typedef struct {
  // constructor
  int (*construct)(
      void
  );
  // destructor
  int (*destruct)(
      void
  );
  // save information to files
  int (*output)(
      const char dirname[]
  );
  // getters
  const config_getter_t get;
} config_t;

extern const config_t config;

#endif // CONFIG_H

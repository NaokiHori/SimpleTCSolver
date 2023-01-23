#include "common.h"
#include "sdecomp.h"
#include "linear_system.h"



int linear_system_finalise(linear_system_t * restrict linear_system){
  // this function is so general that
  //   we need to worry that this might be called
  //   even if linear_system is not initialised
  if(linear_system != NULL){
    common_free(linear_system->x1pncl);
    common_free(linear_system->y1pncl);
    common_free(linear_system->z1pncl);
    common_free(linear_system->tdm_l);
    common_free(linear_system->tdm_c);
    common_free(linear_system->tdm_u);
    sdecomp_transpose_finalise(linear_system->transposer_x1_to_y1);
    sdecomp_transpose_finalise(linear_system->transposer_y1_to_x1);
    sdecomp_transpose_finalise(linear_system->transposer_y1_to_z1);
    sdecomp_transpose_finalise(linear_system->transposer_z1_to_y1);
  }
  // free(NULL) is a defined behaviour
  common_free(linear_system);
  return 0;
}


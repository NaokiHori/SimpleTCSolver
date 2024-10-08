#if !defined(INCLUDE_ARRAY_MACROS_FLUID_T_H)
#define INCLUDE_ARRAY_MACROS_FLUID_T_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define T(I, J, K) (t[(I  ) + (isize+2) * ((J  ) + (jsize+2) * (K  ))])
#define T_NADDS (int [NDIMS][2]){ {1, 1}, {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_T_H

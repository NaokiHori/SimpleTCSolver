#if !defined(INCLUDE_ARRAY_MACROS_FLUID_LYX1_H)
#define INCLUDE_ARRAY_MACROS_FLUID_LYX1_H

// This file is generated by tools/define_arrays.py

// [1 : isize+1], [0 : jsize+1], [0 : ksize+1]
#define LYX1(I, J, K) (lyx1[(I-1) + (isize+1) * ((J  ) + (jsize+2) * (K  ))])
#define LYX1_NADDS (int [NDIMS][2]){ {0, 1}, {1, 1}, {1, 1}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_LYX1_H

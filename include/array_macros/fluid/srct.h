#if !defined(INCLUDE_ARRAY_MACROS_FLUID_SRCT_H)
#define INCLUDE_ARRAY_MACROS_FLUID_SRCT_H

// This file is generated by tools/define_arrays.py

// [1 : isize+0], [1 : jsize+0], [1 : ksize+0]
#define SRCT(I, J, K) (srct[(I-1) + (isize+0) * ((J-1) + (jsize+0) * (K-1))])
#define SRCT_NADDS (int [NDIMS][2]){ {0, 0}, {0, 0}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_FLUID_SRCT_H

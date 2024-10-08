#if !defined(INCLUDE_ARRAY_MACROS_STATISTICS_UZ1_H)
#define INCLUDE_ARRAY_MACROS_STATISTICS_UZ1_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [1 : jsize+0], [1 : ksize+0]
#define UZ1(I, J, K) (uz1[(I  ) + (isize+2) * ((J-1) + (jsize+0) * (K-1))])
#define UZ1_NADDS (int [NDIMS][2]){ {1, 1}, {0, 0}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_STATISTICS_UZ1_H

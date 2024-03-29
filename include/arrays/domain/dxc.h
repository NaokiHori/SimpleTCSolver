#if !defined(INCLUDE_ARRAYS_DOMAIN_DXC_H)
#define INCLUDE_ARRAYS_DOMAIN_DXC_H

/*
 * This file is generated by
 *   tools/define_arrays.py
 */

// [1:isize+1]
// number of additional cells w.r.t. 1, isize
#define DXC_LNADD_0 0
#define DXC_UNADD_0 1
// number of cells
#define DXC_NITEMS_0(I) (I+1)
#define DXC(I) (dxc[(I-1)])

#endif // INCLUDE_ARRAYS_DOMAIN_DXC_H

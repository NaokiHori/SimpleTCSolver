#if !defined(FILEIO_H)
#define FILEIO_H

#include <stdio.h>  // FILE, size_t
#include <mpi.h>    // MPI_Datatype

/* general file opener / closer */
extern FILE *fileio_fopen(const char *path, const char *mode);
extern int fileio_fclose(FILE *stream);

/* prepare directory to be stored */
extern char *fileio_create_dirname(const char prefix[], const int step);
extern int fileio_mkdir(const char dirname[]);

/* for snpyio lib */
// datatypes, which are embedded in NPY files (argument: dtype)
// declared here, defined in src/fileio/entrypoint.c
// 1-byte boolean
extern const char NPY_BOL[];
// 4-byte little-endian integer
extern const char NPY_INT[];
// 8-byte little-endian floating point
extern const char NPY_DBL[];
// wrapper functions of snpyio header writer / reader
// serial io (called by one process)
extern int fileio_r_serial(const char dirname[], const char dsetname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize,       void *data);
extern int fileio_w_serial(const char dirname[], const char dsetname[], const size_t ndims, const size_t *shape, const char dtype[], const size_t totalsize, const void *data);
// multi-dimensional array (called by all processes)
extern int fileio_r_nd_parallel(const MPI_Comm comm, const char dirname[], const char dsetname[], const size_t ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts, const char dtype[], const MPI_Datatype mpi_datatype, void *data);
extern int fileio_w_nd_parallel(const MPI_Comm comm, const char dirname[], const char dsetname[], const size_t ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts, const char dtype[], const MPI_Datatype mpi_datatype, const void *data);

#endif // FILEIO_H

Set of Fortran interfaces to CUDA libraries for GCC/OpenACC.
Modules include cublas, cublas_v2, cublasxt, openacc_cublas, cufft.

cublas_core is an internal module of shared definitions used by some of the
above listed modules.

To compile and .mod file: gfortran -c -O2 -g <each-module-name>.f90

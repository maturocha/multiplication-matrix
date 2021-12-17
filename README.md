# Matrix Multiplication

Multiplication of 2 matrices of N*N weighted by a weight vector W of size N

# Execution

## Secuential
Compile: ``
   gcc secuential_multiplication.c -std=c99 -lm -o secuential
``

Execute: ``
   ./secuential
``

## MPI

Compile: ``
   mpicc mpi_multiplication.c -std=c99 -lm -o mpi
``
Execute: ``
   mpirun -np [nro_process] ./mpi
``

# OpenMP
Compile: ``
   mpicc openmp_multiplication.c -lm -fopenmp -o openmp
``

Execute: ``
   mpirun -np [nro_process] openmp [nro_threads]
``

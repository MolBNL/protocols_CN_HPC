#!/bin/bash

export MPI=8
export OMP=1

export OMP_NUM_THREADS=$OMP
export CUDA_VISIBLE_DEVICES=0

export LMPRUN="mpirun -np $MPI --bind-to none lmp -sf omp -pk omp $OMP"

$LMPRUN -in 01_thermo_GK.lmp

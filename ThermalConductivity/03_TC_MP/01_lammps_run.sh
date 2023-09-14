#!/bin/bash

export MPI=8
export OMP=1

export OMP_NUM_THREADS=$OMP
export CUDA_VISIBLE_DEVICES=0

export LMPRUN="mpirun -np $MPI --bind-to none lmp -sf omp -pk omp $OMP"

$LMPRUN -in 01_TC_MP.lmp -v dir x
$LMPRUN -in 01_TC_MP.lmp -v dir z
$LMPRUN -in 01_TC_MP.lmp -v dir y


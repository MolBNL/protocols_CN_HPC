#!/bin/bash

export MPI=8
export OMP=1

export OMP_NUM_THREADS=$OMP
export CUDA_VISIBLE_DEVICES=0

export LMPRUN="mpirun -np $MPI --bind-to none lmp -sf omp -pk omp $OMP"

$LMPRUN -in 01_NVTeq1.lmp  
$LMPRUN -in 02_NPTeq1.lmp  
$LMPRUN -in 03_NPTeq2.lmp  
$LMPRUN -in 04_NPTeq3.lmp  
$LMPRUN -in 05_NPTeq4.lmp  

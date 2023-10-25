#!/bin/bash

export MPI=8

export RUNPARALLEL="mpirun -np $MPI --bind-to none lmp -sf opt"

mkdir -p Surf
mkdir -p Sys
mkdir -p Scr

lmp -in 01_create_box_and_surf.lmpin  
$RUNPARALLEL -in 02_min_surf.lmpin   
$RUNPARALLEL -in 03_eq_surf_np1.lmpin   
$RUNPARALLEL -in 04_eq_surf_np2.lmpin   
$RUNPARALLEL -in 05_eq_surf_nv.lmpin   
lmp -in 06_create_tip.lmpin  
$RUNPARALLEL -in 07_eq_all.lmpin   
$RUNPARALLEL -in 08_ind.lmpin   

$RUNPARALLEL -in 09_scr.lmpin -var ind -6 
$RUNPARALLEL -in 09_scr.lmpin -var ind -5 
$RUNPARALLEL -in 09_scr.lmpin -var ind -4 
$RUNPARALLEL -in 09_scr.lmpin -var ind -3 
$RUNPARALLEL -in 09_scr.lmpin -var ind -2 
$RUNPARALLEL -in 09_scr.lmpin -var ind -1 
$RUNPARALLEL -in 09_scr.lmpin -var ind 0 
$RUNPARALLEL -in 09_scr.lmpin -var ind 1 
$RUNPARALLEL -in 09_scr.lmpin -var ind 2





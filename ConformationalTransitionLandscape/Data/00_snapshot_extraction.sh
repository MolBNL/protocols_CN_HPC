#!/bin/bash

prepare () { 
mpirun -np 8 pmemd.MPI -O -p "$1".solvated.prmtop \
-i min_all.mdin \
-c "$1".equilibrated_restart.ncrst \
-ref "$1".equilibrated_restart.ncrst \
-o "$1".ene_minimized.out -r "$1".ene_minimized.ncrst

cpptraj << EOF
  parm "$1".solvated.prmtop
  trajin "$1".ene_minimized.ncrst
  strip :WAT,Na+,Cl-
  autoimage 
  trajout "$1".pdb
  go
  quit
EOF
mv "$1".pdb ./reference_pdbs/
}

prepare "ref1"


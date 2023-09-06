#!/usr/bin/env python3

import parmed as pmd
import pytraj as pt


# INPUT FILES (you can find the files uesd here in the provided compressed file)
prmtop = "HSA.TRM.solv.prmtop"
ncrst = "HSA.TRM.solv.equil.ncrst"

# OUTPUT FILES
autoimg_nc = "HSA.TRM.solv.equil.autoimg.nc"
gro_top = "HSA.TRM.top"
gro_coords = "HSA.TRM.solv.equil.autoimg.gro"


# LOAD THE RESTART FILE WITH PYTRAJ, AND REIMAGE ALL THE SOLUTE ATOMS INSIDE THE ORIGINAL BOX
traj = pt.load(filename=ncrst, top=prmtop)
traj = pt.autoimage(traj)

# SAVE THE AUTOIMAGED RESTART FILE
pt.write_traj(filename=autoimg_nc, traj=traj, overwrite=True)

# IDENTIFY THE INDECES OF THE CA ATOMS OF THE RESIDUES WITHIN 0.5 nm OF THE TRM RESIDUE
traj.top.set_reference(traj[0])
bs_atom_indeces = traj.top.select(':TRM<:5.0&@CA')+1
bs_atom_indeces = ','.join(map(str, bs_atom_indeces))

# IDENTIFY THE INDECES OF THE TRM RESIDUE HEAVY ATOMS
drug_atom_indeces = traj.top.select(':TRM&!@H=')+1
drug_atom_indeces = ','.join(map(str, drug_atom_indeces))


print("Now 'parmed' will be used to generate the GROMACS topology and starting coordinates")
amber = pmd.load_file(prmtop, autoimg_nc)
amber.save(gro_top, overwrite=True)
amber.save(gro_coords, overwrite=True)
print("Done! You'll find this new files:")
print(f"- gromacs topology: {gro_top}")
print(f"- gromacs coordinates: {gro_coords}")

print("The string of CA atoms of residues that are within 0.5 nm of the ligand is:")
print(bs_atom_indeces)
print("\n")
print("The string of the ligand heavy atoms is:")
print(drug_atom_indeces)

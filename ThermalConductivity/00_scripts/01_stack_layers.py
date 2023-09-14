#!/usr/bin/env python3

from math import floor as math_floor
from copy import deepcopy
from lmp_classes import LmpDataFile

N_GRA_SHEETS = 3   # How many (stacked) sheets of graphene do you want?

D0_GRA_PAR = 1.5     # initial distance between graphene sheet and paraffin
D0_GRA_GRA = 3     # initial distance between the graphene sheets

GRA_INI_FILE = "../01_prepare_lmpdata/charmm_gui_graphene.lmpdata"
PAR_INI_FILE = "../01_prepare_lmpdata/charmm_gui_step3_input.lmpdata"


# NOTE 1 for paraffin:
# The paraffin data file has a weird concept of molecule, each fragment is considered a molecule.
# With `adjust_molecules` we take care of that and identify the correct number of molecules.

# NOTE 2 for paraffin:
# The pair coeffs are setted in the lmpinput file, and are from (a simplified version of) GAFF.
# To simplify the description here, atom types 1-2 (H) & 3-4 (C) are merged, as well as some bonded types

graphene = LmpDataFile(GRA_INI_FILE)
paraffin = LmpDataFile(PAR_INI_FILE, adjust_molecules=True)

print("Read graphene data file:", graphene)
print("Read paraffin data file:", paraffin)

################################################################################
# adjusting paraffin - # see NOTE 2 above...
################################################################################
for atom in paraffin.atoms:
    if atom.type in [1, 2]:
        atom.type = 2
    elif atom.type in [3, 4]:
        atom.type = 1
for bond in paraffin.bonds:
    if bond.type in [1, 2]:
        bond.type = 1
    elif bond.type in [3, 4]:
        bond.type = 2
for angle in paraffin.angles:
    if angle.type in [1, 2]:
        angle.type = 1
    elif angle.type in [3, 4, 5]:
        angle.type = 3
    elif angle.type in [6, 7]:
        angle.type = 2
for dihedral in paraffin.dihedrals:
    if dihedral.type in range(1, 9):
        dihedral.type = 1
    elif dihedral.type in range(9, 12):
        dihedral.type = 2
    elif dihedral.type in [12, 13]:
        dihedral.type = 3
for atom_idx, atom in enumerate(paraffin.atoms):
    atom_idx_in_mol = atom_idx % paraffin.n_atm_per_mol
    if atom_idx_in_mol in [0, 58]:
        atom.charge = -0.0951
    elif atom_idx_in_mol in range(1, 4) or atom_idx_in_mol in range(59, 62):
        atom.charge = 0.0327
    elif atom_idx_in_mol in [4, 55]:
        atom.charge = -0.0804
    elif atom_idx_in_mol in range(5, 7) or atom_idx_in_mol in range(56, 58):
        atom.charge = 0.0387

################################################################################
# expand paraffin box to exactly the graphene box in x and y directions
# and adjust the image flaxs
################################################################################
stretch_factor_x = graphene.xhi / paraffin.xhi
stretch_factor_y = graphene.yhi / paraffin.yhi
print(f"To match the graphene sheet x direction the stretchis {stretch_factor_x:.4f}")
print(f"To match the graphene sheet y direction the stretchis {stretch_factor_y:.4f}")

paraffin_lenX = paraffin.xhi - paraffin.xlo
paraffin_lenY = paraffin.yhi - paraffin.ylo
for atom in paraffin.atoms:
    if atom.x > paraffin.xhi:
        atom.x -= paraffin_lenX
        atom.nx = 1
    elif atom.x < paraffin.xlo:
        atom.x += paraffin_lenX
        atom.nx = -1
    if atom.y > paraffin.yhi:
        atom.y -= paraffin_lenY
        atom.ny = 1
    elif atom.y < paraffin.ylo:
        atom.y += paraffin_lenY
        atom.ny = -1
    atom.x *= stretch_factor_x
    atom.y *= stretch_factor_y

################################################################################
# adjusting the first graphene sheet,
# which in order is under all the paraffin chains
################################################################################
for atom in graphene.atoms:
    atom.id += paraffin.n_atms
    atom.molecule_id = paraffin.n_mols + 1
    atom.type = 3

################################################################################
# create stacked graphene
# (just the atoms, as there are no bonded interactions in tersoff FF)
################################################################################
stacked_graphene_atoms = []
for n in range(N_GRA_SHEETS):
    for atom in graphene.atoms:
        stacked_graphene_atoms.append(deepcopy(atom))
        if n:
            stacked_graphene_atoms[-1].id += n*graphene.n_atms
            stacked_graphene_atoms[-1].molecule_id += n
            stacked_graphene_atoms[-1].z += n*D0_GRA_GRA 

################################################################################
# place paraffin on top of the stacked graphene
################################################################################
max_gra_z = max(stacked_graphene_atoms, key=lambda k: k.z).z
min_gra_z = min(stacked_graphene_atoms, key=lambda k: k.z).z

min_par_z = min(paraffin.atoms, key=lambda k: k.z).z
for atom in paraffin.atoms:
    atom.z += (max_gra_z - min_par_z + D0_GRA_PAR)
max_par_z = max(paraffin.atoms, key=lambda k: k.z).z

################################################################################
# finalize and write to file
################################################################################
tot_num_atoms = paraffin.n_atms + graphene.n_atms*N_GRA_SHEETS
zlo = min_gra_z - D0_GRA_PAR/2
zhi = max_par_z + D0_GRA_PAR/2

with open('../01_prepare_lmpdata/graphene_paraffin.lmpdata', 'w') as f:
    f.write(f'# {N_GRA_SHEETS} graphene sheets + paraffin \n')
    f.write('\n')
    f.write(f'   {len(stacked_graphene_atoms) + paraffin.n_atms} atoms \n')
    f.write(f'   {paraffin.n_bnds} bonds \n')
    f.write(f'   {paraffin.n_angs} angles \n')
    f.write(f'   {paraffin.n_dihs} dihedrals \n')
    f.write(f'   0 impropers \n')
    f.write('\n')
    f.write('   3 atom types \n')
    f.write('   2 bond types \n')
    f.write('   3 angle types \n')
    f.write('   3 dihedral types \n')
    f.write('   0 improper types \n')
    f.write('\n')
    f.write(f'  {graphene.xlo:.5f} {graphene.xhi:.5f} xlo xhi \n')
    f.write(f'  {graphene.ylo:.5f} {graphene.yhi:.5f} ylo yhi \n')
    f.write(f'  {zlo:.5f} {zhi:.5f} zlo zhi \n')
    f.write('\n')
    f.write('Masses \n')
    f.write('\n')
    f.write('1  12.01 \n')
    f.write('2  1.008 \n')
    f.write('3  12.01 \n')
    f.write('\n')
    f.write('Atoms \n')
    f.write('\n')
    for atom in paraffin.atoms:
        atom.write_to(f)
    for atom in stacked_graphene_atoms:
        atom.write_to(f)
    f.write('\n')
    f.write('Bonds \n')
    f.write('\n')
    for bond in paraffin.bonds:
        bond.write_to(f)
    f.write('\n')
    f.write('Angles \n')
    f.write('\n')
    for angle in paraffin.angles:
        angle.write_to(f)
    f.write('\n')
    f.write('Dihedrals \n')
    f.write('\n')
    for dihedral in paraffin.dihedrals:
        dihedral.write_to(f)
    f.write('\n')

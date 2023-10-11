#!/usr/bin/env python3

from os import makedirs, rename
from subprocess import run
from utils import cd, ReferenceStructure

def filter_contacts_and_write_reference(list_of_structures):
    """selects the contacts and write a draft of a plumed input file"""
    for structure in list_of_structures:
        structure.select_unique_contacts(list_of_others=[s for s in list_of_structures if s != structure])
    with open(f"general_reference_cmap_plumed_input.dat", 'w') as f_plm:
        f_plm.write("UNITS LENGTH=A TIME=ns ENERGY=kcal/mol\n\n")
        f_plm.write("c: CONTACTMAP ...\n")
        starting_idx = 1
        for structure in list_of_structures:
            starting_idx = structure.write_reference_cmap_plumed_creation(f_plm, start_idx=starting_idx)
        f_plm.write("...\n\n")
        f_plm.write(f"PRINT STRIDE=1 ARG=* FILE=general_reference_cmap_plumed_output.dat\n")

def run_plumed_to_create_references(list_of_structures):
    """ Runs the plumed driver to compute the reference distance for every contact"""
    for structure in list_of_structures:
        cmd = "export PLUMED_MAXBACKUP=0"
        cmd += " ; plumed driver --plumed general_reference_cmap_plumed_input.dat"
        cmd += f" --mf_pdb ../reference_pdbs/{structure.name}.pdb"
        run(cmd, shell=True)
        rename(
            "general_reference_cmap_plumed_output.dat",
            f"general_reference_cmap_plumed_output_{structure.name}.dat")

def prepare_plumed_inputs(list_of_structures):
    """ prepare the plumed input files"""
    with open(f"general_reference_cmap_plumed_input.dat") as f:
        lines = [l for l in f if l.startswith("  ATOMS")]
    for structure in list_of_structures:
        with open(f"general_reference_cmap_plumed_output_{structure.name}.dat") as f:
            references = f.readlines()[1].split()
        with open(f"cmap_{structure.name}.dat", 'w') as f_out:
            f_out.write(f"c_{structure.name}: CONTACTMAP ...\n")
            for idx, line in enumerate(lines, start=1):
                f_out.write(f"  {line.split('#')[0].strip()} REFERENCE{idx}={references[idx]}")
                dist_kind = line.split(",")[-1].split()[-1]
                if dist_kind in ["standard", "BB"]:
                    weight = 1.00
                elif dist_kind in ["salt_bridge"]:
                    weight = 3.00
                elif dist_kind in ["polar", "pi_charged", "pi_pi"]:
                    weight = 1.75
                else:
                    weight = 1.00
                f_out.write(f" WEIGHT{idx}={weight}\n")
            f_out.write("  CMDIST\n...\n\n")

    with open(f"metad.dat", 'w') as f_out:
        f_out.write("UNITS LENGTH=A TIME=ns ENERGY=kcal/mol\n\n")
        for structure in list_of_structures:
            f_out.write(f"INCLUDE FILE=cmap_{structure.name}.dat\n")
        final_names = [f"c_{structure.name}" for structure in list_of_structures]
        f_out.write(f"\nPRINT STRIDE=1 ARG={','.join(final_names)} FILE=metad_results.dat\n")

    for structure in list_of_structures:
        cmd = "export PLUMED_MAXBACKUP=0"
        cmd += f" ; plumed driver --plumed metad.dat"
        cmd += f" --mf_pdb ../reference_pdbs/{structure.name}.pdb"
        run(cmd, shell=True)
        rename(
            "metad_results.dat",
            f"metad_results_{structure.name}.dat")


################################################################################################################################
# In this example, we have three different reference structures,
# and the contact map will be computed for all three of them.
# The script can be easily modified to work with a different number of reference structures
################################################################################################################################
def do_things(list_of_names):
    wd = "_".join(list_of_names)
    makedirs(wd, exist_ok=True)
    with cd(wd):
        selected_structures = [s for s in structures if s.name in list_of_names]
        filter_contacts_and_write_reference(selected_structures)
        run_plumed_to_create_references(selected_structures)
        prepare_plumed_inputs(selected_structures)

structure_names_list = ["ref1", "ref2", "ref3"]
structures = [ReferenceStructure(name=name) for name in structure_names_list]
for structure in structures:
    structure.load_atoms()
do_things(structure_names_list)

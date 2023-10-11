# Protocol to characterize the conformational landscape of a protein, starting from two or more known structures.

In this protocol we are going to construct a free energy surface that aims to describe the conformational landscape of a protein, starting from two (or more) equilibrium structures of the protein.


## 1. Select and equilibrate the starting structures
The starting structures of the protein can be retrived, as an example, from the [Protein Data Bank](https://www.rcsb.org). Once the structures have been selected, we should assess that they contains the same aminoacids, starting and ending from the same protein residue. To be able to use the provided scripts, also the same protonation state of each aminoacid must be selected for each structure. Otherwise, the scripts should be adjusted accordingly.

Each system should be equilibrated, to obtain a starting reference structure for each of the protein states.

The equilibrated and energy-minimized structure for each protein should finally be extracted, removing ions and solvent molecules, and saved as a _pdb_ file in the _reference_pdbs_ directory, accessible to the 01_cmap_preparation.py script.

A script ([00_snapshot_extraction.sh](Data/00_snapshot_extraction.sh)) to extract the reference structure from simulations performed in [AMBER](https://ambermd.org/) is provided in the [Data](Data) directory. Some trivial modifications are needed to adapt the script to the user's files.


## 2. Unique contacts identification and contact map creation
The contact map distances will be used as collective variables (CVs) for the construction of the free energy surface. Each CV will be the cotact map distance from the reference structure to each of the others structures (in this example, we'll have 3 structures, called ref1, ref2 and ref3).


To determine the unique contacts in each structure, the distance between each couple of atoms in each structure has to be computed. If such distance is less than a defined cutoff in a system, while greater than the cutoff in the others, the contact is considered unique. Only non-adjacent residues should be considered.
Once the unique contacts, that caracterize a reference structure, are identified, the contact map can be calculated. The contact map for the $i^{th}$ reference structure can be defined as:
$$CM_i=\frac{1}{N}\sum_{\gamma}^N [D_\gamma(R) - D_\gamma(R_i)]^2$$
where $N$ is a normalization constant and $D_\gamma(R)$ is a sigmoidal function defining the formation of the contact, as
$$D_\gamma(R)=\frac{1-(\frac{r_\gamma}{r_{\gamma,0}})^6}{1-(\frac{r_\gamma}{r_{\gamma,0}})^10}$$
, where $r_{\gamma,0}$ is the contact distance in the reference structure in which the contact is unique.

A python script that create such contact maps for a list of reference structures, can be found in the [Data](Data) directory ([01_cmap_preparation.py](Data/01_cmap_preparation.py)). The script is also able to identify different types of interactions between atoms, with different cutoff values for each interaction. Different interaction types are also weighted differently in the creation of the CV. 
Such cutoffs and weights are hard-coded inside the script, but can be easily changed or setted to the same value.
The script have as a final output one file for each reference structure, called _cmap_SOMETHING.dat_, where _SOMETHING_ is the name of the reference structure. These files can be loaded in a plumed metadynamic scipt to define the proper contact map. 


## 3. Use of steered MD to check the effectiveness of the CVs
To check if the created contact map CVs are effectivelly able to drive the protein to explore its different conformations, it can be useful to run steered molecular dynamic (SMD) simulations starting from the reference structure and forcing the structure to morph to each of the other structures, using the appropriate CV. This can be done with the input plumed file [02_smd_to_ref2.dat](Data/02_smd_to_ref2.dat), supposing we are using _ref1_ as the starting structure and pushing the conformation to reach _ref2_.
The _AT0_ value is setted as 1.0 in the example, but should be replaced with the correct initial value, that can be identified with a short (1-10 steps) simulation. Also the force of the pushing can be adjusted as needed.

To check that the SMD is effectivelly pushing the structure to the corresponding reference, a _rmsd_ calculation can be performed, using the target structure as the reference structure. The _rmsd_ value should decrease during the simulation, eventually approacing 0.


## 4. Run metadyniamc simulations with the selected CVs
The script [03_run_metad.dat](Data/03_run_metad.dat) can be used as a template to run the metadyniamc simulations to collect the data nedded to reconstruct the free energy surface. The script should be used in combination with the replica exchange capability of GROMACS, and it's tailored to perform on an HPC cluster with at least 32 GPUS available. As it is, 32 temperature replica exchange simulations should ran in parallel, each starting from one of the reference structures. The relevant data, needed to reconstruct the free energy surface, will be extracted from the replica ran at 300.0 K.
#!/usr/bin/env python3

from contextlib import suppress
import pickle
import lzma
from os.path import isfile, expanduser
from os import getcwd, chdir
from tqdm import tqdm
from numpy import array, zeros, sqrt as np_sqrt

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, destination, selfreturn=False):
        self.destination = expanduser(destination)
        self.selfreturn = selfreturn

    def __enter__(self):
        self.savedPath = getcwd()
        chdir(self.destination)
        if self.selfreturn:
            return self

    def __exit__(self, etype, value, traceback):
        chdir(self.savedPath)


class Atom:
    def __init__(self, res_name='', atom_name='', element='',
                 res_num=1, atom_num=1, x=0.0, y=0.0, z=0.0, kind='', chain_id='',
                 pdb_line=''):
        self.kind = kind
        self.atom_num = atom_num
        self.atom_name = atom_name
        self.altloc = ''
        self.res_name = res_name
        self.chainID = chain_id
        self.res_num = res_num
        self.x = x
        self.y = y
        self.z = z
        self.xyz_array = array([x, y, z])
        self.occupancy = 0.0
        self.tempfactor = 0.0
        self.element = element
        if pdb_line:
            self.pdb_line_parser(pdb_line=pdb_line)
        self.characterize()

    def __str__(self):
        return f"{self.res_name}.{self.res_num}_{self.atom_name}.{self.atom_num}"

    def _is_valid_operand(self, other):
        return hasattr(other, "atom_num") & hasattr(other, "atom_name")

    def __eq__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.atom_num == other.atom_num & self.atom_name == other.atom_name

    def __ne__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return False if self.__eq__(other) else True

    def __lt__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.atom_num < other.atom_num

    def __le__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.atom_num <= other.atom_num

    def __gt__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.atom_num > other.atom_num

    def __ge__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.atom_num >= other.atom_num

    def characterize(self):
        self.is_polar = False
        self.is_aromatic = False
        self.is_backbone = False
        self.formal_charge = None
        if self.atom_name.startswith("O") or self.atom_name.startswith("N") or self.atom_name.startswith("S"):
            self.is_polar = True
        if self.res_name in _negative_residues and self.atom_name in _negative_atoms[self.res_name]:
            self.formal_charge = -1
        if self.res_name in _positive_residues and self.atom_name in _positive_atoms:
            self.formal_charge = 1
        if self.atom_name in _backbone_atoms:
            self.is_backbone = True
        if self.res_name in _aromatic_residues and self.atom_name in _aromatic_atoms:
            self.is_aromatic = True

    def pdb_line_parser(self, pdb_line, prev_at_idx=0):
        self.kind = pdb_line[0:6].strip()
        try:
            self.atom_num = int(pdb_line[6:11].strip())
        except ValueError:
            if pdb_line[6:11] == "*****" and prev_at_idx:
                self.atom_num = prev_at_idx + 1
            else:
                self.atom_num = pdb_line[6:11].strip()
        self.atom_name = pdb_line[12:16].strip()
        self.altloc = pdb_line[16].strip()
        self.res_name = pdb_line[17:20].strip()
        self.chainID = pdb_line[21].strip()
        self.res_num = int(pdb_line[22:28].strip())  # 22:26 standard...
        self.x = float(pdb_line[30:38].strip())
        self.y = float(pdb_line[38:46].strip())
        self.z = float(pdb_line[46:54].strip())
        self.xyz_array = array([self.x, self.y, self.z])
        with suppress(ValueError):
            self.occupancy = float(pdb_line[54:60].strip())
        with suppress(ValueError):
            self.tempfactor = float(pdb_line[60:66].strip())
        with suppress(ValueError):
            self.element = pdb_line[76:78].strip()
        with suppress(ValueError):
            self.charge = float(pdb_line[78:80].strip())
        if prev_at_idx:
            return self.atom_num

    def pdb_line_writer(self, out_file=None):
        pdb_line = f'{self.kind:<6s}'
        try:
            if float(self.atom_num) > 99999 or self.atom_num == "*****":
                pdb_line += f'***** {self.atom_name:>4s}'
            else:
                pdb_line += f'{self.atom_num:5d} {self.atom_name:<4s}'
        except ValueError:
            pdb_line += f'***** {self.atom_name:>4s}'
        pdb_line += f'{self.altloc:1s}' if self.altloc else ' '
        pdb_line += f'{self.res_name:<4s}'
        pdb_line += f'{self.chainID:1s}' if self.chainID else ' '

        def format_coordinate_string(coord_value):
            coord = f'{coord_value:<7.4f}'.rstrip()
            if len(coord) >= 6 and coord[-1] == '0':
                coord = coord[:-1]
            if len(coord) >= 7 and "." in coord:
                coord = f'{float(coord):<6.3f}'
            return coord
        x = format_coordinate_string(self.x)
        y = format_coordinate_string(self.y)
        z = format_coordinate_string(self.z)
        if self.res_num <= 9999:
            pdb_line += f'{self.res_num:4d}    {x:>8s}{y:>8s}{z:>8s}'
        elif self.res_num <= 99999:
            pdb_line += f'{self.res_num:5d}   {x:>8s}{y:>8s}{z:>8s}'
        elif self.res_num <= 999999:
            pdb_line += f'{self.res_num:6d}  {x:>8s}{y:>8s}{z:>8s}'
        elif self.res_num <= 9999999:
            pdb_line += f'{self.res_num:7d} {x:>8s}{y:>8s}{z:>8s}'

        pdb_line += f'{self.occupancy:6.2f}' if self.occupancy else '      '
        pdb_line += f'{self.tempfactor:6.2f}' if self.tempfactor else '      '
        pdb_line += f'{self.element:2s}' if self.element else '  '
        pdb_line += f'{self.charge:8.5f}' if self.charge else '        '
        pdb_line += '\n'
        if not out_file:
            return pdb_line
        out_file.write(pdb_line)

    def calc_distances(self, ref_coords):
        ref_coords = array(ref_coords)
        mycoord = zeros(ref_coords.shape)
        mycoord[:, :] = self.xyz_array
        distances = np_sqrt((mycoord[:, 0] - ref_coords[:, 0]) ** 2 +
                            (mycoord[:, 1] - ref_coords[:, 1]) ** 2 +
                            (mycoord[:, 2] - ref_coords[:, 2]) ** 2)
        self.distances = distances
        self.sorted_dist_idx = array(
            [idx for idx, _d in sorted(enumerate(self.distances), key=lambda k: k[1])])

class Molecule:
    def __init__(self, molecule_name="molecule"):
        self.residues = []
        self.mol_name = molecule_name
        self.previous_residue_num = None

    def __str__(self):
        return self.mol_name

    def __len__(self):
        return len(self.residues)

    def add_residue(self, res_to_add=None):
        self.residues.append(res_to_add if isinstance(res_to_add, Residue) else Residue())

    def add_atom(self, atom_to_add):
        if atom_to_add.res_num != self.previous_residue_num:
            self.previous_residue_num = atom_to_add.res_num
            self.add_residue()
        self.residues[-1].add_atom(atom_to_add=atom_to_add)

    def write_atoms_to_pdb(self, pdb_file, add_ter_after_residues=False, add_ter_after_molecules=True):
        for self.res in self.residues:
            self.res.write_atoms_to_pdb(pdb_file=pdb_file)
            if add_ter_after_residues:
                pdb_file.write("TER\n")
        if add_ter_after_molecules and self.residues:
            pdb_file.write("TER\n")


class Residue:
    def __init__(self, res_name='', res_num=0):
        self.atoms = []
        self.res_name = res_name
        self.res_num = res_num

    def add_atom(self, atom_to_add):
        if not self.res_name:
            self.res_name = atom_to_add.res_name
        if not self.res_num:
            self.res_num = atom_to_add.res_num
        self.atoms.append(atom_to_add)

    @property
    def number_of_atoms(self):
        return len(self.atoms)

    def write_atoms_to_pdb(self, pdb_file):
        for self.atom in self.atoms:
            self.atom.pdb_line_writer(out_file=pdb_file)


class PdbFile:
    def __init__(self, file_path):
        self.comments = []
        self.other_lines = []
        self.molecules = []
        self.reading_new_molecule = True
        self.cryst_dimensions, self.cryst_line = [], ""
        self.file_path = file_path
        if self.file_path: self.read_pdb_file()

    def read_pdb_file(self):
        for line in open(self.file_path, 'r'):
            if line.startswith("CRYST1"): self.cryst_line = line
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                if self.reading_new_molecule:
                    self.reading_new_molecule = False
                    self.add_molecule()
                self.molecules[-1].add_atom(atom_to_add=Atom(pdb_line=line))
            elif line.startswith("TER"):
                self.reading_new_molecule = True

    def add_molecule(self):
        self.molecules.append(Molecule())

    def write_pdb(self, out_file_path, add_ter_after_residues=False, add_ter_after_molecules=True):
        with open(out_file_path, 'w') as self.out_file:
            if self.cryst_dimensions:
                self.out_file.write(
                    'CRYST1{x[0]:>9.3f}{x[1]:>9.3f}{x[2]:>9.3f}{a:>7.2f}{a:>7.2f}{a:>7.2f}\n'.format(
                    x=self.cryst_dimensions, a=90.0))
            elif self.cryst_line:
                self.out_file.write(self.cryst_line)
            for molecule in self.molecules:
                molecule.write_atoms_to_pdb(
                    pdb_file=self.out_file, add_ter_after_residues=add_ter_after_residues,
                    add_ter_after_molecules=add_ter_after_molecules)
            self.out_file.write("END\n")


class Interaction:
    def __init__(self, obj1, obj2, distance, tag=None, kind=None, comment=None, idx=None):
        self.obj1 = obj1
        self.obj2 = obj2
        self.distance = distance
        self.kind = kind
        self.comment = comment
        self.idx = idx

        if isinstance(self.obj1, Atom) and isinstance(self.obj2, Atom):
            self.obj_min = min([obj1, obj2])
            self.obj_max = max([obj1, obj2])
        elif isinstance(self.obj1, list) and isinstance(self.obj2, list):
            self.obj_min = min([obj1[0], obj2[0]])
            self.obj_max = max([obj1[0], obj2[0]])
        else:
            print("> WARNING: obj1 and obj2 are not both of type 'Atom' or 'list', I'm confused, but I try to proceed...")
            self.obj_min = min([obj1, obj2])
            self.obj_max = max([obj1, obj2])

        self.unique_tag = tag if tag is not None else f"{self.obj_min} to {self.obj_max}"
        if isinstance(self.obj_min, Atom) and isinstance(self.obj_max, Atom):
            self.res_involved = f"{self.obj_min.res_name}_{self.obj_min.res_num} to "
            self.res_involved += f"{self.obj_max.res_name}_{self.obj_max.res_num}"

    def characterize(
            self, cutoff_standard=3.5, cutoff_SB=5.5, cutoff_pi_charged=3.5, cutoff_BB=3.75, 
            cutoff_polar=4.0, cutoff_PiPi=3.75):
        if isinstance(self.obj1, Atom) and isinstance(self.obj2, Atom):
            if self.kind is None:
                self.kind = "standard"
            else:
                raise Exception(
                    f"trying to characterize an interaction which has already a kind '{self.kind}'")
            self.cutoff = cutoff_standard
            if self.obj1.formal_charge is not None:
                if self.obj2.formal_charge is not None:
                    if self.obj1.formal_charge * self.obj2.formal_charge < 0:
                        self.kind = "salt_bridge"
                        self.cutoff = cutoff_SB
                    else:
                        self.kind = "salt_repulsion"
                        self.cutoff = 0.0
                elif self.obj2.is_aromatic and self.obj1.formal_charge > 0:
                    self.kind = "pi_charged"
                    self.cutoff = cutoff_pi_charged
                elif self.obj2.is_polar:
                    self.kind = "polar"
                    self.cutoff = cutoff_polar
            elif self.obj1.is_backbone and self.obj2.is_backbone:
                self.kind = "BB"
                self.cutoff = cutoff_BB
            elif self.obj1.is_polar and self.obj2.is_polar:
                self.kind = "polar"
                self.cutoff = cutoff_polar
            elif self.obj1.is_aromatic and self.obj2.is_aromatic:
                self.kind = "pi_pi"
                self.cutoff = cutoff_PiPi
            elif self.obj1.is_aromatic and self.obj2.formal_charge is not None:
                if self.obj2.formal_charge > 0:
                    self.kind = "pi_charged"
                    self.cutoff = cutoff_pi_charged

    def __str__(self):
        repr_str = f"interaction {self.unique_tag} has distance {self.distance:.2f}"
        if self.kind: repr_str += f", is of kind {self.kind}"
        return repr_str

    @staticmethod
    def _is_valid_operand(other):
        return hasattr(other, "obj1") and hasattr(other, "obj2") and hasattr(other, "distance")

    def __eq__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.unique_tag == other.unique_tag

    def __ne__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return False if self.__eq__(other) else True

    def __lt__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.distance < other.distance

    def __le__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.distance <= other.distance

    def __gt__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.distance > other.distance

    def __ge__(self, other):
        if not self._is_valid_operand(other): return NotImplemented
        return self.distance >= other.distance


class ReferenceStructure:
    def __init__(self, name):
        self.name = name
        self.input_pdb = f"./reference_pdbs/{name}.pdb"
        self.atoms_to_skip = "H"  # can be "H" or "H,C" for example

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def load_atoms(self):
        print(f">>>> Working on: {self.name}")
        picke_filename = f"list_of_atoms_for_{self.name}_{self.atoms_to_skip}.pickle.xz"
        if isfile(picke_filename):
            print(">> Loading precomputed list of atoms...")
            with lzma.open(picke_filename, 'rb') as f:
                self.atoms = pickle.load(f)
        else:
            self.atoms = self._collect_atoms()
            self._identify_the_contacts()
            with lzma.open(picke_filename, 'wb') as f:
                pickle.dump(self.atoms, f)

    def _collect_atoms(self):
        print(f">> Loading coordinates: {self.input_pdb}")
        atoms = []
        for res in PdbFile(file_path=self.input_pdb).molecules[0].residues:
            for atom in res.atoms:
                if not any(atom.atom_name.startswith(x) for x in self.atoms_to_skip.split(",")):
                    atoms.append(atom)
        return atoms

    def _identify_the_contacts(self):
        print(">> Identifying the contacts")
        all_atoms_coords = array([atom.xyz_array for atom in self.atoms])
        for atom in tqdm(self.atoms):
            atom.calc_distances(ref_coords=all_atoms_coords)
            atom.contacts = []
            for idx, distance_value in enumerate(atom.distances):
                other_atom = self.atoms[idx]
                if other_atom.res_num - 1 <= atom.res_num <= other_atom.res_num + 1:
                    continue
                d = Interaction(obj1=atom, obj2=other_atom, distance=distance_value)
                d.characterize()
                if d.distance <= d.cutoff:
                    atom.contacts.append(d)

    def select_unique_contacts(self, list_of_others):
        others_names = ' | '.join([other.name for other in list_of_others])
        print(f">> Selecting unique contacts of {self.name} with respect to ({others_names})")
        self.unique_contacts = []
        for idx, atom in enumerate(tqdm(self.atoms)):
            for contact in atom.contacts:
                if contact not in self.unique_contacts:
                    if not any(contact in other.atoms[idx].contacts for other in list_of_others):
                        self.unique_contacts.append(contact)

    def _summarize_contacts_by_res(self, verbose=False, keep_all=True):
        dict_of_contacts_by_res_involved = {}
        for contact in self.unique_contacts:
            if contact.res_involved not in dict_of_contacts_by_res_involved:
                dict_of_contacts_by_res_involved[contact.res_involved] = []
            dict_of_contacts_by_res_involved[contact.res_involved].append(contact)

        characterized_dict = {}
        for res_involved, contacts in dict_of_contacts_by_res_involved.items():
            for c in contacts:
                if res_involved not in characterized_dict:
                    characterized_dict[res_involved] = {"base": [], "advanced": []}
                if c.kind in ["standard", "BB"]:
                    characterized_dict[res_involved]["base"].append(c)
                elif c.kind == "salt_bridge":
                    characterized_dict[res_involved]["advanced"].append(c)
                elif c.kind == "polar":
                    characterized_dict[res_involved]["advanced"].append(c)
                elif c.kind == "pi_charged":
                    characterized_dict[res_involved]["advanced"].append(c)
                elif c.kind == "pi_pi":
                    characterized_dict[res_involved]["advanced"].append(c)
                else:
                    raise ValueError(f"contact {c} doesn't have a valid 'kind', it's '{c.kind}'")

        purged_duplicates_dict = {}
        for res_involved, contacts in characterized_dict.items():
            if res_involved not in purged_duplicates_dict:
                purged_duplicates_dict[res_involved] = {"base": [], "advanced": []}
            if keep_all:
                for c in contacts["advanced"]:
                    purged_duplicates_dict[res_involved]["advanced"].append(c)
                for c in contacts["base"]:
                    purged_duplicates_dict[res_involved]["base"].append(c)
            else:
                if contacts["advanced"]:
                    if len(contacts['base']) > 0 and verbose:
                        print(f"Skipped {len(contacts['base'])} contacts for {res_involved}")
                    for c in contacts["advanced"]:
                        purged_duplicates_dict[res_involved]["advanced"].append(c)
                else:
                    if verbose: print(f"Only {len(contacts['base'])} base contacts for {res_involved}")
                    for c in contacts["base"]:
                        purged_duplicates_dict[res_involved]["base"].append(c)
        return purged_duplicates_dict

    def write_reference_cmap_plumed_creation(self, f, start_idx):
        f.write(f"# Unique contacts of {self.name}\n")
        for idx, contact in enumerate(self.unique_contacts, start=start_idx):
            start_idx += 1
            f.write(f"  ATOMS{idx}={contact.obj1.atom_num},{contact.obj2.atom_num}")
            f.write(f"    SWITCH{idx}={{RATIONAL R_0={contact.distance*100/60:.4f} D_0=0.0 NN=6 MM=10}}")
            f.write(f"    # {contact}")
            f.write("\n")
        return start_idx


_negative_residues = ["ASP", "GLU"]
_positive_residues = ["LYS", "ARG", "HIP"]
_negative_atoms = {
    "ASP": ["OD1", "OD2", "CG"],
    "GLU": ["OE1", "OE2", "CD"]}
_positive_atoms = ["NZ", "NE", "NH1", "NH2", "CZ", "NE1", "ND1"]
_backbone_atoms = ["N", "C", "O", "CA"]
_aromatic_residues = ["PHE", "TYR", "TRP"]
_aromatic_atoms = ["CZ", "CG", "CE1", "CE2", "CD1", "CD2", "CZ2", "CZ3", "CE3", "CH2"]

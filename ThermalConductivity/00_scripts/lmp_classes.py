#!/usr/bin/env python3

from pandas import DataFrame
from typing import IO


class PushBackIterator():
    def __init__(self, iterator):
        self.iterator = iter(iterator)
        self.pushed_back = []
        self.ok_to_push_back = True

    def __iter__(self):
        return self

    def __next__(self):
        if self.pushed_back:
            return self.pushed_back.pop()
        try:
            return next(self.iterator)
        except StopIteration:
            self.ok_to_push_back = False
            raise

    def push_back(self, element):
        if not self.ok_to_push_back:
            raise StopIteration
        self.pushed_back.append(element)


def type_aware_appender(line_to_parse, list_to_append_to):
    for idx, item in enumerate(line_to_parse.split()):
        if '.' in item:
            list_to_append_to[idx][1].append(float(item))
        else:
            try:
                list_to_append_to[idx][1].append(int(item))
            except ValueError:
                try:
                    list_to_append_to[idx][1].append(float(item))
                except ValueError as exc:
                    raise ValueError(
                        f"Problem parsing '{item}' as float/int from line:\n{line_to_parse}\n"
                        ) from exc


class LmpLogFile():
    def __init__(self, log_file):
        self.log_file = log_file
        self.df = self.parse_log_file()

    def parse_log_file(self):
        with open(self.log_file) as f:
            line = f.readline()
            data = []
            while line:
                if line.strip().startswith("Step"):
                    for item in line.split():
                        data.append((item.strip(), []))
                    line = f.readline()
                    while not line.strip().startswith("Loop"):
                        if len(line.split()) != len(data):
                            line = f.readline()
                            break
                        type_aware_appender(line_to_parse=line, list_to_append_to=data)
                        line = f.readline()
                else:
                    line = f.readline()
        return DataFrame.from_dict(dict(data))

    @property
    def n_frames(self):
        return len(self.df)


class LmpDumpFile():
    def __init__(self, dump_file):
        self.dump_file = dump_file
        self.frames = []
        self.parse_dump_file()

    def parse_dump_file(self):
        with open(self.dump_file) as f:
            line = f.readline()
            while line:
                if line.strip().startswith("ITEM: TIMESTEP"):
                    timestep = int(f.readline())
                    line = f.readline()
                elif line.strip().startswith("ITEM: NUMBER OF ATOMS"):
                    num_atoms = int(f.readline())
                    line = f.readline()
                elif line.strip().startswith("ITEM: BOX BOUNDS"):
                    xlo, xhi = map(float, f.readline().split())
                    ylo, yhi = map(float, f.readline().split())
                    zlo, zhi = map(float, f.readline().split())
                    line = f.readline()
                elif line.strip().startswith("ITEM: ATOMS"):
                    current_frame = []
                    for item in line.split()[2:]:
                        current_frame.append((item.strip(), []))
                    line = f.readline()
                    while line and not line.strip().startswith("ITEM: TIMESTEP"):
                        type_aware_appender(line_to_parse=line, list_to_append_to=current_frame)
                        line = f.readline()
                    df = DataFrame.from_dict(dict(current_frame))
                    df.attrs['timestep'] = timestep
                    df.attrs['num_atoms'] = num_atoms
                    df.attrs['box_dims'] = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
                    self.frames.append(df.sort_values(by=['id']).reset_index(drop=True))
                    # timestep = int(f.readline())
                else:
                    line = f.readline()

    @property
    def n_frames(self) -> int:
        return len(self.frames)

def extract_data(l_iter: PushBackIterator, defined_object=None, object_style: str | None = None) -> list:
    parsed_section = []
    items = next(l_iter).split()
    while not items or items[0].isdigit():
        if items:
            if defined_object is None:
                parsed_section.append([])
                for item_idx, item in enumerate(items):
                    try:
                        parsed_section[-1].append(int(item))
                    except ValueError:
                        try:
                            parsed_section[-1].append(float(item))
                        except ValueError:
                            if item.startswith("#"):
                                parsed_section[-1].append(" ".join(items[item_idx:]))
                                break
            else:
                parsed_section.append(defined_object(items, style=object_style))
        try:
            line = next(l_iter)
        except StopIteration:
            return parsed_section
        items = line.split()
    l_iter.push_back(line)
    return parsed_section


class LmpBond():
    def __init__(self, items: list, style: str | None = None) -> None:
        self.id = int(items[0])
        self.type = int(items[1])
        self.at1 = int(items[2])
        self.at2 = int(items[3])

    def write_to(self, out_f: IO) -> None:
        out_f.write(f'{self.id:10d}')
        out_f.write(f' {self.type:5}')
        out_f.write(f' {self.at1:10d} {self.at2:10d}')
        out_f.write('\n')

class LmpAngle():
    def __init__(self, items: list, style: str | None = None) -> None:
        self.id = int(items[0])
        self.type = int(items[1])
        self.at1 = int(items[2])
        self.at2 = int(items[3])
        self.at3 = int(items[4])

    def write_to(self, out_f: IO) -> None:
        out_f.write(f'{self.id:10d}')
        out_f.write(f' {self.type:5}')
        out_f.write(f' {self.at1:10d} {self.at2:10d} {self.at3:10d}')
        out_f.write('\n')

class LmpTorsion():
    def __init__(self, items: list, style: str | None = None) -> None:
        self.id = int(items[0])
        self.type = int(items[1])
        self.at1 = int(items[2])
        self.at2 = int(items[3])
        self.at3 = int(items[4])
        self.at4 = int(items[5])

    def write_to(self, out_f: IO) -> None:
        out_f.write(f'{self.id:10d}')
        out_f.write(f' {self.type:5}')
        out_f.write(f' {self.at1:10d} {self.at2:10d} {self.at3:10d} {self.at4:10d}')
        out_f.write('\n')

class LmpAtom:
    def __init__(self, items: list, style: str | None = "full") -> None:
        self.id = int(items[0])
        self.style = style
        self.nx = self.ny = self.nz = 0

        if self.style == "full":
            self.molecule_id = int(items[1])
            self.type = int(items[2])
            self.charge = float(items[3])
            self.x = float(items[4])
            self.y = float(items[5])
            self.z = float(items[6])
            self.comment = " ".join(items[7:])

    def write_to(self, out_f: IO) -> None:
        out_f.write(f'{self.id:10d}')
        if self.style == "full":
            out_f.write(f' {self.molecule_id:10d} {self.type:5d} {self.charge:10.5f}')
            out_f.write(f' {self.x:10.5f} {self.y:10.5f} {self.z:10.5f}')
            out_f.write(f' {self.nx} {self.ny} {self.nz}')

        out_f.write('\n')

    def __repr__(self) -> str:
        return f"atom {self.id}, type {self.type} at ({self.x:.3f}, {self.y:.3f}, {self.z:.3f}) [{self.comment}]"


class LmpDataFile():
    def __init__(self, data_file: str, adjust_molecules=False) -> None:
        self.data_file = data_file
        self.read_data_file()
        self.n_mols = 1
        self.n_atm_per_mol = self.adjust_molecules_definition() if adjust_molecules else self.n_atms

    def __repr__(self) -> str:
        return(f"{self.data_file}, with {self.n_mols} molecules and {self.n_atm_per_mol} atoms per molecule")
        
    def read_data_file(self) -> None:
        with open(self.data_file) as f:
            lines = f.readlines()
        self.header = lines[0]
        self.n_atms = int(lines[2].split()[0])
        self.n_bnds = int(lines[3].split()[0])
        self.n_angs = int(lines[4].split()[0])
        self.n_dihs = int(lines[5].split()[0])
        self.n_imps = int(lines[6].split()[0])

        iterator = PushBackIterator(lines[7:])
        for l in iterator:
            items = l.split()
            if not items:
                continue
            if "crossterms" in items:
                self.n_crs = int(items[0])
            elif "types" in items:
                if "atom" in items:
                    self.n_atm_types = int(items[0])
                elif "bond" in items:
                    self.n_bnd_types = int(items[0])
                elif "angle" in items:
                    self.n_ang_types = int(items[0])
                elif "dihedral" in items:
                    self.n_dih_types = int(items[0])
                elif "improper" in items:
                    self.n_imp_types = int(items[0])
            elif "xlo" in items:
                self.xlo = float(items[0])
                self.xhi = float(items[1])
            elif "ylo" in items:
                self.ylo = float(items[0])
                self.yhi = float(items[1])
            elif "zlo" in items:
                self.zlo = float(items[0])
                self.zhi = float(items[1])
                break

        for l in iterator:
            if l.startswith("Atoms"):
                self.atoms = extract_data(iterator, defined_object=LmpAtom, object_style="full")
            if l.startswith("Bonds"):
                self.bonds = extract_data(iterator, defined_object=LmpBond)
            if l.startswith("Angles"):
                self.angles = extract_data(iterator, defined_object=LmpAngle)
            if l.startswith("Dihedrals"):
                self.dihedrals = extract_data(iterator, defined_object=LmpTorsion)
            if l.startswith("Pair Coeffs"):
                self.pair_coeffs = extract_data(iterator)
            if l.startswith("Bond Coeffs"):
                self.bond_coeffs = extract_data(iterator)
            if l.startswith("Angle Coeffs"):
                self.angle_coeffs = extract_data(iterator)
            if l.startswith("Dihedral Coeffs"):
                self.dihedral_coeffs = extract_data(iterator)

    def adjust_molecules_definition(self, starting_mol_id: int = 1) -> int:
        real_mol_id = starting_mol_id
        previous_mol_id = n_atm_per_mol = 0
        for atom in self.atoms:
            if atom.molecule_id < previous_mol_id:
                real_mol_id += 1
                n_atm_per_mol = 0
            n_atm_per_mol += 1
            previous_mol_id = atom.molecule_id
            atom.molecule_id = real_mol_id
        self.n_mols = int(self.n_atms/n_atm_per_mol)
        modulus =  self.n_atms % n_atm_per_mol
        if modulus != 0:
            raise ValueError(f"number of atoms ({self.n_atms}) % atoms per chains ({n_atm_per_mol} shoud be 0 (is {modulus})")
        return n_atm_per_mol


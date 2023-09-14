#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from lmp_classes import LmpDumpFile, LmpLogFile

ANG3_TO_CM3 = (1E-8)**3
N_AVOGADRO = 6.0221 * 10**23

GRAPHENE_ATOM_TYPE = 3
GRAPHENE_SHEET_LAST = [10552, 11432, 12312] # A list of the last graphene atom index of each graphene sheet

DELTA_Z = 0.2
OUTPUT_BASENAME = "../02_equilibrate/Density_profile"
input_log_file = '../02_equilibrate/s03_NPTeq2.lmplog'
input_dump_file = '../02_equilibrate/s03_NPTeq2.lammpstrj'


lammps_log = LmpLogFile(log_file=input_log_file)
lammps_dump = LmpDumpFile(dump_file=input_dump_file)

# we logged l_x/y/z at moments of dump, so now we unwarp the coordinates
def unwarp(coord: float, image_flag: int, box_length: float) -> float:
    return coord + (image_flag * box_length)

log_to_dump_frame_ratio = int((lammps_log.n_frames-1)/(lammps_dump.n_frames-1))
for i, frame in enumerate(lammps_dump.frames):
    frame.attrs["lx"] = lammps_log.df.Lx[i*log_to_dump_frame_ratio]
    frame.attrs["ly"] = lammps_log.df.Ly[i*log_to_dump_frame_ratio]
    frame.attrs["lz"] = lammps_log.df.Lz[i*log_to_dump_frame_ratio]
    frame["xu"] = unwarp(coord=frame.x, image_flag=frame.ix, box_length=frame.attrs["lx"])
    frame["yu"] = unwarp(coord=frame.y, image_flag=frame.iy, box_length=frame.attrs["ly"])
    frame["zu"] = unwarp(coord=frame.z, image_flag=frame.iz, box_length=frame.attrs["lz"])


########## Determine average z coordinate of all graphene atoms (type 3) for all sheets
for frame_idx, frame in enumerate(lammps_dump.frames, start=1):
    graphene_sheet_zs = [[] for _ in GRAPHENE_SHEET_LAST]
    for atom_idx, atom_type, atom_z in zip(frame.id, frame.type, frame.zu):
        if atom_type == GRAPHENE_ATOM_TYPE:
            for sheet_idx, cut_idx in enumerate(GRAPHENE_SHEET_LAST):
                if atom_idx <= cut_idx:
                    graphene_sheet_zs[sheet_idx].append(atom_z)
                    break

    graphene_sheet_zs = np.array(graphene_sheet_zs)

    print(f'Frame {frame_idx:10d}:', end="")
    frame.attrs["avg_z"] = []
    for idx, (avg, std) in enumerate(zip(graphene_sheet_zs.mean(axis=1),
                                         graphene_sheet_zs.std(axis=1)), start=1):
        frame.attrs["avg_z"].append(avg)
        print(f' z{idx} {avg:.3f} ({std:.3f}) |', end="")
    print("")


########## compute the density (C and H atom types are hard coded here for simplicity)
length_z = lammps_dump.frames[-1].attrs["lz"]/2
z_steps = int(length_z/DELTA_Z)

plt_z = np.zeros((lammps_dump.n_frames, z_steps))
N_layer = np.zeros((lammps_dump.n_frames, z_steps))
M_layer = np.zeros((lammps_dump.n_frames, z_steps))
dens_layer = np.zeros((lammps_dump.n_frames, z_steps))

for frame_idx, frame in enumerate(lammps_dump.frames):
    dist_zp = max(frame.attrs["avg_z"])
    dist_zm = min(frame.attrs["avg_z"])
    for z in range(z_steps):
        dist_zp += DELTA_Z
        dist_zm -= DELTA_Z
        plt_z[frame_idx, z] = z*DELTA_Z
        for atom_z, atom_type in zip(frame.zu, frame.type):
            if atom_type != GRAPHENE_ATOM_TYPE:
                if (dist_zp - DELTA_Z < atom_z < dist_zp) or (dist_zm + DELTA_Z > atom_z > dist_zm):
                    N_layer[frame_idx, z] += 1
                    if atom_type == 1:
                        M_layer[frame_idx, z] += 12.01
                    elif atom_type == 2:
                        M_layer[frame_idx, z] += 1.008

for frame_idx, frame in enumerate(lammps_dump.frames):
    volume = frame.attrs['lx']*frame.attrs['ly']*DELTA_Z*ANG3_TO_CM3
    # mass is in grams/mole for units metal, hence the avogadro number
    dens_layer[frame_idx, :] = M_layer[frame_idx, :] / volume / N_AVOGADRO

yfin = dens_layer.mean(axis=0)

with open(f'{OUTPUT_BASENAME}.txt', 'w') as f_out:
    f_out.write('#z coordinates [ang]  |  density [g/cm^3]\n')
    for x, y in zip(plt_z[0, :], yfin):
        f_out.write(f'{x:.10e}\t\t{y:.10e}\n')

plt.plot(plt_z[0, :], yfin)
plt.xlabel('distance along z $\\left[\\AA\\right]$', fontsize = 8)
plt.ylabel('$\\rho$ $\\left[\\frac{g}{cm^3}\\right]$', fontsize = 8)
plt.tight_layout()
plt.savefig(f'{OUTPUT_BASENAME}.pdf')

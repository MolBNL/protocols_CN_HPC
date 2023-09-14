#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from MyLammps.lmp_classes import type_aware_appender
import numpy as np
from scipy.stats import linregress
from lmp_classes import LmpLogFile

sns.set()
sns.set_context("paper", font_scale=0.5)

################################################################
# units from formula: eV/(psec*ang*ang) / (K/ang)     -->    eV / (psec*ang*K)
# W = J/s
# kappa units: W / (m*K)        --> J / (s*m*K)

eV_to_J = 1.60218e-19
psec_to_s = 1e-12
ang_to_m = 1e-10

converter = eV_to_J / psec_to_s / ang_to_m
################################################################

class TprofileFile:
    def __init__(self, profile_file):
        self.profile_file = profile_file
        self.df = self.load_profile_file()
        self.slopes1, self.slopes2, self.avg, self.coord = self.compute_slopes()

        self.avg_slope1 = np.mean(self.slopes1)
        self.std_slope1 = np.std(self.slopes1)       
        self.avg_slope2 = np.mean(self.slopes2)
        self.std_slope2 = np.std(self.slopes2)

    def load_profile_file(self):
        with open(self.profile_file) as f:
            data = [("timestep", [])]
            ts = 0
            for line in f:
                if line.startswith("# Chunk "):
                    for item in line.strip('#').split():
                        data.append((item.strip(), []))
                    continue
                if line.startswith("#") or line.startswith("\n"):
                    continue
                if line.startswith(" "):
                    line_with_timestep = f"{ts} {line}"
                    type_aware_appender(line_to_parse=line_with_timestep, list_to_append_to=data)
                else:
                    ts = int(line.split(" ")[0])
        return pd.DataFrame.from_dict(dict(data))

    def compute_slopes(self):
        by_ts = self.df.groupby(["timestep"])
        timesteps = list(by_ts.groups.keys())
        coord = by_ts.get_group(timesteps[0])["Coord1"]
        middle_chunks = int(len(coord)/2)+1
        def get_average_temps(beg, end):
            avg = np.zeros(by_ts.get_group(timesteps[0])["temp"].shape)
            for timestep in timesteps[beg:end]:
                avg += by_ts.get_group(timestep)["temp"].to_numpy()
            return avg / (end-beg)
        slopes1 = []
        slopes2 = []
        beg = 10
        interval = 100
        end = beg + interval
        while end <= len(timesteps):
            avg = get_average_temps(beg, end)
            slope1, _intercept1, _r_value1, _p_value1, _std_err1 = linregress(coord[:middle_chunks], avg[:middle_chunks])
            slope2, _intercept2, _r_value2, _p_value2, _std_err2 = linregress(coord[middle_chunks:], avg[middle_chunks:])
            end += interval
            beg += interval
            slopes1.append(slope1)
            slopes2.append(-slope2)
        return slopes1, slopes2, avg, coord


def check(direction="x", step="eq"):
    if step == "eq":
        log_file = LmpLogFile(f"../03_TC_MP/s01_eq_heat_flux.dir_{direction}.lmplog")
        pro_file = TprofileFile(f"../03_TC_MP/s01_eq_heat_flux.T_profile_{direction}.dat")
    elif step == "k":
        log_file = LmpLogFile(f"../03_TC_MP/s02_compute_kappa.dir_{direction}.lmplog")
        pro_file = TprofileFile(f"../03_TC_MP/s02_compute_kappa.T_profile_{direction}.dat")

    fig, axes = plt.subplots(nrows=3, ncols=3)

    log_file.df['kappa_num'] = log_file.df["f_transf_ek"] / (2 * log_file.df["Time"] * log_file.df["v_area"])
    log_file.df['rough_slope'] = log_file.df["f_aveTdiff"] / (log_file.df["v_Lb"]/2.0)
    log_file.df['kappa_Tdiff'] = log_file.df['kappa_num'] / log_file.df['rough_slope']
    log_file.df['kappa_Tdiff'] *= converter
    log_file.df['kappa_Tslope1'] = log_file.df['kappa_num'] / pro_file.avg_slope1
    log_file.df['kappa_Tslope2'] = log_file.df['kappa_num'] / pro_file.avg_slope2
    log_file.df['kappa_Tslope1'] *= converter
    log_file.df['kappa_Tslope2'] *= converter

    LW = 0.3
    BBOX = dict(boxstyle='round', fc='blanchedalmond', ec='orange', alpha=0.5)

    ax_avg_profile = axes[0, 0]
    ax_t_diff = axes[0,1]
    ax_transf = axes[0,2]
    ax_ro_slo = axes[1,0]
    ax_slopes1 = axes[1,1]
    ax_slopes2 = axes[1,2]
    ax_k_tslo1 = axes[2,1]
    ax_k_tslo2 = axes[2,2]
    ax_k_tdif = axes[2,0]

    ax_avg_profile.scatter(pro_file.coord, pro_file.avg, marker='.')
    ax_transf.plot(log_file.df["f_transf_ek"], linewidth=LW)
    ax_slopes1.plot(pro_file.slopes1, linewidth=LW)
    ax_slopes2.plot(pro_file.slopes2, linewidth=LW)
    ax_t_diff.plot(log_file.df["v_T_diff"], linewidth=LW)
    ax_k_tslo1.plot(log_file.df["kappa_Tslope1"], linewidth=LW)
    ax_k_tslo1.text(0.50, 0.05, 
           f'final = {log_file.df["kappa_Tslope1"].iloc[-1]:.2e}',
           bbox=BBOX, transform=ax_k_tslo1.transAxes)
    ax_k_tslo2.plot(log_file.df["kappa_Tslope2"], linewidth=LW)
    ax_k_tslo2.text(0.50, 0.05, 
           f'final = {log_file.df["kappa_Tslope2"].iloc[-1]:.2e}',
           bbox=BBOX, transform=ax_k_tslo2.transAxes)
    ax_k_tdif.plot(log_file.df["kappa_Tdiff"], linewidth=LW)
    ax_k_tdif.text(0.50, 0.90, 
           f'final = {log_file.df["kappa_Tdiff"].iloc[-1]:.2e}',
           bbox=BBOX, transform=ax_k_tdif.transAxes)
    ax_ro_slo.plot(log_file.df["rough_slope"], linewidth=LW)
    ax_ro_slo.text(0.50, 0.05, 
                   f'final = {log_file.df["rough_slope"].iloc[-1]:.2e}',
                   bbox=BBOX, transform=ax_ro_slo.transAxes)

    ax_avg_profile.set_xlabel(f"{direction} coord [ang]")
    ax_avg_profile.set_ylabel("T [K]")
    ax_t_diff.set_ylabel("T diff [K]")
    ax_t_diff.set_xlabel("frame #")
    ax_transf.set_ylabel("Cumulative transferred KE [eV]")
    ax_transf.set_xlabel("frame #")
    ax_ro_slo.set_ylabel("T gradient from beg-half T")
    ax_ro_slo.set_xlabel("frame #")
    ax_slopes1.set_ylabel("T gradient from fitting (L)")
    ax_slopes1.text(0.05, 0.90, 
                   f'average = {pro_file.avg_slope1:.2e} ({pro_file.std_slope1:.2e})',
                    bbox=BBOX, transform=ax_slopes1.transAxes)
    ax_slopes1.set_xlabel("chunk #")
    ax_slopes2.set_ylabel("T gradient from fitting (R)")
    ax_slopes2.text(0.05, 0.90, 
                   f'average = {pro_file.avg_slope2:.2e} ({pro_file.std_slope2:.2e})',
                    bbox=BBOX, transform=ax_slopes2.transAxes)
    ax_slopes2.set_xlabel("chunk #")
    ax_k_tslo1.set_ylabel("kappa from fitting (L)")
    ax_k_tslo1.set_xlabel("frame #")
    ax_k_tslo2.set_ylabel("kappa from fitting (R)")
    ax_k_tslo2.set_xlabel("frame #")
    ax_k_tdif.set_ylabel("kappa from T_diff")
    ax_k_tdif.set_xlabel("frame #")

    fig.tight_layout()
    fig.savefig(f"check_{direction}.{step}.pdf")


check(direction="x", step="eq")
check(direction="x", step="k")
check(direction="y", step="eq")
check(direction="y", step="k")
check(direction="z", step="eq")
check(direction="z", step="k")

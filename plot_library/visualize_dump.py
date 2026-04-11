#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import iread
import matplotlib
matplotlib.use("Agg")

# 1 amu / Å^3 = 1.66053906660 g / cm^3
AMU_A3_TO_G_CM3 = 1.66053906660


def get_time_from_info(info, fallback):
    for key in ("Time", "time", "t", "md_time", "step_time"):
        if key in info:
            try:
                return float(info[key])
            except Exception:
                pass
    return float(fallback)


def main():
    infile = "dump.xyz"
    if len(sys.argv) > 1:
        infile = sys.argv[1]

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Cannot find input file: {infile}")

    times = []
    a_list = []
    b_list = []
    c_list = []
    density_list = []
    vol_per_atom_list = []
    dvdt_list = []
    frame_idx = 0

    # 获取初始体积
    first_atoms = next(iread(infile, index=0, format="extxyz"))
    init_Volume = first_atoms.get_volume()

    for atoms in iread(infile, index=":", format="extxyz"):
        info = atoms.info
        n_atoms = len(atoms)
        if n_atoms == 0:
            frame_idx += 1
            continue

        volume = atoms.get_volume()
        if volume <= 0:
            frame_idx += 1
            continue

        t = get_time_from_info(info, frame_idx)
        a, b, c = atoms.cell.lengths()

        # 质量密度
        total_mass = atoms.get_masses().sum()
        density = total_mass / volume * AMU_A3_TO_G_CM3

        # 单位原子体积
        vol_per_atom = volume / n_atoms

        # 体积变化率
        dv_rate = volume / init_Volume - 1 


        a_list.append(a)
        b_list.append(b)
        c_list.append(c)
        density_list.append(density)
        vol_per_atom_list.append(vol_per_atom)
        dvdt_list.append(dv_rate)

        frame_idx += 1

    if len(a_list) == 0:
        raise RuntimeError("No valid frames found in input file.")

    
    a_arr = np.asarray(a_list, dtype=float)
    b_arr = np.asarray(b_list, dtype=float)
    c_arr = np.asarray(c_list, dtype=float)
    density_arr = np.asarray(density_list, dtype=float)
    vol_atom_arr = np.asarray(vol_per_atom_list, dtype=float)
    dvdt_arr = np.asarray(dvdt_list, dtype=float)
    times = np.arange(0, len(times), 1)

    # 密度变化：相对于第一帧
    density_change_arr = density_arr - density_arr[0]

    # 导出数据
    df = pd.DataFrame({
        "time": times,
        "a_A": a_arr,
        "b_A": b_arr,
        "c_A": c_arr,
        "density_g_cm3": density_arr,
        "density_change_g_cm3": density_change_arr,
        "volume_per_atom_A3": vol_atom_arr,
        "dVrate": dvdt_arr,
    })
    df.to_csv("analysis_data.csv", index=False, encoding="utf-8-sig")

    # ---------- 作图 ----------
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=180)
    ax1, ax2, ax3, ax4 = axes.ravel()

    color_a = "#4C78A8"
    color_b = "#F58518"
    color_c = "#54A24B"
    color_density = "#E45756"
    color_dvdt = "#72B7B2"
    color_vpa = "#B279A2"

    ax1.plot(times, a_arr, lw=1.8, color=color_a, label="a")
    ax1.plot(times, b_arr, lw=1.8, color=color_b, label="b")
    ax1.plot(times, c_arr, lw=1.8, color=color_c, label="c")
    #ax1.set_title("Cell lengths")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Length (Å)")
    ax1.grid(True, alpha=0.25, linestyle="--")
    ax1.legend(frameon=False)

    ax2.plot(times, density_change_arr, lw=1.8, color=color_density)
    #ax2.set_title("Density change")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("ΔDensity (g/cm³)")
    ax2.grid(True, alpha=0.25, linestyle="--")

    ax3.plot(times, dvdt_arr, lw=1.8, color=color_dvdt)
    #ax3.set_title("Volume change rate")
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Volume change rate")
    ax3.grid(True, alpha=0.25, linestyle="--")

    ax4.plot(times, vol_atom_arr, lw=1.8, color=color_vpa)
    #ax4.set_title("Volume per atom")
    ax4.set_xlabel("Time")
    ax4.set_ylabel("Volume / atom (Å³)")
    ax4.grid(True, alpha=0.25, linestyle="--")

    fig.suptitle(f"Trajectory analysis: {os.path.basename(infile)}", fontsize=15)
    fig.tight_layout()
    fig.savefig("dump_overview.png", bbox_inches="tight")
    plt.close(fig)

    print("Done.")
    print(f"Input : {infile}")
    print("Output: analyze_dump_data.csv")
    print("Output: dump_overview.png")


if __name__ == "__main__":
    main()
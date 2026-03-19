#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 自动使用gpumd计算声子谱
# 使用需求：
# 1. 输入文件：POSCAR 或其他格式的结构文件，run.in 文件
# 2. 输出目录：指定存放计算结果的目录
# 3. 可选参数：k-path 点数、计算波数范围等
# 自动生成 k-path, 并写入 kpoints.in 文件，
# 同时打印出自动识别的晶格类型、空间群号、国际符号、霍尔符号等信息。


import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib
matplotlib.use("Agg")

from ase.io import read, write
from ase.spacegroup.symmetrize import check_symmetry
from ase.dft.kpoints import special_paths, sc_special_points
from gpyumd.atoms import GpumdAtoms
from ase.dft.kpoints import special_paths, sc_special_points
import spglib


def _spg_get(ds, key, default=None):
    # spglib 新版本推荐属性接口；老版本还能 ds[key]
    if hasattr(ds, key):
        return getattr(ds, key)
    try:
        return ds[key]
    except Exception:
        return default

def _extract_centering_from_hall(hall: str) -> str:
    """
    Extract centering letter (P/F/I/R/A/B/C) from Hall symbol.
    Handles leading '-' like '-F 4 2 3'.
    """
    if not hall:
        return ""
    token0 = hall.split()[0]  # e.g. '-F'
    # find first alphabetic character in token0
    for ch in token0:
        if ch.isalpha():
            return ch.upper()
    return ""

def auto_kpath_and_points(atoms, symprec=1e-3):
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)

    spg_no = int(_spg_get(ds, "number", 0))
    hall = _spg_get(ds, "hall", "") or ""
    international = _spg_get(ds, "international", "") or ""
    centering = _extract_centering_from_hall(hall)

    # crystal system by space-group number
    if 195 <= spg_no <= 230:
        system = "cubic"
    elif 168 <= spg_no <= 194:
        system = "hexagonal"
    elif 143 <= spg_no <= 167:
        system = "trigonal"
    elif 75 <= spg_no <= 142:
        system = "tetragonal"
    elif 16 <= spg_no <= 74:
        system = "orthorhombic"
    elif 3 <= spg_no <= 15:
        system = "monoclinic"
    else:
        system = "triclinic"

    # map to ASE kpoint-table keys available in YOUR ASE
    # your ASE keys include: cubic, fcc, bcc, tetragonal, orthorhombic, hexagonal, monoclinic, rhombohedral type 1/2
    if system == "cubic":
        if centering == "F":
            key = "fcc"
        elif centering == "I":
            key = "bcc"
        else:
            key = "cubic"
    elif system == "hexagonal":
        key = "hexagonal"
    elif system == "tetragonal":
        key = "tetragonal"
    elif system == "orthorhombic":
        key = "orthorhombic"
    elif system == "monoclinic":
        key = "monoclinic"
    elif system == "trigonal":
        key = "rhombohedral type 1"
    else:
        raise ValueError(f"Low symmetry (SG={spg_no}, hall='{hall}'), no standard ASE path key.")

    if key not in special_paths or key not in sc_special_points:
        raise ValueError(f"ASE kpoint tables do not contain key='{key}'. "
                         f"Available keys: {list(special_paths.keys())}")

    meta = {
        "spacegroup_no": spg_no,
        "international": international,
        "hall": hall,
        "centering": centering,
        "system": system,
        "ase_key": key,
    }

    return special_paths[key], sc_special_points[key], meta


def gamma_label(s: str) -> str:
    """Convert 'G' to Gamma symbol for plotting."""
    return r'$\Gamma$' if s == "G" else s

def pretty_print_meta(meta):
    keys = ["ase_key", "system", "centering", "spacegroup_no", "international", "hall"]
    print("[INFO] Symmetry identification")
    for k in keys:
        if k in meta:
            print(f"  - {k:14s}: {meta[k]}")

def prepare_inputs(poscar: str, supercell=(12, 12, 12),
                   kpath=None, npoints=4000,
                   outdir=".", overwrite=False):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read structure and convert to GPUMD atoms
    uc = read(poscar)
    uc = GpumdAtoms(uc)
    uc.add_basis()

    # Build supercell for GPUMD
    cfg = uc.repeat([1, 1, 1])
    cfg.wrap()
    cfg = cfg.repeat(list(supercell))

    # Write structure + basis
    model_xyz = outdir / "model.xyz"
    if model_xyz.exists() and not overwrite:
        print(f"[SKIP] {model_xyz} exists (use --overwrite to regenerate)")
    else:
        write(str(model_xyz), cfg)
        cfg.write_basis(directory=str(outdir))
        print(f"[OK] Wrote {model_xyz} and basis files to {outdir}")

    kpath, sp, meta = auto_kpath_and_points(cfg, symprec=1e-2)
    pretty_print_meta(meta)
    print("[INFO] Auto k-path:", kpath)
    print("[INFO] Bravais (via centering+system):", f"{meta['system']} + {meta['centering']}",
      "SG:", meta["spacegroup_no"], meta["international"])


   
    # Write kpoints.in via your wrapper (uses ASE cell.bandpath internally)
    # Note: your write_kpoints returns path.get_linear_kpoint_axis()
    # In ASE this is typically (x, X, labels).
    axis = uc.write_kpoints(path=kpath, npoints=npoints, special_points=sp,
                            filename="kpoints.in", directory=str(outdir))

    # Try to robustly unpack
    if isinstance(axis, tuple) and len(axis) == 3:
        x, X, labels = axis
    else:
        # fallback: some ASE versions return only x
        x, X, labels = axis, None, None

    print(f"[OK] Wrote kpoints.in (npoints={npoints}, path='{kpath}') to {outdir}")
    return outdir, kpath, x, X, labels


def run_gpumd(cmd: str, workdir=".", logfile="gpumd.log"):
    workdir = Path(workdir)
    logpath = workdir / logfile

    print(f"[RUN] {cmd} (cwd={workdir})")
    with open(logpath, "w", encoding="utf-8") as f:
        # shell=True 方便你写 "mpirun -np 4 gpumd" 这类命令
        p = subprocess.run(cmd, cwd=str(workdir), shell=True, stdout=f, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        raise RuntimeError(f"GPUMD failed (return code {p.returncode}). See {logpath}")
    print(f"[OK] GPUMD finished. Log: {logpath}")


def load_omega2_to_thz(omega2_path: Path):
    # omega2.out: usually omega^2 values (can be negative for imaginary modes)
    data = np.loadtxt(str(omega2_path))
    # Convert sqrt(|omega^2|) to frequency with sign, then to THz:
    # Your old code used: sqrt(abs)/ (2*pi) * sign
    nu = np.sqrt(np.abs(data)) / (2.0 * np.pi) * np.sign(data)
    return nu


def plot_phonon(x, X, labels, nu, outpng="phonon.png",
                ylabel=r'$\nu$ (THz)', dpi=300):
    fig, ax = plt.subplots(figsize=(9.5, 7.0), dpi=dpi)

    # Plot all branches
    # nu shape: (nk, nbands)
    ax.plot(x, nu, lw=1.2, color="tab:blue")

    # Symmetry points / vertical lines / xticks
    if X is not None and labels is not None:
        ax.set_xticks(X)
        ax.set_xticklabels([gamma_label(l) for l in labels])
        for xv in X:
            ax.axvline(xv, lw=0.8, alpha=0.6)
    else:
        # fallback: keep default ticks
        pass

    ax.set_xlim(x[0], x[-1])

    # Y range: automatic but clean
    y_min = np.min(nu)
    y_max = np.max(nu)
    # Usually phonon plots show from 0 upward; keep negative if imaginary exists
    pad = 0.05 * (y_max - y_min + 1e-12)
    ax.set_ylim(y_min - pad, y_max + pad)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))

    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel("Wave vector", fontsize=18)

    ax.tick_params(direction="in", top=True, right=True, width=1.5, length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    fig.tight_layout()
    fig.savefig(outpng)
    print(f"[OK] Saved figure: {outpng}")


def main():
    ap = argparse.ArgumentParser(
        description="Prepare GPUMD phonon inputs, run GPUMD, and plot phonon dispersion."
    )
    ap.add_argument("--poscar", default="POSCAR", help="Structure file (POSCAR/vasp/etc.)")
    ap.add_argument("--outdir", default=".", help="Working directory for inputs/outputs")
    ap.add_argument("--supercell", nargs=3, type=int, default=[12, 12, 12], help="Supercell repeat (a b c)")
    ap.add_argument("--kpath", default=None,
                    help="K-path string. Default uses ASE fcc path, e.g. 'GXWKGLUWLK,UX'")
    ap.add_argument("--npoints", type=int, default=4000, help="Total k-points")
    ap.add_argument("--gpumd-cmd", default="gpumd",
                    help="Command to run GPUMD, e.g. 'gpumd' or 'mpirun -np 4 gpumd'")
    ap.add_argument("--skip-run", action="store_true", help="Only prepare + plot (assumes omega2.out exists)")
    ap.add_argument("--skip-prepare", action="store_true", help="Only run + plot (assumes inputs exist)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing model/basis/kpoints")
    ap.add_argument("--omega2", default="omega2.out", help="Omega2 output file name")
    ap.add_argument("--outpng", default="phonon.png", help="Output png name")
    args = ap.parse_args()

    outdir = Path(args.outdir)

    # Prepare
    if not args.skip_prepare:
        outdir, kpath, x, X, labels = prepare_inputs(
            poscar=args.poscar,
            supercell=tuple(args.supercell),
            kpath=args.kpath,
            npoints=args.npoints,
            outdir=str(outdir),
            overwrite=args.overwrite
        )
    else:
        # If skipping prepare, we still need x/X/labels for plotting
        # You can regenerate them quickly from POSCAR (cheap) without writing files
        _, kpath, x, X, labels = prepare_inputs(
            poscar=args.poscar,
            supercell=tuple(args.supercell),
            kpath=args.kpath,
            npoints=args.npoints,
            outdir=str(outdir),
            overwrite=False
        )

    # Run GPUMD
    if not args.skip_run:
        run_gpumd(args.gpumd_cmd, workdir=str(outdir))

    # Check omega2.out exists
    omega2_path = outdir / args.omega2
    if not omega2_path.exists():
        raise FileNotFoundError(f"{omega2_path} not found. Check gpumd.log / your phonon settings.")

    # Load and plot
    nu = load_omega2_to_thz(omega2_path)
    outpng_path = outdir / args.outpng
    plot_phonon(x, X, labels, nu, outpng=str(outpng_path))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dataset_quality_diagnosis.py

一个用于机器学习势训练集的“体检脚本”。

功能：
  1. 读取 extxyz / xyz 数据集中的能量、原子数、元素组成、力。
  2. 用 WLS/OLS 拟合元素参考能基线，诊断成分基线是否合理。
  3. 检查 residual energy 是否仍然强烈依赖原子数、成分空间。
  4. 可选读取 descriptor.out，对局域环境空间做 PCA 诊断。
  5. 输出图像、csv、文本报告，帮助判断数据集是否适合直接训练。

示例：
  python dataset_quality_diagnosis.py train.xyz --descriptor descriptor.out
  python dataset_quality_diagnosis.py train.xyz --fit-mode wls --outdir quality_report

依赖：
  numpy matplotlib ase scikit-learn
  scipy 可选；没有 scipy 时会自动跳过部分相关性统计。
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import ase.io
except Exception as exc:
    raise SystemExit("[ERROR] 需要安装 ase: pip install ase") from exc

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
except Exception as exc:
    raise SystemExit("[ERROR] 需要安装 scikit-learn: pip install scikit-learn") from exc

try:
    from scipy.stats import spearmanr, pearsonr
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


ENERGY_KEYS = [
    "energy", "Energy", "ENERGY", "free_energy", "Free_energy",
    "dft_energy", "DFT_energy", "E", "e",
]
FORCE_KEYS = ["forces", "force", "Forces", "FORCES", "dft_forces", "DFT_forces"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#E0E0E0",
    "grid.linewidth": 0.6,
    "grid.linestyle": "--",
    "figure.dpi": 140,
    "savefig.dpi": 220,
    "savefig.bbox": "tight",
})


def read_energies_from_info(atoms_list):
    vals = []
    key_used = None
    for atoms in atoms_list:
        found = False
        for key in ENERGY_KEYS:
            if key in atoms.info:
                vals.append(float(atoms.info[key]))
                key_used = key if key_used is None else key_used
                found = True
                break
        if not found:
            vals.append(np.nan)
    arr = np.asarray(vals, dtype=float)
    if np.isfinite(arr).all():
        return arr, key_used or "atoms.info"
    return None, None


def read_energies_by_regex(xyz_path: str | Path):
    """兼容普通 xyz 第二行中的 energy = ... / Energy = ... / nergy = ...。"""
    pattern = re.compile(
        r"(?:energy|Energy|ENERGY|nergy)\s*=\s*"
        r"(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    )
    energies = []
    with open(xyz_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                energies.append(float(m.group(1)))
    return np.asarray(energies, dtype=float)


def extract_forces(atoms_list):
    all_forces = []
    has_any = False
    key_used = None
    for atoms in atoms_list:
        arr = None
        for key in FORCE_KEYS:
            if key in atoms.arrays:
                arr = np.asarray(atoms.arrays[key], dtype=float)
                key_used = key if key_used is None else key_used
                break
        if arr is None:
            try:
                arr = atoms.get_forces()
                key_used = key_used or "calculator"
            except Exception:
                pass
        if arr is not None and arr.shape == (len(atoms), 3):
            all_forces.append(arr)
            has_any = True
        else:
            all_forces.append(None)
    return all_forces if has_any else None, key_used


def composition_matrix(atoms_list):
    elements = sorted({s for atoms in atoms_list for s in atoms.get_chemical_symbols()})
    rows = []
    systems = []
    for atoms in atoms_list:
        syms = atoms.get_chemical_symbols()
        c = Counter(syms)
        rows.append([c[e] for e in elements])
        systems.append("-".join(sorted(c.keys())))
    counts = np.asarray(rows, dtype=float)
    natoms = counts.sum(axis=1)
    comp = counts / np.maximum(natoms[:, None], 1.0)
    return elements, counts, comp, natoms, systems


def fit_atomic_reference(counts, energies, mode="wls"):
    natoms = counts.sum(axis=1)
    valid = np.isfinite(energies) & (natoms > 0)
    N = counts[valid]
    E = energies[valid]
    n = natoms[valid]

    if mode == "wls":
        w = 1.0 / np.maximum(n, 1.0)
        A = N * w[:, None]
        b = E * w
    elif mode == "ols":
        A = N
        b = E
    else:
        raise ValueError("mode must be 'wls' or 'ols'")

    mu, *_ = np.linalg.lstsq(A, b, rcond=None)
    Eref = counts @ mu
    raw_pa = energies / np.maximum(natoms, 1.0)
    ref_pa = Eref / np.maximum(natoms, 1.0)
    residual_pa = (energies - Eref) / np.maximum(natoms, 1.0)

    ss_res = float(np.nansum((energies - Eref) ** 2))
    mean_E = float(np.nanmean(energies))
    ss_tot = float(np.nansum((energies - mean_E) ** 2))
    r2_total = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    mask = np.isfinite(raw_pa) & np.isfinite(ref_pa)
    ss_res_pa = float(np.sum((raw_pa[mask] - ref_pa[mask]) ** 2))
    ss_tot_pa = float(np.sum((raw_pa[mask] - np.mean(raw_pa[mask])) ** 2))
    r2_pa = 1.0 - ss_res_pa / ss_tot_pa if ss_tot_pa > 0 else np.nan

    rank = int(np.linalg.matrix_rank(counts))
    if rank < counts.shape[1]:
        cond = float("inf")
    else:
        try:
            cond = float(np.linalg.cond(counts))
        except Exception:
            cond = float("inf")

    info = {
        "fit_mode": mode,
        "rank": rank,
        "n_elements": int(counts.shape[1]),
        "rank_deficient": bool(rank < counts.shape[1]),
        "condition_number": cond,
        "r2_total": float(r2_total),
        "r2_per_atom": float(r2_pa),
        "rmse_residual_eV_per_atom": float(np.sqrt(np.nanmean(residual_pa ** 2))),
        "mae_residual_eV_per_atom": float(np.nanmean(np.abs(residual_pa))),
        "residual_min": float(np.nanmin(residual_pa)),
        "residual_max": float(np.nanmax(residual_pa)),
        "raw_pa_min": float(np.nanmin(raw_pa)),
        "raw_pa_max": float(np.nanmax(raw_pa)),
    }
    return mu, Eref, raw_pa, ref_pa, residual_pa, info


def safe_corr(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return {"pearson": np.nan, "spearman": np.nan}
    x = np.asarray(x)[mask]
    y = np.asarray(y)[mask]
    if np.std(x) == 0 or np.std(y) == 0:
        return {"pearson": np.nan, "spearman": np.nan}
    if SCIPY_OK:
        return {
            "pearson": float(pearsonr(x, y)[0]),
            "spearman": float(spearmanr(x, y)[0]),
        }
    return {"pearson": float(np.corrcoef(x, y)[0, 1]), "spearman": np.nan}


def robust_z(x):
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad < 1e-12:
        std = np.nanstd(x)
        return (x - np.nanmean(x)) / max(std, 1e-12)
    return 0.6745 * (x - med) / mad


def pca_2d(X, scale=True):
    X = np.asarray(X, dtype=float)
    if scale:
        X2 = StandardScaler().fit_transform(X)
    else:
        X2 = X
    n_comp = min(2, X2.shape[1], X2.shape[0])
    pca = PCA(n_components=n_comp)
    Y = pca.fit_transform(X2)
    if n_comp == 1:
        Y = np.column_stack([Y[:, 0], np.zeros(len(Y))])
    return Y, pca


def save_csv(path, header, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            out = []
            for v in row:
                if isinstance(v, str):
                    out.append(v)
                elif isinstance(v, (int, np.integer)):
                    out.append(str(int(v)))
                elif v is None or (isinstance(v, float) and not np.isfinite(v)):
                    out.append("")
                else:
                    out.append(f"{float(v):.10g}")
            f.write(",".join(out) + "\n")


def plot_baseline(outdir, raw_pa, ref_pa, residual, natoms, comp_pc, species_count, info):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    sc = ax.scatter(ref_pa, raw_pa, c=species_count, s=14, alpha=0.75, edgecolors="none")
    lo = np.nanmin([np.nanmin(ref_pa), np.nanmin(raw_pa)])
    hi = np.nanmax([np.nanmax(ref_pa), np.nanmax(raw_pa)])
    pad = 0.05 * max(hi - lo, 1e-12)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="#555555", lw=1.1)
    ax.set_xlabel(r"Fitted composition baseline  $E_{ref}/N$  (eV/atom)")
    ax.set_ylabel(r"Raw energy  $E/N$  (eV/atom)")
    ax.set_title("Baseline fit: raw energy vs fitted baseline", fontweight="bold")
    text = (
        f"fit mode = {info['fit_mode']}\n"
        f"R² total = {info['r2_total']:.5f}\n"
        f"R² per atom = {info['r2_per_atom']:.5f}\n"
        f"RMSE/atom = {info['rmse_residual_eV_per_atom']*1000:.2f} meV\n"
        f"rank = {info['rank']} / {info['n_elements']}\n"
        f"cond(N) = {info['condition_number']:.2e}" if np.isfinite(info['condition_number']) else
        f"fit mode = {info['fit_mode']}\nR² total = {info['r2_total']:.5f}\nR² per atom = {info['r2_per_atom']:.5f}\nRMSE/atom = {info['rmse_residual_eV_per_atom']*1000:.2f} meV\nrank = {info['rank']} / {info['n_elements']}\ncond(N) = inf"
    )
    ax.text(0.04, 0.96, text, transform=ax.transAxes, va="top", ha="left",
            fontsize=8, family="monospace",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.9))
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Number of species")

    ax = axes[0, 1]
    ax.scatter(ref_pa, residual, c=species_count, s=14, alpha=0.75, edgecolors="none")
    ax.axhline(0, ls="--", color="#555555", lw=1.1)
    ax.set_xlabel(r"Fitted composition baseline  $E_{ref}/N$  (eV/atom)")
    ax.set_ylabel(r"Shifted energy  $(E-E_{ref})/N$  (eV/atom)")
    ax.set_title("Residual after composition baseline removal", fontweight="bold")

    ax = axes[1, 0]
    sc = ax.scatter(natoms, residual, c=ref_pa, s=14, alpha=0.75, edgecolors="none")
    ax.axhline(0, ls="--", color="#555555", lw=1.1)
    ax.set_xlabel(r"Number of atoms  $N$")
    ax.set_ylabel(r"Shifted energy  $(E-E_{ref})/N$  (eV/atom)")
    ax.set_title("Residual vs system size", fontweight="bold")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label(r"$E_{ref}/N$ (eV/atom)")

    ax = axes[1, 1]
    sc = ax.scatter(comp_pc[:, 0], residual, c=natoms, s=14, alpha=0.75, edgecolors="none")
    ax.axhline(0, ls="--", color="#555555", lw=1.1)
    ax.set_xlabel("Composition PC1")
    ax.set_ylabel(r"Shifted energy  $(E-E_{ref})/N$  (eV/atom)")
    ax.set_title("Residual in composition space", fontweight="bold")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Number of atoms")

    fig.suptitle("Atomic Reference Energy Baseline Fit Diagnostics", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(outdir / "01_atomic_reference_baseline.png")
    plt.close(fig)


def plot_distributions(outdir, raw_pa, residual, natoms, species_count, force_mags=None):
    nplots = 4 if force_mags is not None else 3
    fig, axes = plt.subplots(1, nplots, figsize=(4 * nplots, 3.6))
    if nplots == 1:
        axes = [axes]
    axes[0].hist(raw_pa[np.isfinite(raw_pa)], bins="auto", alpha=0.85, edgecolor="white")
    axes[0].set_title("Raw E/N")
    axes[0].set_xlabel("eV/atom")
    axes[0].set_ylabel("Count")

    axes[1].hist(residual[np.isfinite(residual)], bins="auto", alpha=0.85, edgecolor="white")
    axes[1].set_title("Shifted residual")
    axes[1].set_xlabel("eV/atom")

    axes[2].hist(natoms[np.isfinite(natoms)], bins="auto", alpha=0.85, edgecolor="white")
    axes[2].set_title("Atom count")
    axes[2].set_xlabel("N")

    if force_mags is not None:
        fm = force_mags[np.isfinite(force_mags)]
        axes[3].hist(fm, bins="auto", alpha=0.85, edgecolor="white")
        axes[3].set_title("Force magnitude")
        axes[3].set_xlabel(r"|F| (eV/Å)")

    fig.suptitle("Dataset Basic Distributions", fontweight="bold")
    plt.tight_layout()
    fig.savefig(outdir / "02_basic_distributions.png")
    plt.close(fig)


def plot_composition_pca(outdir, comp_pc, residual, raw_pa, natoms, systems):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc = axes[0].scatter(comp_pc[:, 0], comp_pc[:, 1], c=residual, s=14, alpha=0.8, edgecolors="none")
    axes[0].set_xlabel("Composition PC1")
    axes[0].set_ylabel("Composition PC2")
    axes[0].set_title("Composition PCA colored by residual", fontweight="bold")
    cb = fig.colorbar(sc, ax=axes[0])
    cb.set_label("Shifted energy (eV/atom)")

    sc = axes[1].scatter(comp_pc[:, 0], comp_pc[:, 1], c=natoms, s=14, alpha=0.8, edgecolors="none")
    axes[1].set_xlabel("Composition PC1")
    axes[1].set_ylabel("Composition PC2")
    axes[1].set_title("Composition PCA colored by N", fontweight="bold")
    cb = fig.colorbar(sc, ax=axes[1])
    cb.set_label("Number of atoms")

    plt.tight_layout()
    fig.savefig(outdir / "03_composition_pca.png")
    plt.close(fig)


def plot_descriptor_pca(outdir, des_pc, pca, residual, natoms, raw_pa):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    labels = ["residual E", "raw E/N", "N"]
    colors = [residual, raw_pa, natoms]
    titles = [
        "Descriptor PCA colored by residual",
        "Descriptor PCA colored by raw E/N",
        "Descriptor PCA colored by N",
    ]
    for ax, c, title, label in zip(axes, colors, titles, labels):
        sc = ax.scatter(des_pc[:, 0], des_pc[:, 1], c=c, s=14, alpha=0.8, edgecolors="none")
        ax.set_xlabel("Descriptor PC1")
        ax.set_ylabel("Descriptor PC2")
        ax.set_title(title, fontweight="bold")
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label(label)
    evr = pca.explained_variance_ratio_
    fig.suptitle(f"Descriptor PCA  (PC1={evr[0]*100:.1f}%, PC2={(evr[1] if len(evr)>1 else 0)*100:.1f}%)", fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(outdir / "04_descriptor_pca.png")
    plt.close(fig)


def evaluate_quality(info, residual, natoms, comp_pc, force_mags=None, descriptor_warning=None):
    flags = []
    score = 100

    rmse = info["rmse_residual_eV_per_atom"]
    r2pa = info["r2_per_atom"]
    res_range = info["residual_max"] - info["residual_min"]
    corr_n = safe_corr(natoms, residual)
    corr_c = safe_corr(comp_pc[:, 0], residual)

    if info["rank_deficient"]:
        flags.append(("SEVERE", "组成矩阵秩亏：某些元素参考能不唯一；需要检查元素覆盖或按体系分组。"))
        score -= 18
    elif info["condition_number"] > 1e8:
        flags.append(("WARN", "组成矩阵条件数很大：元素参考能拟合可能不稳定。"))
        score -= 10

    if rmse > 0.5:
        flags.append(("SEVERE", f"WLS residual RMSE 很大：{rmse:.3f} eV/atom，数据内部差异或异常构型很强。"))
        score -= 25
    elif rmse > 0.2:
        flags.append(("WARN", f"WLS residual RMSE 偏大：{rmse:.3f} eV/atom，建议分组诊断。"))
        score -= 12
    elif rmse > 0.1:
        flags.append(("NOTE", f"WLS residual RMSE 中等：{rmse:.3f} eV/atom。"))
        score -= 5

    if np.isfinite(r2pa) and r2pa < 0.2:
        flags.append(("WARN", f"R² per atom 很低：{r2pa:.3f}，元素线性基线几乎不能解释每原子能量变化。"))
        score -= 10

    if res_range > 1.0:
        flags.append(("WARN", f"平移后能量跨度很大：{res_range:.3f} eV/atom，可能有高能异常结构或多类结构混合。"))
        score -= 10

    if np.isfinite(corr_n["spearman"]) and abs(corr_n["spearman"]) > 0.45:
        flags.append(("WARN", f"residual 与原子数 N 相关性明显：Spearman={corr_n['spearman']:.3f}，存在尺寸/来源偏置。"))
        score -= 10

    if np.isfinite(corr_c["spearman"]) and abs(corr_c["spearman"]) > 0.45:
        flags.append(("WARN", f"residual 与 composition PC1 相关性明显：Spearman={corr_c['spearman']:.3f}，成分偏置未完全去除。"))
        score -= 10

    rz = robust_z(residual)
    outlier_frac = float(np.mean(np.abs(rz) > 3.5))
    if outlier_frac > 0.05:
        flags.append(("WARN", f"robust residual outlier 比例较高：{outlier_frac*100:.1f}%。"))
        score -= 8

    if force_mags is not None and len(force_mags) > 0:
        p99 = float(np.nanpercentile(force_mags, 99))
        p999 = float(np.nanpercentile(force_mags, 99.9))
        if p999 > 50:
            flags.append(("SEVERE", f"存在极大力值：|F| 99.9%={p999:.2f} eV/A，需检查原子重叠或异常构型。"))
            score -= 20
        elif p99 > 20:
            flags.append(("WARN", f"高力尾部较强：|F| 99%={p99:.2f} eV/A，建议单独检查高力构型。"))
            score -= 8

    if descriptor_warning:
        flags.append(("WARN", descriptor_warning))
        score -= 6

    score = int(max(0, min(100, score)))
    if score >= 80:
        grade = "GOOD"
    elif score >= 60:
        grade = "OK_WITH_WARNINGS"
    elif score >= 40:
        grade = "RISKY"
    else:
        grade = "POOR"

    return {
        "score": score,
        "grade": grade,
        "flags": flags,
        "corr_residual_vs_N": corr_n,
        "corr_residual_vs_composition_PC1": corr_c,
        "residual_outlier_fraction": outlier_frac,
    }


def write_report(outdir, args, elements, mu, info, quality, systems, natoms, residual, force_mags, descriptor_info):
    system_counts = Counter(systems)
    n_species = np.asarray([len(s.split("-")) if s else 0 for s in systems])

    lines = []
    lines.append("# Dataset Quality Diagnosis Report")
    lines.append("")
    lines.append(f"Input file: {args.xyz}")
    lines.append(f"Fit mode: {args.fit_mode}")
    lines.append(f"Number of configurations: {len(natoms)}")
    lines.append(f"Elements: {' '.join(elements)}")
    lines.append(f"Quality score: {quality['score']} / 100")
    lines.append(f"Quality grade: {quality['grade']}")
    lines.append("")

    lines.append("## Atomic reference baseline")
    lines.append("")
    for e, v in zip(elements, mu):
        lines.append(f"- mu({e}) = {v:.8f} eV/atom")
    lines.append(f"- rank = {info['rank']} / {info['n_elements']}")
    lines.append(f"- condition number = {info['condition_number']:.6e}" if np.isfinite(info['condition_number']) else "- condition number = inf")
    lines.append(f"- R2 total = {info['r2_total']:.6f}")
    lines.append(f"- R2 per atom = {info['r2_per_atom']:.6f}")
    lines.append(f"- residual RMSE = {info['rmse_residual_eV_per_atom']*1000:.3f} meV/atom")
    lines.append(f"- residual MAE = {info['mae_residual_eV_per_atom']*1000:.3f} meV/atom")
    lines.append(f"- residual range = [{info['residual_min']:.6f}, {info['residual_max']:.6f}] eV/atom")
    lines.append("")

    lines.append("## Dataset structure")
    lines.append("")
    lines.append(f"- atom count range = {int(np.min(natoms))} to {int(np.max(natoms))}")
    lines.append(f"- atom count median = {float(np.median(natoms)):.1f}")
    lines.append(f"- species count distribution = {dict(Counter(n_species))}")
    lines.append("- top chemical systems:")
    for sys, cnt in system_counts.most_common(12):
        lines.append(f"  - {sys}: {cnt}")
    lines.append("")

    if force_mags is not None and len(force_mags) > 0:
        lines.append("## Force statistics")
        lines.append("")
        for q in [50, 90, 95, 99, 99.9, 100]:
            lines.append(f"- |F| p{q} = {np.nanpercentile(force_mags, q):.6f} eV/A")
        lines.append("")

    if descriptor_info:
        lines.append("## Descriptor PCA")
        lines.append("")
        lines.append(f"- descriptor shape = {descriptor_info['shape']}")
        lines.append(f"- PC1 explained variance = {descriptor_info['pc1_var']*100:.3f}%")
        lines.append(f"- PC2 explained variance = {descriptor_info['pc2_var']*100:.3f}%")
        lines.append("")

    lines.append("## Correlations")
    lines.append("")
    cn = quality["corr_residual_vs_N"]
    cc = quality["corr_residual_vs_composition_PC1"]
    lines.append(f"- residual vs N: Pearson={cn['pearson']:.4f}, Spearman={cn['spearman']:.4f}")
    lines.append(f"- residual vs composition PC1: Pearson={cc['pearson']:.4f}, Spearman={cc['spearman']:.4f}")
    lines.append(f"- robust residual outlier fraction = {quality['residual_outlier_fraction']*100:.3f}%")
    lines.append("")

    lines.append("## Flags")
    lines.append("")
    if quality["flags"]:
        for level, msg in quality["flags"]:
            lines.append(f"- [{level}] {msg}")
    else:
        lines.append("- No major warning flags.")
    lines.append("")

    lines.append("## Suggested next checks")
    lines.append("")
    lines.append("1. Inspect 01_atomic_reference_baseline.png: if residual forms vertical stripes or large tails, split by chemical system / atom count / source.")
    lines.append("2. Inspect outlier_configs.csv: export those structures and visualize them in OVITO/VESTA/ASE.")
    lines.append("3. If descriptor.out is provided, inspect 04_descriptor_pca.png: train/test or selected/unselected sets should cover the same descriptor space.")
    lines.append("4. For sampling, prefer group-wise sampling: chemical system -> shifted-energy quantile bins -> descriptor PCA/FPS.")
    lines.append("")

    (outdir / "quality_report.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="ML potential dataset quality diagnosis tool")
    parser.add_argument("xyz", help="input train.xyz / extxyz file")
    parser.add_argument("--descriptor", default=None, help="optional descriptor.out file")
    parser.add_argument("--outdir", default="dataset_quality_report", help="output directory")
    parser.add_argument("--fit-mode", choices=["wls", "ols"], default="wls", help="atomic reference fit mode")
    parser.add_argument("--top-outliers", type=int, default=80, help="number of residual outliers to save")
    args = parser.parse_args()

    xyz_path = Path(args.xyz)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    atoms_list = ase.io.read(str(xyz_path), index=":")
    energies, energy_source = read_energies_from_info(atoms_list)
    if energies is None:
        energies = read_energies_by_regex(xyz_path)
        energy_source = "regex"
    if len(energies) != len(atoms_list):
        raise SystemExit(
            f"[ERROR] 能量数量({len(energies)})与构型数量({len(atoms_list)})不一致。"
            "请确认 xyz 第二行含 energy=...，或 extxyz info 中有 energy 字段。"
        )

    elements, counts, comp, natoms, systems = composition_matrix(atoms_list)
    species_count = np.sum(counts > 0, axis=1)
    mu, Eref, raw_pa, ref_pa, residual, info = fit_atomic_reference(counts, energies, mode=args.fit_mode)

    comp_pc, comp_pca = pca_2d(comp, scale=False)

    forces, force_source = extract_forces(atoms_list)
    force_mags = None
    if forces is not None:
        mags = []
        for arr in forces:
            if arr is not None:
                mags.extend(np.linalg.norm(arr, axis=1).tolist())
        force_mags = np.asarray(mags, dtype=float) if mags else None

    descriptor_info = None
    descriptor_warning = None
    if args.descriptor:
        des = np.loadtxt(args.descriptor)
        if des.ndim == 1:
            des = des.reshape(1, -1)
        if des.shape[0] != len(atoms_list):
            descriptor_warning = f"descriptor 行数({des.shape[0]})与构型数({len(atoms_list)})不一致，跳过 descriptor PCA。"
        else:
            des_pc, des_pca = pca_2d(des, scale=True)
            plot_descriptor_pca(outdir, des_pc, des_pca, residual, natoms, raw_pa)
            evr = des_pca.explained_variance_ratio_
            descriptor_info = {
                "shape": f"{des.shape[0]} x {des.shape[1]}",
                "pc1_var": float(evr[0]),
                "pc2_var": float(evr[1]) if len(evr) > 1 else 0.0,
            }
            save_csv(
                outdir / "descriptor_pca.csv",
                ["index", "descriptor_PC1", "descriptor_PC2", "raw_E_per_atom", "shifted_E_per_atom", "N", "chemical_system"],
                [[i, des_pc[i, 0], des_pc[i, 1], raw_pa[i], residual[i], natoms[i], systems[i]] for i in range(len(natoms))],
            )

    quality = evaluate_quality(info, residual, natoms, comp_pc, force_mags, descriptor_warning)

    plot_baseline(outdir, raw_pa, ref_pa, residual, natoms, comp_pc, species_count, info)
    plot_distributions(outdir, raw_pa, residual, natoms, species_count, force_mags)
    plot_composition_pca(outdir, comp_pc, residual, raw_pa, natoms, systems)

    rz = robust_z(residual)
    order = np.argsort(-np.abs(rz))[:args.top_outliers]
    save_csv(
        outdir / "outlier_configs.csv",
        ["index", "robust_z_residual", "raw_E_per_atom", "shifted_E_per_atom", "Eref_per_atom", "N", "n_species", "chemical_system"],
        [[int(i), rz[i], raw_pa[i], residual[i], ref_pa[i], natoms[i], species_count[i], systems[i]] for i in order],
    )

    save_csv(
        outdir / "config_summary.csv",
        ["index", "raw_E_per_atom", "shifted_E_per_atom", "Eref_per_atom", "N", "n_species", "composition_PC1", "composition_PC2", "chemical_system"],
        [[i, raw_pa[i], residual[i], ref_pa[i], natoms[i], species_count[i], comp_pc[i, 0], comp_pc[i, 1], systems[i]] for i in range(len(natoms))],
    )

    summary_json = {
        "input": str(xyz_path),
        "energy_source": energy_source,
        "force_source": force_source,
        "n_configs": len(atoms_list),
        "elements": elements,
        "fit_info": info,
        "quality": {
            "score": quality["score"],
            "grade": quality["grade"],
            "flags": quality["flags"],
            "corr_residual_vs_N": quality["corr_residual_vs_N"],
            "corr_residual_vs_composition_PC1": quality["corr_residual_vs_composition_PC1"],
            "residual_outlier_fraction": quality["residual_outlier_fraction"],
        },
        "descriptor_info": descriptor_info,
    }
    (outdir / "summary.json").write_text(json.dumps(summary_json, indent=2, ensure_ascii=False), encoding="utf-8")

    write_report(outdir, args, elements, mu, info, quality, systems, natoms, residual, force_mags, descriptor_info)

    print("\n=== Dataset quality diagnosis finished ===")
    print(f"Output directory : {outdir}")
    print(f"Quality score    : {quality['score']} / 100")
    print(f"Quality grade    : {quality['grade']}")
    print("Key files:")
    print(f"  - {outdir / 'quality_report.md'}")
    print(f"  - {outdir / '01_atomic_reference_baseline.png'}")
    print(f"  - {outdir / '02_basic_distributions.png'}")
    print(f"  - {outdir / '03_composition_pca.png'}")
    if descriptor_info:
        print(f"  - {outdir / '04_descriptor_pca.png'}")
    print(f"  - {outdir / 'outlier_configs.csv'}")


if __name__ == "__main__":
    main()

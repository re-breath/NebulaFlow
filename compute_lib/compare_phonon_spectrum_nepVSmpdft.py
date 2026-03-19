#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 该脚本用来比较nep和mpdft的phonon 谱图是否一致

# 使用方式
#python3 comparebeauty.py --poscar POSCAR --omega2 omega2.out --mpjson mp-661_phonon_bs_dfpt.json 

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from ase.io import read
from ase.dft.kpoints import special_paths, sc_special_points
import spglib
import ase.io

THZ_TO_CM1 = 33.35641


def _spg_get(ds, key, default=None):
    if hasattr(ds, key):
        return getattr(ds, key)
    try:
        return ds[key]
    except Exception:
        return default


def _extract_centering_from_hall(hall: str) -> str:
    if not hall:
        return ""
    token0 = hall.split()[0]
    for ch in token0:
        if ch.isalpha():
            return ch.upper()
    return ""


def auto_ase_key(atoms, symprec=1e-2):
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
    spg_no = int(_spg_get(ds, "number", 0))
    hall = _spg_get(ds, "hall", "") or ""
    centering = _extract_centering_from_hall(hall)

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

    if system == "cubic":
        if centering == "F":
            return "fcc"
        elif centering == "I":
            return "bcc"
        else:
            return "cubic"
    elif system == "hexagonal":
        return "hexagonal"
    elif system == "tetragonal":
        return "tetragonal"
    elif system == "orthorhombic":
        return "orthorhombic"
    elif system == "monoclinic":
        return "monoclinic"
    elif system == "trigonal":
        return "rhombohedral type 1"
    else:
        raise ValueError("Low symmetry; no default ASE path key.")


def parse_kpath_labels(kpath: str):
    """
    把 'GXWKGLUWLK,UX' -> ['G','X','W','K','G','L','U','W','L','K','U','X']
    注意：这里保留 K 和 U 两个点，后面再合并成 K|U（显示用）
    """
    out = []
    for part in kpath.split(","):
        part = part.strip()
        for ch in part:
            out.append(ch)
    return out


def gamma_label(lab: str) -> str:
    if lab in ("G", "\\Gamma", "Γ"):
        return r"$\Gamma$"
    return lab


def merge_KU_display(labels, X):
    """
    如果最后一段出现 ... K, U, X 且 K/U 对应位置相同或非常接近 -> 显示为 K|U
    """
    labels = list(labels)
    X = np.array(X, dtype=float)

    merged_labels = [labels[0]]
    merged_X = [X[0]]
    eps = 1e-8

    for lab, x in zip(labels[1:], X[1:]):
        if abs(x - merged_X[-1]) < eps:
            merged_labels[-1] = merged_labels[-1] + "|" + lab
        else:
            merged_labels.append(lab)
            merged_X.append(x)

    return merged_labels, np.array(merged_X)


def gpumd_axis(poscar, kpath=None, npoints=4000, ase_key=None):
    atoms = read(poscar)
    if ase_key is None:
        ase_key = auto_ase_key(atoms)
    sp = sc_special_points[ase_key]
    path = kpath or special_paths[ase_key]
    bp = atoms.cell.bandpath(path, npoints=npoints, special_points=sp)
    x, X, labels = bp.get_linear_kpoint_axis()
    # labels 在 ASE 里通常用 'G' 表示 Γ
    labels = [("G" if l in ("\\Gamma", "Γ") else l) for l in labels]
    # 合并显示用 K|U（若重合）
    labels_disp, X_disp = merge_KU_display(labels, X)
    return path, x, X, labels, X_disp, labels_disp


def load_gpumd_omega2(omega2_path: Path):
    data = np.loadtxt(str(omega2_path))
    nu = np.sqrt(np.abs(data)) / (2.0 * np.pi) * np.sign(data)  # THz
    if nu.ndim == 1:
        nu = nu[:, None]


    return nu


def load_mp_json(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def get_reciprocal_matrix(mp):
    rl = mp["reciprocal_lattice"]
    if isinstance(rl, dict) and "matrix" in rl:
        return np.array(rl["matrix"], dtype=float)
    return np.array(rl, dtype=float)


def cumulative_k_distance(qpoints_frac, rec_mat):
    q = np.array(qpoints_frac, dtype=float)
    k = q @ rec_mat
    dk = np.linalg.norm(k[1:] - k[:-1], axis=1)
    return np.concatenate([[0.0], np.cumsum(dk)])


def find_indices_for_label_sequence(qpoints, labels_dict, label_seq, tol=1e-6):
    """
    labels_dict 是 {label: [qx,qy,qz]}（你的 JSON 就是这个）
    在 qpoints 里按 label_seq 顺序找出现位置（允许重复点）
    """
    q = np.array(qpoints, dtype=float)
    lab2coord = {}
    for lab, coord in labels_dict.items():
        lab = lab.replace("\\Gamma", "G").replace("Γ", "G")
        lab2coord[lab] = np.array(coord, dtype=float)

    def candidates_for(lab):
        if lab not in lab2coord:
            return []
        c = lab2coord[lab]
        d = np.linalg.norm(q - c[None, :], axis=1)
        idxs = np.where(d < tol)[0]
        if idxs.size == 0:
            idxs = np.where(d < 1e-4)[0]
        return idxs.tolist()

    idx_list = []
    prev = -1
    for lab in label_seq:
        cand = candidates_for(lab)
        chosen = None
        for i in cand:
            if i > prev:
                chosen = i
                break
        if chosen is None and cand:
            chosen = cand[0]
        if chosen is None:
            return None  # 匹配失败
        idx_list.append(int(chosen))
        prev = chosen

    # 去重（有时同一点重复）
    out = []
    for i in idx_list:
        if not out or i != out[-1]:
            out.append(i)
    return out


def piecewise_map(x_src, src_ticks, dst_ticks):
    x_src = np.asarray(x_src, float)
    src_ticks = np.asarray(src_ticks, float)
    dst_ticks = np.asarray(dst_ticks, float)

    x_dst = np.empty_like(x_src)
    x_dst[:] = np.nan

    for i in range(len(src_ticks) - 1):
        a0, a1 = src_ticks[i], src_ticks[i + 1]
        b0, b1 = dst_ticks[i], dst_ticks[i + 1]
        if abs(a1 - a0) < 1e-12:
            mask = (x_src >= a0 - 1e-12) & (x_src <= a1 + 1e-12)
            x_dst[mask] = b0
            continue
        lo, hi = (a0, a1) if a0 <= a1 else (a1, a0)
        mask = (x_src >= lo - 1e-12) & (x_src <= hi + 1e-12)
        t = (x_src[mask] - a0) / (a1 - a0)
        x_dst[mask] = b0 + t * (b1 - b0)

    # 兜底填充
    nan = np.isnan(x_dst)
    if np.any(nan):
        good = np.where(~nan)[0]
        for j in np.where(nan)[0]:
            k = good[np.argmin(np.abs(good - j))]
            x_dst[j] = x_dst[k]
    return x_dst

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

def beautify_phonon_ax(
    ax,
    X_ticks,
    labels,
    xlim,
    ylabel=r"Frequencies (cm$^{-1}$)",
    title=None,
    legend_loc="upper right",
    show_grid=False,
    spine_width=1.2,
    tick_width=1.2,
    tick_length=6,
    minor_tick_length=3.5,
    vline_color="0.85",
    vline_lw=0.8,
    vline_alpha=1.0,
    y_nbins=7,
    y_pad_frac=0.03,
):
    """
    美化声子谱坐标轴到论文风格。
    你在画完曲线后调用一次：
        beautify_phonon_ax(ax, X_disp, labels_disp, (x_gp[0], x_gp[-1]), ...)
    """

    # --- x ticks / labels ---
    ax.set_xticks(X_ticks)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_xlim(*xlim)

    # --- vertical lines at high-symmetry points ---
    for xv in X_ticks:
        ax.axvline(xv, color=vline_color, lw=vline_lw, alpha=vline_alpha, zorder=0)

    # --- labels & title ---
    ax.set_xlabel("Wave vector", fontsize=13,fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=13,fontweight="bold")
    if title:
        ax.set_title(title, fontsize=14, pad=10)

    # --- y ticks: clean & consistent ---
    ax.yaxis.set_major_locator(MaxNLocator(nbins=y_nbins))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # pad y-limits a bit (so curves don't touch the frame)
    y0, y1 = ax.get_ylim()
    yr = y1 - y0 if (y1 - y0) > 0 else 1.0
    ax.set_ylim(0, y1 + y_pad_frac * yr)

    # --- ticks style (paper-like) ---
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        width=tick_width,
        length=tick_length,
        labelsize=12,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        width=tick_width * 0.9,
        length=minor_tick_length,
    )

    # --- spines ---
    for s in ax.spines.values():
        s.set_linewidth(spine_width)

    # --- grid (usually OFF for phonons; optional) ---
    if show_grid:
        ax.grid(True, which="major", axis="y", color="0.9", lw=0.6)
        ax.grid(True, which="minor", axis="y", color="0.95", lw=0.5)

    # --- legend: only keep unique labels & make it unobtrusive ---
    handles, leg_labels = ax.get_legend_handles_labels()
    uniq = {}
    for h, lab in zip(handles, leg_labels):
        if lab and lab not in uniq:
            uniq[lab] = h
    if uniq:
        leg = ax.legend(
            uniq.values(),
            uniq.keys(),
            loc=legend_loc,
            frameon=True,
            framealpha=0.85,
            borderpad=0.4,
            handlelength=2.0,
            fontsize=11,
        )
        leg.get_frame().set_linewidth(0.8)

    return ax

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poscar", default="POSCAR", help="VASP格式文件")
    ap.add_argument("--omega2", default="omega2.out", help="GPUMD计算得到的声子频率文件")
    ap.add_argument("--mpjson", help="MP-DFT计算得到的声子谱文件")
    ap.add_argument("--kpath", help="高对称点路径,例如：GMKGALHA,LM,KH")
    ap.add_argument("--npoints", type=int, default=4000, help="GPUMD计算得到的声子频率文件中的点数量")
    ap.add_argument("--unit", choices=["cm1", "thz"], default="cm1", help="单位，默认cm1")
    ap.add_argument("--outpng", default="overlay_compare.png", help="输出图片文件名")
    ap.add_argument("--ase-key", default=None, help="ASE读取POSCAR时的键名")
    args = ap.parse_args()

    # GPUMD axis
    path, x_gp, X_gp, labels_gp, X_disp, labels_disp = gpumd_axis(
        args.poscar, kpath=args.kpath, npoints=args.npoints, ase_key=args.ase_key
    )
    

    # GPUMD data
    nu_gp = load_gpumd_omega2(Path(args.omega2))
    if nu_gp.shape[0] != len(x_gp):
        raise ValueError(f"omega2 nk={nu_gp.shape[0]} but axis nk={len(x_gp)}; check npoints.")

    # MP-DFT data
    mp = load_mp_json(Path(args.mpjson))
    freqs = np.array(mp["frequencies"], dtype=float)
    qpoints = mp["qpoints"]
    if freqs.shape[0] != len(qpoints) and freqs.shape[1] == len(qpoints):
        freqs = freqs.T
    if freqs.shape[0] != len(qpoints):
        raise ValueError(f"MP qpoints={len(qpoints)} but frequencies shape={freqs.shape}")

    rec = get_reciprocal_matrix(mp)
    x_dft = cumulative_k_distance(qpoints, rec)

    # 用 kpath 序列（G X W K G L U W L K U X）在 MP qpoints 里找高对称点 index
    label_seq = parse_kpath_labels(path)
    idx_ticks = find_indices_for_label_sequence(qpoints, mp.get("labels_dict", {}), label_seq, tol=1e-6)
    

    # 生成 DFT->GPUMD 的映射
    if idx_ticks is not None and len(idx_ticks) >= 2 and len(X_gp) >= 2:
        m = min(len(idx_ticks), len(X_gp))
        src_ticks = x_dft[idx_ticks[:m]]
        dst_ticks = X_gp[:m]
        x_dft_mapped = piecewise_map(x_dft, src_ticks, dst_ticks)
    else:
        # 匹配失败：退化为全局线性映射（不会出现你那种“尖峰扭曲”）
        x_dft_mapped = (x_dft / x_dft[-1]) * x_gp[-1]

    # 单位
    if args.unit == "cm1":
        y_gp = nu_gp * THZ_TO_CM1
        y_dft = freqs * THZ_TO_CM1
        ylabel = r"Frequencies (cm$^{-1}$)"
    else:
        y_gp = nu_gp
        y_dft = freqs
        ylabel = "Frequency (THz)"

    fig, ax = plt.subplots(figsize=(10, 7), dpi=300)

    COLOR_GPUMD = "#26292B"   # blue-ish
    COLOR_DFT   = "#9D1919"   # magenta
    LW_GP, LW_DFT = 1.6, 1.4
    A_GP,  A_DFT  = 0.95, 0.80

    lines_gpumd = ax.plot(x_gp, y_gp, color=COLOR_GPUMD, lw=1.6, alpha=0.95, zorder=3)
    lines_gpumd[0].set_label("GPUMD")

    lines_dft = ax.plot(x_dft_mapped, y_dft, color=COLOR_DFT, lw=1.4, alpha=0.80, zorder=2)
    lines_dft[0].set_label("DFT (MP)")


    # 竖线 + 标签（显示用 K|U）
    for xv in X_disp:
        ax.axvline(xv, color="0.7", lw=0.8, alpha=0.7)

    ax.set_xticks(X_disp)
    ax.set_xticklabels([gamma_label(l) for l in labels_disp])

    ax.set_xlim(x_gp[0], x_gp[-1])
    ax.set_ylim(0, None)
    ax.set_xlabel("Wave vector", fontsize=13, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
    #ax.set_title(f"CaF2 phonon dispersion: GPUMD vs DFT(MP)   path={path}")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
    
    ax.legend()

    beautify_phonon_ax(
        ax,
        X_ticks=X_disp,
        labels=[gamma_label(l) for l in labels_disp],
        xlim=(x_gp[0], x_gp[-1]),
        ylabel=r"Frequencies (cm$^{-1}$)",
        legend_loc="upper right",
        show_grid=False,
    )

    fig.tight_layout()
    fig.savefig(args.outpng, bbox_inches="tight")
    print("Saved:", args.outpng)


if __name__ == "__main__":
    main()
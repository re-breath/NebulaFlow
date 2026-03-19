#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NEP/GPUMD vs DFT: 6-panel Nature-style figure
(a) Energy parity
(b) Force parity (hexbin density)
(c) |ΔE| vs |E_DFT| (binned stats on quantile bins)
(d) |ΔF| vs |F_DFT|
(e) relative |ΔE|/(|E|+eps) vs |E_DFT|
(f) relative |ΔF|/(|F|+eps) vs |F_DFT|

Key changes vs your old version:
- Remove the bottom "Fraction/N" subplots (quantile binning makes them ~flat and uninformative).
- Encode bin sample size by marker size (more informative, cleaner, journal-friendly).
- Blue–green palette (MAE = blue, P95 = bluish-green).
- Save both PNG and PDF (vector) for submission.

该版本为新设计的nep画图方案，可以内部使用，更加详细的查看nep在各个区域内的训练效果
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter, MaxNLocator
from sklearn.metrics import r2_score

# =========================
# Files / output
# =========================
ENERGY_FILE = "energy_merge.out"
FORCE_FILE  = "force_merge.out"
OUT_PNG     = "nep_6panel_nature.png"
OUT_PDF     = "nep_6panel_nature.pdf"

# =========================
# Settings
# =========================
NBINS_E = 18
NBINS_F = 20
EPS_E = 1e-6
EPS_F = 1e-6

# Drop left bins for relative-error plots (|DFT| too small -> ratio blows up)
DROP_LEFT_REL_BINS = 2
REL_YMAX = 10.0

# X-axis capping quantile for profile panels / parity force limits
# Use 98/99 to avoid few outliers stretching axes; use 100 to show all.
CAP_Q = 99.98

# =========================
# Nature-like blue–green palette
# =========================
C_MAE  = "#0072B2"   # blue
C_P95  = "#009E73"   # bluish green
C_REF  = "0.15"      # parity line
C_GRID = "0.92"      # very light grid
C_TXT  = "0.20"

# Bands use same hue as MAE but very light alpha
BAND_ALPHA_WIDE = 0.08   # 5–95
BAND_ALPHA_IQR  = 0.16   # 25–75

# =========================
# Matplotlib rc
# =========================
plt.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 1.15,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "savefig.dpi": 600,
    # PDF font embedding (good for journals)
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# =========================
# Utils
# =========================
def rmse(y_true, y_pred) -> float:
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))

def style_axes(ax, grid=True, grid_axis="y"):
    if grid:
        ax.grid(True, axis=grid_axis, color=C_GRID, lw=0.8)
        ax.set_axisbelow(True)
    ax.tick_params(which="major", top=True, right=True, length=6, width=1.15)
    ax.tick_params(which="minor", top=True, right=True, length=3.2, width=1.0)
    if ax.get_xscale() == "linear":
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    if ax.get_yscale() == "linear":
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    for s in ax.spines.values():
        s.set_linewidth(1.15)

def style_logy(ax):
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=12))
    ax.yaxis.set_minor_formatter(NullFormatter())

def panel_tag(ax, letter: str):
    ax.text(0.02, 0.98, f"({letter})",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=12, fontweight="bold", color=C_TXT)

def load_energy(fp: str):
    E = np.loadtxt(fp)
    if E.ndim == 1:
        E = E.reshape(1, -1)
    # allow (nep, dft) or (..., nep, dft)
    if E.shape[1] == 2:
        nep, dft = E[:, 0], E[:, 1]
    else:
        nep, dft = E[:, -2], E[:, -1]
    return nep, dft

def load_force(fp: str):
    F = np.loadtxt(fp)
    if F.ndim == 1:
        F = F.reshape(1, -1)
    if F.shape[1] < 6:
        raise ValueError(f"force file has too few columns: {F.shape[1]}")
    if F.shape[1] > 6:
        F = F[:, -6:]  # take last 6 columns: Fx,Fy,Fz (nep) then Fx,Fy,Fz (dft)
    return F

def quantile_edges(x, nbins: int):
    x = np.asarray(x, float)
    x = np.clip(x, 1e-12, None)
    qs = np.linspace(0, 1, nbins + 1)
    edges = np.quantile(x, qs)
    edges = np.unique(edges)
    if len(edges) < 3:
        raise ValueError("Too few unique bin edges; distribution degenerate.")
    return edges

def binned_stats(x_scale, err_abs, nbins: int):
    x = np.asarray(x_scale, float)
    e = np.asarray(err_abs, float)
    x = np.clip(x, 1e-12, None)

    edges = quantile_edges(x, nbins)
    nb = len(edges) - 1
    idx = np.digitize(x, edges[1:-1], right=False)

    xc=[]; mae=[]; p95=[]; p25=[]; p75=[]; p05=[]; p95b=[]; cnt=[]
    for b in range(nb):
        m = (idx == b)
        if not np.any(m):
            continue
        xb = x[m]; eb = e[m]
        xc.append(np.median(xb))
        mae.append(np.mean(eb))
        p95v = np.percentile(eb, 95)
        p95.append(p95v)
        p25.append(np.percentile(eb, 25))
        p75.append(np.percentile(eb, 75))
        p05.append(np.percentile(eb, 5))
        p95b.append(p95v)
        cnt.append(int(m.sum()))

    return (edges, np.array(xc), np.array(mae), np.array(p95),
            np.array(p25), np.array(p75), np.array(p05), np.array(p95b), np.array(cnt))

def drop_left_bins(edges, xc, mae, p95, p25, p75, p05, p95b, cnt, n_drop: int):
    if n_drop <= 0:
        return edges, xc, mae, p95, p25, p75, p05, p95b, cnt
    n_drop = min(n_drop, len(cnt)-1)
    return (edges[n_drop:], xc[n_drop:], mae[n_drop:], p95[n_drop:],
            p25[n_drop:], p75[n_drop:], p05[n_drop:], p95b[n_drop:], cnt[n_drop:])

def marker_sizes_from_cnt(cnt, smin=18, smax=60):
    cnt = np.asarray(cnt, float)
    if cnt.size == 0:
        return cnt
    cmin, cmax = float(cnt.min()), float(cnt.max())
    if abs(cmax - cmin) < 1e-12:
        return np.full_like(cnt, (smin + smax)/2.0)
    return smin + (smax - smin) * (cnt - cmin) / (cmax - cmin)

# =========================
# Plotters
# =========================
def plot_parity_energy(ax, x_nep, y_dft):
    # Use scatter; if huge points, consider hexbin like force panel.
    ax.scatter(x_nep, y_dft, s=10, color=C_MAE, alpha=0.25, linewidths=0, rasterized=True)

    lo = min(x_nep.min(), y_dft.min())
    hi = max(x_nep.max(), y_dft.max())
    pad = 0.03 * (hi - lo + 1e-12)
    lo -= pad; hi += pad

    ax.plot([lo, hi], [lo, hi], color=C_REF, lw=1.8)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    ax.set_xlabel("NEP energy (eV/atom)")
    ax.set_ylabel("DFT energy (eV/atom)")
    style_axes(ax, grid=True, grid_axis="both")

    txt = f"RMSE = {rmse(y_dft, x_nep):.3g} eV/atom\nR$^2$ = {r2_score(y_dft, x_nep):.3f}"
    ax.text(0.04, 0.95, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="0.35", alpha=0.90))

def plot_parity_force(ax, x_nep, y_dft, lim):
    # Hexbin with log binning looks more "publication-grade" than raw grayscale
    hb = ax.hexbin(
        x_nep, y_dft,
        gridsize=90,
        cmap="Greys",
        bins="log",     # density contrast
        mincnt=1,
        linewidths=0,
        extent=(-lim, lim, -lim, lim),
        rasterized=True,
    )
    ax.plot([-lim, lim], [-lim, lim], color=C_REF, lw=1.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax.set_xlabel(r"NEP force (eV/$\AA$) [x,y,z merged]")
    ax.set_ylabel(r"DFT force (eV/$\AA$) [x,y,z merged]")
    style_axes(ax, grid=True, grid_axis="both")

    txt = f"RMSE = {rmse(y_dft, x_nep):.3g} eV/$\\AA$\nR$^2$ = {r2_score(y_dft, x_nep):.3f}"
    ax.text(0.04, 0.95, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="0.35", alpha=0.90))

def plot_profile(ax, xc, mae, p95, p25, p75, p05, p95b, cnt,
                 xlabel, ylabel, ycap=None, band_color=C_MAE, cap_note=None):
    # bands
    ax.fill_between(xc, p05, p95b, color=band_color, alpha=BAND_ALPHA_WIDE, lw=0, zorder=1)
    ax.fill_between(xc, p25, p75, color=band_color, alpha=BAND_ALPHA_IQR,  lw=0, zorder=2)

    # marker size encodes bin sample size
    s = marker_sizes_from_cnt(cnt, smin=18, smax=58)

    # MAE
    ax.plot(xc, mae, color=C_MAE, lw=2.2, zorder=3)
    ax.scatter(xc, mae, s=s, color=C_MAE, edgecolor="white", linewidth=0.6, zorder=4)

    # P95
    ax.plot(xc, p95, color=C_P95, lw=2.2, zorder=3)
    ax.scatter(xc, p95, s=s, color=C_P95, marker="s", edgecolor="white", linewidth=0.6, zorder=4)

    style_logy(ax)
    if ycap is not None:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(bottom=max(ymin, 1e-8), top=ycap)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # nature-like: only y-grid, very light
    style_axes(ax, grid=True, grid_axis="y")
    ax.grid(False, axis="x")

    # optional "cap note" (helps reviewers understand axis clipping)
    if cap_note:
        ax.text(0.98, 0.06, cap_note,
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="0.35")

# =========================
# Main
# =========================
def main():
    # ----- load -----
    E_nep, E_dft = load_energy(ENERGY_FILE)
    dE = E_dft - E_nep

    F = load_force(FORCE_FILE)
    Fx_nep, Fy_nep, Fz_nep = F[:, 0], F[:, 1], F[:, 2]
    Fx_dft, Fy_dft, Fz_dft = F[:, 3], F[:, 4], F[:, 5]

    F_nep = np.concatenate([Fx_nep, Fy_nep, Fz_nep])
    F_dft = np.concatenate([Fx_dft, Fy_dft, Fz_dft])
    dF = F_dft - F_nep

    # ----- axis caps -----
    capE = float(np.percentile(np.abs(E_dft), CAP_Q))
    capF = float(np.percentile(np.abs(F_dft), CAP_Q))

    # for profiles, clip x only (errors not clipped)
    xE = np.clip(np.abs(E_dft), 1e-12, capE)
    xF = np.clip(np.abs(F_dft), 1e-12, capF)

    errE_abs = np.abs(dE)
    errF_abs = np.abs(dF)
    errE_rel = errE_abs / (np.abs(E_dft) + EPS_E)
    errF_rel = errF_abs / (np.abs(F_dft) + EPS_F)

    # ----- binned stats (quantile bins) -----
    e_edges, e_xc, e_mae, e_p95, e_p25, e_p75, e_p05, e_p95b, e_cnt = binned_stats(xE, errE_abs, NBINS_E)
    f_edges, f_xc, f_mae, f_p95, f_p25, f_p75, f_p05, f_p95b, f_cnt = binned_stats(xF, errF_abs, NBINS_F)

    er_edges, er_xc, er_mae, er_p95, er_p25, er_p75, er_p05, er_p95b, er_cnt = binned_stats(xE, errE_rel, NBINS_E)
    fr_edges, fr_xc, fr_mae, fr_p95, fr_p25, fr_p75, fr_p05, fr_p95b, fr_cnt = binned_stats(xF, errF_rel, NBINS_F)

    # drop left bins for relative error
    er_edges, er_xc, er_mae, er_p95, er_p25, er_p75, er_p05, er_p95b, er_cnt = drop_left_bins(
        er_edges, er_xc, er_mae, er_p95, er_p25, er_p75, er_p05, er_p95b, er_cnt, DROP_LEFT_REL_BINS
    )
    fr_edges, fr_xc, fr_mae, fr_p95, fr_p25, fr_p75, fr_p05, fr_p95b, fr_cnt = drop_left_bins(
        fr_edges, fr_xc, fr_mae, fr_p95, fr_p25, fr_p75, fr_p05, fr_p95b, fr_cnt, DROP_LEFT_REL_BINS
    )

    cap_note = None
    if CAP_Q < 100:
        cap_note = f"x-axis capped at P{CAP_Q}(|DFT|)"

    # ----- layout: clean 3×2 (no extra subpanels) -----
    fig = plt.figure(figsize=(10.0, 9.0), dpi=600)
    gs = fig.add_gridspec(3, 2, hspace=0.34, wspace=0.25)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    panel_tag(ax_a, "a"); panel_tag(ax_b, "b")
    panel_tag(ax_c, "c"); panel_tag(ax_d, "d")
    panel_tag(ax_e, "e"); panel_tag(ax_f, "f")

    # parity panels
    plot_parity_energy(ax_a, E_nep, E_dft)
    plot_parity_force(ax_b, F_nep, F_dft, lim=capF if CAP_Q < 100 else float(np.max(np.abs(F_dft))))

    # profile panels (bands + MAE/P95 with marker sizes = cnt)
    plot_profile(
        ax_c, e_xc, e_mae, e_p95, e_p25, e_p75, e_p05, e_p95b, e_cnt,
        xlabel=r"$|E_{\mathrm{DFT}}|$ (eV/atom)",
        ylabel=r"$|\Delta E|$ (eV/atom)",
        ycap=None,
        band_color=C_MAE,
        cap_note=cap_note
    )

    plot_profile(
        ax_d, f_xc, f_mae, f_p95, f_p25, f_p75, f_p05, f_p95b, f_cnt,
        xlabel=r"$|F_{\mathrm{DFT}}|$ (eV/$\AA$)",
        ylabel=r"$|\Delta F|$ (eV/$\AA$)",
        ycap=None,
        band_color=C_MAE,
        cap_note=cap_note
    )

    plot_profile(
        ax_e, er_xc, er_mae, er_p95, er_p25, er_p75, er_p05, er_p95b, er_cnt,
        xlabel=r"$|E_{\mathrm{DFT}}|$ (eV/atom)",
        ylabel=r"$|\Delta E|/(|E|+\epsilon)$",
        ycap=REL_YMAX,
        band_color=C_MAE,
        cap_note=cap_note
    )

    plot_profile(
        ax_f, fr_xc, fr_mae, fr_p95, fr_p25, fr_p75, fr_p05, fr_p95b, fr_cnt,
        xlabel=r"$|F_{\mathrm{DFT}}|$ (eV/$\AA$)",
        ylabel=r"$|\Delta F|/(|F|+\epsilon)$",
        ycap=REL_YMAX,
        band_color=C_MAE,
        cap_note=cap_note
    )

    # global legend (top center)
    h1 = plt.Line2D([0], [0], color=C_MAE, lw=2.2, marker="o", ms=5,
                    markerfacecolor=C_MAE, markeredgecolor="white", markeredgewidth=0.6)
    h2 = plt.Line2D([0], [0], color=C_P95, lw=2.2, marker="s", ms=5,
                    markerfacecolor=C_P95, markeredgecolor="white", markeredgewidth=0.6)
    fig.legend([h1, h2], ["MAE", "P95(|Δ|)"], loc="upper center",
               ncol=2, frameon=False, bbox_to_anchor=(0.52, 0.995))

    fig.tight_layout(rect=[0, 0, 1, 0.975])

    fig.savefig(OUT_PNG, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", OUT_PNG)
    print("Saved:", OUT_PDF)

if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from gpyumd.math import running_ave
from typing import Dict, Optional, Tuple
mpl.use("Agg")
import gpyumd.util as util
import pandas as pd
import re



def load_kappa(filename: str = "kappa.out", directory: Optional[str] = None) -> Dict[str, np.ndarray]:
    """
    本地版本 load_kappa：避免 gpyumd 接口变化/弃用警告。
    读取 GPUMD HNEMD 输出 kappa.out（5列：kxi kxo kyi kyo kz）。
    """
    kappa_path = util.get_path(directory, filename)
    data = pd.read_csv(kappa_path, sep=r"\s+", header=None)

    labels = ["kxi", "kxo", "kyi", "kyo", "kz"]
    out: Dict[str, np.ndarray] = {}
    for i, key in enumerate(labels):
        out[key] = data[i].to_numpy(dtype=float)
    return out




AW = 2         
FS = 16        
DT_NS = 0.001  


COLOR_OTHERS = "#23395B" 
COLOR_MEAN   = "#C14953" 

mpl.rcParams.update({
    "font.size": FS,
    "axes.linewidth": AW,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": False,    # 关键：不使用 top ticks
    "ytick.right": False,  # 关键：不使用 right ticks
})


def set_axes(ax: plt.Axes) -> None:
    tl_major, tl_minor, tw = 8, 4, 2
    ax.tick_params(which="major", length=tl_major, width=tw, top=False, right=False)
    ax.tick_params(which="minor", length=tl_minor, width=tw, top=False, right=False)


def last_fraction_mean(x: np.ndarray, frac: float = 0.10) -> float:
    """取数组最后 frac 部分的均值。"""
    n = x.size
    if n == 0:
        return float("nan")
    k = max(1, int(np.ceil(n * frac)))
    return float(np.mean(x[-k:]))


def load_kz_running_average(kappa_file: Path, dt_ns: float) -> Tuple[np.ndarray, np.ndarray]:
    """读取 kappa 文件，返回 t(ns) 与 kz 的 running average 曲线。"""
    kappa = load_kappa(str(kappa_file))
    nsteps = int(kappa["kxi"].shape[0])
    t = np.arange(1, nsteps + 1, dtype=float) * dt_ns
    kz_ra = running_ave(kappa["kz"], t)
    return t, kz_ra


def stack_mean_curve(curves: List[np.ndarray]) -> np.ndarray:
    """
    对多条曲线做逐点平均（要求长度一致；不一致时会裁剪到最短长度）。
    """
    min_len = min(c.size for c in curves)
    trimmed = np.stack([c[:min_len] for c in curves], axis=0)
    return trimmed.mean(axis=0)

def find_kappa_index_files(folder: Path = Path(".")) -> List[Path]:
    """
    自动查找 folder 下所有形如 kappa_%d.out 的文件，并按 %d 数字排序。
    例如：kappa_0.out, kappa_1.out, ...
    """
    rx = re.compile(r"^kappa_(\d+)\.out$")
    hits: List[Tuple[int, Path]] = []

    for p in folder.iterdir():
        if not p.is_file():
            continue
        m = rx.match(p.name)
        if m:
            hits.append((int(m.group(1)), p))

    hits.sort(key=lambda x: x[0])
    return [p for _, p in hits]


def plot_hnemd_kz_with_tail10_mean(
    out_png: str = "hnemd_mul_z.png",
    title: str = "kappa z direction",
    tail_frac: float = 0.10,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=300)
    set_axes(ax)

    # 读取并绘制“其他”多条曲线
    run_files = find_kappa_index_files(Path("."))
    t_ref: Optional[np.ndarray] = None
    kz_curves: List[np.ndarray] = []
    run_tail_means: List[Tuple[str, float]] = []

    for f in run_files:
        if not f.exists():
            print(f"[跳过] 找不到文件: {f}")
            continue

        t, kz_ra = load_kz_running_average(f, DT_NS)
        if t_ref is None:
            t_ref = t

        kz_curves.append(kz_ra)

        tail_mean = last_fraction_mean(kz_ra, tail_frac)
        run_tail_means.append((f.stem, tail_mean))
        print(f"{f.name} | <kz>_last{int(tail_frac*100)}% = {tail_mean:.6g}")

        ax.plot(t, kz_ra, color=COLOR_OTHERS, lw=1.4, alpha=0.35)

    if not kz_curves:
        raise RuntimeError("没有成功读取任何 kappa_i.out 文件，无法作图。")

    kz_mean = stack_mean_curve(kz_curves)

    min_len = kz_mean.size
    t_plot = (t_ref[:min_len] if t_ref is not None else np.arange(1, min_len + 1) * DT_NS)

    ax.plot(t_plot, kz_mean, color=COLOR_MEAN, lw=2.8, alpha=0.95, label="Mean of others")

    mean_tail_mean = last_fraction_mean(kz_mean, tail_frac)
    print(f"[平均曲线]  | <kz>_last{int(tail_frac*100)}% = {mean_tail_mean:.6g}")


    lines = [rf"$\langle \kappa_z \rangle_{{last\,{int(tail_frac*100)}\%}}$ = {mean_tail_mean:.4g}"]

    arr = np.array(run_tail_means, dtype=object)  # 每行 [name, value]
    np.savetxt("Everyhnemd.txt", arr, fmt="%s,%.6f", header="run,mean_kz_last10%", comments="")
    with open("Everyhnemd.txt", "a", encoding="utf-8") as f:
        f.write(f"MeanCurve,{mean_tail_mean:.6f}\n")


    text = "\n".join(lines)
    text +=" W/m/K"
    ax.text(
        0.98, 0.98, text,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=FS - 2,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="none", alpha=0.85),
    )

    ax.set_xlabel("time (ns)")
    ax.set_ylabel(r"$\kappa$ (W/m/K)")
    ax.set_title(title)

    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[完成] 已保存: {out_png}")


if __name__ == "__main__":
    plot_hnemd_kz_with_tail10_mean(
        out_png="hnemd_mul_z.png",
        title="kappa z direction",
        tail_frac=0.10,
    )

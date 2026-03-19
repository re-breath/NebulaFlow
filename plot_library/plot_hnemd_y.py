#!/usr/bin/env python3
import numpy as np
import pandas as pd
from typing import Dict, Optional, Tuple

import matplotlib as mpl
mpl.use("Agg")  # 必须在 pyplot 之前
import matplotlib.pyplot as plt

from gpyumd.math import running_ave
import gpyumd.util as util


def load_kappa_local(filename: str = "kappa.out", directory: Optional[str] = None) -> Dict[str, np.ndarray]:
    """本地版本 load_kappa：避免 gpyumd 接口变化/弃用警告。"""
    kappa_path = util.get_path(directory, filename)
    data = pd.read_csv(kappa_path, sep=r"\s+", header=None)

    labels = ["kxi", "kxo", "kyi", "kyo", "kz"]
    out: Dict[str, np.ndarray] = {}
    for i, key in enumerate(labels):
        out[key] = data[i].to_numpy(dtype=float)
    return out


def style_matplotlib(fontsize: int = 10, axes_linewidth: float = 1.0) -> None:
    """小图（5x4）适配的风格设置。"""
    mpl.rcParams.update({
        "font.size": fontsize,
        "axes.linewidth": axes_linewidth,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })


def estimate_plateau(
    t: np.ndarray,
    y: np.ndarray,
    tail_ns: Optional[float] = None,
    tail_frac: float = 0.2,
    min_points: int = 200,
) -> Tuple[float, float, Tuple[int, int]]:
    """平台段估计：取末尾一段做 mean ± std。"""
    t = np.asarray(t)
    y = np.asarray(y)

    n = y.size
    if n < 5:
        raise ValueError("数据点太少，无法估计平台段。")

    if tail_ns is not None:
        t_start = t[-1] - float(tail_ns)
        i0 = int(np.searchsorted(t, t_start, side="left"))
    else:
        if not (0 < tail_frac <= 1):
            raise ValueError("tail_frac 必须在 (0,1] 之间。")
        i0 = int(max(0, np.floor(n * (1 - tail_frac))))

    i1 = n
    if (i1 - i0) < min_points:
        i0 = max(0, i1 - min_points)

    seg = y[i0:i1]
    mean = float(np.mean(seg))
    std = float(np.std(seg, ddof=1)) if seg.size > 1 else 0.0
    return mean, std, (i0, i1)


def main(
    dt_ns: float = 0.001,
    kappa_file: str = "kappa.out",
    out_png: str = "hnemd.png",
    out_txt: str = "kappa_dealed.txt",
    tail_ns: Optional[float] = None,  # None -> 用 tail_frac
    tail_frac: float = 0.1,
):
    style_matplotlib(fontsize=10, axes_linewidth=1.0)

    kappa = load_kappa_local(filename=kappa_file)

    n = kappa["kyi"].shape[0]
    t = np.arange(1, n + 1) * dt_ns  # ns

    # running average：分别RA后再相加（与你原来 y 脚本一致）
    kyi_ra = running_ave(kappa["kyi"], t)
    kyo_ra = running_ave(kappa["kyo"], t)
    ky_ra = kyi_ra + kyo_ra

    # 如果你希望“先相加再RA”，改成：
    # ky = kappa["kyi"] + kappa["kyo"]
    # ky_ra = running_ave(ky, t)

    ky_mean, ky_std, (i0, i1) = estimate_plateau(t, ky_ra, tail_ns=tail_ns, tail_frac=tail_frac)
    t0, t1 = float(t[i0]), float(t[i1 - 1])

    print(f"ky（平台段平均）: {ky_mean:.2f} ± {ky_std:.2f} W/m/K (std)")
    print(f"平台段: {t0:.3f} ~ {t1:.3f} ns, points={i1 - i0}")

    header = (
        "time(ns)  kappa_y_running_average(W/m/K)\n"
        f"plateau_mean={ky_mean:.6f}  plateau_std={ky_std:.6f}  plateau_range_ns=[{t0:.6f}, {t1:.6f}]"
    )
    np.savetxt(out_txt, np.column_stack([t, ky_ra]), fmt="%.6f", header=header)

    fig, ax = plt.subplots(figsize=(5, 4), dpi=220)

    # 主曲线
    ax.plot(t, ky_ra, linewidth=1.6, color="#6C81A2")

    # 平台阴影 + 平均虚线
    ax.axvspan(t0, t1, alpha=0.12, color="#E79D28")   # 柔和一点的阴影
    ax.axhline(ky_mean, linestyle="--", linewidth=2.0, alpha=0.55, color="#E07A5F")

    ax.set_xlabel("time (ns)", fontweight="bold", fontsize=10)
    ax.set_ylabel(r"$\mathbf{\kappa_y}$ (W/m/K)", fontweight="bold", fontsize=10)


    ax.text(
        0.97, 0.97,
        r"$\mathbf{\kappa_y}$ = " + f"{ky_mean:.2f} $\pm$ {ky_std:.2f} W/m/K",
        transform=ax.transAxes,
        ha="right", va="top",
        fontweight="bold",
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", alpha=0.0),
        color="#544497",
    )


    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

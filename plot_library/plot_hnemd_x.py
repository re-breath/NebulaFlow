#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict, Optional, Tuple
mpl.use("Agg")
from gpyumd.math import running_ave
import gpyumd.util as util


def load_kappa_local(filename: str = "kappa.out", directory: Optional[str] = None) -> Dict[str, np.ndarray]:
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


def style_matplotlib(fontsize: int = 12, axes_linewidth: float = 2.5) -> None:
    """统一画图风格：更接近你给的那张图。"""
    mpl.rcParams.update({
        "font.size": fontsize,
        "axes.linewidth": axes_linewidth,
        "xtick.direction": "in",
        "ytick.direction": "in",
    })


def estimate_plateau(
    t: np.ndarray,
    y: np.ndarray,
    tail_ns: Optional[float] = 2.0,
    tail_frac: float = 0.2,
    min_points: int = 200,
) -> Tuple[float, float, Tuple[int, int]]:
    """
    平台段估计：取末尾一段做 mean ± std。
    - 优先用 tail_ns（末尾多少 ns），若为 None 则用 tail_frac（末尾比例）
    """
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
    out_png: str = "hnemd_x.png",
    out_txt: str = "kappa_dealed.txt",
    tail_ns: Optional[float] = None,   # 末尾平台段长度（ns）；不想用固定时长就设为 None
    tail_frac: float = 0.1,           # tail_ns=None 时生效：取末尾比例
):
    style_matplotlib(fontsize=10, axes_linewidth=1.0)

    kappa = load_kappa_local(filename=kappa_file)

    # 时间轴
    n = kappa["kxi"].shape[0]
    t = np.arange(1, n + 1) * dt_ns  # ns

    # running average（按你当前做法：分别RA后再相加）
    kxi_ra = running_ave(kappa["kxi"], t)
    kxo_ra = running_ave(kappa["kxo"], t)
    kx_ra = kxi_ra + kxo_ra

    # 如果你想“先相加再RA”，用下面两行替代上面三行：
    # kx = kappa["kxi"] + kappa["kxo"]
    # kx_ra = running_ave(kx, t)

    # 平台估计
    kx_mean, kx_std, (i0, i1) = estimate_plateau(t, kx_ra, tail_ns=tail_ns, tail_frac=tail_frac)
    t0, t1 = float(t[i0]), float(t[i1 - 1])

    print(f"kx（平台段平均）: {kx_mean:.2f} ± {kx_std:.2f} W/m/K (std)")
    print(f"平台段: {t0:.3f} ~ {t1:.3f} ns, points={i1 - i0}")

    # 保存 txt：两列 time, kappa_x_running_average
    header = (
        "time(ns)  kappa_x_running_average(W/m/K)\n"
        f"plateau_mean={kx_mean:.6f}  plateau_std={kx_std:.6f}  plateau_range_ns=[{t0:.6f}, {t1:.6f}]"
    )
    np.savetxt(out_txt, np.column_stack([t, kx_ra]), fmt="%.6f", header=header)

    # 画图（风格尽量贴你那张图）
    fig, ax = plt.subplots(figsize=(5, 4), dpi=220)

    ax.plot(t, kx_ra, linewidth=1.6, color="#6C81A2")

    # 平台阴影 + 平均虚线（虚线偏红一点更醒目，但不指定颜色也可以；你要完全默认就删 color）
    ax.axvspan(t0, t1, alpha=0.15, color="#12F4A1")
    ax.axhline(kx_mean, linestyle="--", linewidth=2, alpha=0.6, color="#737CF6")

    ax.set_xlabel("time (ns)", fontweight="bold", fontsize=10)
    ax.set_ylabel(r"$\mathbf{\kappa_x}$ (W/m/K)", fontweight="bold", fontsize=10)

    ax.text(
        0.97, 0.97,
        r"$\mathbf{\kappa_x}$ = "
        f"{kx_mean:.2f} $\pm$ {kx_std:.2f} W/m/K",
        transform=ax.transAxes,
        ha="right", va="top",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", alpha=0.0)  # 不要框就 alpha=0
    )

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

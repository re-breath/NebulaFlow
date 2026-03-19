#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("Agg")

from gpyumd.load import load_kappa
from gpyumd.math import running_ave
import gpyumd.util as util
import pandas as pd
from typing import Dict

def style_matplotlib(fontsize=2, axes_linewidth=1, tick_len=3, tick_w=1.5, minor_len=2):
    """统一设置绘图风格（更像论文图）。"""
    mpl.rcParams.update({
        "font.size": fontsize,
        "axes.linewidth": axes_linewidth,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": tick_len,
        "ytick.major.size": tick_len,
        "xtick.major.width": tick_w,
        "ytick.major.width": tick_w,
        "xtick.minor.size": minor_len,
        "ytick.minor.size": minor_len,
        "xtick.minor.width": tick_w,
        "ytick.minor.width": tick_w,
    })


def compute_running_averages(kappa, t, keys):
    """对指定 keys 计算 running average，并写回到 kappa 字典中。"""
    for k in keys:
        if k in kappa:
            kappa[f"{k}_ra"] = running_ave(kappa[k], t)


def estimate_plateau_kappa(t, kz_ra, tail_frac=0.2, tail_ns=None, min_points=50):
    """
    用 running average 的末尾“平台段”估计 kappa（更合理：不是取最后一个点）。
    - tail_ns: 取最后多少 ns（优先级更高）
    - tail_frac: 否则取最后多少比例（0~1）
    返回：mean, std, (i0, i1)
    """
    t = np.asarray(t)
    kz_ra = np.asarray(kz_ra)

    n = kz_ra.size
    if n < 5:
        raise ValueError("数据点太少，无法估计稳态段。")

    if tail_ns is not None:
        t_start = t[-1] - float(tail_ns)
        i0 = int(np.searchsorted(t, t_start, side="left"))
    else:
        tail_frac = float(tail_frac)
        if not (0 < tail_frac <= 1):
            raise ValueError("tail_frac 必须在 (0, 1] 之间")
        i0 = int(max(0, np.floor(n * (1 - tail_frac))))

    i1 = n

    # 保底：末段点数太少就向前扩一点，避免 std 无意义
    if (i1 - i0) < min_points:
        i0 = max(0, i1 - min_points)

    seg = kz_ra[i0:i1]
    mean = float(np.mean(seg))
    std = float(np.std(seg, ddof=1)) if seg.size > 1 else 0.0
    return mean, std, (i0, i1)

def load_kappa(filename: str = "kappa.out", directory: str = None) -> Dict[str, np.ndarray]:
    """
    Loads data from kappa.out GPUMD output file which contains HNEMD kappa.

    Args:
        filename: The kappa data file
        directory: Directory containing kappa data file

    Returns:
        Dictionary with keys corresponding to the columns in 'kappa.out'.
         Units are [kxi, kxo, kyi, kyo, kz -> W(m^-1)(K^-1)]
    """
    kappa_path = util.get_path(directory, filename)
    data = pd.read_csv(kappa_path, sep=r"\s+", header=None)
    labels = ['kxi', 'kxo', 'kyi', 'kyo', 'kz']
    out = dict()
    for i, key in enumerate(labels):
        out[key] = data[i].to_numpy(dtype='float')
    return out


def main(
    dt_ns=0.001,
    out_png="hnemd.png",
    out_txt="kappa_dealed.txt",
    tail_ns=None,      # 平台段取最后多少 ns（改成 None 则用 tail_frac）
    tail_frac=0.1,    # 平台段取最后多少比例（当 tail_ns=None 时生效）
):
    style_matplotlib(fontsize=12, axes_linewidth=2.0)

    kappa = load_kappa()

    if "kz" not in kappa:
        raise KeyError("kappa 数据中找不到 'kz'。请确认 load_kappa() 结果包含 kz。")

    n = kappa["kz"].shape[0]
    t = np.arange(1, n + 1) * dt_ns  # ns

    # 需要的话可加更多分量
    compute_running_averages(kappa, t, keys=["kxi", "kxo", "kyi", "kyo", "kz"])
    kz_ra = kappa.get("kz_ra", kappa["kz"])

    kappa_mean, kappa_std, (i0, i1) = estimate_plateau_kappa(
        t, kz_ra, tail_frac=tail_frac, tail_ns=tail_ns
    )
    t0, t1 = float(t[i0]), float(t[i1 - 1])

    print(f"z方向热导率（平台段平均）: {kappa_mean:.2f} ± {kappa_std:.2f} W/m/K  (std)")
    print(f"平台段区间: {t0:.3f} ~ {t1:.3f} ns, 点数: {i1 - i0}")

    # 保存两列：t, kz_ra，并把平台信息写进 header
    header = (
        "time(ns)  kappa_z_running_average(W/m/K)\n"
        f"plateau_mean={kappa_mean:.6f}  plateau_std={kappa_std:.6f}  "
        f"plateau_range_ns=[{t0:.6f}, {t1:.6f}]"
    )
    np.savetxt(out_txt, np.column_stack([t, kz_ra]), fmt="%.6f", header=header)

    # 作图
    fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
    ax.plot(t, kz_ra, linewidth=1.2, color="#587798")

    # 标出平台段（阴影） + 平均值虚线
    ax.axvspan(t0, t1, alpha=0.15, color="#DD901C")
    ax.axhline(kappa_mean, linewidth=2.0, linestyle="--", color="#F4AA9C")

    ax.set_xlabel("time (ns)",fontweight="bold",fontsize=10)
    ax.set_ylabel(r"$\mathbf{{\kappa}_z}$ (W/m/K)",fontweight="bold",fontsize=10)
    #ax.set_title(r"HNEMD: $\mathbf{{\kappa}_z}$ (running average)",fontweight="bold")
    ax.tick_params(axis='both', labelsize=8)

    # 右上角标注
    text = (
        r"$\mathbf{{\kappa}_z}$ = " 
        f"{kappa_mean:.2f} ± {kappa_std:.2f} W/m/K\n"
    )
    ax.text(
        0.97, 0.97, text,
        transform=ax.transAxes,
        ha="right", va="top",
        #bbox=dict(boxstyle="round,pad=0.35", alpha=0.85, color="#B8EFBD",edgecolor="#5779E7"),
        color="#865D8A",
        fontsize=6,
        fontweight="bold"
    )


    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

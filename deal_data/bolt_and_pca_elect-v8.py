"""
bolt_and_pca_elect-v7.py

变更说明（相对 v6）：
  - 代码整理：删除重复导入、未使用的模块和函数、注释掉的旧代码
  - 功能不变：多峰自动检测 + 分峰 Boltzmann + PCA 最远点联合采样
"""

# ── 标准库 ────────────────────────────────────────────────────────────
import re
import os
import subprocess
import datetime
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

# ── 第三方库 ──────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")                        # 非交互后端，必须在 plt 之前设置

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from sklearn.decomposition import PCA
import ase.io

# ── rich 导入（带降级回退） ───────────────────────────────────────────
try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.text import Text
    from rich import box
    from rich.progress import (
        Progress, SpinnerColumn, TextColumn,
        BarColumn, TimeElapsedColumn,
    )
    _RICH = True
    console = Console()
except ImportError:
    _RICH = False
    console = None

# ── 全局绘图风格 ──────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.9,
    "axes.grid":         True,
    "grid.color":        "#E0E0E0",
    "grid.linewidth":    0.6,
    "grid.linestyle":    "--",
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "figure.dpi":        150,
    "savefig.dpi":       200,
    "savefig.bbox":      "tight",
})

# ══════════════════════════════════════════════════════════════════════
#  NEP 头部关键字集合（全局常量，供 from_nep_txt 和 payload 统计使用）
# ══════════════════════════════════════════════════════════════════════

_NEP_HEADER_KEYS: frozenset[str] = frozenset({
    # 结构参数（nep.txt 中常见）
    "zbl",            # ZBL 势参数（可选行，如 "zbl 1 2"）
    "cutoff",         # 截断半径
    "n_max",          # 径向/角向最大展开阶数
    "basis_size",     # 基组大小
    "l_max",          # 角动量最大值
    "ann",            # 神经网络结构（旧格式）
    "neuron",         # 神经网络结构（新格式，与 ann 等价）
    # 训练控制参数（nep.in 中出现；nep.txt 中通常不出现，保留作容错）
    "lambda_1", "lambda_e", "lambda_f", "lambda_v", "lambda_shear",
    "batch", "epoch", "population", "prediction", "output_descriptor",
    "force_delta", "version", "type",
})


# ══════════════════════════════════════════════════════════════════════
#  终端输出辅助
# ══════════════════════════════════════════════════════════════════════

def _print(msg="", **kw):
    if _RICH:
        console.print(msg, **kw)
    else:
        print(msg)


def _rule(title=""):
    if _RICH:
        console.print()
        console.print(
            f"[bold #39FB90][[/bold #39FB90]"
            f"[bold #ABE2F9]{title}[/bold #ABE2F9]"
            f"[bold #39FB90]][/bold #39FB90]"
        )
    else:
        print(f"\n——>  {title}")


# ══════════════════════════════════════════════════════════════════════
#  结构化日志
# ══════════════════════════════════════════════════════════════════════

class StructuredLogger:
    """将运行过程写入带时间戳的日志文件。"""

    def __init__(self, path: str = "bolt_pca_run.log"):
        self.path = path
        self._f = open(path, "w", encoding="utf-8")
        self._write_header()

    def _ts(self) -> str:
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def _write_header(self):
        self._f.write("=" * 70 + "\n")
        self._f.write("  bolt_and_pca_elect 运行日志\n")
        self._f.write(f"  开始时间：{self._ts()}\n")
        self._f.write("=" * 70 + "\n\n")
        self._f.flush()

    def section(self, title: str):
        self._f.write(f"\n{'─'*70}\n")
        self._f.write(f"[{self._ts()}]  {title}\n")
        self._f.write(f"{'─'*70}\n")
        self._f.flush()

    def info(self, msg: str):
        self._f.write(f"[{self._ts()}]  INFO  {msg}\n")
        self._f.flush()

    def data(self, label: str, value):
        self._f.write(f"  {label:<35} {value}\n")
        self._f.flush()

    def write_raw(self, text: str):
        self._f.write(text + "\n")
        self._f.flush()

    def close(self, success: bool = True):
        self._f.write(f"\n{'='*70}\n")
        status = "成功" if success else "异常退出"
        self._f.write(f"  结束时间：{self._ts()}   状态：{status}\n")
        self._f.write("=" * 70 + "\n")
        self._f.close()


# ══════════════════════════════════════════════════════════════════════
#  NEP 势函数信息解析
# ══════════════════════════════════════════════════════════════════════

@dataclass
class NEPPotentialInfo:
    path: Path
    version: int
    n_types: int
    elements: List[str]
    cutoff_radial: float
    cutoff_angular: float
    max_neighbors_radial: Optional[int]
    max_neighbors_angular: Optional[int]
    n_max_radial: int
    n_max_angular: int
    basis_size_radial: int
    basis_size_angular: int
    l_max_3body: int
    l_max_4body: int
    l_max_5body: int
    n_neurons: int
    ann_extra: List[str]
    n_descriptor: int
    n_descriptor_parameters: int
    n_nn_parameters: int
    n_payload_values: int
    expected_payload_values: int

    @classmethod
    def from_nep_txt(cls, filename: str | Path) -> "NEPPotentialInfo":
        path  = Path(filename)
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

        # 去掉行内注释与空行
        clean = [l.split("#", 1)[0].strip() for l in lines]
        clean = [l for l in clean if l]

        if len(clean) < 5:
            raise ValueError("文件行数不足，不像有效的 NEP 势函数文件。")

        # ── 首行：版本 + 元素 ──────────────────────────────────────────
        first = clean[0].split()
        m = re.search(r"\d+", first[0])
        if not m:
            raise ValueError(f"无法从首行解析 NEP 版本：{clean[0]}")
        version  = int(m.group())
        n_types  = int(first[1])
        elements = first[2 : 2 + n_types]
        if len(elements) != n_types:
            raise ValueError(
                f"期望 {n_types} 个元素，实际得到 {len(elements)}：{clean[0]}"
            )

        # ── 按关键字建立索引（跳过首行；未知关键字静默跳过） ───────────
        kw_map: dict[str, list[str]] = {}
        zbl_tok: list[str] | None = None

        for line in clean[1:]:
            tok = line.split()
            if not tok:
                continue
            key = tok[0].lower()

            # 新增：单独记录 zbl，避免把 nep4_zbl 误识别成 zbl
            if key == "zbl" and zbl_tok is None:
                zbl_tok = tok
                continue

            # 只索引已知头部关键字，且每个关键字只取第一次出现
            if key in _NEP_HEADER_KEYS and key not in kw_map:
                kw_map[key] = tok

        def _require(key: str) -> list[str]:
            if key not in kw_map:
                raise ValueError(
                    f"nep.txt 中未找到必需关键字行：'{key}'  "
                    f"（已检测到的关键字：{sorted(kw_map.keys())}）"
                )
            return kw_map[key]

        tok            = _require("cutoff")
        cutoff_radial  = float(tok[1])
        cutoff_angular = float(tok[2])

        max_nb_r = int(float(tok[3])) if len(tok) >= 5 else None
        max_nb_a = int(float(tok[4])) if len(tok) >= 5 else None

        tok      = _require("n_max")
        n_max_r  = int(tok[1])
        n_max_a  = int(tok[2])

        tok     = _require("basis_size")
        basis_r = int(tok[1])
        basis_a = int(tok[2])

        tok    = _require("l_max")
        l_vals = [int(x) for x in tok[1:]]
        while len(l_vals) < 3:
            l_vals.append(0)          # 4-body / 5-body 默认为 0
        l3, l4, l5 = l_vals[:3]

        ann_tok = kw_map.get("ann") or kw_map.get("neuron")
        if ann_tok is None:
            raise ValueError("nep.txt 中未找到 'ann' 或 'neuron' 行。")
        n_neurons = int(ann_tok[1])
        ann_extra = ann_tok[2:]

        zbl = None
        if zbl_tok is not None:
            if len(zbl_tok) < 2:
                raise ValueError(f"zbl 行格式不正确：{' '.join(zbl_tok)}")
            try:
                zbl = int(zbl_tok[-1])   # 例如 zbl 1 2 -> 取 2
            except ValueError as e:
                raise ValueError(f"zbl 行最后一列不是整数：{' '.join(zbl_tok)}") from e

        n_descriptor = (
            (n_max_r + 1)
            + (n_max_a + 1) * l3
            + ((n_max_a + 1) if l4 > 0 else 0)
            + ((n_max_a + 1) if l5 > 0 else 0)
        )
        n_desc_params = n_types ** 2 * (
            (n_max_r + 1) * (basis_r + 1)
            + (n_max_a + 1) * (basis_a + 1)
        )
        if version == 3:
            n_nn_params = (n_descriptor + 2) * n_neurons + 1
        elif version >= 4:
            n_nn_params = (n_descriptor + 2) * n_neurons * n_types + 1
        else:
            raise ValueError(f"仅支持 NEP3/NEP4，当前版本 NEP{version}")

        payload: list[float] = []
        for line in clean[1:]:
            tok0 = line.split()
            if not tok0:
                continue
            key0 = tok0[0].lower()
            if key0 in _NEP_HEADER_KEYS or key0 == "zbl":
                continue
            for t in tok0:
                try:
                    payload.append(float(t))
                except ValueError:
                    pass

        expected = n_desc_params + n_nn_params + n_descriptor

        obj = cls(
            path=path, version=version, n_types=n_types, elements=elements,
            cutoff_radial=cutoff_radial, cutoff_angular=cutoff_angular,
            max_neighbors_radial=max_nb_r, max_neighbors_angular=max_nb_a,
            n_max_radial=n_max_r, n_max_angular=n_max_a,
            basis_size_radial=basis_r, basis_size_angular=basis_a,
            l_max_3body=l3, l_max_4body=l4, l_max_5body=l5,
            n_neurons=n_neurons, ann_extra=ann_extra,
            n_descriptor=n_descriptor,
            n_descriptor_parameters=n_desc_params,
            n_nn_parameters=n_nn_params,
            n_payload_values=len(payload),
            expected_payload_values=expected,
        )

        obj.zbl = zbl
        return obj

    def to_nep_in(self, output_descriptor: int = 1, prediction: int = 1) -> str:
        lines = [
            f"version         {self.version}",
            f"type            {self.n_types} {' '.join(self.elements)}",
            f"cutoff          {self.cutoff_radial:g}  {self.cutoff_angular:g}",
            f"n_max           {self.n_max_radial}  {self.n_max_angular}",
            f"basis_size      {self.basis_size_radial}  {self.basis_size_angular}",
            f"l_max           {self.l_max_3body}  {self.l_max_4body}  {self.l_max_5body}",
            f"neuron          {self.n_neurons}",
        ]

        # 如果存在 zbl，就写入 nep.in
        zbl = getattr(self, "zbl", None)
        if zbl is not None:
            lines.append(f"zbl             {zbl}")

        lines.extend([
            "",
            f"prediction      {prediction}",
            f"output_descriptor {output_descriptor}",
        ])
        return "\n".join(lines) + "\n"


    def write_nep_in(self, filename: str | Path = "nep.in",
                     output_descriptor: int = 1) -> None:
        Path(filename).write_text(
            self.to_nep_in(output_descriptor=output_descriptor),
            encoding="utf-8",
        )

    def summary(self) -> str:
        return (
            f"NEP{self.version} potential: {' '.join(self.elements)}\n"
            f"cutoff radial/angular: {self.cutoff_radial}, {self.cutoff_angular}\n"
            f"max neighbors radial/angular: "
            f"{self.max_neighbors_radial}, {self.max_neighbors_angular}\n"
            f"n_max radial/angular: {self.n_max_radial}, {self.n_max_angular}\n"
            f"basis_size radial/angular: "
            f"{self.basis_size_radial}, {self.basis_size_angular}\n"
            f"l_max 3/4/5-body: "
            f"{self.l_max_3body}, {self.l_max_4body}, {self.l_max_5body}\n"
            f"hidden neurons: {self.n_neurons}\n"
            f"descriptor dimension: {self.n_descriptor}\n"
            f"descriptor parameters: {self.n_descriptor_parameters}\n"
            f"NN parameters: {self.n_nn_parameters}\n"
            f"payload values in file: {self.n_payload_values}\n"
            f"expected payload values: {self.expected_payload_values}"
        )

    def print_rich(self):
        if not _RICH:
            print(self.summary())
            return

        title_text = Text()
        title_text.append(f"NEP{self.version}", style="bold bright_cyan")
        title_text.append("   元素：", style="bold white")
        title_text.append(" ".join(self.elements), style="bold yellow")
        console.print(title_text)

        tbl = Table(box=box.ROUNDED, show_header=True,
                    header_style="bold white", min_width=52)
        tbl.add_column("参数", style="white", no_wrap=True)
        tbl.add_column("径向 (Radial)",  justify="right", style="white")
        tbl.add_column("角向 (Angular)", justify="right", style="white")
        tbl.add_row("截断半径 cutoff  (Å)",
                    f"{self.cutoff_radial:g}", f"{self.cutoff_angular:g}")
        tbl.add_row("最大近邻数 max_neighbors",
                    str(self.max_neighbors_radial), str(self.max_neighbors_angular))
        tbl.add_row("n_max",       str(self.n_max_radial),    str(self.n_max_angular))
        tbl.add_row("basis_size",  str(self.basis_size_radial), str(self.basis_size_angular))
        console.print(tbl)

        tbl2 = Table(box=box.SIMPLE, show_header=False, min_width=52)
        tbl2.add_column("key", style="white")
        tbl2.add_column("val", style="bold white")
        tbl2.add_row("l_max  3/4/5-body",
                     f"{self.l_max_3body}  /  {self.l_max_4body}  /  {self.l_max_5body}")
        tbl2.add_row("隐藏神经元 neurons",    str(self.n_neurons))
        tbl2.add_row("描述符维度 n_descriptor", str(self.n_descriptor))
        tbl2.add_row("描述符参数数",           str(self.n_descriptor_parameters))
        tbl2.add_row("NN 参数数",             str(self.n_nn_parameters))

        ok            = self.n_payload_values == self.expected_payload_values
        payload_style = "bold bright_green" if ok else "bold bright_red"
        payload_icon  = "✔" if ok else "✘"
        tbl2.add_row(
            "文件参数值 / 预期",
            Text(
                f"{self.n_payload_values}  /  {self.expected_payload_values}"
                f"  {payload_icon}",
                style=payload_style,
            ),
        )
        console.print(tbl2)


# ══════════════════════════════════════════════════════════════════════
#  峰检测数据类
# ══════════════════════════════════════════════════════════════════════

@dataclass
class EnergyPeakResult:
    """
    峰检测与分配结果。

    属性
    ────
    n_peaks       : 峰数量
    peak_centers  : 各峰能量位置 (eV/atom)
    boundaries    : 各峰左右边界，boundaries[k]~boundaries[k+1] 为第 k 峰区间
    labels        : 每个构型的峰编号（0-indexed）
    peak_counts   : 各峰构型数
    peak_indices  : peak_indices[k] = 第 k 峰所有构型在原数组中的索引
    kde_x, kde_y  : KDE 曲线数据（供可视化使用）
    """
    n_peaks:      int
    peak_centers: np.ndarray
    boundaries:   np.ndarray
    labels:       np.ndarray
    peak_counts:  List[int]
    peak_indices: List[np.ndarray]
    kde_x:        np.ndarray
    kde_y:        np.ndarray


# ══════════════════════════════════════════════════════════════════════
#  峰检测
# ══════════════════════════════════════════════════════════════════════

def find_energy_peaks(
    ave_energy: np.ndarray,
    sensitivity: int = 5,
    bw_method: str | float = "silverman",
    bw_factor: Optional[float] = None,
    n_grid: int = 2000,
    min_prominence: Optional[float] = None,
    min_peak_distance_eV: Optional[float] = None,
    merge_small_peak_ratio: float = 0.05,
) -> EnergyPeakResult:
    """
    对平均能量分布做 KDE，自动检测多峰并将每个构型归入对应峰区间。

    Parameters
    ----------
    sensitivity : 1~10，峰检测灵敏度（越大越灵敏，默认 5）
    merge_small_peak_ratio : 构型数占比低于该值的峰并入最近相邻峰
    """
    ave_energy = np.asarray(ave_energy, dtype=float)
    n_configs  = len(ave_energy)

    sensitivity = int(np.clip(sensitivity, 1, 10))
    if bw_factor is None:
        bw_factor = 10 ** (-(sensitivity - 1) / 9.0)
    _prom_coef = max(0.01, 0.08 - sensitivity * 0.006)
    _dist_coef = max(0.01, 0.06 - sensitivity * 0.004)

    lo_pct = np.percentile(ave_energy, 0.2)
    hi_pct = np.percentile(ave_energy, 99.8)
    e_span = hi_pct - lo_pct
    margin = e_span * 0.06
    x_grid = np.linspace(lo_pct - margin, hi_pct + margin, n_grid)

    if isinstance(bw_method, str):
        kde_ref   = gaussian_kde(ave_energy, bw_method=bw_method)
        bw_scaled = kde_ref.factor * bw_factor
        kde       = gaussian_kde(ave_energy, bw_method=bw_scaled)
    else:
        kde = gaussian_kde(ave_energy, bw_method=bw_method)
    y_grid = kde(x_grid)

    if min_prominence is None:
        min_prominence = y_grid.max() * _prom_coef
    if min_peak_distance_eV is None:
        min_peak_distance_eV = e_span * _dist_coef

    grid_spacing      = (x_grid[-1] - x_grid[0]) / (n_grid - 1)
    min_peak_dist_pts = max(1, int(min_peak_distance_eV / grid_spacing))

    peak_pts, _ = find_peaks(y_grid, prominence=min_prominence,
                              distance=min_peak_dist_pts)
    if len(peak_pts) == 0:
        peak_pts = np.array([np.argmax(y_grid)])

    peak_centers = x_grid[peak_pts]

    boundaries = [ave_energy.min()]
    for i in range(len(peak_pts) - 1):
        lo_pt, hi_pt = peak_pts[i], peak_pts[i + 1]
        valley_pt    = lo_pt + np.argmin(y_grid[lo_pt: hi_pt + 1])
        boundaries.append(float(x_grid[valley_pt]))
    boundaries.append(ave_energy.max())
    boundaries = np.array(boundaries)

    labels = np.full(n_configs, -1, dtype=int)
    for k in range(len(peak_centers)):
        lo, hi = boundaries[k], boundaries[k + 1]
        mask = (ave_energy >= lo) & (
            (ave_energy <= hi) if k == len(peak_centers) - 1 else (ave_energy < hi)
        )
        labels[mask] = k
    for idx in np.where(labels == -1)[0]:
        labels[idx] = int(np.argmin(np.abs(peak_centers - ave_energy[idx])))

    peak_counts = [int((labels == k).sum()) for k in range(len(peak_centers))]

    if len(peak_centers) > 1 and merge_small_peak_ratio > 0:
        small       = {k for k, cnt in enumerate(peak_counts)
                       if cnt / n_configs < merge_small_peak_ratio}
        valid_peaks = set(range(len(peak_centers))) - small
        if small:
            for k in small:
                neighbors = [n for n in [k - 1, k + 1] if n in valid_peaks]
                if not neighbors:
                    valid_peaks.add(k)
                    continue
                nearest = min(neighbors,
                              key=lambda n: abs(peak_centers[n] - peak_centers[k]))
                labels[labels == k] = nearest

            unique_peaks = sorted(valid_peaks)
            remap        = {old: new for new, old in enumerate(unique_peaks)}
            labels       = np.array([remap[l] for l in labels])
            peak_centers = peak_centers[unique_peaks]

            bnds_new = [boundaries[0]]
            for i in range(len(unique_peaks) - 1):
                bnds_new.append(boundaries[unique_peaks[i] + 1])
            bnds_new.append(boundaries[-1])
            boundaries  = np.array(bnds_new)
            peak_counts = [int((labels == k).sum())
                           for k in range(len(peak_centers))]

    peak_indices = [np.where(labels == k)[0] for k in range(len(peak_centers))]

    return EnergyPeakResult(
        n_peaks=len(peak_centers), peak_centers=peak_centers,
        boundaries=boundaries, labels=labels,
        peak_counts=peak_counts, peak_indices=peak_indices,
        kde_x=x_grid, kde_y=y_grid,
    )


def _peak_custom_snr(ave_energy: np.ndarray, result: EnergyPeakResult) -> float:
    """SNR = 相邻峰平均间距 / 构型数加权平均峰内标准差。"""
    if result.n_peaks < 2:
        return 0.0
    c = result.peak_centers
    n = np.array(result.peak_counts, dtype=float)
    mean_sep = float(np.mean(np.diff(c)))
    stds = np.array([
        np.std(ave_energy[result.peak_indices[k]], ddof=0)
        if len(result.peak_indices[k]) > 1 else 0.0
        for k in range(result.n_peaks)
    ])
    weighted_std = float(np.dot(n, stds) / n.sum())
    return float(mean_sep / weighted_std) if weighted_std > 0 else float("inf")

def auto_find_n_peaks(
    ave_energy: np.ndarray,
    metric: str = "snr",
    merge_small_peak_ratio: float = 0.05,
    plot_diagnosis: bool = True,
    save_diagnosis: Optional[str] = "peak_auto_diagnosis.png",
    show_diagnosis: bool = False,
    verbose: bool = True,
    peak_penalty_rate: float = 0.05,
    **peak_kwargs,
) -> EnergyPeakResult:
    """
    自动扫描 sensitivity=1~10，用惩罚后 SNR（或 CH）选出最优分峰方案。

    Parameters
    ----------
    metric : "snr"（默认）| "calinski_harabasz"
    peak_penalty_rate : float
        峰数惩罚系数，仅在 metric="snr" 时生效。
        惩罚公式：SNR_adj(K) = SNR(K) - peak_penalty_rate × (K-1) × SNR_max
        默认 0.05，即每多一个峰扣除 SNR_max 的 5%。
        典型取值建议：
          0.02 → 轻惩罚，允许较多峰
          0.05 → 中等惩罚（默认），SNR 相近时倾向更少峰
          0.10 → 重惩罚，显著偏向少峰
    """
    try:
        from sklearn.metrics import calinski_harabasz_score
        _SKLEARN = True
    except ImportError:
        _SKLEARN = False
        if metric == "calinski_harabasz":
            print("[警告] scikit-learn 未安装，自动切换为 metric='snr'")
            metric = "snr"

    ave_energy = np.asarray(ave_energy, dtype=float)
    X = ave_energy.reshape(-1, 1)

    results_by_k: dict = {}
    scan_log: list[tuple] = []

    # ── 扫描所有灵敏度 ────────────────────────────────────────────────
    for sens in range(1, 11):
        result = find_energy_peaks(
            ave_energy, sensitivity=sens,
            merge_small_peak_ratio=merge_small_peak_ratio,
            **peak_kwargs,
        )
        k   = result.n_peaks
        snr = _peak_custom_snr(ave_energy, result)
        ch  = (float(calinski_harabasz_score(X, result.labels))
               if k >= 2 and _SKLEARN else 0.0)
        scan_log.append((sens, k, snr, ch))
        score = snr if metric == "snr" else ch
        if k not in results_by_k or score > results_by_k[k][1]:
            results_by_k[k] = (result, score, snr, ch)

    ks   = sorted(results_by_k.keys())
    snrs = {k: results_by_k[k][2] for k in ks}
    chs  = {k: results_by_k[k][3] for k in ks}

    # ── 计算惩罚后 SNR ────────────────────────────────────────────────
    # SNR_adj(K) = SNR(K) - peak_penalty_rate × (K-1) × SNR_max
    # 物理含义：每多一个峰，需要用 peak_penalty_rate × SNR_max 的 SNR "代价" 来买
    snr_max = max(snrs.values()) if snrs else 1.0
    penalty_per_peak = peak_penalty_rate * snr_max
    adjusted_snrs = {
        k: snrs[k] - penalty_per_peak * max(0, k - 1)
        for k in ks
    }

    # ── 终端扫描摘要 ──────────────────────────────────────────────────
    if verbose:
        print("\n── sensitivity 扫描摘要 " + "─" * 48)
        if metric == "snr":
            print(f"  {'sens':>5}  {'K':>4}  {'SNR':>8}  {'SNR_adj':>10}  {'CH':>10}")
            print("  " + "-" * 46)
            prev_k = None
            for sens, k, snr, ch in scan_log:
                marker = " ←" if k != prev_k else ""
                ch_str = f"{ch:>10.1f}" if ch > 0 else f"{'—':>10}"
                print(
                    f"  {sens:>5}  {k:>4}  {snr:>8.3f}"
                    f"  {adjusted_snrs[k]:>10.3f}  {ch_str}{marker}"
                )
                prev_k = k
        else:
            print(f"  {'sens':>5}  {'K':>4}  {'SNR':>8}  {'CH':>10}")
            print("  " + "-" * 34)
            prev_k = None
            for sens, k, snr, ch in scan_log:
                marker = " ←" if k != prev_k else ""
                ch_str = f"{ch:>10.1f}" if ch > 0 else f"{'—':>10}"
                print(f"  {sens:>5}  {k:>4}  {snr:>8.3f}  {ch_str}{marker}")
                prev_k = k
        print("  " + "─" * 52)

    # ── 选择最优 K ────────────────────────────────────────────────────
    k_candidates = [k for k in ks if k >= 2]
    if not k_candidates:
        best_k = 1
        print("[提示] 所有灵敏度均只检测到 1 个峰，返回单峰结果。")
    else:
        if metric == "snr":
            best_k = max(k_candidates, key=lambda k: adjusted_snrs[k])
            if adjusted_snrs[best_k] < 0:
                # 极端情况：所有惩罚后分数为负，退回最小峰数
                best_k = min(k_candidates)
                print("[警告] 所有 K 的惩罚后 SNR 均为负，已选最小峰数作为保底。")
        else:
            best_k = max(k_candidates, key=lambda k: chs[k])

        if metric == "snr" and snrs.get(best_k, 0) < 1.0:
            print(
                f"[警告] 最优 K={best_k} 的原始 SNR={snrs[best_k]:.3f} < 1.0，"
                f"峰间距小于峰内宽度，数据可能近似单峰。"
            )

    best_result = results_by_k[best_k][0]

    if verbose:
        print(
            f"\n  峰数惩罚：α={peak_penalty_rate:.3f}，"
            f"SNR_max={snr_max:.3f}，每多一峰扣 {penalty_per_peak:.4f}"
        )
        print(
            f"  最优分峰方案：K = {best_k}  "
            f"（原始 SNR = {snrs.get(best_k, 0):.3f}，"
            f"惩罚后 SNR_adj = {adjusted_snrs.get(best_k, 0):.3f}）"
        )
        for k in range(best_result.n_peaks):
            cnt   = best_result.peak_counts[k]
            ratio = cnt / len(ave_energy) * 100
            print(
                f"    Peak {k+1}: 中心 = {best_result.peak_centers[k]:.4f} eV/atom  "
                f"构型数 = {cnt} ({ratio:.1f}%)  "
                f"区间 [{best_result.boundaries[k]:.4f}, "
                f"{best_result.boundaries[k+1]:.4f}]"
            )

    if plot_diagnosis:
        _plot_peak_diagnosis(
            ave_energy=ave_energy, results_by_k=results_by_k,
            ks=ks, best_k=best_k, metric=metric,
            snrs=snrs, chs=chs,
            adjusted_snrs=adjusted_snrs,
            peak_penalty_rate=peak_penalty_rate,
            save_path=save_diagnosis, show=show_diagnosis,
        )

    return best_result


def _plot_peak_diagnosis(
    ave_energy, results_by_k, ks, best_k, metric, snrs, chs,
    adjusted_snrs=None, peak_penalty_rate=0.0,
    save_path=None, show=False,
):

    n_configs = len(ave_energy)
    fig = plt.figure(figsize=(13, 5))
    gs  = fig.add_gridspec(1, 2, width_ratios=[1, 1.8], wspace=0.35)
    ax_left  = fig.add_subplot(gs[0])
    ax_right = fig.add_subplot(gs[1])
    ax_r2    = ax_right.twinx()

    raw_vals = [snrs.get(k, 0) for k in ks]
    adj_vals = ([adjusted_snrs.get(k, 0) for k in ks]
                if adjusted_snrs is not None else None)

    bar_colors = ["#D79294" if k == best_k else "#9EC6E0" for k in ks]
    x_pos      = np.arange(len(ks))

    bars = ax_left.bar(
        x_pos, raw_vals, width=0.55,
        color=bar_colors, alpha=0.50,
        edgecolor="white", linewidth=0.7, zorder=3,
        label="SNR (raw)",
    )
    y_max_bar = max(max(raw_vals), 0.01)

    if adj_vals is not None:
        adj_colors = ["#A33040" if k == best_k else "#2A6BA0" for k in ks]
        ax_left.bar(
            x_pos, [max(v, 0) for v in adj_vals], width=0.28,
            color=adj_colors, alpha=0.88,
            edgecolor="white", linewidth=0.5, zorder=4,
            label="SNR_adj (penalized)",
        )

        for xi, (k, av) in enumerate(zip(ks, adj_vals)):
            if av < 0:
                ax_left.annotate(
                    f"{av:.2f}", xy=(xi, 0.01), ha="center",
                    va="bottom", fontsize=7.5, color="#A33040",
                )

    for bar, k, val in zip(bars, ks, raw_vals):
        adj = adjusted_snrs.get(k, val) if adjusted_snrs else val
        label_str = (f"SNR={val:.2f}\nadj={adj:.2f}"
                     if adjusted_snrs is not None else f"{val:.2f}")
        ax_left.text(
            bar.get_x() + bar.get_width() / 2,
            max(val, adj if adj > 0 else 0) + y_max_bar * 0.04,
            label_str,
            ha="center", va="bottom", fontsize=7.5,
            color="#CB181D" if k == best_k else "#333333",
            fontweight="bold" if k == best_k else "normal",
        )

    if metric == "snr":
        ax_left.axhline(1.0, color="#F6882B", linestyle="--",
                        linewidth=1.1, alpha=0.9, zorder=5)
        ax_left.text(len(ks) - 0.45, 1.02, "SNR = 1.0",
                     fontsize=7.5, color="#F6882B", va="bottom", ha="right")

    if best_k in ks:
        ax_left.text(
            ks.index(best_k), -y_max_bar * 0.12,
            "★ optimal", ha="center",
            fontsize=8.5, color="#CB181D", fontweight="bold",
        )

    penalty_note = (
        f"α={peak_penalty_rate:.3f}" if peak_penalty_rate > 0 else ""
    )
    y_label = (
        f"SNR  (penalty {penalty_note})" if metric == "snr"
        else "Calinski-Harabasz Score"
    )
    ax_left.set_xticks(x_pos)
    ax_left.set_xticklabels([str(k) for k in ks])
    ax_left.set_xlabel("Number of Peaks  K", fontsize=10, labelpad=6)
    ax_left.set_ylabel(y_label, fontsize=9, labelpad=6)
    ax_left.set_title(
        f"Peak Count Auto-Selection\nBest K = {best_k}  ({metric.upper()})",
        fontsize=11, fontweight="bold", pad=8,
    )
    ax_left.set_ylim(bottom=0, top=y_max_bar * 1.40)
    ax_left.tick_params(labelsize=9)
    ax_left.set_axisbelow(True)
    ax_left.legend(fontsize=8, loc="upper right", framealpha=0.88)

    best_result = results_by_k[best_k][0]
    n_peaks     = best_result.n_peaks
    colors      = _gen_colors(n_peaks)

    x_lo     = np.percentile(ave_energy, 0.1)
    x_hi     = np.percentile(ave_energy, 99.9)
    x_margin = (x_hi - x_lo) * 0.12
    view_lo, view_hi = x_lo - x_margin, x_hi + x_margin

    _, bin_edges = np.histogram(ave_energy, bins="auto")
    bin_edges = bin_edges[
        (bin_edges >= view_lo - 1e-9) & (bin_edges <= view_hi + 1e-9)
    ]
    if len(bin_edges) < 3:
        _, bin_edges = np.histogram(
            ave_energy[(ave_energy >= view_lo) & (ave_energy <= view_hi)], bins=50
        )
    bin_width = np.mean(np.diff(bin_edges))

    ax_right.bar(
        (bin_edges[:-1] + bin_edges[1:]) / 2,
        np.histogram(ave_energy, bins=bin_edges)[0],
        width=bin_width * 0.95, color="#BDBDBD", alpha=0.40, zorder=1,
    )

    y_max_right    = 0
    legend_handles = []
    for k in range(n_peaks):
        energies_k = ave_energy[best_result.peak_indices[k]]
        cnt        = best_result.peak_counts[k]
        ratio      = cnt / n_configs * 100
        color      = colors[k]
        hist_k, _ = np.histogram(energies_k, bins=bin_edges)
        ax_right.bar(
            (bin_edges[:-1] + bin_edges[1:]) / 2, hist_k,
            width=bin_width * 0.95, color=color, alpha=0.65, zorder=2,
        )
        y_max_right = max(y_max_right, hist_k.max())
        legend_handles.append(
            mpatches.Patch(color=color, alpha=0.75,
                           label=f"Peak {k+1}  ({cnt}, {ratio:.1f}%)"),
        )

    y_max_right = max(y_max_right, 1)
    ax_right.set_ylim(0, y_max_right * 1.38)

    for k in range(n_peaks):
        energies_k = ave_energy[best_result.peak_indices[k]]
        color      = colors[k]
        peak_e     = best_result.peak_centers[k]
        hist_k, _ = np.histogram(energies_k, bins=bin_edges)
        y_top     = hist_k.max() if hist_k.max() > 0 else 1
        ax_right.text(
            peak_e, y_top + y_max_right * 0.04,
            f"Peak {k+1}\n{peak_e:.3f}",
            ha="center", va="bottom", fontsize=8,
            fontweight="bold", color=color, zorder=6,
        )
        ax_right.plot([peak_e, peak_e],
                      [y_top * 1.005, y_top + y_max_right * 0.033],
                      color=color, linewidth=1.0, alpha=0.8, zorder=5)

    for k in range(1, n_peaks):
        bnd = best_result.boundaries[k]
        ax_right.axvline(bnd, color="#888888", linestyle="--",
                         linewidth=0.9, alpha=0.7, zorder=3)

    kde_mask = (best_result.kde_x >= view_lo) & (best_result.kde_x <= view_hi)
    ax_r2.plot(best_result.kde_x[kde_mask], best_result.kde_y[kde_mask],
               color="#333333", linewidth=2.0, alpha=0.75, zorder=4)
    ax_r2.set_ylabel("KDE density (a.u.)", fontsize=9, color="#444444", labelpad=6)
    ax_r2.tick_params(axis="y", labelcolor="#555555", labelsize=8)
    ax_r2.set_ylim(bottom=0)
    ax_r2.spines["right"].set_visible(True)

    ax_right.set_xlim(view_lo, view_hi)
    ax_r2.set_xlim(view_lo, view_hi)
    ax_right.set_xlabel("Average Energy (eV/atom)", fontsize=10, labelpad=6)
    ax_right.set_ylabel("Count", fontsize=10, labelpad=6)
    ax_right.set_title(
        f"Optimal Decomposition  (K = {best_k})",
        fontsize=11, fontweight="bold", pad=8,
    )
    ax_right.tick_params(labelsize=9)
    ax_right.set_axisbelow(True)

    all_patch  = mpatches.Patch(color="#BDBDBD", alpha=0.5, label="All configs")
    kde_handle = Line2D([0], [0], color="#333333", lw=2.0, label="KDE density")
    ax_right.legend(
        handles=[all_patch] + legend_handles + [kde_handle],
        fontsize=8, loc="upper right",
        framealpha=0.92, edgecolor="#CCCCCC", handlelength=1.4,
    )

    fig.suptitle(
        f"Auto Peak Selection — Diagnosis  "
        f"(penalty α={peak_penalty_rate:.3f})",
        fontsize=13, fontweight="bold", y=1.01,
    )
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"[✔]  诊断图已保存：{save_path}")
    if show:
        plt.show()
    plt.close(fig)


_BASE_PALETTE = [
    "#9FBDD7", "#D9989A", "#8EBD9E", "#8C7EAF", "#BA907B",
    "#4DA6E1", "#D66496", "#6FA299", "#946693", "#4B759F",
]


def _gen_colors(n: int) -> list[str]:
    """返回 n 个视觉区分度高的颜色；超过 10 个时从 HSV 色环均匀采样。"""
    if n <= len(_BASE_PALETTE):
        return _BASE_PALETTE[:n]
    import colorsys
    colors = []
    for i in range(n):
        r, g, b = colorsys.hsv_to_rgb(i / n, 0.72, 0.82)
        colors.append(f"#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}")
    return colors


def plot_energy_peaks(
    ave_energy: np.ndarray,
    result: EnergyPeakResult,
    title: str = "Energy Distribution — Peak Decomposition",
    save_path: Optional[str] = "energy_peaks.png",
    figsize: tuple = (9, 4.5),
    show: bool = False,
    hist_bins: int | str = "auto",
    xlim_percentile: tuple = (0.1, 99.9),
    xlabel: str = "Average Energy (eV/atom)",
    ylabel: str = "Count",
) -> plt.Figure:
    """分峰可视化"""
    ave_energy = np.asarray(ave_energy, dtype=float)
    n_configs  = len(ave_energy)
    colors     = _gen_colors(result.n_peaks)

    fig, ax1 = plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()
    ax2.set_zorder(ax1.get_zorder() - 1)
    ax1.patch.set_visible(False)

    x_lo     = np.percentile(ave_energy, xlim_percentile[0])
    x_hi     = np.percentile(ave_energy, xlim_percentile[1])
    x_margin = (x_hi - x_lo) * 0.12
    view_lo, view_hi = x_lo - x_margin, x_hi + x_margin

    _, all_edges = np.histogram(ave_energy, bins=hist_bins)
    bin_edges = all_edges[
        (all_edges >= view_lo - 1e-9) & (all_edges <= view_hi + 1e-9)
    ]
    if len(bin_edges) < 3:
        _, bin_edges = np.histogram(
            ave_energy[(ave_energy >= view_lo) & (ave_energy <= view_hi)], bins=50)
    bin_width = np.mean(np.diff(bin_edges))

    ax1.bar((bin_edges[:-1] + bin_edges[1:]) / 2,
            np.histogram(ave_energy, bins=bin_edges)[0],
            width=bin_width * 0.95, color="#BDBDBD", alpha=0.45,
            label="All configs", zorder=1)

    legend_handles = []
    y_max_global   = 0
    for k in range(result.n_peaks):
        energies_k = ave_energy[result.peak_indices[k]]
        cnt        = result.peak_counts[k]
        ratio      = cnt / n_configs * 100
        color      = colors[k]
        label      = f"Peak {k + 1}  —  {cnt} configs ({ratio:.1f}%)"

        hist_k, _ = np.histogram(energies_k, bins=bin_edges)
        ax1.bar((bin_edges[:-1] + bin_edges[1:]) / 2, hist_k,
                width=bin_width * 0.95, color=color, alpha=0.65, zorder=2)
        y_max_global = max(y_max_global, hist_k.max())
        legend_handles.append(mpatches.Patch(color=color, alpha=0.75, label=label))

    y_max_global = max(y_max_global, 1)
    ax1.set_ylim(0, y_max_global * 1.38)

    for k in range(result.n_peaks):
        energies_k = ave_energy[result.peak_indices[k]]
        color      = colors[k]
        peak_e     = result.peak_centers[k]
        hist_k, _  = np.histogram(energies_k, bins=bin_edges)
        y_top      = hist_k.max() if hist_k.max() > 0 else 1

        ax1.text(peak_e, y_top + y_max_global * 0.04,
                 f"Peak {k + 1}\n{peak_e:.3f}",
                 ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                 color=color, zorder=6)
        ax1.plot([peak_e, peak_e],
                 [y_top * 1.005, y_top + y_max_global * 0.035],
                 color=color, linewidth=1.0, alpha=0.8, zorder=5)

    for k in range(1, result.n_peaks):
        bnd = result.boundaries[k]
        ax1.axvline(bnd, color="#888888", linestyle="--",
                    linewidth=0.9, alpha=0.7, zorder=3)
        # ax1.text(bnd, y_max_global * 0.01, f" {bnd:.3f}",
        #          fontsize=7, color="#666666",
        #          va="bottom", ha="left", rotation=90, zorder=4)

    kde_mask = (result.kde_x >= view_lo) & (result.kde_x <= view_hi)
    ax2.plot(result.kde_x[kde_mask], result.kde_y[kde_mask],
             color="#333333", linewidth=2.0, alpha=0.80, zorder=4)
    ax2.set_ylabel("KDE density (a.u.)", fontsize=9.5, color="#444444", labelpad=8)
    ax2.tick_params(axis="y", labelcolor="#555555", labelsize=8)
    ax2.set_ylim(bottom=0)
    ax2.spines["right"].set_visible(True)

    ax1.set_xlim(view_lo, view_hi)
    ax2.set_xlim(view_lo, view_hi)
    ax1.set_xlabel(xlabel, fontsize=11, labelpad=6)
    ax1.set_ylabel(ylabel, fontsize=11, labelpad=6)
    ax1.set_title(title, fontsize=12, fontweight="bold", pad=10)
    ax1.tick_params(labelsize=9)

    all_patch  = mpatches.Patch(color="#BDBDBD", alpha=0.6, label="All configs")
    kde_handle = Line2D([0], [0], color="#333333", lw=2.0, label="KDE density")
    ax1.legend(handles=[all_patch] + legend_handles + [kde_handle],
               fontsize=8.5, loc="upper right",
               framealpha=0.92, edgecolor="#CCCCCC",
               handlelength=1.4, handleheight=0.9)
    ax1.grid(axis="y", alpha=0.35, linewidth=0.6, linestyle="--")
    ax1.set_axisbelow(True)

    plt.tight_layout(pad=1.2)
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
        print(f"[✔]  分峰图已保存：{save_path}")
    if show:
        plt.show()
    plt.close(fig)
    return fig


def read_xyz_energy(file_name: str) -> list[float]:
    """从 .xyz 文件中提取所有能量值（匹配 'nergy = ...' 格式）。"""
    #pattern = re.compile(r"nergy\s*=\s*(-?\d+\.\d+)")  # 第一版本
    pattern = re.compile(r"nergy\s*=\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)")

    energies = []
    try:
        with open(file_name, "r") as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    energies.append(float(m.group(1)))
    except (FileNotFoundError, IOError) as e:
        _print(f"[red]错误：{e}[/red]")
    return energies


def get_atom_avg_energy(file_name: str) -> np.ndarray:
    """返回每个构型的平均能量（eV/atom）。"""
    atoms     = ase.io.read(file_name, index=":")
    energy    = np.array(read_xyz_energy(file_name))
    atoms_num = np.array([len(a) for a in atoms])
    valid     = atoms_num != 0
    return energy[valid] / atoms_num[valid]



def calculate_probabilities(energies: np.ndarray) -> np.ndarray:
    """Boltzmann 概率（在给定能量集合内归一化）。"""
    energies = np.asarray(energies)
    w = np.exp(-energies)
    return w / w.sum()


def sample_by_energy_intervals(
    ave_energy: np.ndarray,
    probabilities: np.ndarray,
    des: np.ndarray,
    num_intervals: int,
    num_samples: int,
) -> list[int]:
    """
    在 num_intervals 个等间距能量区间内，按 Boltzmann 概率分配配额，
    各区间内用 PCA + 最远点采样挑选多样性最高的构型。

    Returns
    -------
    list[int]
        在传入数组中的局部索引列表。
    """
    energy_bounds  = np.linspace(ave_energy.min(), ave_energy.max(), num_intervals + 1)
    selected_indices = []

    for i in range(num_intervals):
        interval_indices = np.where(
            (ave_energy >= energy_bounds[i]) & (ave_energy < energy_bounds[i + 1])
        )[0]
        if len(interval_indices) == 0:
            continue

        interval_prob_sum    = probabilities[interval_indices].sum()
        interval_sample_count = int(round(interval_prob_sum * num_samples))

        if interval_sample_count == 0 or interval_sample_count >= len(interval_indices):
            selected_indices.extend(interval_indices)
            continue

        interval_des = des[interval_indices]
        pca_des      = PCA(n_components=min(len(interval_des), 2)).fit_transform(
            interval_des
        )
        distances = squareform(pdist(pca_des))

        picked = []
        for _ in range(interval_sample_count):
            if len(picked) == len(interval_indices):
                break
            idx_max = np.unravel_index(np.argmax(distances), distances.shape)
            picked.append(interval_indices[idx_max[0]])
            distances[idx_max[0], :] = 0
            distances[:, idx_max[0]] = 0

        selected_indices.extend(picked)

    return selected_indices


def get_descriptors(nepfile: str, logger: StructuredLogger = None) -> np.ndarray:
    """解析 NEP 势函数，生成 nep.in，调用 nep 程序计算描述符并返回。"""
    info = NEPPotentialInfo.from_nep_txt(nepfile)

    _rule(" NEP 势函数解析 ")
    info.print_rich()
    if logger:
        logger.section("NEP 势函数信息")
        logger.write_raw(info.summary())

    info.write_nep_in("nep.in", output_descriptor=1)
    if logger:
        logger.info("已写入 nep.in（output_descriptor=1）")

    for fname in ("nep.in", "train.xyz"):
        if not os.path.exists(fname):
            _print(f"[bold red]✘  {fname} 不存在，脚本终止。[/bold red]")
            if logger:
                logger.info(f"错误：{fname} 不存在，退出。")
            exit(1)

    if os.path.exists("descriptor.out"):
        os.remove("descriptor.out")
        if logger:
            logger.info("已删除旧的 descriptor.out")

    _rule(" 运行 NEP 计算 ")
    nep_log_path = "nep_run.log"
    with open(nep_log_path, "w", encoding="utf-8") as log_f:
        if _RICH:
            with Progress(
                SpinnerColumn(spinner_name="dots"),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(bar_width=30),
                TimeElapsedColumn(),
                console=console, transient=True,
            ) as prog:
                task = prog.add_task("[cyan]调用 nep 计算描述符...", total=None)
                subprocess.run(["nep"], stdout=log_f,
                               stderr=subprocess.STDOUT, check=True)
                prog.update(task, description="[green]✔  nep 计算完成")
            _print("[green]✔  descriptor.out 已生成[/green]")
        else:
            print("运行 nep 中……")
            subprocess.run(["nep"], stdout=log_f,
                           stderr=subprocess.STDOUT, check=True)
            print("nep 计算完成。")

    if logger:
        logger.info(f"nep 进程完成，详见 {nep_log_path}")

    des = np.loadtxt("descriptor.out")
    if logger:
        logger.info(f"descriptor.out 已读取，形状：{des.shape}")
    return des


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="多峰 bolt + PCA 联合采样脚本")
    parser.add_argument("file_name",   help="输入 .xyz 文件名（任意路径）")
    parser.add_argument("num_samples", type=int, help="目标采样数量")
    parser.add_argument(
        "--workdir", default=None,
        help="工作目录名称（默认：bolt_pca_<xyz文件名去后缀>）",
    )
    args = parser.parse_args()

    src_xyz  = os.path.abspath(args.file_name)
    src_nep  = os.path.abspath("nep.txt")
    stem     = os.path.splitext(os.path.basename(src_xyz))[0]
    workdir  = os.path.abspath(args.workdir if args.workdir else f"bolt_pca_{stem}")

    # ── 准备工作目录 ────────────────────────────────────────────────────
    os.makedirs(workdir, exist_ok=True)

    # 清理上次遗留的输出文件（避免 ase.io.write append 模式叠加旧数据）
    stale_files = [
        "train.xyz", "nep.txt", "nep.in", "descriptor.out",
        "elect_bolt_pca.xyz", "unelect_bolt_pca.xyz",
        "bolt_pca_run.log", "nep_run.log",
        "energy_peaks.png", "peak_auto_diagnosis.png", "elect_bolt_and_pca.png",
    ]
    for fname in stale_files:
        fpath = os.path.join(workdir, fname)
        if os.path.exists(fpath):
            os.remove(fpath)

    # 将输入文件复制进工作目录
    import shutil
    shutil.copy2(src_xyz, os.path.join(workdir, "train.xyz"))
    shutil.copy2(src_nep, os.path.join(workdir, "nep.txt"))

    # 切换到工作目录（后续所有操作均在此目录下进行）
    os.chdir(workdir)

    file_name   = "train.xyz"
    num_samples = args.num_samples
    num_intervals = 30
    nepfile       = "nep.txt"

    ELECT_FILE   = "elect_bolt_pca.xyz"
    UNELECT_FILE = "unelect_bolt_pca.xyz"

    # ── 初始化日志 ──────────────────────────────────────────────────────
    logger = StructuredLogger("bolt_pca_run.log")
    logger.section("启动参数")
    logger.data("输入文件",      file_name)
    logger.data("采样数量",      num_samples)
    logger.data("能量区间数",    num_intervals)
    logger.data("NEP 势函数文件", nepfile)

    # ── 启动横幅 ────────────────────────────────────────────────────────
    if _RICH:
        console.print()
        console.print(
            Panel.fit(
                "[bold cyan]多峰 bolt + PCA 联合采样[/bold cyan]  [dim]v7[/dim]\n"
                f"  输入文件：[yellow]{file_name}[/yellow]   "
                f"目标采样量：[yellow]{num_samples}[/yellow]",
                border_style="bright_blue",
            )
        )
        console.print()
    else:
        print(f"\n{'='*60}")
        print(f"  多峰 bolt + PCA 联合采样 v7")
        print(f"  输入：{file_name}   采样量：{num_samples}")
        print(f"{'='*60}\n")

    # ── 读取构型数据 ─────────────────────────────────────────────────────
    _rule(" 读取构型数据 ")
    atoms      = ase.io.read(file_name, index=":")
    energy     = read_xyz_energy(file_name)
    ave_energy = get_atom_avg_energy(file_name)

    _print(f"  构型总数：[bold white]{len(atoms)}[/bold white]   "
           f"能量条目：[bold white]{len(energy)}[/bold white]")
    logger.section("数据读取")
    logger.data("构型总数",   len(atoms))
    logger.data("能量条目数", len(energy))
    logger.data("平均能量范围 (eV/atom)",
                f"{ave_energy.min():.6f}  ~  {ave_energy.max():.6f}")

    # ── 采样数量合法性检查 ──────────────────────────────────────────────
    if not (1 <= num_samples <= len(ave_energy)):
        msg = (f"采样数量应在 1 ~ {len(ave_energy)} 范围内，"
               f"当前输入：{num_samples}")
        _print(f"[bold red]✘  {msg}[/bold red]")
        logger.info(f"错误：{msg}")
        logger.close(success=False)
        exit(1)

    # ── 自动分峰 ────────────────────────────────────────────────────────
    _rule(" 自动分峰（SNR 最优化）")
    peak_result = auto_find_n_peaks(
        ave_energy,
        metric                 = "snr",
        merge_small_peak_ratio = 0.05,
        plot_diagnosis         = True,
        save_diagnosis         = "peak_auto_diagnosis.png",
        verbose                = True,
    )
    plot_energy_peaks(ave_energy, peak_result, save_path="energy_peaks.png")

    _print(f"  检测到 [bold white]{peak_result.n_peaks}[/bold white] 个峰")
    logger.section("分峰结果")
    for k in range(peak_result.n_peaks):
        cnt   = peak_result.peak_counts[k]
        ratio = cnt / len(ave_energy) * 100
        logger.data(f"Peak {k+1}",
                    f"{cnt} 个构型 ({ratio:.1f}%)  "
                    f"区间 [{peak_result.boundaries[k]:.4f}, "
                    f"{peak_result.boundaries[k+1]:.4f}]")

    # ── 计算描述符 ──────────────────────────────────────────────────────
    des = get_descriptors(nepfile, logger=logger)

    # ── 分峰 Bolt + PCA 采样 ────────────────────────────────────────────
    _rule(" 分峰 Bolt + PCA 采样 ")
    sampled_indices: list[int] = []

    for k, global_idxs in enumerate(peak_result.peak_indices):
        # 按各峰构型数比例分配配额
        quota = int(round(num_samples * len(global_idxs) / len(ave_energy)))
        quota = max(1, min(quota, len(global_idxs)))

        _print(f"  Peak {k+1}：{len(global_idxs)} 个构型  →  配额 [bold]{quota}[/bold]")
        logger.data(f"Peak {k+1} 配额", f"{quota} / {len(global_idxs)}")

        sub_energy = ave_energy[global_idxs]
        sub_des    = des[global_idxs]
        sub_probs  = calculate_probabilities(sub_energy)

        # 区间数按峰内构型数自动缩放，避免小峰区间过多
        sub_intervals = max(5, min(num_intervals, len(global_idxs) // 20))
        sub_selected  = sample_by_energy_intervals(
            sub_energy, sub_probs, sub_des, sub_intervals, quota
        )
        sampled_indices.extend(global_idxs[sub_selected].tolist())

    unelected_indices = list(set(range(len(ave_energy))) - set(sampled_indices))

    logger.section("采样结果")
    logger.data("已选中数量", len(sampled_indices))
    logger.data("未选中数量", len(unelected_indices))
    logger.data("已选中索引", sampled_indices)
    logger.data("未选中索引", unelected_indices)

    # ── 终端展示采样索引 ────────────────────────────────────────────────
    if _RICH:
        DISPLAY_LIMIT = 30

        def _chunk_lines(idxs, per_row=6, width=4):
            if not idxs:
                return ["(none)"]
            rows = [idxs[i: i + per_row] for i in range(0, len(idxs), per_row)]
            return [" ".join(f"{x:>{width}d}" for x in row) for row in rows]

        sel_sorted   = sorted(sampled_indices)
        unsel_sorted = sorted(unelected_indices)
        sel_lines    = _chunk_lines(sel_sorted[:DISPLAY_LIMIT])
        unsel_lines  = _chunk_lines(unsel_sorted[:DISPLAY_LIMIT])

        if len(sel_sorted) > DISPLAY_LIMIT:
            sel_lines.append(f"  … +{len(sel_sorted) - DISPLAY_LIMIT} more")
        if len(unsel_sorted) > DISPLAY_LIMIT:
            unsel_lines.append(f"  … +{len(unsel_sorted) - DISPLAY_LIMIT} more")

        tbl = Table(box=box.MINIMAL_DOUBLE_HEAD, show_header=True,
                    header_style="bold", expand=True)
        tbl.add_column(f"Selected [{len(sampled_indices)}]",   style="#FFFFFF")
        tbl.add_column(f"Unselected [{len(unelected_indices)}]", style="#FFFFFF")
        for i in range(max(len(sel_lines), len(unsel_lines))):
            tbl.add_row(
                sel_lines[i]   if i < len(sel_lines)   else "",
                unsel_lines[i] if i < len(unsel_lines) else "",
            )
        console.print(tbl)

    # ── 写出构型文件 ─────────────────────────────────────────────────────
    for fname in (ELECT_FILE, UNELECT_FILE):
        if os.path.exists(fname):
            os.remove(fname)

    for i in sampled_indices:
        ase.io.write(ELECT_FILE,   atoms[i], append=True)
    for j in unelected_indices:
        ase.io.write(UNELECT_FILE, atoms[j], append=True)

    logger.info(f"{ELECT_FILE} 写出完成")
    logger.info(f"{UNELECT_FILE} 写出完成")

    # ── 采样前后能量分布对比图 ──────────────────────────────────────────
    _, bin_edges     = np.histogram(ave_energy, bins="auto")
    ave_energy_elect = get_atom_avg_energy(ELECT_FILE)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(ave_energy,       bins=bin_edges, alpha=0.45,
            color="#a56cc1", label="All")
    ax.hist(ave_energy_elect, bins=bin_edges, alpha=0.55,
            color="#d62728", label="Selected")
    ax.set_title("Histogram of Average Energy Values")
    ax.set_xlabel("Average Energy (eV/atom)")
    ax.set_ylabel("Count")
    ax.legend()
    ax.grid(axis="y", alpha=0.5)
    plt.tight_layout()
    plt.savefig("elect_bolt_and_pca.png", dpi=150)
    plt.close(fig)

    # ── 运行摘要 ────────────────────────────────────────────────────────
    _rule(" 采样结果 ")
    if _RICH:
        summary_tbl = Table(
            box=box.SIMPLE_HEAVY, show_header=True,
            header_style="bold bright_white",
            border_style="bright_blue", pad_edge=True, expand=False,
        )
        summary_tbl.add_column("项目", style="white", no_wrap=True, justify="right")
        summary_tbl.add_column("内容", style="white", overflow="fold")

        summary_tbl.add_row("输入文件", file_name)
        summary_tbl.add_row("构型总数",  Text(str(len(energy)), style="bold white"))
        summary_tbl.add_row("已选中",   Text(str(len(sampled_indices)),   style="#73F173"))
        summary_tbl.add_row("未选中",   Text(str(len(unelected_indices)), style="#FDC0C0"))
        summary_tbl.add_section()
        summary_tbl.add_row("输出文件（选中）",  ELECT_FILE)
        summary_tbl.add_row("输出文件（未选中）", UNELECT_FILE)
        summary_tbl.add_row("能量分布图",       "elect_bolt_and_pca.png")
        summary_tbl.add_row("峰检测诊断图",     "peak_auto_diagnosis.png")
        summary_tbl.add_row("分峰结果图",       "energy_peaks.png")
        summary_tbl.add_row("运行日志",         "bolt_pca_run.log")
        console.print(summary_tbl)
    else:
        print("\n── 运行完成 ──")
        print(f"  已选中：{len(sampled_indices)} 个 → {ELECT_FILE}")
        print(f"  未选中：{len(unelected_indices)} 个 → {UNELECT_FILE}")
        print(f"  能量图：elect_bolt_and_pca.png")
        print(f"  日志：bolt_pca_run.log")

    logger.close(success=True)


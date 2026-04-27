"""
bolt_and_pca_elect-v4.py
变更说明：在v4 的基础上，对于内部bug进行了修正，主要修正图的展示，以及采样数据的展示
单峰采样完整版本
"""

from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
import ase.io
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import os
from dataclasses import dataclass
from pathlib import Path
import re
from typing import List, Optional
import subprocess
import datetime
import textwrap


try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.text import Text
    from rich import box
    from rich.rule import Rule
    from rich.columns import Columns
    from rich.progress import (
        Progress, SpinnerColumn, TextColumn,
        BarColumn, TimeElapsedColumn
    )
    _RICH = True
    console = Console()
except ImportError:          # 没有 rich 时退回普通 print
    _RICH = False
    console = None

# ── 统一打印辅助 ──────────────────────────────────────────────────────

def _print(msg="", **kw):
    if _RICH:
        console.print(msg, **kw)
    else:
        print(msg)

def _rule(title=""):
    if _RICH:
        console.print()
        console.print(f"[bold #39FB90][[/bold #39FB90][bold #ABE2F9]{title}[/bold #ABE2F9][bold #39FB90]][/bold #39FB90]")
        #console.print()

    else:
        print(f"\n{'——>'}  {title}")

# ── 结构化日志辅助 ────────────────────────────────────────────────────

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
        self._f.write(f"  bolt_and_pca_elect 运行日志\n")
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

# ── NEP Potential 信息类（与 v4 完全相同） ────────────────────────────

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
        path = Path(filename)
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

        clean_lines = []
        for line in lines:
            line = line.split("#", 1)[0].strip()
            if line:
                clean_lines.append(line)

        if len(clean_lines) < 6:
            raise ValueError("This file is too short to be a normal NEP potential file.")

        first = clean_lines[0].split()
        version_match = re.search(r"\d+", first[0])
        if not version_match:
            raise ValueError(f"Cannot parse NEP version from line: {clean_lines[0]}")

        version = int(version_match.group())
        n_types = int(first[1])
        elements = first[2:2 + n_types]

        if len(elements) != n_types:
            raise ValueError(
                f"Expected {n_types} elements, but got {len(elements)}: {elements}"
            )

        cutoff_tokens = clean_lines[1].split()
        if cutoff_tokens[0].lower() != "cutoff":
            raise ValueError(f"Expected cutoff line, got: {clean_lines[1]}")

        cutoff_radial = float(cutoff_tokens[1])
        cutoff_angular = float(cutoff_tokens[2])

        max_neighbors_radial = None
        max_neighbors_angular = None
        if len(cutoff_tokens) >= 5:
            max_neighbors_radial = int(float(cutoff_tokens[3]))
            max_neighbors_angular = int(float(cutoff_tokens[4]))

        nmax_tokens = clean_lines[2].split()
        if nmax_tokens[0].lower() != "n_max":
            raise ValueError(f"Expected n_max line, got: {clean_lines[2]}")

        n_max_radial = int(nmax_tokens[1])
        n_max_angular = int(nmax_tokens[2])

        basis_tokens = clean_lines[3].split()
        if basis_tokens[0].lower() != "basis_size":
            raise ValueError(f"Expected basis_size line, got: {clean_lines[3]}")

        basis_size_radial = int(basis_tokens[1])
        basis_size_angular = int(basis_tokens[2])

        lmax_tokens = clean_lines[4].split()
        if lmax_tokens[0].lower() != "l_max":
            raise ValueError(f"Expected l_max line, got: {clean_lines[4]}")

        l_values = [int(x) for x in lmax_tokens[1:]]
        while len(l_values) < 3:
            l_values.append(0)

        l_max_3body, l_max_4body, l_max_5body = l_values[:3]

        ann_tokens = clean_lines[5].split()
        if ann_tokens[0].lower() != "ann":
            raise ValueError(f"Expected ANN line, got: {clean_lines[5]}")

        n_neurons = int(ann_tokens[1])
        ann_extra = ann_tokens[2:]

        n_radial = n_max_radial + 1
        n_3body = (n_max_angular + 1) * l_max_3body
        n_4body = (n_max_angular + 1) if l_max_4body > 0 else 0
        n_5body = (n_max_angular + 1) if l_max_5body > 0 else 0
        n_descriptor = n_radial + n_3body + n_4body + n_5body

        n_descriptor_parameters = n_types ** 2 * (
            (n_max_radial + 1) * (basis_size_radial + 1)
            + (n_max_angular + 1) * (basis_size_angular + 1)
        )

        if version == 3:
            n_nn_parameters = (n_descriptor + 2) * n_neurons + 1
        elif version >= 4:
            n_nn_parameters = (n_descriptor + 2) * n_neurons * n_types + 1
        else:
            raise ValueError(f"This parser currently expects NEP3/NEP4, got NEP{version}.")

        payload_values = []
        for line in clean_lines[6:]:
            for token in line.split():
                try:
                    payload_values.append(float(token))
                except ValueError:
                    pass

        n_payload_values = len(payload_values)

        expected_payload_values = (
            n_descriptor_parameters
            + n_nn_parameters
            + n_descriptor
        )

        return cls(
            path=path,
            version=version,
            n_types=n_types,
            elements=elements,
            cutoff_radial=cutoff_radial,
            cutoff_angular=cutoff_angular,
            max_neighbors_radial=max_neighbors_radial,
            max_neighbors_angular=max_neighbors_angular,
            n_max_radial=n_max_radial,
            n_max_angular=n_max_angular,
            basis_size_radial=basis_size_radial,
            basis_size_angular=basis_size_angular,
            l_max_3body=l_max_3body,
            l_max_4body=l_max_4body,
            l_max_5body=l_max_5body,
            n_neurons=n_neurons,
            ann_extra=ann_extra,
            n_descriptor=n_descriptor,
            n_descriptor_parameters=n_descriptor_parameters,
            n_nn_parameters=n_nn_parameters,
            n_payload_values=n_payload_values,
            expected_payload_values=expected_payload_values,
        )

    def to_nep_in(
        self,
        output_descriptor: int = 1,
        prediction: int = 1,
        include_neuron: bool = True,
    ) -> str:
        if output_descriptor not in (0, 1, 2):
            raise ValueError("output_descriptor must be 0, 1, or 2.")

        lines = [
            f"version         {self.version}",
            f"type            {self.n_types} {' '.join(self.elements)}",
            f"cutoff          {self.cutoff_radial:g}  {self.cutoff_angular:g}",
            f"n_max           {self.n_max_radial}  {self.n_max_angular}",
            f"basis_size      {self.basis_size_radial}  {self.basis_size_angular}",
            f"l_max           {self.l_max_3body}  {self.l_max_4body}  {self.l_max_5body}",
        ]

        if include_neuron:
            lines.append(f"neuron          {self.n_neurons}")

        lines += [
            "",
            f"prediction      {prediction}",
            f"output_descriptor {output_descriptor}",
        ]

        return "\n".join(lines) + "\n"

    def write_nep_in(
        self,
        filename: str | Path = "nep.in",
        output_descriptor: int = 1,
    ) -> None:
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

        # ── 标题 Panel ──
        title_text = Text()
        title_text.append(f"NEP{self.version}", style="bold bright_cyan")
        title_text.append("   元素：", style="bold white")
        title_text.append(" ".join(self.elements), style="bold yellow")
        console.print(title_text)

        # ── 参数总览表 ──
        tbl = Table(box=box.ROUNDED, show_header=True,
                    header_style="bold cyan", min_width=52)
        tbl.add_column("参数", style="white", no_wrap=True)
        tbl.add_column("径向 (Radial)", justify="right", style="bright_green")
        tbl.add_column("角向 (Angular)", justify="right", style="bright_magenta")

        tbl.add_row("截断半径 cutoff  (Å)",
                    f"{self.cutoff_radial:g}", f"{self.cutoff_angular:g}")
        tbl.add_row("最大近邻数 max_neighbors",
                    str(self.max_neighbors_radial),
                    str(self.max_neighbors_angular))
        tbl.add_row("n_max",
                    str(self.n_max_radial), str(self.n_max_angular))
        tbl.add_row("basis_size",
                    str(self.basis_size_radial), str(self.basis_size_angular))
        console.print(tbl)

        # ── l_max + ANN ──
        tbl2 = Table(box=box.SIMPLE, show_header=False, min_width=52)
        tbl2.add_column("key", style="white")
        tbl2.add_column("val", style="bold white")
        tbl2.add_row("l_max  3/4/5-body",
                     f"{self.l_max_3body}  /  {self.l_max_4body}  /  {self.l_max_5body}")
        tbl2.add_row("隐藏神经元 neurons", str(self.n_neurons))
        tbl2.add_row("描述符维度 n_descriptor", str(self.n_descriptor))
        tbl2.add_row("描述符参数数", str(self.n_descriptor_parameters))
        tbl2.add_row("NN 参数数", str(self.n_nn_parameters))

        # payload 校验着色
        ok = self.n_payload_values == self.expected_payload_values
        payload_style = "bold bright_green" if ok else "bold bright_red"
        payload_icon  = "✔" if ok else "✘"
        tbl2.add_row(
            "文件参数值 / 预期",
            Text(f"{self.n_payload_values}  /  {self.expected_payload_values}"
                 f"  {payload_icon}", style=payload_style)
        )
        console.print(tbl2)



def sample_by_energy_intervals(ave_energy, probabilities, des, num_intervals, num_samples):
    energy_bounds = np.linspace(min(ave_energy), max(ave_energy), num_intervals+1)
    selected_indices = []

    for i in range(num_intervals):
        interval_indices = np.where(
            (ave_energy >= energy_bounds[i]) & (ave_energy < energy_bounds[i+1])
        )[0]
        if len(interval_indices) == 0:
            continue

        interval_prob_sum = probabilities[interval_indices].sum()
        interval_sample_count = int(round(interval_prob_sum * num_samples))

        if interval_sample_count >= len(interval_indices) or interval_sample_count == 0:
            selected_indices.extend(interval_indices)
            continue

        interval_des = des[interval_indices]
        pca = PCA(n_components=min(len(interval_des), 2))
        pca_des = pca.fit_transform(interval_des)
        distances = squareform(pdist(pca_des))

        interval_selected_indices = []
        for _ in range(interval_sample_count):
            if len(interval_selected_indices) == len(interval_indices):
                break
            idx_max = np.unravel_index(np.argmax(distances), distances.shape)
            interval_selected_indices.append(interval_indices[idx_max[0]])
            distances[idx_max[0], :] = 0
            distances[:, idx_max[0]] = 0

        selected_indices.extend(interval_selected_indices)

    return selected_indices


def read_atoms_num(atoms):
    nums = []
    for i in range(len(atoms)):
        nums.append(len(atoms[i]))
    return nums


def read_energy_line(line):
    pattern = re.compile(r"nergy\s*=\s*(-?\d+\.\d+)")
    match = pattern.search(line)
    if match:
        return float(match.group(1))
    return None


def read_xyz_energy(file_name):
    energies = []
    try:
        with open(file_name, "r") as file:
            for line in file:
                energy = read_energy_line(line)
                if energy is not None:
                    energies.append(energy)
    except FileNotFoundError:
        _print(f"[red]错误：文件 {file_name} 不存在。[/red]")
    except IOError:
        _print(f"[red]错误：无法读取文件 {file_name}。[/red]")
    return energies


# def get_atom_avg_energy(file_name):
#     energy    = np.array(read_xyz_energy(file_name))
#     atoms_num = np.array(read_atoms_num(atoms))
#     print(f"atoms_num: {atoms_num}")
#     valid_indices = atoms_num != 0
#     energy    = energy[valid_indices]
#     atoms_num = atoms_num[valid_indices]
#     return energy / atoms_num

def get_atom_avg_energy(file_name):
    _atoms    = ase.io.read(file_name, index=":")  
    energy    = np.array(read_xyz_energy(file_name))
    atoms_num = np.array(read_atoms_num(_atoms))
    valid     = atoms_num != 0
    return energy[valid] / atoms_num[valid]



def calculate_probabilities(energies):
    energies = np.array(energies)
    probabilities = np.exp(-energies) / np.sum(np.exp(-energies))
    return probabilities


def sample_configurations(energies, probabilities, num_samples=600):
    """根据概率采样构型（本脚本中暂不使用）。"""
    indices = np.arange(len(energies))
    sampled_indices = np.random.choice(indices, size=num_samples,
                                       p=probabilities, replace=False)
    return sampled_indices


def get_descriptors(nepfile, trainfilename="train.xyz", logger: StructuredLogger = None):
    info = NEPPotentialInfo.from_nep_txt("nep.txt")

    # ── 屏幕美化输出 ──
    _rule(" NEP 势函数解析 ")
    info.print_rich()

    # ── 写入日志 ──
    if logger:
        logger.section("NEP 势函数信息")
        logger.write_raw(info.summary())

    info.write_nep_in("nep.in", output_descriptor=1)
    if logger:
        logger.info("已写入 nep.in（output_descriptor=1）")

    # ── 文件检查 ──
    for fname in ("nep.in", "train.xyz"):
        if not os.path.exists(fname):
            _print(f"[bold red]✘  {fname} 不存在，脚本终止。[/bold red]")
            if logger:
                logger.info(f"错误：{fname} 不存在，退出。")
            exit()

    if os.path.exists("descriptor.out"):
        os.remove("descriptor.out")
        if logger:
            logger.info("已删除旧的 descriptor.out")

    # ── 运行 nep ──
    _rule(" 运行 NEP 计算 ")
    nep_log_path = "nep_run.log"

    if _RICH:
        with Progress(
            SpinnerColumn(spinner_name="dots"),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=30),
            TimeElapsedColumn(),
            console=console,
            transient=True,
        ) as prog:
            task = prog.add_task("[cyan]调用 nep 计算描述符...", total=None)
            with open(nep_log_path, "w", encoding="utf-8") as log:
                subprocess.run(
                    ["nep"],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True
                )
            prog.update(task, description="[green]✔  nep 计算完成")
        _print("[green]✔  descriptor.out 已生成[/green]")
    else:
        print("运行 nep 中……")
        with open(nep_log_path, "w", encoding="utf-8") as log:
            subprocess.run(
                ["nep"],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True
            )
        print("nep 计算完成。")

    if logger:
        logger.info(f"nep 进程完成，详见 {nep_log_path}")

    des = np.loadtxt("descriptor.out")
    if logger:
        logger.info(f"descriptor.out 已读取，形状：{des.shape}")
    return des



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="bolt + PCA 联合采样脚本"
    )
    parser.add_argument("file_name",   help="输入 .xyz 文件名")
    parser.add_argument("num_samples", type=int, help="采样数量")
    args = parser.parse_args()
    file_name   = args.file_name
    num_samples = args.num_samples

    num_intervals = 30
    nepfile       = "nep.txt"

    elect_file_name = "elect_bolt_pca.xyz"
    unelect_file_name = "unelect_bolt_pca.xyz"

    # ── 初始化日志 ──
    logger = StructuredLogger("bolt_pca_run.log")
    logger.section("启动参数")
    logger.data("输入文件", file_name)
    logger.data("采样数量", num_samples)
    logger.data("能量区间数", num_intervals)
    logger.data("NEP 势函数文件", nepfile)

    # ── 启动横幅 ──
    if _RICH:
        console.print()
        console.print(
            Panel.fit(
                "[bold cyan]bolt + PCA 联合采样[/bold cyan]  [dim]v4[/dim]\n"
                f"  输入文件：[yellow]{file_name}[/yellow]   "
                f"目标采样量：[yellow]{num_samples}[/yellow]",
                border_style="bright_blue",
            )
        )
        console.print()
    else:
        print(f"\n{'='*60}")
        print(f"  bolt + PCA 联合采样 v4")
        print(f"  输入：{file_name}   采样量：{num_samples}")
        print(f"{'='*60}\n")


    _rule(" 读取构型数据 ")
    atoms       = ase.io.read(file_name, index=":")
    energy      = read_xyz_energy(file_name)
    atoms_num   = read_atoms_num(atoms)
    ave_energy  = get_atom_avg_energy(file_name)

    _print(f"  构型总数：[bold white]{len(atoms)}[/bold white]   "
           f"能量条目：[bold white]{len(energy)}[/bold white]")
    logger.section("数据读取")
    logger.data("构型总数", len(atoms))
    logger.data("能量条目数", len(energy))
    logger.data("平均能量范围 (eV/atom)",
                f"{ave_energy.min():.6f}  ~  {ave_energy.max():.6f}")

    # ── 获取描述符 ──
    des = get_descriptors(nepfile, file_name, logger=logger)

    # ── 采样数量校验 ──
    if num_samples > len(ave_energy) or num_samples <= 0:
        msg = (f"采样数量不符合要求，应在 1 ~ {len(ave_energy)} 范围内。"
               f"当前输入：{num_samples}")
        _print(f"[bold red]✘  {msg}[/bold red]")
        logger.info(f"错误：{msg}")
        logger.close(success=False)
        exit()

    # ── 采样 ──
    _rule(" PCA + Boltzmann 采样 ")
    probabilities    = calculate_probabilities(ave_energy)
    sampled_indices  = sample_by_energy_intervals(
        ave_energy, probabilities, des, num_intervals, num_samples
    )
    unelected_indices = list(set(range(len(ave_energy))) - set(sampled_indices))

    logger.section("采样结果")
    logger.data("已选中数量", len(sampled_indices))
    logger.data("未选中数量", len(unelected_indices))
    logger.data("已选中索引", sampled_indices)
    logger.data("未选中索引", unelected_indices)
    from rich.align import Align

    # ── 结果展示 ──
    if _RICH:
        console = Console()

        def _chunk_lines(idxs, per_row=10, width=4):
            if not idxs:
                return ["(none)"]
            rows = [idxs[i:i + per_row] for i in range(0, len(idxs), per_row)]
            return [
                " ".join(f"{x:>{width}d}" for x in row)
                for row in rows
            ]

        sel_lines = _chunk_lines(sorted(sampled_indices), per_row=6)
        unsel_lines = _chunk_lines(sorted(unelected_indices), per_row=6)
        nrows = max(len(sel_lines), len(unsel_lines))

        tbl = Table(
            box=box.MINIMAL_DOUBLE_HEAD,
            show_header=True,
            header_style="bold",
            expand=True,
        )

        tbl.add_column(f"Selected [{len(sampled_indices)}]", style="#74FCF8")
        tbl.add_column(f"Unselected [{len(unelected_indices)}]", style="#FF9E9E")

        for i in range(nrows):
            tbl.add_row(
                sel_lines[i] if i < len(sel_lines) else "",
                unsel_lines[i] if i < len(unsel_lines) else "",
            )

        console.print(tbl)
    else:
        print(f"选中的构型索引：{sampled_indices}")
        print(f"未选中的构型索引：{unelected_indices}")
        print(f"选中数量：{len(sampled_indices)}")
        print(f"未选中数量：{len(unelected_indices)}")

    # ── 写出文件 ──
    #_rule(" 写出构型文件 ")
    for fname in (elect_file_name, unelect_file_name):
        if os.path.exists(fname):
            os.remove(fname)

    for i in sampled_indices:
        ase.io.write(elect_file_name, atoms[i], append=True)
    for j in unelected_indices:
        ase.io.write(unelect_file_name, atoms[j], append=True)

    # _print("  [green]✔[/green]  elect_bolt_pca.xyz      已写出")
    # _print("  [red]✘[/red]  bolt_pca_unelect.xyz   已写出")
    logger.info(f"{elect_file_name} 写出完成")
    logger.info(f"{unelect_file_name} 写出完成")

    # ── 能量分布图 ──
    #_rule(" 生成能量分布图 ")
    favorite_colors = [
        '#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd',
        '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff'
    ]

    print(f"[DEBUG] All ave_energy: n={len(ave_energy)}, min={ave_energy.min():.4f}, max={ave_energy.max():.4f}")

    plt.switch_backend('Agg')

    _, bin_edges = np.histogram(ave_energy, bins='auto')
    plt.hist(ave_energy, bins=bin_edges, alpha=0.4,
             color=favorite_colors[1], label='All')

    atoms      = ase.io.read('elect_bolt_pca.xyz', index=":")
    ave_energy = get_atom_avg_energy('elect_bolt_pca.xyz')
    plt.hist(ave_energy, bins=bin_edges, alpha=0.4,
             color=favorite_colors[3], label='Elect')
    plt.title('Histogram of Average Energy Values')
    plt.legend()
    #plt.xlim(-7.5, -5)
    plt.xlabel('Average Energy')
    plt.ylabel('Count')
    plt.grid(axis='y', alpha=0.75)
    plt.savefig('elect_bolt_and_pca.png')

    # _print("  [green]✔[/green]  elect_bolt_and_pca.png  已保存")
    # logger.info("elect_bolt_and_pca.png 已保存")

    _rule(" 采样结果 ")

    if _RICH:
        console = Console()

        summary_tbl = Table(
            #title="运行摘要",
            title_style="bold cyan",
            box=box.SIMPLE_HEAVY,
            show_header=True,
            header_style="bold bright_white",
            border_style="bright_blue",
            pad_edge=True,
            expand=False,
        )

        summary_tbl.add_column("项目", style="bold cyan", no_wrap=True, justify="right")
        summary_tbl.add_column("内容", style="white", overflow="fold")

        summary_tbl.add_row("输入文件", file_name)
        summary_tbl.add_row("构型总数", Text(str(len(energy)), style="bold white"))
        summary_tbl.add_row("已选中", Text(str(len(sampled_indices)), style="bold bright_green"))
        summary_tbl.add_row("未选中", Text(str(len(unelected_indices)), style="bold bright_red"))
        summary_tbl.add_section()
        summary_tbl.add_row("输出文件（选中）", "elect_bolt_pca.xyz")
        summary_tbl.add_row("输出文件（未选中）", "bolt_pca_unelect.xyz")
        summary_tbl.add_row("能量分布图", "elect_bolt_and_pca.png")
        summary_tbl.add_row("运行日志", "bolt_pca_run.log")

        console.print(summary_tbl)
    else:
        print("\n── 运行完成 ──")
        print(f"  已选中：{len(sampled_indices)} 个 → elect_bolt_pca.xyz")
        print(f"  未选中：{len(unelected_indices)} 个 → bolt_pca_unelect.xyz")
        print(f"  能量图：elect_bolt_and_pca.png")
        print(f"  日志：bolt_pca_run.log")

    logger.close(success=True)


<div align="center">

![NebulaFlow Logo](logo/NebulaFlow-logo.png)

**🌌 一行命令，完成复杂科学计算 / One Command, Complex Science Done.**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL-green.svg)](#)
[![Shell](https://img.shields.io/badge/shell-bash%20%7C%20zsh-orange.svg)](#)

</div>

---

## 目录 / Table of Contents

- [What is NebulaFlow?](#what-is-nebulaflow)
- [项目结构 / Project Structure](#项目结构--project-structure)
- [快速开始 / Quick Start](#快速开始--quick-start)
- [安装 / Installation](#安装--installation)
- [环境要求 / Requirements](#环境要求--requirements)
- [目前支持的计算 / Supported Calculations](#目前支持的计算--supported-calculations)
- [使用手册 / User Manual](#使用手册--user-manual)
- [贡献 / Contributing](#贡献--contributing)

---

## What is NebulaFlow?

**NebulaFlow** 是一款辅助**分子动力学（MD）**和**密度泛函理论（DFT）**计算的命令行工具。核心理念是用 **Shell 编排 Python/C++/其他语言**，让用户**只需一条命令**完成一系列复杂任务。

**NebulaFlow** is a command-line tool that assists **molecular dynamics (MD)** and **density functional theory (DFT)** calculations. Its core philosophy is to **use Shell to orchestrate Python, C++, and other languages**, allowing users to complete complex workflows with a **single command**.

### 为什么选择 NebulaFlow？ / Why NebulaFlow?

- 🔰 **零门槛 / Low barrier to entry** — 不需要会写 Python 脚本，一行命令搞定数据处理
- ⚡ **高效 / High efficiency** — 管道式处理，多 GPU 自动调度，自动化工作流
- 🔧 **灵活 / Flexible** — Shell 作为胶水语言，组合现有工具和脚本，可轻松扩展
- 📦 **可复用 / Reusable** — 所有函数模块化，开箱即用

---

## 项目结构 / Project Structure

```
NebulaFlow/
├── rebreath-env-function        # 🧠 核心函数库 / Core function library (entry point)
├── gpuq                         # 📊 GPU 队列管理器 / GPU queue management
├── NebulaFlowinstaller.sh       # 📦 一键安装脚本 / One-click installer
├── .config                      # ⚙️  计算软件配置 / Executable paths config
│
├── envsrc/                      # 🗂️  按功能分类的环境函数 / Environment functions by domain
│   ├── gpumdEnvFunction.sh      #    GPUMD / NEP 训练 / HNEMD 热导率
│   ├── vaspEnvFunction.sh       #    VASP 计算 / 声子谱 / 弹性模量
│   ├── dealDataEnvFunction.sh   #    数据处理 / 晶粒分析 / 配位数
│   ├── ioEnvFunction.sh         #    文件格式转换 (xyz↔POSCAR↔data↔cif)
│   ├── cp2kEnvFunction.sh       #    CP2K 输入文件处理
│   ├── lammpsEnvFunction.sh     #    Lammps 相关
│   └── ...
│
├── nebula_pylib/                # 🐍 Python 核心库 / Python core library
│   └── nebula.py                #    xyz 构型读写 / Config 类 / thermo 解析
│
├── deal_data/                   # 📐 数据处理脚本 / Data processing scripts
│   ├── dataset_quality_diagnosis.py  # 数据集质量诊断工具
│   ├── elect_rely_*.py          #    训练集筛选 (力/能量/位力)
│   ├── bolt_and_pca_elect*.py   #    PCA 结构筛选
│   └── ...
│
├── auto/                        # 🤖 自动化工作流 / Automated workflows
│   ├── auto_vasp_nep1/          #    VASP → NEP 全自动训练流程
│   ├── modify_nep/              #    NEP 势迭代优化
│   └── phono_specturm_vasp/     #    VASP+phonopy 声子谱
│
├── plot_library/                # 📈 可视化与绘图 / Visualization & plotting
├── compute_lib/                 # 🧮 计算库 / Computation library
├── sh_lib/                      # 📜 Shell 辅助脚本 / Shell helper scripts
├── inp_lib/                     # 📋 输入文件模板 / Input templates (VASP/CP2K/Lammps)
├── Manual/                      # 📖 中英文使用手册 / User manuals
└── expample/                    # 🧪 示例文件与测试用例 / Examples & test cases
```

---

## 快速开始 / Quick Start

```bash
# 1. 安装 / Install
git clone https://github.com/re-breath/NebulaFlow.git
cd NebulaFlow
bash NebulaFlowinstaller.sh
source ~/.bashrc
# 看到 "NebulaFlow library loaded O.<" 即安装成功
# If you see "NebulaFlow library loaded O.<", you're all set!

# 2. 配置计算软件路径 / Configure executable paths
vim ~/.rebreath/.config
# 设置 gpumd_exe, vasp_exe 等路径 / Set paths for gpumd_exe, vasp_exe, etc.

# 3. 开始使用 / Start using
# 示例：从 VASP OUTCAR 文件构建 NEP 训练集
deal_outcar_to_train   # 将当前目录所有 OUTCAR 转为 train.xyz

# 示例：检查数据集质量
check_dataset_quality  # 对 train.xyz 进行多维度诊断

# 示例：在空闲 GPU 上运行计算
gpuq submit 'gpumd'
gpuq list              # 查看任务队列
```

---

## 安装 / Installation

### 一键安装 / One-click Install

```bash
git clone https://github.com/re-breath/NebulaFlow.git
cd NebulaFlow
bash NebulaFlowinstaller.sh
source ~/.bashrc
```

或者从 release 包安装 / Or install from release:
```bash
tar -zxvf NebulaFlow-1.0.tar.gz
cd NebulaFlow-1.0
bash NebulaFlowinstaller.sh
source ~/.bashrc
```

### 更新 / Update

```bash
update_NebulaFlow   # 自动从 GitHub 拉取最新代码并重新安装
```

### Windows 用户注意事项 / Notes for Windows Users

在 Windows 转 Linux 环境中（如 Git Bash → WSL），可能出现格式问题导致报错。使用以下命令修复：

On Windows-to-Linux environments (e.g. Git Bash → WSL), format issues may cause errors. Fix with:

```bash
dos2unix ~/.rebreath/*
```

> **建议 / Recommendation**: 在 **WSL** 中使用 NebulaFlow 以获得最佳体验。
> Use **WSL** for the best experience with NebulaFlow.

---

## 环境要求 / Requirements

### 系统 / System

| 项目 | 要求 / Requirement |
|------|------|
| 操作系统 / OS | Linux (推荐/Recommended) 或 WSL |
| Shell | Bash ≥ 4.0 |
| Python | Python 3.8+ (默认使用 `python3`) |

### Python 依赖 / Python Dependencies

NebulaFlow 使用 Shell 调用 Python 库，许多命令**无需额外依赖**即可使用。部分特殊功能需要以下库：

NebulaFlow uses Shell to call Python libraries. Many commands work **without extra dependencies**. Some specialized functions require the following:

| 库 / Library | 用途 / Purpose | 优先级 / Priority |
|-------------|---------------|-----------------|
| `numpy` | 数值计算 / Numerical computation | 必需 / Required |
| `matplotlib` | 画图 / Plotting | 推荐 / Recommended |
| `ase` | 原子结构读写 / Atomic structure I/O | 推荐 / Recommended |
| `scikit-learn` | PCA、聚类 / PCA, clustering | 可选 / Optional |
| `ovito` | 晶粒分析 / Grain analysis | 可选 / Optional |
| `calorine` | GPUMD 弹性模量 / Elastic moduli | 可选 / Optional |
| `pynep` | NEP 势分析 / NEP analysis | 可选 / Optional |
| `gpumd-wizard` | GPUMD 辅助 / GPUMD utility | 可选 / Optional |
| `pymatgen` | 吸附位点 / Adsorption sites | 可选 / Optional |

> 💡 用到时再安装即可 / Install when needed: `pip install <library>`

### 外部软件 / External Tools

| 软件 / Tool | 用途 / Purpose |
|------------|---------------|
| GPUMD | MD 模拟 / Molecular dynamics |
| VASP | DFT 计算 / DFT calculations |
| vaspkit | VASP 辅助 / VASP utility |
| phonopy | 声子谱计算 / Phonon spectrum |

---

## 目前支持的计算 / Supported Calculations

### GPUMD 相关 / GPUMD Related

| 功能 / Feature | 命令 / Command | 状态 / Status |
|------|---------|------|
| NEP 训练集构建 / Training set construction | `deal_outcar_to_train` | ✅ |
| NEP 势训练 / NEP training | `nep_train_auto` | ✅ |
| NEP 结果可视化 / NEP result visualization | `plot_nep`, `plot_ultimate_nep` | ✅ |
| 数据集质量诊断 / Dataset quality diagnosis | `check_dataset_quality` | ✅ |
| HNEMD 热导率 / HNEMD thermal conductivity | `start_mul_hnemd`, `deal_hnemd_data` | ✅ |
| 弹性模量 / Elastic moduli | `compute_elastic_moduli` | ✅ |
| 声子谱 / Phonon spectrum | `compute_phonon_spectrum` | ✅ |
| 晶胞扩胞 / Cell expansion | `cell_expansion` | ✅ |
| GPU 自动调度 / GPU auto-scheduling | `free_gpu_run`, `gpuq` | ✅ |

### VASP 相关 / VASP Related

| 功能 / Feature | 命令 / Command | 状态 / Status |
|------|---------|------|
| 单点能批量计算 / Batch single-point | `run_all_vasp_job` | ✅ |
| AIMD 模拟 / AIMD simulation | `vasp_aimd_auto` | ✅ |
| OUTCAR → xyz 转换 | `deal_outcar_to_train` | ✅ |
| SLURM 任务提交 / SLURM job submission | `vasprun`, `vasprun_dcu` | ✅ |
| 应力-应变曲线 / Stress-strain curve | `plot_stress_strain_curve` | ✅ |
| 声子谱 (VASP+phonopy) | `ctrl_phono_specturm` | ✅ |

---

## 使用手册 / User Manual

详细的使用说明请参阅 Manual 目录 / For detailed usage instructions, see the Manual directory:

- 🇨🇳 [中文版使用手册](https://github.com/re-breath/NebulaFlow/blob/main/Manual/use_detail_Chinese.md)
- 🇺🇸 [English User Manual](https://github.com/re-breath/NebulaFlow/blob/main/Manual/use_detail_English.md)

---

## 贡献 / Contributing

本项目目前由 [rebreath](https://github.com/re-breath) 独立维护。

This project is currently maintained by [rebreath](https://github.com/re-breath).

- 发现 Bug 或有功能建议？欢迎提交 [Issue](https://github.com/re-breath/NebulaFlow/issues)
- Found a bug or have a feature request? Open an [Issue](https://github.com/re-breath/NebulaFlow/issues)
- 代码风格请遵循项目中的注释规范
- Please follow the project's comment conventions when contributing code

### 开发路线 / Roadmap

- [ ] 添加更多 CP2K 和 Lammps 支持
- [ ] 完善测试用例 / Add test cases
- [ ] 添加 CI/CD / Set up CI/CD
- [ ] 支持更多 ML 势函数格式 / Support more ML potential formats

---

<div align="center">

**🌌 NebulaFlow — 让科学计算更简单 / Making Scientific Computing Simpler**

</div>

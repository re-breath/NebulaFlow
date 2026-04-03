import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")


def read_thermo_file(filepath: str) -> pd.DataFrame:
    """
    读取 GPUMD / thermo.out 文件
    quantity T K U Pxx Pyy Pzz Pyz Pxz Pxy ax ay az bx by bz cx cy cz
    """
    filepath = Path(filepath)

    columns18 = ["T", "K", "U", "Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy", "ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]
    columns12 = ["T", "K", "U", "Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy", "ax", "by", "cz"]

    temp_df = pd.read_csv(
        filepath,
        sep=r'\s+',
        nrows=1,  # 只读第一行数据
        header=None,
        engine="python"
    )
    col_count = temp_df.shape[1]  # 获取列数


    columns = columns18 if col_count == 18 else columns12


    df = pd.read_csv(
        filepath,
        sep=r'\s+',
        names=columns,
        engine="python"
    )

    return df


def compute_derived_quantities(
    df: pd.DataFrame,
    strict: bool = True,
    volume_abs: bool = True
) -> pd.DataFrame:
    """
    从 thermo DataFrame 中提取并计算：step, T, U, Pavg, V
    ✅ 自动适配 18列 / 12列 两种 GPUMD thermo.out 格式
    
    18列：T K U Pxx~Pxy ax ay az bx by bz cx cy cz
    12列：T K U Pxx~Pxy ax by cz
    """
    # ===================== 1. 基础列校验（两种格式通用）=====================
    base_cols = ["T", "U", "Pxx", "Pyy", "Pzz"]
    missing_base = [c for c in base_cols if c not in df.columns]
    if missing_base:
        raise ValueError(f"缺少必需列：{missing_base}")

    # 复制数据并转为数值类型
    work = df.copy()
    for col in work.columns:
        work[col] = pd.to_numeric(work[col], errors="coerce")

    # 检查缺失值
    if strict and work.isna().any().any():
        raise ValueError("数据存在非法值/缺失值，请检查 thermo.out 文件")


    out = pd.DataFrame(index=work.index)
    out["step"] = np.arange(len(work), dtype=int)  # 步数
    out["T"] = work["T"]                           # 温度
    out["U"] = work["U"]                           # 势能
    out["Pavg"] = (work["Pxx"] + work["Pyy"] + work["Pzz"]) / 3  # 平均压强


    if all(c in df.columns for c in ["ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]):

        a = work[["ax", "ay", "az"]].values
        b = work[["bx", "by", "bz"]].values
        c = work[["cx", "cy", "cz"]].values
        cell = np.stack([a, b, c], axis=1)
        V = np.linalg.det(cell)

    elif all(c in df.columns for c in ["ax", "by", "cz"]):

        V = work["ax"] * work["by"] * work["cz"]

    else:
        raise ValueError("无法识别晶格列！仅支持 18列 / 12列 标准 thermo.out")


    if volume_abs:
        V = np.abs(V)
    out["V"] = V


    if strict and (not np.isfinite(out["V"]).all() or (out["V"] <= 0).any()):
        raise ValueError("体积计算异常：存在 NaN/Inf/非正体积，请检查数据")

    return out



def save_summary(data: pd.DataFrame, summary_file: str) -> None:
    stats = data[["T", "U", "Pavg", "V"]].agg(["count", "mean", "std", "min", "max"])
    with open(summary_file, "w", encoding="utf-8") as f:
        f.write("thermo.out 结果统计\n")
        f.write("=" * 60 + "\n\n")
        f.write("字段说明：\n")
        f.write("T    : 温度\n")
        f.write("U    : 能量\n")
        f.write("Pavg : 平均压强 = (Pxx + Pyy + Pzz) / 3\n")
        f.write("V    : 晶胞体积 = det([a, b, c])\n\n")
        f.write("统计结果：\n")
        f.write(stats.to_string())
        f.write("\n")


def plot_thermo(data: pd.DataFrame, fig_file: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=150)
    axes = axes.ravel()

    plots = [
        ("T", "Temperature (T)"),
        ("U", "Energy (U)"),
        ("Pavg", "Average Pressure (Pavg)"),
        ("V", "Volume (V)")
    ]

    for ax, (col, title) in zip(axes, plots):
        ax.plot(data["step"], data[col], linewidth=1.0)
        ax.set_title(title)
        ax.set_xlabel("Step index")
        ax.set_ylabel(col)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(fig_file, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Visualize thermo.out and save processed results.")
    parser.add_argument("--input", default="thermo.out", help="输入 thermo.out 文件路径")
    parser.add_argument("--csv", default="thermo_processed.csv", help="保存提取结果的 CSV 文件")
    parser.add_argument("--summary", default="thermo_summary.txt", help="保存统计结果的文本文件")
    parser.add_argument("--fig", default="thermo_plot.png", help="输出图像文件名")
    args = parser.parse_args()

    df = read_thermo_file(args.input)
    data = compute_derived_quantities(df)

    # 保存提取后的结果
    data.to_csv(args.csv, index=False)

    # 保存统计信息
    save_summary(data, args.summary)

    # 绘图
    plot_thermo(data, args.fig)

    # 屏幕打印前10行和后10行
    print("\n" + "=" * 60)
    print("提取后的主要物理量（前 10 行）")
    print("=" * 60)
    print(data.head(10).to_string(index=False))

    print("\n" + "=" * 60)
    print("提取后的主要物理量（后 10 行）")
    print("=" * 60)
    print(data.tail(10).to_string(index=False))

    print("\n" + "=" * 60)
    print("统计摘要")
    print("=" * 60)
    print(data[["T", "U", "Pavg", "V"]].agg(["mean", "std", "min", "max"]).to_string())

    print("\n输出文件：")
    print(f"1. 处理后数据: {args.csv}")
    print(f"2. 统计摘要  : {args.summary}")
    print(f"3. 可视化图片: {args.fig}")


if __name__ == "__main__":
    main()

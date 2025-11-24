#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hybridization analysis for carbon fiber systems using ASE.
Reads an .xyz file and classifies each carbon atom as sp, sp2, or sp3
based on nearest neighbors and bond angles.
Also computes orientation factor <cos^2 φ> relative to the fiber axis.
"""

import sys
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
import matplotlib.pyplot as plt


# ==========================
# 控制台参数区（可调整）
# ==========================
MAX_NEIGHBORS = 4               # 只取最近的前4个邻居
ANGLE_TOL = 15.0                # 角度容差 (度)
PRINT_DETAIL = False             # 是否打印每个原子的详细信息
# ==========================


def angle_between(v1, v2):
    """计算两个向量的夹角（度）"""
    cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cosang = np.clip(cosang, -1.0, 1.0)
    return np.degrees(np.arccos(cosang))


def classify_atom(neighbors, positions, idx):
    """根据邻居数和角度判断杂化类型"""
    n_neighbors = len(neighbors)
    if n_neighbors == 2:
        v1 = positions[neighbors[0]] - positions[idx]
        v2 = positions[neighbors[1]] - positions[idx]
        ang = angle_between(v1, v2)
        if abs(ang - 180.0) < ANGLE_TOL:
            return "sp"
    elif n_neighbors == 3:
        angles = []
        for i in range(3):
            for j in range(i+1, 3):
                v1 = positions[neighbors[i]] - positions[idx]
                v2 = positions[neighbors[j]] - positions[idx]
                angles.append(angle_between(v1, v2))
        avg_ang = np.mean(angles)
        if abs(avg_ang - 120.0) < ANGLE_TOL:
            return "sp2"
    elif n_neighbors == 4:
        angles = []
        for i in range(4):
            for j in range(i+1, 4):
                v1 = positions[neighbors[i]] - positions[idx]
                v2 = positions[neighbors[j]] - positions[idx]
                angles.append(angle_between(v1, v2))
        avg_ang = np.mean(angles)
        if abs(avg_ang - 109.5) < ANGLE_TOL:
            return "sp3"
    return "unknown"


def find_fiber_axis(positions):
    """通过PCA自动识别纤维轴方向"""
    # 去中心化
    centered = positions - np.mean(positions, axis=0)
    cov = np.cov(centered.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # 最大特征值对应的方向
    axis = eigvecs[:, np.argmax(eigvals)]
    return axis / np.linalg.norm(axis)


def main():
    if len(sys.argv) < 2:
        print("用法: python hybridization_carbon.py input.xyz")
        sys.exit(1)

    xyz_file = sys.argv[1]
    atoms = read(xyz_file)
    cutoffs = natural_cutoffs(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    positions = atoms.get_positions()
    results = []

    # 自动识别纤维轴
    fiber_axis = find_fiber_axis(positions)
    print(f"Detected fiber axis: {fiber_axis}")

    cos2_list = []

    for i, atom in enumerate(atoms):
        if atom.symbol != "C":
            continue
        indices, offsets = nl.get_neighbors(i)
        dists = [atoms.get_distance(i, j) for j in indices]
        sorted_neighbors = [j for _, j in sorted(zip(dists, indices))[:MAX_NEIGHBORS]]

        hybrid = classify_atom(sorted_neighbors, positions, i)
        results.append((i, hybrid))

        # 取向因子计算：对每个键向量计算 cos²φ
        for j in sorted_neighbors:
            v = positions[j] - positions[i]
            v /= np.linalg.norm(v)
            cos_phi = np.dot(v, fiber_axis)
            cos2_list.append(cos_phi**2)

    # 总结统计
    counts = {"sp":0, "sp2":0, "sp3":0, "unknown":0}
    for _, h in results:
        counts[h] += 1

    print("\n=== Hybridization Summary ===")
    for k,v in counts.items():
        print(f"{k}: {v}")

    # 结晶率
    total_carbon = counts["sp"] + counts["sp2"] + counts["sp3"]
    crystal_rate = total_carbon / len(atoms) * 100
    print(f"crystallization ratio: {crystal_rate:.2f}%")

    # 取向程度
    cos2_avg = np.mean(cos2_list)
    print(f"<cos²φ> = {cos2_avg:.3f}")
    # 绘制直方图
    plt.figure(figsize=(6,4))
    plt.hist(cos2_list, bins=50, range=(0,1), color="steelblue", edgecolor="black")
    plt.axvline(cos2_avg, color="red", linestyle="--", label=f"average = {cos2_avg:.3f}")
    plt.xlabel("cos²φ")
    plt.ylabel("frequency")
    plt.title("The cos²φ distribution of the carbon fiber system")
    plt.legend()
    plt.tight_layout()
    plt.savefig("cos²φ_distribute.png")

if __name__ == "__main__":
    main()

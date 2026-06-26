"""
NEP训练集 — 按平均位力筛选构型 / Filter training set by average virial
==========================================================================
该脚本使用来去除平均位力（对角线分量）不在指定范围内的构型。
从xyz文件中提取每个构型的位力张量对角元(xx, yy, zz)，除以原子数得到
每原子平均位力，筛选出三个分量均在 [min_virial, max_virial] 范围内的构型。

功能 / What it does:
  - 读取extxyz训练集文件 (使用ASE)
  - 从注释行正则提取位力对角分量
  - 检查 xx/yy/zz 是否均在范围内
  - 输出 virial_fit.xyz 和 vir_unfit.xyz
  - 生成位力分布直方图 virial_distribution.png (含 x/y/z 分量)

使用场景 / When to use:
  训练NEP前剔除极端应力状态的构型（如过度压缩/拉伸的构型位力会异常）。
  Before NEP training, remove configurations with extreme stress states.

使用方法 / Usage:
  python elect_rely_virial.py <train.xyz> <min_virial> <max_virial>

示例 / Example:
  python elect_rely_virial.py train.xyz -50 50
  # 筛选每原子位力对角分量在[-50, 50] eV/atom范围内的构型
  # Output: virial_fit.xyz, vir_unfit.xyz, virial_distribution.png

依赖 / Dependencies: numpy, ase, matplotlib
Author: rebreath
"""
# 该脚本使用来去除位力不在指定范围内的构型
import os
import numpy as np
import ase.io
import re
import matplotlib.pyplot as plt
import sys
def read_atoms_num(atoms):
    nums=[]
    for i in range(len(atoms)):
        nums.append(len(atoms[i]))
    return nums

def read_tot_virial_line(line):
    """
    从文件的一行中提取 virial tensor 的对角线分量。
    """
    pattern = re.compile(r"irial=\"([-?\d+\.\d+]+)\s+[-?\d+\.\d+]+\s+[-?\d+\.\d+]+\s+[-?\d+\.\d+]+\s+([-?\d+\.\d+]+)\s+[-?\d+\.\d+]+\s+[-?\d+\.\d+]+\s+[-?\d+\.\d+]+\s+([-?\d+\.\d+]+)\"")
    match = pattern.search(line)
    if match:
        # 提取三个对角线分量
        virial_diagonal = [float(match.group(i)) for i in range(1, 4)]
        return virial_diagonal
    return None

def read_virial(file_name):
    """计算每个构型的平均 virial。"""
    virials = []
    try:
        with open(file_name, "r") as file:
            for line in file:
                virial = read_tot_virial_line(line)
                if virial is not None:
                    virials.append(virial)
    except FileNotFoundError:
        print(f"文件 {file_name} 不存在。")
    except IOError:
        print(f"无法读取文件 {file_name}。")
    return virials  

def get_atom_avg_virial(file_name,atoms):
    """计算每个构型的平均 virial。"""
    virials=np.array(read_virial(file_name))
    if virials is None or len(virials) == 0:
        print("警告：未能读取到有效的位力数据。")
        return None
    atoms_num=np.array(read_atoms_num(atoms))
    #print("原子的数量分布",len(atoms_num))
    valid_indices = atoms_num != 0
    if len(virials) != len(atoms_num):
        print("警告：virials 和 atoms_num 的长度不匹配。")
        return None
    virials = np.array(virials)[valid_indices]
    #print(virials)
    atoms_num = np.array(atoms_num)[valid_indices]

    virials = virials[valid_indices]  # 确保这一步之后 virials 形状仍是 (N, 3)
    atoms_num = atoms_num[valid_indices][:, np.newaxis]  # 将 atoms_num 重塑为 (N, 1) 形状
    avg_virial = virials / atoms_num
    return avg_virial

def filter_virials(file_name, atoms, min_virial, max_virial):
    """
    根据对角线位力分量的平均值筛选出任何超出指定阈值的构型。
    """
    ave_virials = get_atom_avg_virial(file_name, atoms)  
    if ave_virials is None:
        print("警告：无法获取平均位力数据。")
        return []
        
    unfit_indices = []

    for i, virial_avg in enumerate(ave_virials):
        # 如果平均位力的任何一个分量超出阈值，则记录该构型索引
        if any(v < min_virial or v > max_virial for v in virial_avg):
            unfit_indices.append(i)
            
    return unfit_indices 


electfile=sys.argv[1]

min_virial=float(sys.argv[2])
max_virial=float(sys.argv[3])

atoms = ase.io.read(electfile, index=":")
atoms_num = read_atoms_num(atoms)


#检查位力不再在指定范围内的构型
vir=read_virial(electfile)
#print("位力的分布",len(vir))
av_virial = get_atom_avg_virial(electfile,atoms)
#print("位力的分布",av_virial)

Vir_unfit_index = filter_virials(electfile, atoms, min_virial,max_virial)
if Vir_unfit_index:  # 如果不为空
    Vir_unfit_index_avg = [av_virial[i] for i in Vir_unfit_index]
    # print("位力不在指定范围内的构型索引：", Vir_unfit_index)
    # print("位力不在指定范围内的构型平均位力：", Vir_unfit_index_avg)


if os.path.exists("virial_fit.xyz"):
    os.remove("virial_fit.xyz")

if os.path.exists("vir_unfit.xyz"):
    os.remove("vir_unfit.xyz")

fit_Vir_rational = [atoms[i] for i in range(len(atoms)) if i not in Vir_unfit_index]
ase.io.write("virial_fit.xyz", fit_Vir_rational)


unfit_Vir = [atoms[i] for i in  Vir_unfit_index]
ase.io.write("vir_unfit.xyz", unfit_Vir)
# print("Number of configs : ",len(atoms))
# print(f"Average enery range : {min(av_energy)}  {max(av_energy)}")
# print("Enery fit number : ",len(fit_E))
# print("Enery unfit number : ",len(unfit_E))

all_av_virs = np.array(av_virial).flatten()
print("Number of configs : ",len(atoms))
print(f"Average virial range : {min(all_av_virs)}  {max(all_av_virs)}")
print("Number of configs fit virial : ",len(fit_Vir_rational))
print("Number of configs unfit virial : ",len(unfit_Vir))


alpha = 0.5
my_favorite_colors = ['#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd', '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff']
x_axie = np.arange(len(av_virial))
plt.hist(av_virial[:,1], bins=50, color=my_favorite_colors[4],label="virial_y",alpha=alpha)
plt.hist(av_virial[:,0], bins=50, color=my_favorite_colors[0],label="virial_x",alpha=alpha)
plt.hist(av_virial[:,2], bins=50, color=my_favorite_colors[7],label="virial_z",alpha=alpha)

plt.xlabel("Virial (ev/atom)")
plt.ylabel("Frequency")
plt.grid(True)
plt.legend()
plt.savefig("virial_distribution.png")



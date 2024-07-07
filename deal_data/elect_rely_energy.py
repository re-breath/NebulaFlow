#该脚本使用来去除位力，平均能量大于一定值的构型
import os
import numpy as np
import ase.io
import re
import sys
import matplotlib.pyplot as plt

def read_atoms_num(atoms):
    nums=[]
    for i in range(len(atoms)):
        nums.append(len(atoms[i]))
    return nums
def read_energy_line(line):
    """
    从文件的一行中提取能量值。
    """
    pattern = re.compile(r"nergy\s*=\s*(-?\d+\.\d+)")
    match = pattern.search(line)
    if match:
        return float(match.group(1))
    return None        

def read_xyz_energy(file_name):
    """
    从给定的 .xyz 文件中读取并返回所有能量值的列表。
    """
    energies = []
    try:
        with open(file_name, "r") as file:
            for line in file:
                energy = read_energy_line(line)
                if energy is not None:
                    energies.append(energy)
    except FileNotFoundError:
        print(f"文件 {file_name} 不存在。")
    except IOError:
        print(f"无法读取文件 {file_name}。")
    return energies
def get_atom_avg_energy(file_name,atoms):
    """计算每个构型的平均能量。"""
    energy=np.array(read_xyz_energy(file_name))
    atoms_num=np.array(read_atoms_num(atoms))

    valid_indices = atoms_num != 0
    energy = energy[valid_indices]
    atoms_num = atoms_num[valid_indices]
    return energy/atoms_num
def elect_from_av_E(av_energy,min,max):
    """检查能量的平均值不在指定范围内的构型。"""
    unfit_index=[]
    for i in range(len(atoms)):
        if min>av_energy[i] or av_energy[i] > max:
            unfit_index.append(i)
    return unfit_index

electfile = sys.argv[1]
min_E = float(sys.argv[2])
max_E = float(sys.argv[3])


# electfile="train.xyz"


atoms = ase.io.read(electfile, index=":")
atoms_num = read_atoms_num(atoms)

#检查能量不再在指定范围内的构型

av_energy = get_atom_avg_energy(electfile,atoms)
E_unfit_index = elect_from_av_E(av_energy, min_E, max_E)
E_unfit_index_avg = [av_energy[i] for i in E_unfit_index]
# print("能量不在指定范围内的构型索引：", E_unfit_index)
# print("能量不在指定范围内的构型平均能量：", E_unfit_index_avg)


#输出筛选后的文件train_E_Vir_rational.xyz 能量不合格的文件E_unfit.xyz 位力不合格的文件Vir_unfit.xyz
if os.path.exists("energy_fited.xyz"):
    os.remove("energy_fited.xyz")
if os.path.exists("energy_unfited.xyz"):
    os.remove("energy_unfited.xyz")

fit_E = [atoms[i] for i in range(len(atoms)) if i not in E_unfit_index]
ase.io.write("energy_fited.xyz", fit_E)

unfit_E = [atoms[i] for i in  E_unfit_index]
ase.io.write("energy_unfited.xyz", unfit_E)

# fit_energy = [av_energy[i] for i in range(len(atoms)) if i not in E_unfit_index]

alpha = 0.5
my_favorite_colors = ['#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd', '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff']

plt.hist(av_energy, bins=100, color=my_favorite_colors[3],label="all_energy",alpha=alpha)
# plt.hist(fit_energy, bins=50, color=my_favorite_colors[2],label="fit_energy",alpha=alpha)
print("Number of configs : ",len(atoms))
print(f"Average enery range : {min(av_energy)}  {max(av_energy)}")
print("Enery fit number : ",len(fit_E))
print("Enery unfit number : ",len(unfit_E))

plt.xlabel("Energy (ev/atom)")
plt.ylabel("Frequency")
plt.grid(True)
plt.legend()
plt.savefig("energy_distribution.png")



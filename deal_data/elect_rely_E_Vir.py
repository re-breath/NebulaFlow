#该脚本使用来去除位力，平均能量大于一定值的构型
import os
import numpy as np
import ase.io
import re
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
    print("原子的数量分布",len(atoms_num))
    valid_indices = atoms_num != 0
    if len(virials) != len(atoms_num):
        print("警告：virials 和 atoms_num 的长度不匹配。")
        return None
    virials = np.array(virials)[valid_indices]
    print(virials)
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
def elect_from_av_E(av_energy,min,max):
    """检查能量的平均值不在指定范围内的构型。"""
    unfit_index=[]
    for i in range(len(atoms)):
        if min>av_energy[i] or av_energy[i] > max:
            unfit_index.append(i)
    return unfit_index


electfile="all_wz.xyz"


atoms = ase.io.read(electfile, index=":")
atoms_num = read_atoms_num(atoms)

#检查能量不再在指定范围内的构型

av_energy = get_atom_avg_energy(electfile,atoms)
E_unfit_index = elect_from_av_E(av_energy, -8, -6)
E_unfit_index_avg = [av_energy[i] for i in E_unfit_index]
print("能量不在指定范围内的构型索引：", E_unfit_index)
print("能量不在指定范围内的构型平均能量：", E_unfit_index_avg)

#检查位力不再在指定范围内的构型
vir=read_virial(electfile)
print("位力的分布",len(vir))
av_virial = get_atom_avg_virial(electfile,atoms)
print("位力的分布",av_virial)

Vir_unfit_index = filter_virials(electfile, atoms, -5,5)
if Vir_unfit_index:  # 如果不为空
    Vir_unfit_index_avg = [av_virial[i] for i in Vir_unfit_index]
    print("位力不在指定范围内的构型索引：", Vir_unfit_index)
    print("位力不在指定范围内的构型平均位力：", Vir_unfit_index_avg)

#输出筛选后的文件train_E_Vir_rational.xyz 能量不合格的文件E_unfit.xyz 位力不合格的文件Vir_unfit.xyz
if os.path.exists("train_E_Vir_rational.xyz"):
    os.remove("train_E_Vir_rational.xyz")
if os.path.exists("E_unfit.xyz"):
    os.remove("E_unfit.xyz")
if os.path.exists("Vir_unfit.xyz"):
    os.remove("Vir_unfit.xyz")

fit_E_Vir_rational = [atoms[i] for i in range(len(atoms)) if i not in E_unfit_index and i not in Vir_unfit_index]
ase.io.write("train_E_Vir_rational.xyz", fit_E_Vir_rational)

unfit_E = [atoms[i] for i in  E_unfit_index]
ase.io.write("E_unfit.xyz", unfit_E)

unfit_Vir = [atoms[i] for i in  Vir_unfit_index]
ase.io.write("Vir_unfit.xyz", unfit_Vir)



import matplotlib.pyplot as plt
#画出位力与能量的分布图
#控制透明度
alpha = 0.5
my_favorite_colors = ['#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd', '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff']
x_axie = np.arange(len(av_virial))
plt.hist(av_virial[:,1], bins=50, color=my_favorite_colors[4],label="virial_y",alpha=alpha)
plt.hist(av_virial[:,0], bins=50, color=my_favorite_colors[0],label="virial_x",alpha=alpha)
plt.hist(av_virial[:,2], bins=50, color=my_favorite_colors[7],label="virial_z",alpha=alpha)
plt.hist(av_energy, bins=50, color=my_favorite_colors[3],label="energy",alpha=alpha)
plt.xlabel("Virial and Energy (ev/atom)")
plt.ylabel("Frequency")
plt.grid(True)
plt.legend()
plt.savefig("virial_energy_distribution.png")



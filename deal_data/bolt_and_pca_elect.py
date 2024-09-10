from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
import ase.io
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import os
def sample_by_energy_intervals(ave_energy, probabilities, des, num_intervals, num_samples):
    """
    按照能量区间采样构型，并确保描述符间的距离尽量大。
    """
    # 将平均能量排序并分为num_intervals个区间
    energy_bounds = np.linspace(min(ave_energy), max(ave_energy), num_intervals+1)
    selected_indices = []
    
    for i in range(num_intervals):
        # 确定当前区间内的构型
        interval_indices = np.where((ave_energy >= energy_bounds[i]) & (ave_energy < energy_bounds[i+1]))[0]
        if len(interval_indices) == 0: continue
        
        # 计算区间内的构型应该采样的数量
        interval_prob_sum = probabilities[interval_indices].sum()
        interval_sample_count = int(round(interval_prob_sum * num_samples))
        
        # 如果区间内的构型不够采样，则全部选择
        if interval_sample_count >= len(interval_indices) or interval_sample_count == 0:
            selected_indices.extend(interval_indices)
            continue

        # 对区间内的构型进行PCA和最远点采样
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
def get_atom_avg_energy(file_name):
    """计算每个构型的平均能量。"""
    energy=np.array(read_xyz_energy(file_name))
    atoms_num=np.array(read_atoms_num(atoms))

    valid_indices = atoms_num != 0
    energy = energy[valid_indices]
    atoms_num = atoms_num[valid_indices]
    return energy/atoms_num
def calculate_probabilities(energies):
    """计算每个构型的概率."""
    energies = np.array(energies)
    probabilities = np.exp(-energies) / np.sum(np.exp(-energies))  # 使用Boltzmann分布
    return probabilities

def sample_configurations(energies, probabilities, num_samples=600):
    """根据概率采样构型.在本例子中不采用该方法"""
    indices = np.arange(len(energies))
    sampled_indices = np.random.choice(indices, size=num_samples, p=probabilities,replace=False)
    return sampled_indices

def get_descriptors(nepfile,trainfilename="train.xyz"):
    a = ase.io.read(file_name, ':')
    calc = NEP(nepfile)
    print(calc)
    des = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in a])
    return des

if __name__=="__main__":
    
    # file_name=sys.argv[1]   #指定进行挑选的文件
    # num_samples=int(sys.argv[2])  #指定采样数量

    file_name = "train_all.xyz"
    num_samples=900
    num_intervals=30
    nepfile="nep_all6.txt"

    atoms=ase.io.read(file_name,index=":")
    energy=read_xyz_energy(file_name)
    atoms_num=read_atoms_num(file_name)
    ave_energy=get_atom_avg_energy(file_name)
    des=get_descriptors(nepfile,file_name)
    if (num_samples>len(ave_energy) or num_samples <= 0 ):
        print(f"采样数量不符合要求，应该在1到{len(ave_energy)}范围内采样")
        exit()
    probabilities = calculate_probabilities(ave_energy)
    #print(probabilities)

    sampled_indices = sample_by_energy_intervals(ave_energy, probabilities, des, num_intervals, num_samples)

    unelected_indices=list(set(range(len(ave_energy)))-set(sampled_indices))
    print(f"选中的构型索引为：{sampled_indices}")
    print(f"未选中的构型索引为：{unelected_indices}")

    #检查是否存在bolt_elect.xyz，bolt_unelect.xyz文件，如果存在则删除
    if os.path.exists("bolt_elect.xyz"):
        os.remove("bolt_elect.xyz")
    if os.path.exists("bolt_unelect.xyz"):
        os.remove("bolt_unelect.xyz")

    for i in sampled_indices:
        ase.io.write(f"elect_bolt_pca.xyz",atoms[i],append=True)
    for j in unelected_indices:
        ase.io.write(f"bolt_pca_unelect.xyz",atoms[j],append=True)
    print(f"选中的构型数量为{len(sampled_indices)}")
    print(f"未选中的构型数量为{len(unelected_indices)}")

    #输出选择的构型的能量分布图
    
    favorite_colors=['#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd', '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff']

    plt.switch_backend('Agg')
    plt.hist(ave_energy,bins=100,alpha=0.4,color=favorite_colors[1],label='All')

    atoms=ase.io.read('elect_bolt_pca.xyz',index=":")
    ave_energy=get_atom_avg_energy('elect_bolt_pca.xyz')
    plt.hist(ave_energy,bins=100,alpha=0.4,color=favorite_colors[3],label='Elect')
    plt.title('Histogram of Average Energy Values')  
    plt.legend()
    plt.xlim(-7.5, -5)
    plt.xlabel('Average Energy')  
    plt.ylabel('Probability Density') 
    plt.grid(axis='y', alpha=0.75) 
    plt.savefig('elect_bolt_and_pca.png')
"""
NebulaFlow Python 核心库 — xyz 构型数据读写
==============================================
专门处理 GPUMD 风格的 extended xyz 格式构型数据，建立原子构型对象（Config）。
提供高效的读取、写入、统计和信息提取功能。

Core library for reading/writing GPUMD-style extended xyz trajectory files.
Provides the Config class for representing atomic configurations and
the read_xyz / write_xyz functions for I/O.

主要功能 / Key Features:
  - read_xyz:  读取xyz文件，返回Config对象列表 / Read xyz into list of Config objects
  - write_xyz: 将Config对象写回xyz文件 / Write Config objects back to xyz
  - read_thermo: 解析GPUMD的thermo.out热力学输出 / Parse GPUMD thermo output
  - Config类: 提供atom_num, lattice, energy, virial, position, force, atom_type属性

典型使用场景 / Typical Usage:
  1. 读取NEP训练集并筛选构型: read train.xyz → filter configs → write back
  2. 分析MD轨迹的能量演化: read dump.xyz → extract energies → plot
  3. 提取特定构型的受力信息: read xyz → config.force → analysis

依赖 / Dependencies:
  - numpy, re (标准库)

Author: rebreath
"""
# 该库为专门处理xyz类型的构型数据，建立原子构型的对象使用
# 主要应用为引入read_xyz函数，读取xyz文件，得到原子构型的对象列表，并提供相应的属性和方法。
import numpy as np
import re


compile_num = re.compile(r'^\d+\n') 
compile_energy = re.compile(r'Energy\s*=\s*(-?\d+\.\d+)[eE][Vv]') 
compile_virial = re.compile(r'irial\s*=\s*\"\s*([-?(\d+)?\.\d+\s+]+)\"')
compile_lattice = re.compile(r'attice\s*=\s*\"\s*([-\d.\s]+)\"')

class Config:
    """
    原子构型对象 / Atomic Configuration Object
    ============================================
    该类为原子构型的对象，包含了原子数，晶格常数，能量，晶胞，原子坐标，原子力，原子类型。
    Represents one configuration (frame) from a GPUMD extended xyz file.

    属性 / Attributes:
        atom_num (int):    原子总数 / Total number of atoms
        lattice (list):    晶格常数 (9元素) / Lattice constants in order [ax,ay,az,bx,by,bz,cx,cy,cz]
        energy (float):    总能量 (eV) / Total energy in eV
        virial (ndarray):  位力张量 (9分量) / Virial tensor [xx,xy,xz,yx,yy,yz,zx,zy,zz]
        position (ndarray): 原子坐标 (Nx3) / Atomic positions in Angstrom
        force (ndarray):   原子受力 (Nx3) / Atomic forces in eV/Angstrom
        atom_type (list):  原子元素符号列表 / List of atomic species e.g. ['Si','Si','O',...]

    使用场景 / When to use:
        读取NEP训练集后需要按能量/受力筛选构型、提取特定构型分析、修改晶格参数等。

    使用示例 / Usage:
        >>> configs = read_xyz('train.xyz')
        >>> config = configs[0]
        >>> print(config.atom_num)   # 64
        >>> print(config.energy)     # -1234.56
        >>> print(config.lattice)    # [5.43, 0.0, 0.0, 0.0, 5.43, 0.0, 0.0, 0.0, 5.43]
    """
    def __init__(self,lines,index1,index2):
        self.lines = lines
        self.index1 = index1
        self.index2 = index2
        self.atom_num = int(self.lines[self.index1].strip())
        self.second_row = self.lines[self.index1 + 1].strip()
        self.lattice,self.virial,self.energy = self.analy_second_row()
        self.position = self.get_atom_position()
        self.force = self.get_atom_force()
        self.atom_type = self.get_atom_type()

    # 
    def analy_second_row(self):
        """分析 second_row 数据，并提取 lattice, virial 和 energy 参数"""
        lattice = None  
        virial = None  
        energy = None  
        match_lattice = compile_lattice.search(self.second_row)
        match_virial = compile_virial.search(self.second_row)
        match_energy = compile_energy.search(self.second_row)
        if match_lattice:
            lattice_str = match_lattice.group(1)
            lattice = [float(i) for i in lattice_str.split()]
        if match_virial:
            virial_str = match_virial.group(1)
            virial = np.array([float(i) for i in virial_str.split()])
        if match_energy:
            energy = float(match_energy.group(1))  
        return lattice,virial,energy
    def get_atom_type(self):
        """获得原子类型"""
        atom_type = []
        for i in range(self.index1+2,self.index2):
            columns = self.lines[i].strip().split()
            if len(columns) > 0:
                atom_type.append(columns[0])
            
        return atom_type
    
    def get_atom_position(self):
        """获得原子坐标"""
        position = []
        if self.index2 == -1:
            limit = len(self.lines)
        else:
            limit = self.index2 + 1

        for i in range(self.index1+2,limit):  
            data = self.lines[i].strip().split()[1:4]
            if len(data) == 3:  
                position.append(data)              
        position = np.array(position,dtype=float)
        return position
    
    def get_atom_force(self):
        """获得原子力"""
        force = []
        if self.index2 == -1:
            limit = len(self.lines)
        else:
            limit = self.index2 + 1
        for i in range(self.index1+2,limit):
            if len(self.lines[i].strip().split()) <7:
                continue
            else:
                force.append(self.lines[i].strip().split()[4:7])
        force = np.array(force,dtype=float)
        return force
    

def get_index(lines):
    """获得原子的索引"""
    config_row = []

    for line_num,line in enumerate(lines):
        if compile_num.match(line):
            config_row.append(line_num)
    config_row.append(-1)
    return config_row
def read_xyz(filename:str='train.xyz') -> list :
    """
    读取xyz文件，得到原子构型Config对象的列表
    Read an extended xyz file and return a list of Config objects.

    参数 / Args:
        filename: GPUMD风格的xyz文件路径 / Path to GPUMD-style xyz file.

    返回 / Returns:
        list: Config对象列表 / List of Config objects, one per frame.

    使用场景 / When to use:
        NEP训练集加载、MD轨迹分析、构型筛选等的第一步操作。
        The first step in any xyz-based analysis pipeline.

    使用示例 / Usage:
        >>> configs = read_xyz('train.xyz')
        >>> len(configs)        # 构型总数
        2500
        >>> configs[0].energy   # 第一个构型的能量
        -1234.567
        >>> configs[-1].atom_num  # 最后一个构型的原子数
        64
    """
    global read_xyz_filename
    read_xyz_filename = filename
    with open(filename, 'r') as file:
        lines = file.readlines()
    config_row = get_index(lines)
    configs = [Config(lines,config_row[i],config_row[i+1]) for i in range(len(config_row)-1)]
    return configs
    
def write_xyz(filename:str,configs):
    """ase.io中write函数有bug会导致一些构型力没有写入,该函数专门用于写入xyz格式的构型数据"""
    if read_xyz_filename:
        pass
    else:
        print("Warning: The xyz file has been read, the write function will not work correctly.")
        return
    if isinstance(configs,list):
        write_xyz_list(filename,configs)
    elif isinstance(configs,Config):
        write_xyz_config(filename,configs)
    else:
        print("Warning: The input is not a Config object or a list of Config objects.")
    

def write_xyz_list(filename:str,configlist):
    with open(read_xyz_filename, 'r') as f:
            lines = f.readlines()
            with open(filename, 'w') as w:
                for i in configlist:
                    w.writelines(lines[i.index1:i.index2])

def write_xyz_config(filename:str,config):
    with open(read_xyz_filename, 'r') as f:
            lines = f.readlines()
            with open(filename, 'w') as w:
                w.writelines(lines[config.index1:config.index2])

def read_thermo(filename):
    """
    读取GPUMD的输出文件thermo.out / Read GPUMD thermo output file.

    参数 / Args:
        filename: thermo.out文件路径 / Path to thermo.out.

    返回 / Returns:
        dict: 包含T, K, U, Px...等物理量的字典，每个值都是numpy数组。
              Dictionary with keys T, K, U, Px, Py, Pz... mapping to numpy arrays.

    列宽说明 / Column layout:
        12列: T, K, U, Px, Py, Pz, Pyz, Pxz, Pxy, Lx, Ly, Lz  (正交盒子)
        18列: 同上 + ax,ay,az, bx,by,bz, cx,cy,cz              (三斜盒子)

    使用场景 / When to use:
        MD模拟完成后分析温度/能量/压力/晶格演化。
        Analyze temperature, energy, pressure, and lattice evolution after MD.

    使用示例 / Usage:
        >>> thermo = read_thermo('thermo.out')
        >>> T = thermo['T']    # 温度数组
        >>> P = thermo['Px']   # x方向压强
        >>> print(f"Average T: {np.mean(T):.1f} K")
    """
    data = np.loadtxt(filename)
    thermo = dict()
    if data.shape[1] == 12:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','Lx','Ly','Lz']
        for i in range(12):
            thermo[labels[i]] = data[:,i]
    if data.shape[1] == 18:
        labels = ['T','K','U','Px','Py','Pz','Pyz','Pxz','Pxy','ax','ay','az','bx','by','bz','cx','cy','cz']
        for i in range(18):
            thermo[labels[i]] = data[:,i]
    return thermo


if __name__ == '__main__':
    configs = read_xyz('train.xyz')
    print(len(configs))
    print(configs[1861].atom_num)
    print(configs[1861].lattice)
    print(configs[1861].virial)
    print(configs[1861].energy)
    print(type(configs[1861].second_row))
    print(configs[1861].position)
    print(configs[1861].force)
    list_configs = [configs[0],configs[1861]]
    write_xyz('train_w.xyz',list_configs)
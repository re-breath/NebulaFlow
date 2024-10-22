#该库为专门处理xyz类型的构型数据，建立原子构型的对象使用
#主要应用为引入read_xyz函数，读取xyz文件，得到原子构型的对象列表，并提供相应的属性和方法。
import numpy as np
import re


compile_num = re.compile(r'^\d+\n') 
compile_energy = re.compile(r'Energy\s*=\s*(-?\d+\.\d+)[eE][Vv]') 
compile_virial = re.compile(r'irial\s*=\s*\"\s*([-?(\d+)?\.\d+\s+]+)\"')
compile_lattice = re.compile(r'attice\s*=\s*\"\s*([-\d.\s]+)\"')

class Config:
    """该类为原子构型的对象，包含了原子数，晶格常数，能量，晶胞，原子坐标，原子力，原子类型"""
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
    """读取xyz文件，得到原子构型的对象列表"""
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
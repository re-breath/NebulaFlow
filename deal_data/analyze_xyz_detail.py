import numpy as np
import re


def read_xyz_with_lattice(filename):
    """
    读取包含晶格信息的XYZ文件
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 读取原子总数
    try:
        num_atoms = int(lines[0].strip())
    except:
        print("错误: 第一行应该包含原子总数")
        return None
    
    # 解析第二行的晶格信息
    lattice_info = lines[1].strip()
    
    # 提取晶格参数 (格式: Lattice="a1 a2 a3 a4 a5 a6 a7 a8 a9")
    lattice_match = re.search(r'Lattice="([^"]+)"', lattice_info)
    if lattice_match:
        lattice_str = lattice_match.group(1)
        lattice_values = list(map(float, lattice_str.split()))
    else:
        print("警告: 未找到晶格参数")
        lattice_values = None
    
    # 解析周期性边界条件
    pbc_match = re.search(r'pbc="([^"]+)"', lattice_info)
    if pbc_match:
        pbc_str = pbc_match.group(1)
        pbc = [x == 'T' for x in pbc_str.split()]
    else:
        pbc = [False, False, False]
    
    # 解析原子数据
    atoms = []
    elements = []
    positions = []
    
    for i in range(2, min(2 + num_atoms, len(lines))):
        line = lines[i].strip()
        if not line:
            continue
        
        parts = line.split()
        if len(parts) >= 4:
            element = parts[0]
            pos = list(map(float, parts[1:4]))
            
            atoms.append({
                'element': element,
                'position': pos
            })
            elements.append(element)
            positions.append(pos)
    
    return {
        'num_atoms': num_atoms,
        'lattice_params': lattice_values,
        'pbc': pbc,
        'atoms': atoms,
        'elements': elements,
        'positions': np.array(positions)
    }

def analyze_structure(data):
    """
    分析结构信息
    """
    if data is None:
        return
    
    print("=" * 50)
    print("结构分析报告")
    print("=" * 50)
    
    # 1. 基本信息
    print(f"\n1. 基本信息:")
    print(f"   - 原子总数: {data['num_atoms']}")
    
    # 统计元素种类和数量
    elements_count = {}
    for atom in data['atoms']:
        elem = atom['element']
        elements_count[elem] = elements_count.get(elem, 0) + 1
    
    print(f"   - 元素种类: {len(elements_count)} 种")
    for elem, count in elements_count.items():
        print(f"     * {elem}: {count} 个原子 ({count/data['num_atoms']*100:.1f}%)")
    
    # 2. 晶格信息
    print(f"\n2. 晶格信息:")
    if data['lattice_params']:
        # 将9个参数转换为3x3矩阵
        lattice_matrix = np.array(data['lattice_params']).reshape(3, 3)
        print(f"   - 晶格矩阵 (Å):")
        for i in range(3):
            print(f"     [{lattice_matrix[i,0]:10.6f} {lattice_matrix[i,1]:10.6f} {lattice_matrix[i,2]:10.6f}]")
        
        # 计算晶格长度
        a = np.linalg.norm(lattice_matrix[0])
        b = np.linalg.norm(lattice_matrix[1])
        c = np.linalg.norm(lattice_matrix[2])
        
        print(f"\n   - 晶格长度:")
        print(f"     a = {a:.4f} Å")
        print(f"     b = {b:.4f} Å")
        print(f"     c = {c:.4f} Å")
        
        # 计算晶格角度
        def angle(v1, v2):
            return np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
        
        alpha = angle(lattice_matrix[1], lattice_matrix[2])
        beta = angle(lattice_matrix[0], lattice_matrix[2])
        gamma = angle(lattice_matrix[0], lattice_matrix[1])
        
        print(f"\n   - 晶格角度:")
        print(f"     α (b-c夹角) = {alpha:.2f}°")
        print(f"     β (a-c夹角) = {beta:.2f}°")
        print(f"     γ (a-b夹角) = {gamma:.2f}°")
        
        # 计算晶胞体积
        volume = np.abs(np.linalg.det(lattice_matrix))
        print(f"\n   - 晶胞体积: {volume:.4f} Å³")
    else:
        print("   - 无晶格参数信息")
    
    # 3. 周期性边界条件
    print(f"\n3. 周期性边界条件:")
    pbc_names = ['a方向', 'b方向', 'c方向']
    for i, (is_pbc, name) in enumerate(zip(data['pbc'], pbc_names)):
        status = "是" if is_pbc else "否"
        print(f"   - {name}周期性: {status}")
    
    # 4. 密度计算
    print(f"\n4. 密度计算:")
    if data['lattice_params'] and data['atoms']:
        # 原子质量 (单位: g/mol)
        atomic_masses = {
            'C': 12.011,   # 碳
            'H': 1.008,    # 氢
            'O': 15.999,   # 氧
            'N': 14.007,   # 氮
            'Si': 28.085,  # 硅
            'Ge': 72.63,   # 锗
            # 可以根据需要添加更多元素
        }
        
        # 计算总质量
        total_mass = 0.0
        for atom in data['atoms']:
            elem = atom['element']
            if elem in atomic_masses:
                total_mass += atomic_masses[elem]
            else:
                print(f"   警告: 未知元素 '{elem}'，使用默认质量 12.0 g/mol")
                total_mass += 12.0
        
        # 阿伏伽德罗常数
        avogadro = 6.02214076e23
        
        # 单位转换: g/mol -> g (每个晶胞)
        mass_per_cell = total_mass / avogadro  # 克
        
        # 体积转换: Å³ -> cm³
        volume_cm3 = volume * 1e-24
        global density
        
        # 计算密度
        if volume_cm3 > 0:
            density = mass_per_cell / volume_cm3  # g/cm³
            print(f"   - 总原子质量: {total_mass:.4f} g/mol")
            print(f"   - 晶胞质量: {mass_per_cell:.4e} g")
            print(f"   - 晶胞体积: {volume_cm3:.4e} cm³")
            print(f"   - 计算密度: {density:.4f} g/cm³")

    
    # 5. 坐标范围
    print(f"\n5. 坐标范围:")
    positions = data['positions']
    if len(positions) > 0:
        min_coords = positions.min(axis=0)
        max_coords = positions.max(axis=0)
        ranges = max_coords - min_coords
        
        print(f"   - X方向: {min_coords[0]:.2f} 到 {max_coords[0]:.2f} Å (跨度: {ranges[0]:.2f} Å)")
        print(f"   - Y方向: {min_coords[1]:.2f} 到 {max_coords[1]:.2f} Å (跨度: {ranges[1]:.2f} Å)")
        print(f"   - Z方向: {min_coords[2]:.2f} 到 {max_coords[2]:.2f} Å (跨度: {ranges[2]:.2f} Å)")
    
    print("\n" + "=" * 50)
    print("分析完成")
    print("=" * 50)

def save_summary(data, output_file="structure_summary.txt"):
    """
    保存分析结果到文件
    """
    with open(output_file, 'w') as f:
        
        f.write(f"atoms {data['num_atoms']}\n")
        
        
        if data['lattice_params']:
            lattice_matrix = np.array(data['lattice_params']).reshape(3, 3)
            a = np.linalg.norm(lattice_matrix[0])
            b = np.linalg.norm(lattice_matrix[1])
            c = np.linalg.norm(lattice_matrix[2])
            volume = np.abs(np.linalg.det(lattice_matrix))
            
            f.write(f"lattice_a  {a:.4f} Å\n")
            f.write(f"lattice_b  {b:.4f} Å\n")
            f.write(f"lattice_c  {c:.4f} Å\n")
            f.write(f"volume  {volume:.4f} Å³\n")
            f.write(f"density  {density:.4f} g/cm³\n")

def main():

    import sys 

    filename = sys.argv[1]  # 替换为你的文件名
    
    print(f"正在分析文件: {filename}")
    print("-" * 50)

    # 读取文件
    data = read_xyz_with_lattice(filename)
    
    if data:
        # 分析结构
        analyze_structure(data)
        
        # 保存摘要
        save_summary(data)
        print(f"\n分析结果已保存到: structure_summary.txt")
        
        return data
    else:
        print("无法读取文件")
        return None

if __name__ == "__main__":
    # 直接运行脚本时，将你的文件内容保存为文件
    # 或者修改filename变量指向你的文件
    main()


#!/usr/bin/env python3

# 排序 LAMMPS data 文件中的原子 ID
# 输入：原始 data 文件路径
# 输出：排序后的 data 文件路径

import sys

def sort_section(lines):
    atoms = []
    for line in lines:
        if not line.strip():
            continue
        parts = line.split()
        try:
            atom_id = int(parts[0])
        except ValueError:
            atoms.append((999999999, line))
            continue
        atoms.append((atom_id, line))
    atoms.sort(key=lambda x: x[0])
    return [line for _, line in atoms]

def process_data(infile, outfile):
    with open(infile, 'r') as f:
        lines = f.readlines()

    output_lines = []
    atoms_section = []
    in_atoms = False

    for line in lines:
        if line.strip().startswith("Atoms"):
            output_lines.append(line)
            output_lines.append('\n') 
            in_atoms = True
            continue
        if in_atoms:
            if line.strip().startswith(("Velocities","Bonds","Angles","Dihedrals","Impropers")):
                # 排序并写出 Atoms 部分
                output_lines.extend(sort_section(atoms_section))
                output_lines.append('\n') 
                atoms_section = []
                in_atoms = False
                output_lines.append(line)
            else:
                atoms_section.append(line)
        else:
            output_lines.append(line)

    # 如果文件以 Atoms 部分结束
    if in_atoms:
        output_lines.extend(sort_section(atoms_section))

    with open(outfile, 'w') as f:
        f.writelines(output_lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python sort_lammps_data.py input.data output.data")
        sys.exit(1)
    process_data(sys.argv[1], sys.argv[2])
    print(f"已生成排序后的 data 文件: {sys.argv[2]}")
# ======================================================================
# File:         ioEnvFunction.sh
# Project:      NebulaFlow
# Description:  文件格式转换函数库 — xyz/POSCAR/data/cif/pdb 互转
#               File format conversion: xyz, POSCAR, LAMMPS data, cif, pdb.
# Author:       rebreath
# Dependencies: python3, ase, ovito (optional)
# ======================================================================


# ---------------------------------------------------------------------------
# Function: tran_xyz2pos  (tran_xyz2poscar)
# 功能: 将xyz文件转换为VASP POSCAR格式
# 场景: 从其他软件导出的xyz结构需要转为VASP输入格式时使用。
# Usage: tran_xyz2pos <file.xyz>
# Example:
#   tran_xyz2pos model.xyz
#   # 输出: POSCAR_convered
# ---------------------------------------------------------------------------
tran_xyz2pos() {
    python3 << EOF
import ase.io
filename='$1'; xyzinfo=ase.io.read(filename)
ase.io.write('POSCAR_convered',xyzinfo,format='vasp')
EOF
}


# ---------------------------------------------------------------------------
# Function: tran_pos2xyz
# 功能: 将VASP POSCAR转换为xyz格式（extxyz，包含Lattice信息）
# 场景: 需要将VASP结构转为GPUMD兼容的xyz格式时使用。
# Usage: tran_pos2xyz <POSCAR> [output.xyz]
# Example:
#   tran_pos2xyz POSCAR model.xyz
#   # 输出: model.xyz (或 model_conversed.xyz)
# ---------------------------------------------------------------------------
tran_pos2xyz() {
    python3 - "$1" "${2:-model_conversed.xyz}" <<'PY'
import sys; from ase.io import read, write
atoms=read(sys.argv[1]); write(sys.argv[2],atoms,format="extxyz")
print(f"Written: {sys.argv[2]}")
PY
}


# ---------------------------------------------------------------------------
# Function: tran_data2xyz
# 功能: 将LAMMPS data文件转换为xyz格式
# 场景: 使用LAMMPS建模后用NebulaFlow/GPUMD进行处理时需要转换格式。
# Usage: tran_data2xyz <file.data> [output.xyz] [atom_style]
# Example:
#   tran_data2xyz system.data model.xyz atomic
#   # 输出: model.xyz (extxyz格式)
# ---------------------------------------------------------------------------
tran_data2xyz() {
    python3 - "$1" "${2:-model_conversed.xyz}" "${3:-atomic}" <<'PY'
import sys; from ase.io import read, write
atoms=read(sys.argv[1],format="lammps-data",atom_style=sys.argv[3])
write(sys.argv[2],atoms,format="extxyz")
print(f"Written: {sys.argv[2]}")
PY
}


# ---------------------------------------------------------------------------
# Function: tran_xyz2data
# 功能: 将xyz文件转换为LAMMPS data格式
#       需指定元素顺序以正确分配原子类型编号。
# 场景: 需要将GPUMD/NEP模型转为LAMMPS输入时使用。
# Usage: tran_xyz2data <model.xyz> <output.data> <elem1> [elem2 ...]
# Example:
#   tran_xyz2data model.xyz conf.lmp C
#   tran_xyz2data model.xyz conf.lmp O H    # 多元素体系
#   # 输出: conf.lmp (LAMMPS data格式)
# ---------------------------------------------------------------------------
tran_xyz2data() {
    python3 - "$1" "$2" "${@:3}" <<'PY'
import re,sys,tempfile; from pathlib import Path; from ase.io import read,write
if len(sys.argv)<4: print("Usage: tran_xyz2data model.xyz conf.lmp C [H O ...]"); sys.exit(1)
xyz_file=Path(sys.argv[1]); lmp_file=sys.argv[2]; elements=sys.argv[3:]
if not xyz_file.exists(): raise FileNotFoundError(f"File not found: {xyz_file}")
with open(xyz_file,"r") as f: lines=f.readlines()
if len(lines)<2: raise RuntimeError("xyz file too short")
lines[1]=re.sub(r"\blattice\s*=","Lattice=",lines[1],flags=re.IGNORECASE)
lines[1]=re.sub(r"\bproperties\s*=","Properties=",lines[1],flags=re.IGNORECASE)
with tempfile.NamedTemporaryFile(mode="w",suffix=".xyz",delete=False,encoding="utf-8") as f:
    f.writelines(lines); tmp_xyz=f.name
try:
    atoms=read(tmp_xyz,format="extxyz",index=0)
    write(lmp_file,atoms,format="lammps-data",atom_style="atomic",specorder=elements,masses=True)
finally: Path(tmp_xyz).unlink(missing_ok=True)
print(f"Written: {lmp_file}  |  Element order: {elements}")
PY
}


# ---------------------------------------------------------------------------
# Function: tran_xyz2pdb
# 功能: 将xyz文件转换为PDB格式（蛋白质/分子可视化常用格式）
# Usage: tran_xyz2pdb <file.xyz>
# Example:
#   tran_xyz2pdb model.xyz
#   # 输出: convered.pdb
# ---------------------------------------------------------------------------
tran_xyz2pdb(){
    python3 << EOF
import ase.io
filename='$1'; xyzinfo=ase.io.read(filename)
ase.io.write('convered.pdb',xyzinfo)
EOF
}


# ---------------------------------------------------------------------------
# Function: tran_cif2xyz
# 功能: 将CIF文件转换为xyz格式（extxyz，保留Lattice信息）
# 场景: 从Materials Project或COD数据库下载CIF后转换为计算用xyz格式。
# Usage: tran_cif2xyz <file.cif>
# Example:
#   tran_cif2xyz mp-2741_Si.cif
#   # 输出: mp-2741_Si.xyz
# ---------------------------------------------------------------------------
tran_cif2xyz(){
    python3 << EOF
import os,ase.io
filename='$1'; atoms=ase.io.read(filename)
atoms.info.clear()
keep={'numbers','positions'}
for k in list(atoms.arrays.keys()):
    if k not in keep: del atoms.arrays[k]
outname=os.path.splitext(os.path.basename(filename))[0]+'.xyz'
ase.io.write(outname,atoms,format='extxyz')
print(outname)
EOF
}


# ---------------------------------------------------------------------------
# Function: tran_xyz2cif / tran_pos2cif
# 功能: 将xyz或其他格式转换为CIF格式
#       通用转换，支持ase能读取的所有格式。
# 场景: 需要提交晶体结构到数据库或使用晶体学软件时转为CIF。
# Usage: tran_xyz2cif <input_file> [output.cif]
# Example:
#   tran_xyz2cif model.xyz output.cif
#   tran_xyz2cif POSCAR si.cif
# ---------------------------------------------------------------------------
tran_xyz2cif (){
    cifname=${2:-'file.cif'}
    python3  <<'EOF'
from ase.io import read, write
import sys
def convert(file_in,file_out):
    atoms=read(file_in); write(file_out,atoms)
infile='$1'; outfile='$cifname'
convert(infile,outfile)
EOF
}

# 该文件用于处理输入输出文件
# 进行文件格式间的转转换

#调用ovito转换文件的格式
tran_xyz2pos() {
    python3 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'POSCAR_convered', 'vasp')"
}

tran_xyz2pos() {
    python3 << EOF
import ase.io
filename = '$1'
xyzinfo = ase.io.read(filename)
ase.io.write('POSCAR_convered', xyzinfo, format='vasp')
EOF
}

tran_pos2xyz() {     
       	python3 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'model_conversed.xyz', 'xyz',columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])"
}

tran_data2xyz(){
    python3 -c "from ovito.io import import_file, export_file; pipeline = import_file('$1'); export_file(pipeline, 'model_conversed.xyz', 'xyz',columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])"
}

tran_xyz2pdb(){
    python3 << EOF
import ase.io
filename = '$1'
xyzinfo = ase.io.read(filename)
ase.io.write('convered.pdb', xyzinfo)
EOF
}

tran_data2pdb(){
    data_to_xyz $1
    xyz_to_pdb model_conversed.xyz
}

tran_pos2xyz(){
    python3 << EOF
import ase.io
filename = '$1'
xyzinfo = ase.io.read(filename)
ase.io.write('model_conversed.xyz', xyzinfo, format='xyz')
EOF
}

tran_cif2xyz(){
    python3 << EOF
import os
import ase.io

filename = '$1'
atoms = ase.io.read(filename)

atoms.info.clear()

keep = {'numbers', 'positions'}
for k in list(atoms.arrays.keys()):
    if k not in keep:
        del atoms.arrays[k]

# extxyz才有 Lattice
outname = os.path.splitext(os.path.basename(filename))[0] + '.xyz'
ase.io.write(outname, atoms, format='extxyz')
print(outname)
EOF
}


tran_xyz2cif (){ 
    #其实不只是xyz到cif，很多其他的格式也可以到达cif的格式
    #不过注意xyz文件如果想要使用ovito来进行转化的时候需要注意，最好只输出xyz方向的坐标，多了其他的东西可能会无法识别。
    cifname={2:-'file.cif'}
    python3  <<'EOF'
from ase.io import read, write

def convert(file_in, file_out):
    # 读取XYZ文件
    atoms = read(file_in)
    # 输出为CIF文件
    write(file_out, atoms)
infile = '$1'
outfile = '$cifname
convert(infile, outfile)
EOF
}
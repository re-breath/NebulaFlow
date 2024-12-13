# 转化xyz.in类型的二维材料，该类型的文件比较古老，使用需要谨慎
import numpy as np
lines = []
file_name = "./xyz10.in"
with open(file_name, 'r') as f:
    lines = f.readlines()
    #print("lines",lines[2:10])

length, width = map(float, lines[1].split()[3:5])
lattice_line = f'Lattice="{length} 0.0 0.0 0.0 {width} 0.0 0.0 0.0 3.35" Properties=species:S:1:pos:R:3:group:I:1 pbc="T T F"\n'
lines[0]=lines[0].split()[0]+'\n'

lines[1] = lattice_line
for i in range(2, len(lines)):
    lines[i] = 'C ' + str(lines[i].split()[3:5][0]) + ' '+str(lines[i].split()[3:5][1]) +' 1.675 ' +str(lines[i].split()[1]) + ' \n'
with open(f'{file_name}.xyz', 'w') as f:
    f.writelines(lines)


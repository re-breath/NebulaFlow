import dpdata
import random
import numpy as np
import os
import shutil

for i in range(1, 51):
    dirname = f'train-{i}'
    if not os.path.exists(dirname):
        os.makedirs(dirname)

source_files = ['INCAR', 'POTCAR', 'KPOINTS', 'POSCAR']
for file in source_files:
    for i in range(1, 51):
        src = file
        dst = f'./train-{i}/{file}'
        shutil.copyfile(src, dst)

directory = os.getcwd()

perturbed_system = dpdata.System('POSCAR').perturb(pert_num=50,
                                                    cell_pert_fraction=0.02,
                                                    atom_pert_distance=0.25,
                                                    atom_pert_style='normal')

for i in range(50):

    poscar_filename = f'POSCAR{i}'

    perturbed_system.to_vasp_poscar(poscar_filename, frame_idx=i)

    shutil.move(os.path.join(directory, poscar_filename), f'./train-{i+1}/POSCAR')

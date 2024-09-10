from ase.io import read, write
import numpy as np
import sys
import os
import shutil
from ovito.io import import_file, export_file
from ovito.modifiers import WrapPeriodicImagesModifier

def create_train_folders():
    if not os.path.exists('train_folders'):
        os.makedirs('train_folders')

    for i in range(1, 21):
        folder_name = f'train-{i}'
        folder_path = os.path.join('train_folders', folder_name)
        os.makedirs(folder_path)

def split_train_xyz():
    with open('dump.xyz', 'r') as file:
        lines = file.readlines()

    group_size = 34
    num_groups = len(lines) // group_size

    for i in range(num_groups):
        start_index = i * group_size
        end_index = start_index + group_size

        group_lines = lines[start_index:end_index]
        group_filename = f'train_group{i + 1}.xyz'
        group_filepath = os.path.join('train_folders', f'train-{i + 1}', group_filename)

        with open(group_filepath, 'w') as group_file:
            group_file.writelines(group_lines)

def convert_train_xyz_to_poscar():
    train_folder_path = os.path.join(os.getcwd(), 'train_folders')
    folders = os.listdir(train_folder_path)

    for folder_name in folders:
        folder_path = os.path.join(train_folder_path, folder_name)
        xyz_files = [f for f in os.listdir(folder_path) if f.endswith('.xyz')]

        for xyz_file in xyz_files:
            xyz_filepath = os.path.join(folder_path, xyz_file)
            poscar_filepath = os.path.join(folder_path, xyz_file.replace('.xyz', '.vasp'))

            # Import the XYZ file using ovito
            pipeline = import_file(xyz_filepath)

            # Use the WrapPeriodicImagesModifier to wrap periodic images
            pipeline.modifiers.append(WrapPeriodicImagesModifier())


            export_file(pipeline, poscar_filepath, 'vasp')


            new_poscar_filepath = os.path.join(folder_path, 'POSCAR')
            os.rename(poscar_filepath, new_poscar_filepath)

create_train_folders()
split_train_xyz()
convert_train_xyz_to_poscar()

source_files = ['INCAR', 'POTCAR', 'KPOINTS']
for file in source_files:
    for j in range(1, 21):
        src = file
        folder_name = f'train-{j}'
        folder_path = os.path.join('train_folders', folder_name)
        dst = os.path.join(folder_path, os.path.basename(file))
        shutil.copyfile(src, dst)

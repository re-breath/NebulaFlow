import dpdata

perturbed_system = dpdata.System('./CONTCAR').perturb(pert_num=100, 
    cell_pert_fraction=0.03, 
    atom_pert_distance=0.1, 
    atom_pert_style='normal')
print(perturbed_system.data)

for i in range(100):
    perturbed_system.to_vasp_poscar('POSCAR%s'%i, frame_idx=i)
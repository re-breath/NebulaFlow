import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import AffineTransformationModifier, WrapPeriodicImagesModifier

in_file  = 'optimized_molecule.data'
out_file = 'shift.data'
extra_shift = np.array([10, 10, 10])

pipeline = import_file(in_file)
data = pipeline.compute()

# 1) 盒子中心
box_center = 0.5 * np.sum(data.cell[:, :3], axis=0)

# 2) 分子质心（假设整张 data 就是一个分子）
com = np.mean(data.particles.positions, axis=0)

# 3) 总平移量 = 盒中心 - 质心 + 用户微调
shift = box_center - com + extra_shift

# 4) 应用平移 + wrap
pipeline.modifiers.append(
    AffineTransformationModifier(translation=shift)
)
pipeline.modifiers.append(
    WrapPeriodicImagesModifier()
)

export_file(pipeline, out_file, 'lammps/data', atom_style='full')
print(f'Moved & wrapped -> {out_file}')
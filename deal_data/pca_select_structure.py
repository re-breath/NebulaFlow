"""
Select structures
===========================

This example shows how to select structures from dataset
"""
from pynep.calculate import NEP
from pynep.select import FarthestPointSample
from ase.io import read, write
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


a = read('train.xyz', ':')
calc = NEP("nep-5.txt")
print(calc)
des = np.array([np.mean(calc.get_property('descriptor', i), axis=0) for i in a])
sampler = FarthestPointSample(min_distance=0.006)
#selected_i = sampler.select(des, [])
#write('test.xyz', [a[i] for  i in selected_i])


selected_i = sampler.select(des, [])
selected_structures = [a[i] for i in selected_i]
write('elect.xyz', selected_structures)

# 新增代码，将未选中的构型写入另一个文件
all_indices = set(range(len(a)))
unselected_i = list(all_indices - set(selected_i))
unselected_structures = [a[i] for i in unselected_i]
write('unselected.xyz', unselected_structures)


reducer = PCA(n_components=2)
reducer.fit(des)
proj = reducer.transform(des)
plt.scatter(proj[:,0], proj[:,1], label='all data')
selected_proj = reducer.transform(np.array([des[i] for i in selected_i]))
plt.scatter(selected_proj[:,0], selected_proj[:,1], label='selected data')
plt.legend()
plt.axis('off')
plt.savefig('select.png')

